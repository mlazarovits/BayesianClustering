#include "BHCJetSkimmer.hh"
#include "BayesCluster.hh"

void BHCJetSkimmer::Skim(){
	cout << "Writing skim to: " << _oname << endl;
	cout << "Using clustering strategy";
	if(_strategy == NlnN)
		cout << " NlnN (Delauney)" << endl;
	else if(_strategy == N2)
		cout << " N2 (naive)" << endl;
	else if(_strategy == gmmOnly)
		cout << " GMM only" << endl;
	else
		cout << " undefined. Please use SetStrategy(i) with i == 0 (NlnN), 1 (N2)" << endl;
	
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	//cout << "oname " << _oname << endl;	
	MakeProcCats(_oname, false);

	//cout << "n procs: " << _procCats.size() << endl;
	//for(auto proc : _procCats) cout << "proc: " << proc.plotName << endl;

	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = 2*acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(0.,2,2); 
	//for time smearing (energy dependent)
	double tres_c = 0.2;
	double tres_n = sqrt(1 - tres_c*tres_c)*_gev;	

	
	map<string, Matrix> params;
	vector<node*> trees;
	vector<Jet> rhs;


	_prod->PrintPreselection();
        
	//for computational time
	vector<double> x_nrhs, y_time;
	vector<double> x_nrhs_subcl, y_time_subcl;
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	int SKIP = 1;
	BayesCluster* algo = nullptr;
	clock_t t;
	for(int i = _evti; i < _evtj; i+=SKIP){
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		if(i % (SKIP) == 0) cout << "evt: " << i << " of " << _nEvts;
		_prod->GetRecoJets(_recojets, i);
		////fill gen jet histograms
		_prod->GetGenJets(_genjets, i);
		if(_genjets.size() < 1 && _recojets.size() < 1){ cout << endl; continue; }

		if(i % SKIP == 0) cout << " with " << _recojets.size() << " reco jets and " << _genjets.size() << " gen jets" << endl;
		///do GMM only option
		for(int j = 0; j < _recojets.size(); j++){
			 _recojets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			cout << "SubClustering reco jet #" << j << "..." << endl;	
			algo = new BayesCluster(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresNoise*_gev, _tresStoch*_gev); 
			if(_smear) algo->SetDataSmear(smear);
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
		
			t = clock();
			GaussianMixture* gmm = algo->SubCluster();
			t = clock() - t;
			y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			//set constituents
			vector<double> norms;
			gmm->GetNorms(norms);
			for(int k = 0; k < gmm->GetNClusters(); k++){
				Jet subcl(gmm->GetModel(k), norms[k]/_gev, gmm->GetPi(k), BayesPoint({_pvx, _pvy, _pvz})); 
				_recojets[j].AddConstituent(subcl);
				auto params = gmm->GetDataStatistics(k);
				Matrix mean = params["mean"];
				Matrix cov = params["cov"];
				for(int p = 0; p < _procCats.size(); p++){
					_procCats[p].hists1D[0][98]->Fill(mean.at(0,0));
					_procCats[p].hists1D[0][99]->Fill(mean.at(1,0));
					_procCats[p].hists1D[0][100]->Fill(mean.at(2,0));
					_procCats[p].hists1D[0][101]->Fill(sqrt(cov.at(0,0)));
					_procCats[p].hists1D[0][102]->Fill(sqrt(cov.at(1,1)));
					_procCats[p].hists1D[0][103]->Fill(sqrt(cov.at(2,2)));
				}
			}
		
			rhs.clear();
		}
		FillRecoJetHists();
		//only does above
		if(_strategy == gmmOnly){
			cout << endl;
			continue;
		}
		_prod->GetRecHits(rhs, i);
		
		//safety
		if(rhs.size() < 1) continue;

		x_nrhs.push_back((double)rhs.size());
		//for(int r = 0; r < rhs.size(); r++){
		//	rhTime->Fill(rhs[r].t());
		//}

		//assume detector radius is constant and equal for all rhs (all rhs in event are recorded in same type of detector)
		//this should be true for all events
		vector<JetPoint> rh = rhs[0].GetJetPoints(); //only 1 rh
		_radius = sqrt(rh[0].x()*rh[0].x() + rh[0].y()*rh[0].y());	
		
		//decayId = 0 -> had
		//decayId = 1 -> lep
		//0 && 0 = 0 -> fully had
		//1 && 1 = 1 -> fully lep
		//0 && 1 = 2 -> semi lep 
		//if(_base->Top_decayId->size() > 1){
		//	if(_base->Top_decayId->at(0) == _base->Top_decayId->at(1))
		//		_topDecayType = _base->Top_decayId->at(0);
		//	else
		//		_topDecayType = 2;
		//}
		//else _topDecayType = -1;
	//cout << "\ntopDecayType " << _topDecayType << endl;
	//gen particles in event
	//for(int g = 0; g < _base->genpart_ngenpart; g++){
	//	cout << "gen particle # " << g << " id " << _base->genpart_id->at(g) << " momidx " << _base->genpart_momIdx->at(g) << " E " << _base->genpart_energy->at(g) << " eta " << _base->genpart_eta->at(g) << " phi " << _base->genpart_phi->at(g) << endl;
	//}
	////gen jets
	//cout << "n gen jets " << _base->Jet_genNJet << " " << _genjets.size() << endl;
	//for(int g = 0; g < _base->Jet_genNJet; g++){
	//	cout << "gen jet # " << g << " E " << _base->Jet_genEnergy->at(g) << " eta " << _base->Jet_genEta->at(g) << " phi " << _base->Jet_genPhi->at(g) << endl;
	//}
	

		FillResolutionHists();
		//get PV info
		_pvx = _base->PV_x;
		_pvy = _base->PV_y;
		_pvz = _base->PV_z;
		if(i % SKIP == 0) cout << " and " << rhs.size() << " total rhs" << endl;
		cout << "Clustering..." << endl;	
		BayesCluster* algo = new BayesCluster(rhs);
		if(_smear) algo->SetDataSmear(smear);
		algo->SetMeasErrParams(_cell, _tresCte, _tresNoise*_gev, _tresStoch*_gev); 
		algo->SetThresh(_thresh);
		algo->SetAlpha(_alpha);
		algo->SetSubclusterAlpha(_emAlpha);
		algo->SetVerbosity(_verb);
		algo->SetPriorParameters(_prior_params);
		//run clustering
		//delauney NlnN version
		if(_strategy == NlnN){
			//start clock
			t = clock();
			trees = algo->NlnNCluster();
		}
		//N^2 version
		else if(_strategy == N2){
			//start clock
			t = clock();
			trees = algo->N2Cluster();
		}
		t = clock() - t;
		y_time.push_back((double)t/CLOCKS_PER_SEC);
		//cout <<  "y time entry " << y_time[y_time.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
		comptime->Fill((double)t/CLOCKS_PER_SEC);
		//clean trees (remove mirror point or nullptrs)
		CleanTrees(trees);
		//transform trees (nodes) to jets
		TreesToJets();
		//fill model histograms with trees
		//for subclusters TODO - move the hists filled in FillModelHists to FillPredJetsHists
		//FillModelHists();	
		//fill pred jet hists with jets
		FillPredJetHists();
		cout << endl;
	}
	graphs[1] = new TGraph(x_nrhs_subcl.size(), &x_nrhs_subcl[0], &y_time_subcl[0]);
	graphs[1]->SetName("nrhs_comptime_subcl");
	graphs[1]->SetTitle("nrhs_comptime_subcl");
	graphs[1]->GetXaxis()->SetTitle("# rhs");
	graphs[1]->GetYaxis()->SetTitle("computational time GMM only (sec)");
	
	graphs[0] = new TGraph(x_nrhs.size(), &x_nrhs[0], &y_time[0]);
	graphs[0]->SetName("nrhs_comptime");
	graphs[0]->SetTitle("nrhs_comptime");
	graphs[0]->GetXaxis()->SetTitle("# rhs");
	graphs[0]->GetYaxis()->SetTitle("computational time (sec)");

	WriteOutput(ofile);

	cout << "Wrote skim to: " << _oname << endl;
}

