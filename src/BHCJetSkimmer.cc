#include "BHCJetSkimmer.hh"
#include "BayesCluster.hh"

void BHCJetSkimmer::Skim(){
	cout << "Writing skim to: " << _oname << endl;
	cout << "Using clustering strategy";
	if(_strategy == NlnN)
		cout << " NlnN (Delauney)" << endl;
	else if(_strategy == N2)
		cout << " N2 (naive)" << endl;
	else
		cout << " undefined. Please use SetStrategy(i) with i == 0 (NlnN), 1 (N2), 2 (MM)" << endl;
	
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");

	cout << "oname " << _oname << endl;	
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

	double alpha = 0.1;
	double emAlpha = 0.5;
	double thresh = 1.;
	
	map<string, Matrix> params;
	vector<node*> trees;
	vector<Jet> rhs;

        
	//for computational time
	vector<double> x_nrhs, y_time;
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	int SKIP = 1;
	for(int i = _evti; i < _evtj; i+=SKIP){
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		if(i % (SKIP) == 0) cout << "evt: " << i << " of " << _nEvts;
		_prod->GetRecHits(rhs, i);
		x_nrhs.push_back((double)rhs.size());
		//for(int r = 0; r < rhs.size(); r++){
		//	rhTime->Fill(rhs[r].t());
		//}
	
	

		////fill gen jet histograms
		vector<Jet> genjets;
		_prod->GetGenJets(genjets, i);
	//	if(jets.size() < 1){ cout << endl; continue; }

	
		//if(i % (SKIP) == 0) cout << " with " << jets.size() << " jets to cluster and " << _phos.size() << " photons";
		if(i % SKIP == 0) cout << " with " << rhs.size() << " rhs" << endl;

		cout << "Clustering..." << endl;	
		clock_t t;
		BayesCluster* algo = new BayesCluster(rhs);
		if(_smear) algo->SetDataSmear(smear);
		if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n);
		algo->SetThresh(thresh);
		algo->SetAlpha(alpha);
		algo->SetSubclusterAlpha(emAlpha);
		algo->SetVerbosity(_verb);
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
	
		comptime->Fill((double)t/CLOCKS_PER_SEC);	
		//fill model histograms with trees
		//transform trees to jets

		FillPredJetHists(trees);
	}
	graphs[0] = new TGraph(x_nrhs.size(), &x_nrhs[0], &y_time[0]);
	graphs[0]->SetName("nrhs_comptime");
	graphs[0]->SetTitle("nrhs_comptime");
	graphs[0]->GetXaxis()->SetTitle("# rhs");
	graphs[0]->GetYaxis()->SetTitle("computational time (sec)");

	WriteOutput(ofile);

	cout << "Wrote skim to: " << _oname << endl;
}

