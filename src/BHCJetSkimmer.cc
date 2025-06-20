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
	else if(_strategy == NlnNonAK4)
		cout << " NlnN (Delauney) with reco AK4 rechits" << endl;
	else
		cout << " undefined. Please use SetStrategy(i) with i == 0 (NlnN), 1 (N2)" << endl;

	//event selection
	cout << "Using event selection";
	if(_sel == boostW)
		cout << " boosted Ws (fully hadronic)" << endl;
	else if(_sel == boostTop)
		cout << " boosted tops (fully hadronic)" << endl;
	else if(_sel == QCDdijets)
		cout << " QCD dijets" << endl;
	else
		cout << " default (2+ gen partons that are not tops, including gluons)" << endl;
	
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	//cout << "oname " << _oname << endl;	
	MakeProcCats(_oname, true);

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
		//event level selection
		//at least 1 gen jet
		_base->GetEntry(i);
		if(i % (SKIP) == 0) cout << "evt: " << i << " of " << _nEvts << endl;
	
		int ngenpart = _base->genpart_ngenpart;

		if(_sel == boostW || _sel == boostTop){
			//reject fully leptonic and semi-leptonic W decays for ttbar - later may want to turn off to look at just displaced b decays
			//only look at hadronic W decays
			int nW_lep = 0;
			for(int g = 0; g < ngenpart; g++){
				int genmomidx = _base->genpart_momIdx->at(g);
				if(genmomidx == -1) continue;
				//skip particles that didn't come from a W
				if(fabs(_base->genpart_id->at(genmomidx)) != 24) continue;
				int genid = fabs(_base->genpart_id->at(g));
				//only check for !neutrinos
				if(genid == 11 || genid == 13 || genid == 15) nW_lep++; 
			}
			if(i % SKIP == 0) cout << " has " << nW_lep << " leptonic Ws" << endl;
			//if any Ws in event decay leptonically, don't count
			if(nW_lep > 0) continue;
			//if all Ws in event decay leptonically, don't count
			//if(nW - nW_lep == 0) continue;
		}

		int nW = 0;
		//do event selection here based on enum
		if(_sel == boostW){
			int totW = 0;
			//at least 1 W quark with pt, E requirements
			for(int g = 0; g < ngenpart; g++){
				if(fabs(_base->genpart_id->at(g)) == 6){
					if(i % SKIP == 0) cout << " has top with pt " << _base->genpart_pt->at(g) << " and energy " << _base->genpart_energy->at(g) << " and id " << _base->genpart_id->at(g) << " and pz " << _base->genpart_pz->at(g) << endl;
				}
				if(fabs(_base->genpart_id->at(g)) == 5){
					if(i % SKIP == 0) cout << " has b quark with pt " << _base->genpart_pt->at(g) << " and energy " << _base->genpart_energy->at(g) << " and id " << _base->genpart_id->at(g) << " and pz " << _base->genpart_pz->at(g) << " eta " << _base->genpart_eta->at(g) << " and phi " << _base->genpart_phi->at(g) << endl;
				}
				if(fabs(_base->genpart_id->at(g)) != 24) continue;
				totW++;
				if(_base->genpart_pt->at(g) < _minWPt) continue;
				if(i % SKIP == 0) cout << " has W with pt " << _base->genpart_pt->at(g) << " and energy " << _base->genpart_energy->at(g) << " and id " << _base->genpart_id->at(g) << " and pz " << _base->genpart_pz->at(g) << endl;
				nW++;
			}	
			//at least 1 W	
			if(nW < 1){ 
				if(i % SKIP == 0) cout << " has no Ws that pass pt > " << _minWPt << endl;
				continue;
			}
			//check # b's - should match # Ws
			int nb = count(_base->genpart_id->begin(), _base->genpart_id->end(), 5);
			nb += count(_base->genpart_id->begin(), _base->genpart_id->end(), -5);
			if(i % SKIP == 0) cout << " has " << nb << " gen b's and " << totW << " total gen W's" << endl; 
			if(nb != nW) continue;


		}
		else if(_sel == boostTop){
			//at least 1 top quark with pt, E requirements
			int ntop = 0;
			for(int g = 0; g < ngenpart; g++){
				if(fabs(_base->genpart_id->at(g)) != 6) continue;
				//if(_base->genpart_energy->at(g) < _minTopE) continue;
				if(_base->genpart_pt->at(g) < _minTopPt) continue;
				if(i % SKIP == 0) cout << " has top with pt " << _base->genpart_pt->at(g) << " and energy " << _base->genpart_energy->at(g) << " and id " << _base->genpart_id->at(g) << " and pz " << _base->genpart_pz->at(g) << endl;
				ntop++;
			}	

			if(ntop < 1){
				if(i % SKIP == 0) cout << " has no tops that pass pt > " << _minTopPt << endl;
				continue;
			}


		}
		else if(_sel == QCDdijets){
			//at least two gen partons to be reconstructed as jets in event (ie saved gen partons)
			int nparton = 0;
			vector<int> p_ids = {1,2,3,4,5,21};
			for(int g = 0; g < ngenpart; g++){
				if(find(p_ids.begin(), p_ids.end(), fabs(_base->genpart_id->at(g))) == p_ids.end()) continue;
				nparton++;
			}	
			if(nparton < 2) continue;

		}
		//default selection
		else{
			//at least two gen partons to be reconstructed as jets in event (ie saved gen partons)
			int nparton = 0;
			vector<int> p_ids = {1,2,3,4,5};
			for(int g = 0; g < ngenpart; g++){
				if(find(p_ids.begin(), p_ids.end(), fabs(_base->genpart_id->at(g))) == p_ids.end()) continue;
				nparton++;
			}	
			if(nparton < 2) continue;

		}
	


		_prod->GetGenJets(_genAK4jets, _genAK8jets, _genAK15jets, i);
		//if(_genjets.size() < 1){ cout << endl; continue; }
		_prod->GetRecoJets(_recoAK4jets, _recoAK8jets, _recoAK15jets, i);
		//if(_recoAK4jets.size() < 1){ cout << endl; continue; }
		_prod->GetGenParticles(_genparts, i);
		if(i % SKIP == 0) cout << " has " << _recoAK4jets.size() << " AK4 reco jets and " << _genAK4jets.size() << " AK4 gen jets" << endl;
		if(i % SKIP == 0) cout << " has " << _recoAK8jets.size() << " AK8 reco jets and " << _genAK8jets.size() << " AK8 gen jets" << endl;
		if(i % SKIP == 0) cout << " has " << _recoAK15jets.size() << " AK15 reco jets and " << _genAK15jets.size() << " AK15 gen jets" << endl;
		
		FillGenHists();
		
		//get PV info
		_pvx = _base->PV_x;
		_pvy = _base->PV_y;
		_pvz = _base->PV_z;

		int nsubcls_tot = 0;
		///do GMM only option
		for(int j = 0; j < _recoAK4jets.size(); j++){
			 _recoAK4jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			cout << "SubClustering reco AK4 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
			algo = new BayesCluster(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			if(_smear) algo->SetDataSmear(smear);
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
			algo->SetNGhosts(_nGhosts);			

			t = clock();
			GaussianMixture* gmm = algo->SubCluster();
			t = clock() - t;
			y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			_recoAK4jets[j].SetModel(gmm, _gev);
			
			nsubcls_tot += gmm->GetNClusters();
			cout << " jet has " << gmm->GetNClusters() << " subclusters" << endl;// with parameters" << endl;
			for(int k = 0; k < gmm->GetNClusters(); k++){
				auto params = gmm->GetDataStatistics(k);
				Matrix mean = params["mean"];
				Matrix cov = params["cov"];
				//cout << "cluster #" << k << endl;
				//cout << "mean" << endl; mean.Print();
				//cout << "cov" << endl; cov.Print();
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
		for(int p = 0; p < _procCats.size(); p++){
			_procCats[p].hists1D[0][141]->Fill(nsubcls_tot);
			_procCats[p].hists2D[0][39]->Fill(nsubcls_tot, (int)_recoAK4jets.size());
		}
		//do GMM for AK8
		for(int j = 0; j < _recoAK8jets.size(); j++){
			 _recoAK8jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			cout << "SubClustering reco AK8 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
			algo = new BayesCluster(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			if(_smear) algo->SetDataSmear(smear);
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
			algo->SetNGhosts(_nGhosts);			
			
			t = clock();
			GaussianMixture* gmm = algo->SubCluster();
			t = clock() - t;
			y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			_recoAK8jets[j].SetModel(gmm, _gev);
			
			cout << " jet has " << gmm->GetNClusters() << " subclusters" << endl;// with parameters" << endl;
			rhs.clear();
		}

		//do GMM for AK15
		for(int j = 0; j < _recoAK15jets.size(); j++){
			 _recoAK15jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			cout << "SubClustering reco AK15 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
			algo = new BayesCluster(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			if(_smear) algo->SetDataSmear(smear);
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
			algo->SetNGhosts(_nGhosts);			
			
			t = clock();
			GaussianMixture* gmm = algo->SubCluster();
			t = clock() - t;
			y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			_recoAK15jets[j].SetModel(gmm, _gev);
			
			cout << " jet has " << gmm->GetNClusters() << " subclusters" << endl;// with parameters" << endl;
			rhs.clear();
		}



		FillRecoJetHists();
		//only does above
		if(_strategy == gmmOnly){
			continue;
		}
		if(_strategy == NlnNonAK4){
			//use rhs from reco ak4 jets
			rhs.clear();
			for(int j = 0; j < _recoAK4jets.size(); j++){
				vector<Jet> jet_rhs; 
				_recoAK4jets[j].GetJets(jet_rhs);
				for(auto rh : jet_rhs) rhs.push_back(rh);
			}
			

		}
		else{
			//get all rhs in event
			_prod->GetRecHits(rhs, i);
		}
		for(auto rh : rhs){
			_procCats[1].hists1D[0][131]->Fill(rh.t());
		}
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
		
		//FillResolutionHists(); - does gen matching...again?
		if(i % SKIP == 0) cout << " and " << rhs.size() << " total rhs" << endl;
		cout << "Clustering..." << endl;	
		BayesCluster* algo = new BayesCluster(rhs);
		if(_smear) algo->SetDataSmear(smear);
		algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
		algo->SetThresh(_thresh);
		algo->SetAlpha(_alpha);
		algo->SetSubclusterAlpha(_emAlpha);
		algo->SetVerbosity(_verb);
		algo->SetPriorParameters(_prior_params);
		algo->CheckMerges(_check_merges);
		algo->SetNGhosts(_nGhosts);		

		//run clustering
		//delauney NlnN version
		if(_strategy == NlnN || _strategy == NlnNonAK4){
			//start clock
			t = clock();
			_trees = algo->NlnNCluster();
		}
		//N^2 version
		else if(_strategy == N2){
			//start clock
			t = clock();
			_trees = algo->N2Cluster();
		}
		t = clock() - t;
		y_time.push_back((double)t/CLOCKS_PER_SEC);
		//cout <<  "y time entry " << y_time[y_time.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
		comptime->Fill((double)t/CLOCKS_PER_SEC);

		cout << _trees.size() << " trees" << endl;
		//transform trees (nodes) to jets
		TreesToJets();
		
		//do jet-based evt selection for BHC jets
		if(_sel == boostW){
			//need at least nW*2 jets (ie 1 nW = 1 W jet and 1 b jet)
			if(_predJets.size() < nW*2){
				cout << "have " << _predJets.size() << " need at least " << nW*2 << " jets to match to " << nW << " W(s) and an equal number of b(s)" << endl;
				continue;
			}
			//match subleading jets to gen b's using dR and Eratio - jets are sorted in TreesToJets
			vector<Jet> subjets, leadjets;
			//if 1 W, expect 1 b; 2 Ws = 2 b;s
			if(nW == 1){
				leadjets.push_back(_predJets[0]);
				for(int j = 1; j < _predJets.size(); j++)
					subjets.push_back(_predJets[j]);
			}
			else{
				leadjets.push_back(_predJets[0]);
				leadjets.push_back(_predJets[1]);
				for(int j = 2; j < _predJets.size(); j++)
					subjets.push_back(_predJets[j]);
			}
			vector<int> genbMatchIdxs(subjets.size(),-1);
			GenericMatchJet(subjets,_genparts, genbMatchIdxs, 5); //dr match BHC jets to gen Ws
			//want subleading jets to be very good matches to their respective gen b quarks (these are the "tags")
			//good match = small dr, Eratio close to 1
			bool drGood = true;
			bool EratioGood = true;
			vector<bool> bMatch(subjets.size(), 0);
			//cout << "BHC jet matches to b-jets are " << endl;
			for(int j = 0; j < subjets.size(); j++){
				int genbidx = genbMatchIdxs[j];
				if(genbidx == -1) continue;
				//check dr
				double dr = dR(subjets[j].eta(), subjets[j].phi(), _genparts[genbidx].eta(), _genparts[genbidx].phi());
				//check Eratio
				double Eratio = subjets[j].E()/_genparts[genbidx].E();
				cout << "sublead jet #" << j << " has eta " << subjets[j].eta() << " phi " << subjets[j].phi() << " and E " << subjets[j].E() << endl;
				cout << " gen-matched to b #" << genbidx << "  with eta " << _genparts[genbidx].eta() << " phi " << _genparts[genbidx].phi() << " and E " << _genparts[genbidx].E() << endl; 
				if(dr > 0.1) drGood = false;
				if(Eratio < 0.5 || Eratio > 1.5) EratioGood = false;
				cout << " gen match dr " << dr << " Eratio " << Eratio << " match criteria " << (drGood && EratioGood) << endl;
				bMatch[j] = drGood && EratioGood;
			}
			//need at least 2 subleading jets to pass
			int nPass = 0;
			for(auto b : bMatch){
				if(b) nPass++;
			}

			//could also add requirement that top mass is made with best W+b jet combinations
			//if ALL good b-jet matches, remove from "predJets" list s.t. predJets are only W candidate jets
			//then gen-match remaining jet(s) to Ws (ie the "probes") - done in FillPredJets
			if(nPass > 1){
				cout << "good b-jet matches with BHC jets - continue with boosted W selection" << endl;
				_predJets = leadjets;
			}
			//else skip
			else{
				cout << "BHC jets did not match gen b's well enough, skip" << endl;
				continue; //else, skip this event
			}
		}
		if(_sel == boostTop){
			if(_predJets.size() < 2) continue; //want at least two boosted tops in an event
			//if there are any jets that are good matches to bs then the top can't be boosted, skip
			vector<int> genbMatchIdxs(_predJets.size(),-1);
			GenericMatchJet(_predJets,_genparts, genbMatchIdxs, 5); //dr match BHC jets to gen Ws
			//want subleading jets to be very good matches to their respective gen b quarks (these are the "tags")
			//good match = small dr, Eratio close to 1
			bool drGood = true;
			bool EratioGood = true;
			vector<bool> bMatch(_predJets.size(), 0);
			for(int j = 0; j < _predJets.size(); j++){
				int genbidx = genbMatchIdxs[j];
				//check dr
				double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genparts[genbidx].eta(), _genparts[genbidx].phi());
				//check Eratio
				double Eratio = _predJets[j].E()/_genparts[genbidx].E();
				if(dr > 0.1) drGood = false;
				if(Eratio < 0.8 || Eratio > 1.2) EratioGood = false;
				bMatch[j] = drGood && EratioGood;
			}
			//if ANY good matches, skip	
			if(std::any_of(bMatch.begin(), bMatch.end(), [](bool v) { return v;})){
				continue;
			}

		}






		//fill pred jet hists with jets
		FillPredJetHists();
		nsubcls_tot = 0;
		for(int j = 0; j < _predJets.size(); j++){
			nsubcls_tot += _predJets[j].GetNConstituents();
		}
		for(int p = 0; p < _procCats.size(); p++){
			_procCats[p].hists1D[0][142]->Fill(nsubcls_tot);
			_procCats[p].hists2D[0][42]->Fill(nsubcls_tot, (int)_predJets.size());
		}
		
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

