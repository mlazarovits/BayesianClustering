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
	if(_sel == singW) //single W production
		cout << " single Ws (fully hadronic)" << endl;
	else if(_sel == boostTop) //ttbar production
		cout << " boosted tops (fully hadronic)" << endl;
	else if(_sel == QCDdijets) //self-explanatory
		cout << " QCD dijets" << endl;
	else
		cout << " default (2+ gen partons that are not tops, including gluons)" << endl;
	
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
	if(_evtj > _nEvts) _evtj = _nEvts;
	int SKIP = 1;
	BayesCluster* algo = nullptr;
	clock_t t;
	double phiwindow = acos(-1)/2;
	//for event display
	map<string, BayesPoint> plot_centers; //plot_centers[i][j] is plot center in dim j for plot i (see above for different plot idxs)
	map<string, BayesPoint> plot_widths;
	for(int i = _evti; i < _evtj; i+=SKIP){
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		//event level selection
		//at least 1 gen jet
		_base->GetEntry(i);
		if(i % (SKIP) == 0) cout << "evt: " << i  << " (ntuple: " << _base->event << ") " << " of " << _nEvts << endl;
	
		int ngenpart = _base->genpart_ngenpart;
		_prod->GetGenParticles(_genparts, i);
		
		map<int,int> topidx_decayidx;
		//do event selection here based on enum
		if(_sel == singW){ //single W production
			//if any leptonic W decays, continue
			int nWs = _base->W_decayId->size();
			if(nWs < 1) continue;
			bool lepW = _base->W_decayId->at(0);
			//cout << "lepTop evt - " << lepTop << endl;
			//at least 1 W with pt, E requirements
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				if(fabs(_base->genpart_id->at(genidx)) != 24) continue;
				//min pt req - set by ptHatMin of sample
				//if(_genparts[g].pt() < _minWPt){
				//	continue;
				//}
				//make sure W has hadronic decay
				bool lep = _base->W_decayId->at(genidx);
				if(lep){ cout << "Fully leptonic W - skipping this W" << endl; continue;}
				if(i % SKIP == 0) cout << " has W with pt " << _base->genpart_pt->at(genidx) << " and energy " << _base->genpart_energy->at(genidx) << " and id " << _base->genpart_id->at(genidx) << " and pz " << _base->genpart_pz->at(genidx) << " eta " << _base->genpart_eta->at(genidx) << " phi " << _base->genpart_phi->at(genidx) << " decay id " << lep << endl;
			


				//get decay products
				for(int gg = 0; gg < _genparts.size(); gg++){
					if(gg == g) continue;
					int ggenidx = _genparts[gg].GetUserIdx();
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
				cout << "saving W daughter - id " << _base->genpart_id->at(ggenidx) << " eta " << _genparts[gg].eta() << " phi " << _genparts[gg].phi() << " energy " << _genparts[gg].e() << endl;
					//phi distribution debugging - only look at events w/ jets (ie quarks) in some window around 0 and 2pi
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
					_genq.push_back(_genparts[gg]);
				}
				//skip W if not enough q's are in phi window
				_genW.push_back(_genparts[g]);
			}
			
			//at least 1 W	
			if(_genW.size() < 1){ 
				if(i % SKIP == 0) cout << " has no hadronic Ws that pass pt > " << _minWPt << endl;
				continue;
			}
			//if W+gluon, save gen gluon as well as gen W (save on same level as W for consistent gen-matching)
			if(_oname.find("Wgluon") != string::npos){
				for(int g = 0; g < _genparts.size(); g++){
					int genidx = _genparts[g].GetUserIdx();
					if(fabs(_base->genpart_id->at(genidx)) == 21){
						if(i % SKIP == 0) cout << " has gluon with pt " << _base->genpart_pt->at(genidx) << " and energy " << _base->genpart_energy->at(genidx) << " and id " << _base->genpart_id->at(genidx) << " and pz " << _base->genpart_pz->at(genidx) << " eta " << _base->genpart_eta->at(genidx) << " phi " << _base->genpart_phi->at(genidx) << endl;
						_genglu.push_back(_genparts[g]);	
					}
				}
			}
		}
		else if(_sel == boostTop){
			//if any leptonic Top decays, continue
			int nTops = _base->Top_decayId->size();
			if(nTops < 1) continue;
			bool lepTop = _base->Top_decayId->at(0);
			//cout << "lepTop evt - " << lepTop << endl;
			//at least 1 W with pt, E requirements
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				if(fabs(_base->genpart_id->at(genidx)) != 6) continue;
				//make sure Top has hadronic decay
				bool lep = _base->Top_decayId->at(genidx);
				if(lep){ cout << "Fully leptonic Top - skipping this Top" << endl; continue;}
				if(i % SKIP == 0) cout << " has Top with pt " << _base->genpart_pt->at(genidx) << " and energy " << _base->genpart_energy->at(genidx) << " and id " << _base->genpart_id->at(genidx) << " and pz " << _base->genpart_pz->at(genidx) << " eta " << _base->genpart_eta->at(genidx) << " phi " << _base->genpart_phi->at(genidx) << " decay id " << lep << endl;
			

				//get decay products
				for(int gg = 0; gg < _genparts.size(); gg++){
					if(gg == g) continue;
					int ggenidx = _genparts[gg].GetUserIdx();
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
				cout << "saving Top daughter - id " << _base->genpart_id->at(ggenidx) << " eta " << _genparts[gg].eta() << " phi " << _genparts[gg].phi() << " energy " << _genparts[gg].e() << endl;
					//phi distribution debugging - only look at events w/ jets (ie quarks) in some window around 0 and 2pi
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
					_genq.push_back(_genparts[gg]);
				}
				_genTop.push_back(_genparts[g]);
			}
			
			//at least 1 Top	
			if(_genTop.size() < 1){ 
				if(i % SKIP == 0) cout << " has no hadronic Tops that pass pt > " << _minTopPt << endl;
				continue;
			}
		}
		else if(_sel == QCDdijets){
			//at least two gen partons to be reconstructed as jets in event (ie saved gen partons)
			vector<int> p_ids = {1,2,3,4,5,21};
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				if(find(p_ids.begin(), p_ids.end(), fabs(_base->genpart_id->at(genidx))) == p_ids.end()) continue;
				cout << "saving q/g - id " << _base->genpart_id->at(genidx) << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " energy " << _genparts[g].e() << endl;
				_genq.push_back(_genparts[g]);
			}
			if(_genq.size() < 2) continue;

		}
		//default selection
		else{
			/*
			//debugging phi selection
			int nphi_parts = 0;
			for(int g = 0; g < _genparts.size(); g++){
				double phi = _genparts[g].phi_02pi();
				int genidx = _genparts[g].GetUserIdx();
				if(phi > phiwindow && phi < 2*acos(-1) - phiwindow) continue;
				nphi_parts++;	
				cout << "counting particle - id " << _base->genpart_id->at(genidx) << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " energy " << _genparts[g].e() << endl;
			}
			if(nphi_parts < 1) continue; //skip if no gen parts in phi areas
			*/
			////at least two gen partons to be reconstructed as jets in event (ie saved gen partons)
			//int nparton = 0;
			//vector<int> p_ids = {1,2,3,4,5};
			//for(int g = 0; g < ngenpart; g++){
			//	if(find(p_ids.begin(), p_ids.end(), fabs(_base->genpart_id->at(g))) == p_ids.end()) continue;
			//	nparton++;
			//}	
			//if(nparton < 2) continue;

		}
	

		_prod->GetGenJets(_genAK4jets, _genAK8jets, _genAK15jets, i);
		//if(_genjets.size() < 1){ cout << endl; continue; }
		_prod->GetRecoJets(_recoAK4jets, _recoAK8jets, _recoAK15jets, i);
		//if(_recoAK4jets.size() < 1){ cout << endl; continue; }
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
			
			if(_strategy == gmmOnly) cout << "SubClustering reco AK4 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
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
			if(_strategy == gmmOnly) cout << " jet has " << gmm->GetNClusters() << " subclusters" << endl;// with parameters" << endl;
			for(int k = 0; k < gmm->GetNClusters(); k++){
				auto params = gmm->GetDataStatistics(k);
				Matrix mean = params["mean"];
				Matrix cov = params["cov"];
				//cout << "cluster #" << k << endl;
				//cout << "mean" << endl; mean.Print();
				//cout << "cov" << endl; cov.Print();
			}
		
			rhs.clear();
		}
		for(int p = 0; p < _procCats.size(); p++){
			_procCats[p].hists1D[0][84]->Fill(nsubcls_tot);
			_procCats[p].hists2D[0][27]->Fill(nsubcls_tot, (int)_recoAK4jets.size());
		}
		//do GMM for AK8
		for(int j = 0; j < _recoAK8jets.size(); j++){
			 _recoAK8jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			if(_strategy == gmmOnly) cout << "SubClustering reco AK8 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
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
			
			if(_strategy == gmmOnly) cout << " jet has " << gmm->GetNClusters() << " subclusters" << endl;// with parameters" << endl;
			rhs.clear();
		}

		//do GMM for AK15
		for(int j = 0; j < _recoAK15jets.size(); j++){
			 _recoAK15jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			if(_strategy == gmmOnly) cout << "SubClustering reco AK15 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
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
			
			if(_strategy == gmmOnly) cout << " jet has " << gmm->GetNClusters() << " subclusters" << endl;// with parameters" << endl;
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
			_procCats[1].hists1D[0][77]->Fill(rh.t());
		}
		//safety
		if(rhs.size() < 1) continue;
		x_nrhs.push_back((double)rhs.size());

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
		if(_sel == singW){
			//at least 1 BHC jet
			if(_predJets.size() < 1) continue;

		}


		if(_sel == boostTop){
			//at least 1 BHC jet
			if(_predJets.size() < 1) continue;
		}


		//fill event display if specified _evt2disp
		if(i == _evt2disp){
			cout << "Displaying event " << i << endl;
			//fill hists for overall event
			for(auto rh : rhs){
				double w;
				if(_evt2disp_z == 0)
					w = rh.E();
				else if(_evt2disp_z == 1){
					w = rh.t();
					_procCats[0].hists2D[0][129]->GetZaxis()->SetTitle("time [ns]");
					for(auto h : _evtdisps_obj){
						h->GetZaxis()->SetTitle("time [ns]");
					}
				}
				else
					w = rh.E(); //default energy weighted
				_procCats[0].hists2D[0][129]->Fill(rh.eta(), rh.phi(), w);

			}
			double eta = -999;
			double phi = -999;
			//save gen particles as tmarkers
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				int id = fabs(_base->genpart_id->at(genidx));
				string matchstring;
				int gidx = -1;
				//TODO: set to pretty colors with hex codes and make colors related ie for quarks (light vs b different but similar)
				if(id == 24){ //W
					
					//check if first w event display is already filled
					if(_evtdisps_obj[0]->GetEntries() == 0){
						//fill if not
						gidx = 0;
						matchstring = "W";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					else{ //fill next one
						gidx = 1;
						matchstring = "W2";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					_plot_particles.push_back(TMarker(eta, phi, kOpenTriangleUp));
					//_plot_particles[_plot_particles.size()-1].SetMarkerColor(kRed-4);
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);

				}
				if(id == 6){ //top
					//check if first top event display is already filled
					if(_evtdisps_obj[7]->GetEntries() == 0){
						//fill if not
						gidx = 7;
						matchstring = "top";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					else{ //fill next one
						gidx = 8;
						matchstring = "top2";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					_plot_particles.push_back(TMarker(eta, phi, kOpenSquare));
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kOrange-4);
				}
				else if(id == 1 || id == 2 || id == 3){ //light quarks
					//check if first q event display is already filled
					if(_evtdisps_obj[3]->GetEntries() == 0){
						//fill if not
						gidx = 3;
						matchstring = "q1";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					else{ //fill next one
						gidx = 4;
						matchstring = "q2";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					_plot_particles.push_back(TMarker(eta, phi, kPlus));
					//_plot_particles[_plot_particles.size()-1].SetMarkerColor(kAzure+4);
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);
				}
				else if(id == 5){ //b quark
					//check if first b event display is already filled
					if(_evtdisps_obj[5]->GetEntries() == 0){
						//fill if not
						gidx = 5;
						matchstring = "b";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					else{ //fill next one
						gidx = 6;
						matchstring = "b2";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					_plot_particles.push_back(TMarker(eta, phi, kMultiply));
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kCyan-4);
				}
				else if(id == 21){ //gluon
					gidx = 2;
					matchstring = "gluon";
					eta = _genparts[g].eta();
					phi = _genparts[g].phi_02pi();
					_plot_particles.push_back(TMarker(eta, phi, kStar));
					//_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlue-4);
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);
				}
				if(gidx == -1) continue; //skip if no gen particles specified
				plot_centers[matchstring] = BayesPoint({eta, phi}); //gen particle in question should be centered at (0,0) in local eta, phi coords for local (not global) evt disp
				//set by deltaR \approx 2*m/pT, taking deltaEta = deltaPhi
				double deta = 4*_genparts[g].m()/(_genparts[g].pt());
				double dphi = deta;
				plot_widths[matchstring] = BayesPoint({deta,dphi});
				if(_evt2disp_z == 1){ //update labels to time
					_evtdisps_obj[gidx]->GetZaxis()->SetTitle("time [ns]");
				}
				//fill hists for this gen particle (hist idx) - this is done in WriteEventDisplay() now
				/*
				PointCollection rh_pts;
				for(auto rh : rhs){
					double w;
					if(_evt2disp_z == 0)
						w = rh.E();
					else if(_evt2disp_z == 1){
						w = rh.t();
					}
					else
						w = rh.E(); //default energy weighted
					//center according to main gen particle
					BayesPoint rh_pt({rh.eta(), rh.phi()}); //save as BayesPoint to do correct circular translation to (0,0)
					rh_pt.SetWeight(w);
					rh_pts += rh_pt;	
				}
cout << "matchstring " << matchstring << " gidx " << gidx << " eta " << eta << " phi " << phi << " global center " << endl; plot_centers[matchstring].Print();
				//translate into local eta, phi coords
				rh_pts.Translate(plot_centers[matchstring].at(0),0);
				rh_pts.CircularTranslate(plot_centers[matchstring].at(1),1);
				for(int r = 0; r < rh_pts.GetNPoints(); r++){
					_evtdisps_obj[gidx]->Fill(rh_pts.at(r).at(0), rh_pts.at(r).at(1), rh_pts.at(r).w());
				}
				*/	
				
			}
			if(eta == -999 && phi == -999) continue; //no gen particles specified
cout << "eta " << eta << " phi " << phi << endl;
			//save BHC jets as ellipses
			for(int j = 0; j < _predJets.size(); j++){
				TEllipse el = PlotEll(_predJets[j]);
				el.SetLineColor(j+2);	
				el.SetFillStyle(0);	
				_ellipses.push_back(el);	
				_jellipses[j] = el;			
	
				TMarker m(_predJets[j].eta(), _predJets[j].phi(), kOpenStar);
				m.SetMarkerColor(j+2);
				//plot jet centers
				//_plot_particles.push_back(m); //30 = open star may need to change if not rendering	
				_jcenters[j] = m;			
	
				int nk = _predJets[j].GetNConstituents();
				for(int k = 0; k < nk; k++){
					Jet subcl = _predJets[j].GetConstituent(k);
					TEllipse sub_el = PlotEll(subcl);
					sub_el.SetLineColor(j+2);	
					sub_el.SetFillStyle(0);	
					sub_el.SetLineStyle(9);	
					//_plot jet centers
					TMarker sub_m(subcl.eta(), subcl.phi(), kOpenDiamond); 
					sub_m.SetMarkerColor(j+2);
					_subclellipses[j][k] = sub_el;
					_subclcenters[j][k] = sub_m;
	
				}
			}
		}

		//fill pred jet hists with jets
		FillPredJetHists();
		nsubcls_tot = 0;
		for(int j = 0; j < _predJets.size(); j++){
			nsubcls_tot += _predJets[j].GetNConstituents();
		}
		_procCats[1].hists2D[0][30]->Fill(nsubcls_tot, (int)_predJets.size());
		
		cout << endl;
	
		_genTop.clear();
		_genW.clear();
		_genb.clear();
		_genq.clear();
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

	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	ofile->cd();
	WriteOutput(ofile, plot_centers, plot_widths);
	ofile->Close();
	_infile->Close();	
	cout << "Wrote skim to: " << _oname << endl;
}

