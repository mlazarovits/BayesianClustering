#include "BHCJetSkimmer.hh"
#include "BayesCluster.hh"


void BHCJetSkimmer::Skim(){
	InitMapTree();
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
	clock_t t;
	double phiwindow = acos(-1)/2;
	//for event display
	map<string, BayesPoint> plot_centers; //plot_centers[i][j] is plot center in dim j for plot i (see above for different plot idxs)
	map<string, BayesPoint> plot_widths;
	vector<TCanvas*> jetcvs;
	for(int i = _evti; i < _evtj; i+=SKIP){
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		//event level selection
		//at least 1 gen jet
		_base->GetEntry(i);
		_obs.at("evt") = (double)i;
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
					_genq.push_back(_genparts[gg]);
				}
				//skip W if not enough q's are in phi window
				cout << "saving W - id " << _base->genpart_id->at(genidx) << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " energy " << _genparts[g].e() << endl;
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
						//if(i % SKIP == 0) cout << " has gluon with pt " << _base->genpart_pt->at(genidx) << " and energy " << _base->genpart_energy->at(genidx) << " and id " << _base->genpart_id->at(genidx) << " and pz " << _base->genpart_pz->at(genidx) << " eta " << _base->genpart_eta->at(genidx) << " phi " << _base->genpart_phi->at(genidx) << endl;
						cout << "saving gluon - id " << _base->genpart_id->at(genidx) << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " energy " << _genparts[g].e() << endl;
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
				bool lep = _base->Top_decayId->at(_genTop.size());
				if(lep){ cout << "Fully leptonic Top - skipping this Top" << endl; continue;}
				if(i % SKIP == 0) cout << " has Top with pt " << _base->genpart_pt->at(genidx) << " and energy " << _base->genpart_energy->at(genidx) << " and id " << _base->genpart_id->at(genidx) << " and pz " << _base->genpart_pz->at(genidx) << " eta " << _base->genpart_eta->at(genidx) << " phi " << _base->genpart_phi->at(genidx) << " decay id " << lep << endl;
				cout << "saving top - id " << _base->genpart_id->at(genidx) << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " energy " << _genparts[g].e() << endl;
			

				//get decay products
				int Widx = -1;
				for(int gg = 0; gg < _genparts.size(); gg++){
					if(gg == g) continue;
					int ggenidx = _genparts[gg].GetUserIdx();
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
				cout << "saving Top daughter - id " << _base->genpart_id->at(ggenidx) << " eta " << _genparts[gg].eta() << " phi " << _genparts[gg].phi() << " energy " << _genparts[gg].e() << endl;
					if(fabs(_base->genpart_id->at(ggenidx)) == 5)
						_genb.push_back(_genparts[gg]);
					if(fabs(_base->genpart_id->at(ggenidx)) == 24){
						Widx = ggenidx;
						_genW.push_back(_genparts[gg]);
					}
						
				}
				//get W daughters
				if(Widx != -1){
					for(int gg = 0; gg < _genparts.size(); gg++){
						if(gg == g) continue;
						if(gg == Widx) continue;
						int ggenidx = _genparts[gg].GetUserIdx();
						if(_base->genpart_momIdx->at(ggenidx) != Widx) continue;
					cout << "saving Top granddaughter - id " << _base->genpart_id->at(ggenidx) << " eta " << _genparts[gg].eta() << " phi " << _genparts[gg].phi() << " energy " << _genparts[gg].e() << endl;
						_genq.push_back(_genparts[gg]);		
					}

				}
				_genTop.push_back(_genparts[g]);
			}
			
			//at least 1 Top	
			if(_genTop.size() < 1){ 
				if(i % SKIP == 0) cout << " has no hadronic Tops" << endl;
				continue;
			}

		}
		else if(_sel == QCDdijets){
			//at least two gen partons to be reconstructed as jets in event (ie saved gen partons)
			vector<int> p_ids = {1,2,3,4};
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				if(find(p_ids.begin(), p_ids.end(), fabs(_base->genpart_id->at(genidx))) == p_ids.end()) continue;
				cout << "saving q - id " << _base->genpart_id->at(genidx) << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " energy " << _genparts[g].e() << " momidx " << _base->genpart_momIdx->at(genidx) << endl;
				_genq.push_back(_genparts[g]);
			}
			if(_genq.size() < 2) continue;

		}
		//default selection
		else{
		}
		
		FillGenParticleHists();

		_prod->GetGenJets(_genAK4jets, _genAK8jets, _genAK15jets, i);
		//if(_genjets.size() < 1){ cout << endl; continue; }
		_prod->GetRecoJets(_recoAK4jets, _recoAK8jets, _recoAK15jets, i);
		//if(_recoAK4jets.size() < 1){ cout << endl; continue; }
		if(i % SKIP == 0) cout << " has " << _recoAK4jets.size() << " AK4 reco jets and " << _genAK4jets.size() << " AK4 gen jets" << endl;
		if(i % SKIP == 0) cout << " has " << _recoAK8jets.size() << " AK8 reco jets and " << _genAK8jets.size() << " AK8 gen jets" << endl;
		if(i % SKIP == 0) cout << " has " << _recoAK15jets.size() << " AK15 reco jets and " << _genAK15jets.size() << " AK15 gen jets" << endl;
		
		FillGenJetHists();
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
			auto algo = make_unique<BayesCluster>(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
			algo->SetNGhosts(_nGhosts);			

			unique_ptr<GaussianMixture> gmm = algo->SubCluster();
			//y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			//comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			_recoAK4jets[j].SetModel(gmm.get(), _gev);
			
			nsubcls_tot += _recoAK4jets[j].GetNConstituents();
			if(_strategy == gmmOnly) cout << " jet has " << _recoAK4jets[j].GetNConstituents() << " subclusters" << endl;// with parameters" << endl;
			rhs.clear();
		}
		//do GMM for AK8
		for(int j = 0; j < _recoAK8jets.size(); j++){
			 _recoAK8jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			if(_strategy == gmmOnly) cout << "SubClustering reco AK8 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
			auto algo = make_unique<BayesCluster>(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
			algo->SetNGhosts(_nGhosts);			
			
			unique_ptr<GaussianMixture> gmm = algo->SubCluster();
			//y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			//comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			_recoAK8jets[j].SetModel(gmm.get(), _gev);
			
			if(_strategy == gmmOnly) cout << " jet has " << _recoAK8jets[j].GetNConstituents() << " subclusters" << endl;// with parameters" << endl;
			rhs.clear();
		}

		//do GMM for AK15
		for(int j = 0; j < _recoAK15jets.size(); j++){
			 _recoAK15jets[j].GetJets(rhs);
			//safety
			if(rhs.size() < 1) continue;
			x_nrhs_subcl.push_back((double)rhs.size());
			
			if(_strategy == gmmOnly) cout << "SubClustering reco AK15 jet #" << j << " with " << rhs.size() << " rec hits..." << endl;	
			auto algo = make_unique<BayesCluster>(rhs);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(_verb);
			algo->SetPriorParameters(_prior_params);
			algo->SetNGhosts(_nGhosts);			
			
			unique_ptr<GaussianMixture> gmm = algo->SubCluster();
			//y_time_subcl.push_back((double)t/CLOCKS_PER_SEC);
			//cout <<  "y time_subcl entry " << y_time_subcl[y_time_subcl.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
			//comptime_subcl->Fill((double)t/CLOCKS_PER_SEC);
			
			_recoAK15jets[j].SetModel(gmm.get(), _gev);
			
			if(_strategy == gmmOnly) cout << " jet has " << _recoAK15jets[j].GetNConstituents() << " subclusters" << endl;// with parameters" << endl;
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
			if(_zoom_window){
				//get relevant gen particles from hard scattering
				PointCollection genpts;
				if(_sel == singW)
					genpts += BayesPoint({_genW[0].eta(), _genW[0].phi()});
				else if(_sel == QCDdijets){
					for(auto part : _genq)
						genpts += BayesPoint({part.eta(), part.phi()});
				}
				else if(_sel == boostTop){
					for(auto part : _genTop)
						genpts += BayesPoint({part.eta(), part.phi()});
				}
				_prod->GetRecHits(rhs, i, genpts);
			}
			else
				_prod->GetRecHits(rhs, i);
		}
		//safety
		if(rhs.size() < 1) continue;
		x_nrhs.push_back((double)rhs.size());

		//assume detector radius is constant and equal for all rhs (all rhs in event are recorded in same type of detector)
		//this should be true for all events
		JetPoint rh;
		rhs[0].GetJetPointAt(0,rh); //only 1 rh
		_radius = sqrt(rh.x()*rh.x() + rh.y()*rh.y());
		if(i % SKIP == 0) cout << " and " << rhs.size() << " total rhs" << endl;
		cout << "Clustering..." << endl;	
		auto algo = std::make_unique<BayesCluster>(rhs);
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
			algo->NlnNCluster(_trees);
		}
		//N^2 version
		else if(_strategy == N2){
			//start clock
			t = clock();
			algo->N2Cluster(_trees);
		}
		t = clock() - t;
cout << "BHC time " << (double)t/CLOCKS_PER_SEC << " secs" << endl;
		//y_time.push_back((double)t/CLOCKS_PER_SEC);
		//cout <<  "y time entry " << y_time[y_time.size()-1] << " " << (double)t/CLOCKS_PER_SEC << endl;	
		//comptime->Fill((double)t/CLOCKS_PER_SEC);

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
					EvtDisplay_etaCell_phiCell->GetZaxis()->SetTitle("time [ns]");
					for(auto h : _evtdisps_obj){
						h->GetZaxis()->SetTitle("time [ns]");
					}
				}
				else
					w = rh.E(); //default energy weighted
				if(w == 0) continue;
				EvtDisplay_etaCell_phiCell->Fill(rh.eta(), rh.phi(), w);

			}
			double eta = -999;
			double phi = -999;
			//save gen particles as tmarkers
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				int id = fabs(_base->genpart_id->at(genidx));
cout << "genpart id " << id << endl;
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
						matchstring = "top1";
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
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);
				}
				else if(id == 1 || id == 2 || id == 3 || id == 4){ //light quarks
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
cout << " saving plot info for light quark" << endl;
					_plot_particles.push_back(TMarker(eta, phi, 28));
					//_plot_particles[_plot_particles.size()-1].SetMarkerColor(kAzure+4);
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);
				}
				else if(id == 5){ //b quark
					//check if first b event display is already filled
					if(_evtdisps_obj[5]->GetEntries() == 0){
						//fill if not
						gidx = 5;
						matchstring = "b1";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					else{ //fill next one
						gidx = 6;
						matchstring = "b2";
						eta = _genparts[g].eta();
						phi = _genparts[g].phi_02pi();
					}
					_plot_particles.push_back(TMarker(eta, phi, 46));
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);
				}
				else if(id == 21){ //gluon
					gidx = 2;
					matchstring = "gluon";
					eta = _genparts[g].eta();
					phi = _genparts[g].phi_02pi();
					_plot_particles.push_back(TMarker(eta, phi, 42));
					//_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlue-4);
					_plot_particles[_plot_particles.size()-1].SetMarkerColor(kBlack);
				}
				if(gidx == -1) continue; //skip if no gen particles specified
				plot_centers[matchstring] = BayesPoint({eta, phi}); //gen particle in question should be centered at (0,0) in local eta, phi coords for local (not global) evt disp
				//set by deltaR \approx 2*m/pT, taking deltaEta = deltaPhi
				double deta = 2*_genparts[g].m()/(_genparts[g].pt());
				double dphi = deta;
				plot_widths[matchstring] = BayesPoint({deta,dphi});
				if(_evt2disp_z == 1){ //update labels to time
					_evtdisps_obj[gidx]->GetZaxis()->SetTitle("time [ns]");
				}
				
			}
			if(eta == -999 && phi == -999) continue; //no gen particles specified
			//save BHC jets as ellipses
			int coloridx;
			int customColorBase = 200;
			//TColor *teaGreen = new TColor(customColorBase, 79./255., 107./255., 56./255.);
			//TColor *fernGreen = new TColor(customColorBase, 222./255., 239./255., 183./255.);
			TColor *deepSaffron = new TColor(customColorBase, 255./255., 149./255., 5./255.);
			for(int j = 0; j < _predJets.size(); j++){
				TEllipse el = PlotEll(_predJets[j]);
				if(j == 0)
					coloridx = customColorBase;//"#DEEFB7");//kGreen +4; 
				else if(j == 1)
					coloridx = 800;
				else
					coloridx = j+4;
			if(j == 0) cout << "coloridx " << coloridx << endl;
				el.SetLineColor(coloridx);	
				el.SetFillStyle(0);	
				//el.SetLineStyle(2);
				el.SetLineStyle(1);
				el.SetLineWidth(2);	
				_ellipses.push_back(el);	
				_jellipses[j] = el;			
	
				TMarker m(_predJets[j].eta(), _predJets[j].phi(), kOpenStar);
				m.SetMarkerColor(coloridx);
				//plot jet centers
				_jcenters[j] = m;			
	
				int nk = _predJets[j].GetNConstituents();
				Jet subcl;
				for(int k = 0; k < nk; k++){
					_predJets[j].GetConstituent(k, subcl);
					TEllipse sub_el = PlotEll(subcl);
					sub_el.SetLineColor(coloridx);	
					sub_el.SetFillStyle(0);	
					//sub_el.SetLineStyle(1);
					sub_el.SetLineStyle(2);
					sub_el.SetLineWidth(3);
					//_plot jet centers
					TMarker sub_m(subcl.eta(), subcl.phi(), kFullCircle);
					sub_m.SetMarkerSize(0.7); 
					sub_m.SetMarkerColor(coloridx);
					_subclellipses[j][k] = sub_el;
					_subclcenters[j][k] = sub_m;
	
				}
			}
			DrawJetEventDisplays(jetcvs);
		}
		//fill pred jet hists with jets
		FillPredJetHists();
		nsubcls_tot = 0;
		for(int j = 0; j < _predJets.size(); j++){
			nsubcls_tot += _predJets[j].GetNConstituents();
		}
		
		cout << endl;
	
		_genTop.clear();
		_genW.clear();
		_genb.clear();
		_genq.clear();


		_tree->Fill();
		_reset();
	}
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	_tree->SetDirectory(ofile);
	ofile->cd();
	_tree->Write();
	for(auto cv : jetcvs) cv->Write();
	WriteOutput(ofile, plot_centers, plot_widths);
	ofile->Close();
	_infile->Close();	
	cout << "Wrote skim to: " << _oname << endl;
}

