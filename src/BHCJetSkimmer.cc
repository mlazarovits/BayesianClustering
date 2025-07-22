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
				//min pt req
				if(_genparts[g].pt() < _minWPt){
					continue;
				}
				//make sure W has hadronic decay
				bool lep = _base->W_decayId->at(genidx);
				if(lep){ cout << "Fully leptonic W - skipping this W" << endl; continue;}
				if(i % SKIP == 0) cout << " has W with pt " << _base->genpart_pt->at(genidx) << " and energy " << _base->genpart_energy->at(genidx) << " and id " << _base->genpart_id->at(genidx) << " and pz " << _base->genpart_pz->at(genidx) << " eta " << _base->genpart_eta->at(genidx) << " phi " << _base->genpart_phi->at(genidx) << " decay id " << lep << endl;
			


				int nqs_phi = 0;
				//get decay products
				for(int gg = 0; gg < _genparts.size(); gg++){
					if(gg == g) continue;
					int ggenidx = _genparts[gg].GetUserIdx();
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
				cout << "saving W daughter - id " << _base->genpart_id->at(ggenidx) << " eta " << _genparts[gg].eta() << " phi " << _genparts[gg].phi() << " energy " << _genparts[gg].e() << endl;
					//phi distribution debugging - only look at events w/ jets (ie quarks) in some window around 0 and 2pi
					if(_base->genpart_momIdx->at(ggenidx) != genidx) continue;
					//double phi = _genparts[gg].phi_02pi();
					//if(phi > phiwindow && phi < 2*acos(-1) - phiwindow) continue;
					//nqs_phi++;	
					_genq.push_back(_genparts[gg]);
				}
				//skip W if not enough q's are in phi window
				//if(nqs_phi < 1) continue;
				_genW.push_back(_genparts[g]);
			}	
			//at least 1 W	
			if(_genW.size() < 1){ 
				if(i % SKIP == 0) cout << " has no hadronic Ws that pass pt > " << _minWPt << endl;
				continue;
			}
		}
		else if(_sel == boostTop){
			//if any leptonic top decays, continue
			int nTops = _base->Top_decayId->size();
			if(nTops < 1) continue;
			bool lepTop = _base->Top_decayId->at(0);
			//cout << "lepTop top 1 - " << lepTop << endl;
			if(nTops > 1){
				lepTop = lepTop && _base->Top_decayId->at(1);
				//cout << "lepTop top 2 - " << _base->Top_decayId->at(1) << endl;
			}
			//cout << "lepTop evt - " << lepTop << endl;
			if(lepTop){ cout << "Fully leptonic top decay - skipping" << endl; continue;}
			//tops are ordered highest to lowest energy, top decay ids follow the same order (ie first top decay id is associated to first top)
			//at least 1 top quark with pt, E requirements
			int ntop = 0;
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				if(fabs(_base->genpart_id->at(genidx)) != 6) continue;
				//cout << "genparts idx " << g << " topidx " << genidx << " e " << _base->genpart_energy->at(genidx) << " decay idx " << ntop << endl;
				topidx_decayidx[genidx] = ntop;
				ntop++;
			}
			ntop = 0;
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				if(fabs(_base->genpart_id->at(genidx)) != 6) continue;
				if(_genparts[g].pt() < _minTopPt) continue;
				//make sure top has hadronic decay (0 - hadronic, 1 - leptonic)
				//topdecayid: 0 = had, 1 = lep
				bool lep = _base->Top_decayId->at(topidx_decayidx[genidx]);
				//cout << "topidx " << g << " E " << _base->genpart_energy->at(genidx) << " decay idx " << topidx_decayidx[g] << " decay id " << _base->Top_decayId->at(topidx_decayidx[g]) << " lep " << lep << endl;	
				if(lep) continue;
				if(i % SKIP == 0) cout << " has hadronic top with pt " << _genparts[g].pt() << " and energy " << _genparts[g].E() << " and id " << _base->genpart_id->at(genidx) << " and pz " << _genparts[g].pz() << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << " decayid " << _base->Top_decayId->at(topidx_decayidx[genidx]) << endl;
				ntop++;
			}	
			if(ntop < 1){
				if(i % SKIP == 0) cout << " has no hadronic tops that pass pt > " << _minTopPt << endl;
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
				for(int p = 0; p < _procCats.size(); p++){
					//_procCats[p].hists1D[0][98]->Fill(mean.at(0,0));
					//_procCats[p].hists1D[0][99]->Fill(mean.at(1,0));
					//_procCats[p].hists1D[0][100]->Fill(mean.at(2,0));
					//_procCats[p].hists1D[0][101]->Fill(sqrt(cov.at(0,0)));
					//_procCats[p].hists1D[0][102]->Fill(sqrt(cov.at(1,1)));
					//_procCats[p].hists1D[0][103]->Fill(sqrt(cov.at(2,2)));
				}
			}
		
			rhs.clear();
		}
		for(int p = 0; p < _procCats.size(); p++){
			_procCats[p].hists1D[0][85]->Fill(nsubcls_tot);
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
			_procCats[1].hists1D[0][78]->Fill(rh.t());
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
			if(_predJets.size() < 2) continue; //want at least two boosted tops in an event
			//if there are any jets that are good matches to bs then the top can't be boosted, skip
			vector<int> genbMatchIdxs(_predJets.size(),-1);
			GenericMatchJet(_predJets,_genparts, genbMatchIdxs, 5); //dr match BHC jets to gen Ws
			//want subleading jets to be very good matches to their respective gen b quarks (these are the "tags")
			//good match = small dr, Eratio close to 1
			vector<bool> bMatch(_predJets.size(), 0);
			vector<Jet> topCand;
			for(int j = 0; j < _predJets.size(); j++){
				bool drGood = true;
				bool EratioGood = true;
				int genbidx = genbMatchIdxs[j];
				if(genbidx == -1){
					topCand.push_back(_predJets[j]);
					continue;
				}
				//check dr
				double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genparts[genbidx].eta(), _genparts[genbidx].phi());
				//check Eratio
				double Eratio = _predJets[j].E()/_genparts[genbidx].E();
				//cout << "_pred jet #" << j << " has eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " and E " << _predJets[j].E() << endl;
				//cout << " gen-matched to b #" << genbidx << "  with eta " << _genparts[genbidx].eta() << " phi " << _genparts[genbidx].phi() << " and E " << _genparts[genbidx].E() << " dr " << dr << " Eratio " << Eratio << endl; 
				if(dr > 0.1) drGood = false;
				if(Eratio < 0.8 || Eratio > 1.2) EratioGood = false;
				bMatch[j] = drGood && EratioGood;
				//if bad b match - want to save associated top
				if(!(drGood && EratioGood)){
					//cout << "bad b match - saving associated top " << endl;
					int genidx = _genparts[genbidx].GetUserIdx();
					int bmom = _base->genpart_momIdx->at(genidx);
					for(int g = 0; g < _genparts.size(); g++){
						int ggenidx = _genparts[g].GetUserIdx();
						if(fabs(_base->genpart_id->at(ggenidx)) != 6) continue;
						//make sure top still passes gen pt req
						if(_base->genpart_pt->at(ggenidx) < _minTopPt) continue;
						//needs to be hadronic
						bool lep = _base->Top_decayId->at(topidx_decayidx[ggenidx]);
						if(lep) continue;	
						if(ggenidx == bmom){ //save top
							cout << "saving top # " << ggenidx << " for b " << genidx << " with energy " << _genparts[g].E() << " eta " << _genparts[g].eta() << " phi " << _genparts[g].phi() << endl;
							_genTop.push_back(_genparts[g]);
							break;
						}	
					}
					topCand.push_back(_predJets[j]);

				}
				else{
					_genb.push_back(_genparts[j]);
		
				}
			}
			cout << "# gen bs " << _genb.size() << " # top candidates " << topCand.size() << " # gen tops " << _genTop.size() << endl;
			//if there's more than 1 good gen b match (ie not a candidate for a boosted top), skip
			if(_genb.size() > 1){
				continue;
			}
			else{
				_predJets = topCand;
			}

		}


		//fill event display if specified _evt2disp
		if(i == _evt2disp){
			cout << "Displaying event " << i << endl;
			for(auto rh : rhs){
				_procCats[0].hists2D[0][86]->Fill(rh.eta(), rh.phi(), rh.E());
				_procCats[1].hists2D[0][86]->Fill(rh.eta(), rh.phi(), rh.E());
			}
			//save gen particles as tmarkers
			for(int g = 0; g < _genparts.size(); g++){
				int genidx = _genparts[g].GetUserIdx();
				double eta = _genparts[g].eta();
				double phi = _genparts[g].phi_02pi();
				int id = fabs(_base->genpart_id->at(genidx));
				//TODO: set to pretty colors with hex codes and make colors related ie for quarks (light vs b different but similar)
				if(id == 24){ //W
					plot_particles.push_back(TMarker(eta, phi, kPlus));
					plot_particles[plot_particles.size()-1].SetMarkerColor(kRed-4);
				}
				else if(id == 1 || id == 2 || id == 3){ //light quarks
					plot_particles.push_back(TMarker(eta, phi, kCircle));
					plot_particles[plot_particles.size()-1].SetMarkerColor(kAzure+4);
				}
				else if(id == 5){ //b quark
					plot_particles.push_back(TMarker(eta, phi, kOpenSquare));
					plot_particles[plot_particles.size()-1].SetMarkerColor(kCyan-4);
				}
				else if(id == 21){ //gluon
					plot_particles.push_back(TMarker(eta, phi, kStar));
					plot_particles[plot_particles.size()-1].SetMarkerColor(kBlue-4);
				}
			}
			//save BHC jets as ellipses
			for(int j = 0; j < _predJets.size(); j++){
				Matrix mu, cov;
				_predJets[j].GetClusterParams(mu, cov);
				//ellipse center
				double eta = mu.at(0,0); 
				double phi = mu.at(1,0);

				//get 2D matrix for jet size
				Matrix cov2D(2,2);
				Get2DMat(cov,cov2D);	
				vector<double> eigvals;
				vector<Matrix> eigvecs;
				cov2D.eigenCalc(eigvals, eigvecs);
				
				//define radii (r1 > r2)
				double r1, r2;
				Matrix leadvec;
				if(eigvals[0] > eigvals[1]){
					r1 = eigvals[0];
					leadvec = eigvecs[0];
					r2 = eigvals[1];
				}
				else{
					r1 = eigvals[1];
					leadvec = eigvecs[1];
					r2 = eigvals[0];

				}
				//define angle of r1 rotation
				double theta = atan2(leadvec.at(1,0), leadvec.at(0,0));
				theta = 180 * theta/(4*atan(1)); //put to degrees
				ellipses.push_back(TEllipse(eta, phi, r1, r2, 0, 360, theta));	
				//plot jet centers
				plot_particles.push_back(TMarker(eta, phi, 30)); //30 = open star may need to change if not rendering	
				ellipses[ellipses.size()-1].SetLineColor(kGreen-4);	
				plot_particles[plot_particles.size()-1].SetMarkerColor(kGreen-4);	

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
	WriteOutput(ofile);
	//write ellispes + tmarkers to root file (if i can?)
	//for(int el = 0; el < ellipses.size(); el++){
	//	ellispes[el].Write();
	//}
	//for(int m = 0; m < plot_particles.size(); m++){
	//	plot_particles[m].Write();
	//}
	ofile->Close();
	_infile->Close();	
	cout << "Wrote skim to: " << _oname << endl;
}

