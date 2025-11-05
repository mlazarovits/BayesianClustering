#include "BaseProducer.hh"
#include <iomanip>

void BaseProducer::GetTrueJets(vector<Jet>& jets, int evt, double gev){
        if(gev == -1) gev = _gev;
        double px, py, pz, pt, phi, eta;
        jets.clear();
	//true = skip
	//false = keep (ok)
	bool hemVeto = false;	
	if(evt > _nEvts) return;

        _base->GetEntry(evt);
        int nJets = (int)_base->Jet_energy->size();
        int nrhs, rhidx;
	bool jetid;
	double dr, deta, dphi, eme, et;
	double timecorr, drh, dpv, calibfactor;

	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;

	BayesPoint vtx(3);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);
	vector<unsigned int> rhs;
	for(int j = 0; j < nJets; j++){
                pt = _base->Jet_pt->at(j);
                phi = _base->Jet_phi->at(j);
                eta = _base->Jet_eta->at(j);

                px = pt*cos(phi);
                py = pt*sin(phi);
                pz = pt*sinh(eta);

		rhs.clear();
		rhs = _base->Jet_drRhIds->at(j);
                nrhs = rhs.size();

		//using transverse energy
		et = sqrt(_base->Jet_mass->at(j)*_base->Jet_mass->at(j) + pt*pt);
		eme = et*(_base->Jet_neEmEF->at(j)+_base->Jet_chEmEF->at(j));
		//Jet selection
                if(_base->Jet_pt->at(j) < _minpt) continue;
                if(fabs(_base->Jet_eta->at(j)) > _minobjeta) continue;
		if(eme < _mineme) continue;

		//create jet id (based on tight 2017 Run II recommendations)
		//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
		
		jetid = (_base->Jet_neEmEF->at(j) < 0.9) && (_base->Jet_neHEF->at(j) < 0.9) && (_base->Jet_nConstituents->at(j) > 1)
                        && (_base->Jet_chHEF->at(j) > 0.0) && (_base->Jet_chHM->at(j) > 0);
                if(!jetid) continue;
		
		//hem veto?
		if(_year == 2018 && _data){
			hemVeto = ( (_base->Evt_run >= 319077) && (eta > -1.58) && (eta < -1.34) && (phi > 4.8) && (phi < 5.4) );
			//skip whole event
			if(hemVeto) return;
		}


                Jet jet(px, py, pz, _base->Jet_energy->at(j));
                jet.SetVertex(vtx);
		jet.SetUserIdx(j);
		//set rec hits in jet
		for(int r = 0; r < rhs.size(); r++){
                        unsigned int rhid = rhs[r];
                        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				//if rh is in endcap, skip
				if(fabs(_base->ECALRecHit_eta->at(rhidx)) > 1.479) continue;
				//remove timing reco (ratio) failed fits
				if(_base->ECALRecHit_time->at(rhidx) == 0.) continue;
				//energy cut
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;
				if(_maxrhE != -999 && _base->ECALRecHit_energy->at(rhidx) > _maxrhE) continue;				
				//spike rejection? - only studied for rhE > 4 GeV
				if(_spikes && _data){
					if(1 - _base->ECALRecHit_swCross->at(rhidx) < 0.02*log10(_base->ECALRecHit_energy->at(rhidx))+0.02)
						continue;
				
					//add tighter sw+ cut
					if(_base->ECALRecHit_swCross->at(rhidx) > 0.4) continue;
				}
				
				//TOF from 0 to rh location
				drh = _base->ECALRecHit_0TOF->at(rhidx);
				//TOF from PV to rh location
				dpv = _base->ECALRecHit_pvTOF->at(rhidx); 
				if(_spatial_corr)
					timecorr = drh - dpv;
				else
					timecorr = 0;
				//redo dr matching tighter - dr = 0.5
				dr = sqrt(deltaR2(_base->Jet_eta->at(j), _base->Jet_phi->at(j), _base->ECALRecHit_eta->at(rhidx), _base->ECALRecHit_phi->at(rhidx)));
				if(dr > 0.5) continue;				


				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh;
				float time = _timecalibTool->getCorrectedTime(_base->ECALRecHit_time->at(rhidx), _base->ECALRecHit_ampres->at(rhidx), _base->ECALRecHit_ID->at(rhidx), _base->Evt_run, _timecalibTag, _mctype);
				_rhIdToRes[_base->ECALRecHit_ID->at(rhidx)] = (double)_timecalibTool->getTimeResoltuion( _base->ECALRecHit_ampres->at(rhidx), _base->ECALRecHit_ID->at(rhidx), _base->Evt_run, _timecalibTag, _mctype);	
				/*
				double time = _base->ECALRecHit_time->at(rhidx);
				if(_calib){
                          		calibfactor = GetTimeCalibrationFactor(_base->ECALRecHit_ID->at(rhidx), (int)_base->Evt_run);
				}
				else{
					calibfactor = 0;	
				}
				if(_timesmear){
					time = SmearRecHitTime(_base->ECALRecHit_ampres->at(rhidx), time);
				}
				*/
				time = time + timecorr;
				rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                _base->ECALRecHit_rhz->at(rhidx), (double)time);
				//cout << "_spatial_corr " << _spatial_corr << " jet raw time " <<  _base->ECALRecHit_time->at(rhidx) << " saved rh time " << rh.t() << endl;
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
	//cout << "adding rh with x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << rh.t() << " calib " << calibfactor << endl;			
                                
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				rh.SetUserIdx(rhidx);
				//cut out mist if specified
				if(_mistcuts){
					if(rh.t() < -0.5 && rh.E() < 500 && rh.E() > 10) continue;
					if(rh.t() > 0.5 && rh.E() < 500 && rh.E() > 10) continue;
				}
				//specify if this rechit as tripped a gain switch - if so, the time will be assigned a large uncertainty
				if(_base->ECALRecHit_hasGS1->at(rhidx) || _base->ECALRecHit_hasGS6->at(rhidx))
					rh.SetInvalidTime();
				jet.AddRecHit(rh);
                        }

                }
//cout << "BaseProducer::GetTrueJets - jet #" << jets.size() << " of " << nJets << " j - " << j << " energy: " << _base->Jet_energy->at(j) << " phi: " << _base->Jet_phi->at(j) << " eta: " << _base->Jet_eta->at(j) << " has " << rhs.size() << " rhs - event #" << _base->Evt_event << " charged EMF: " << _base->Jet_chEmEF->at(j) << " neutral EMF: " << _base->Jet_neEmEF->at(j) << " charged HF: " << _base->Jet_chHEF->at(j) << " neutral HF: " << _base->Jet_neHEF->at(j) <<  endl;
		//jet.Print();
		if(jet.GetNRecHits() < _minnrhs) continue;

//		cout << "\njet #" << jets.size() << " nconstituents " << _base->Jet_nConstituents->at(j) << " total rhs " << _base->Jet_drRhIds->at(j).size() << " ana rhs " << jet.GetNPoints() << " event " << evt << endl;
		jets.push_back(jet);
        }
}

void BaseProducer::GetTruePhotons(vector<Jet>& phos, int evt, double gev){
        if(gev == -1) gev = _gev;
	double px, py, pz, pt, phi, eta;
        phos.clear();
        if(evt > _nEvts) return;
        _base->GetEntry(evt);
        int nPhos = (int)_base->Photon_energy->size();
	//only take leading and subleading (if these exist)
	int selPhoCount = 0; //shouldnt be incremented to >2
	//if(nPhos > 2) nPhos = 2;
        int nrhs, rhidx;
	double timecorr, calibfactor; 
	double drh, dpv;
	int scidx;
	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;

	//true = skip
	//false = keep (ok)
	bool hemVeto = false;	
	BayesPoint vtx(3);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);
//cout << "\n evt " << evt << " n total phos " << nPhos << endl;
	for(int p = 0; p < nPhos; p++){
		//if selected photons # is already 2, return (only want 2 highest pt photons that pass selection)
		if(selPhoCount == 2) return;
		//if excluded flag is up, skip (means this one should be skipped in favor of an OOT in same collection)
		if(_base->Photon_excluded->at(p)) continue;
		
		pt = _base->Photon_pt->at(p);
                phi = _base->Photon_phi->at(p);
                eta = _base->Photon_eta->at(p);

                px = pt*cos(phi);
                py = pt*sin(phi);
                pz = pt*sinh(eta);

		scidx = _base->Photon_scIndex->at(p);

                nrhs = _base->SuperCluster_rhIds->at(scidx).size();
		//hem veto?
		if(_year == 2018 && _data){
			hemVeto = ( (_base->Evt_run >= 319077) && (eta > -1.58) && (eta < -1.34) && (phi > 4.8) && (phi < 5.4) );
			//skip whole event
			if(hemVeto) return;
		}

		//Photon selection
                if(_base->Photon_pt->at(p) < _minpt) continue;
		if(fabs(_base->Photon_eta->at(p)) > _minobjeta) continue;
		//isolation cuts
		bool iso;
		bool trksum;
		bool ecalrhsum;
		bool htowoverem;
		if(_isocut){
                        trksum = _base->Photon_trkSumPtSolidConeDR04->at(p) < 6.0;
                        ecalrhsum = _base->Photon_ecalRHSumEtConeDR04->at(p) < 10.0;
                        htowoverem = _base->Photon_hadTowOverEM->at(p) < 0.02;
                        iso = trksum && ecalrhsum && htowoverem;
                	if(!iso) continue;
                }

                Jet pho(px, py, pz, _base->Photon_energy->at(p));
                pho.SetVertex(vtx);
		pho.SetUserIdx(p);
		//set rec hits in photon
		vector<unsigned int> rhs = _base->SuperCluster_rhIds->at(scidx);
		vector<float> fracs = _base->SuperCluster_rhFracs->at(scidx);
		double rhe;
		int nrhs = 0;
		vector<unsigned int> jrhids;
		//cout << rhs.size() << " in SC " << rhids.size() << " in ECAL" << endl;
                for(int r = 0; r < rhs.size(); r++){
                        unsigned int rhid = rhs[r];
                        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				//skip rhs that have already been looked at - avoids duplicates in SC
				auto jrhit = std::find(jrhids.begin(), jrhids.end(), rhid);
				if(jrhit != jrhids.end()) continue;
				//if rh is in endcap, skip
				if(fabs(_base->ECALRecHit_eta->at(rhidx)) > 1.479) continue;
				//remove timing reco (ratio) failed fits
				if(_base->ECALRecHit_time->at(rhidx) == 0.) continue;
				//energy cut
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;				


				//TOF from 0 to rh location
				drh = _base->ECALRecHit_0TOF->at(rhidx);
				//TOF from PV to rh location - use this to improve time covariance
				dpv = _base->ECALRecHit_pvTOF->at(rhidx); 
				if(_spatial_corr)
					timecorr = drh - dpv;
				else
					timecorr = 0;
				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh;
				float time = _timecalibTool->getCorrectedTime(_base->ECALRecHit_time->at(rhidx), _base->ECALRecHit_ampres->at(rhidx), _base->ECALRecHit_ID->at(rhidx), _base->Evt_run, _timecalibTag, _mctype);	
				_rhIdToRes[_base->ECALRecHit_ID->at(rhidx)] = (double)_timecalibTool->getTimeResoltuion( _base->ECALRecHit_ampres->at(rhidx), _base->ECALRecHit_ID->at(rhidx), _base->Evt_run, _timecalibTag, _mctype);	
				/*
				double time = _base->ECALRecHit_time->at(rhidx);
				if(_calib){
                          		calibfactor = GetTimeCalibrationFactor(_base->ECALRecHit_ID->at(rhidx), (int)_base->Evt_run);
				}
				else{
					calibfactor = 0;	
				}
				//cout << "og time " << time;
				*/
				time = time + timecorr;
				//cout << " calibrated time " << time - timecorr << " calib factor " << calibfactor << endl;
				if(_timesmear){
					time = SmearRecHitTime(_base->ECALRecHit_ampres->at(rhidx), time);
				}
				rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                _base->ECALRecHit_rhz->at(rhidx), (double)time);
                               
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
//	cout << "adding rh with x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) <<  " etajetpoint " << rh.eta() << " phi " << _base->ECALRecHit_phi->at(rhidx) << " phijp " << rh.phi() << " timecorr " << timecorr << " calib " << calibfactor << endl;			
				
				rhe = _base->ECALRecHit_energy->at(rhidx);
				//multiply energy by hitsAndFractions fraction
				//indexed within supercluster
				if(_applyFrac){
					rhe *= fracs[r];
					if(rhe < 0.1) continue;
				}					
				//cout << std::setprecision(10) << "rh e og: " << _base->ECALRecHit_energy->at(rhidx) << " frac: " << fracs[r] << " rhE*frac " << _base->ECALRecHit_energy->at(rhidx)*fracs[r] << " new rh e: " << rhe << endl;
				rh.SetEnergy(rhe);
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				//specify if this rechit as tripped a gain switch - if so, the time will be assigned a large uncertainty
				if(_base->ECALRecHit_hasGS1->at(rhidx) || _base->ECALRecHit_hasGS6->at(rhidx))
					rh.SetInvalidTime();
	//cout << "adding rh with rhidx " << rhidx << " x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) << " phi " << _base->ECALRecHit_phi->at(rhidx) << " nrhs so far " << nrhs << " r " << r << " rhid " << rhid << " counts in SC " << count(rhs.begin(), rhs.end(), rhid) << " counts in ECAL " << count(rhids.begin(), rhids.end(), rhid) << endl;
				nrhs++; 
                                pho.AddRecHit(rh);
				//cout << "_spatial_corr " << _spatial_corr << " timecorr " << timecorr << " photon raw time " <<  _base->ECALRecHit_time->at(rhidx) << " saved rh time " << rh.t() << endl;
                		jrhids.push_back(_base->ECALRecHit_ID->at(rhidx));
		        }
			//else cout << "rechit with ID " << rhid << " not found in ECAL collection" << endl;

                }
		//cout << "pt " << _base->Photon_pt->at(p) << " eta " << _base->Photon_eta->at(p) << " # rechits " << pho.GetNRecHits() << endl;	
		if(pho.GetNRecHits() < 2) continue;
		selPhoCount++;
	//	cout << jrhids.size() << " nrhs in pho" << endl;
	//	for(auto rh : jrhids) cout << "rh id  " << rh << " count " << count(jrhids.begin(), jrhids.end(), rh) << endl;

		phos.push_back(pho);
		jrhids.clear();
        }


}
int BaseProducer::GetTrueSuperClusters(vector<Jet>& supercls, int evt, double gev){
        if(gev == -1) gev = _gev;
	double px, py, pz, pt, phi, eta, theta, et, E;
        supercls.clear();
        if(evt > _nEvts){ cout << "evt " << evt << " of " << _nEvts << endl; return -1;}
        _base->GetEntry(evt);
        int nSCs = (int)_base->SuperCluster_energy->size();
        int nrhs, rhidx;
	double timecorr, calibfactor, drh, dpv;
	int scidx;

	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;

	//true = skip
	//false = keep (ok)
	bool hemVeto = false;	
	BayesPoint vtx(3);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);
//cout << "this event " << evt << " ntuple event " << _base->Evt_event << endl;
	vector<JetPoint> sc_rhs;
//cout << "producer found " << nSCs << " scs from ntuple for event" << endl;
	for(int sc = 0; sc < nSCs; sc++){
		//cout << "sc #" << sc << endl;
                phi = _base->SuperCluster_phi->at(sc);
                eta = _base->SuperCluster_eta->at(sc);
		theta = 2*atan2(1,exp(eta));

                nrhs = _base->SuperCluster_rhIds->at(sc).size();

		//if excluded flag is up, skip (means this one should be skipped in favor of an OOT in same collection)
		//depreciated as of ntuple v30
		//if(_base->SuperCluster_excluded->at(sc)) continue;
	
		//hem veto?
		if(_year == 2018 && _data){
			hemVeto = ( (_base->Evt_run >= 319077) && (eta > -1.58) && (eta < -1.34) && (phi > 4.8) && (phi < 5.4) );
			//skip whole event
			if(hemVeto) return -1;
		}

		//SuperCluster selection
                E = _base->SuperCluster_energy->at(sc);
		et = E*sin(theta);
//cout << "sc Et " << et << " eta " << eta << endl;
                if(et < _minpt) continue;
//cout << "sc passed et cut " << endl;
		if(fabs(eta) > _minobjeta) continue;
//cout << "sc passed eta cut " << eta << endl;
                
		//set rec hits in sc
		vector<unsigned int> rhs = _base->SuperCluster_rhIds->at(sc);
		vector<float> fracs = _base->SuperCluster_rhFracs->at(sc);
		double rhe;
		int nrhs = 0;
		vector<unsigned int> jrhids;
		//cout << rhs.size() <<  " rhs in SC " << rhids.size() << " rhs in ECAL" << endl;
                for(int r = 0; r < rhs.size(); r++){
                        unsigned int rhid = rhs[r];
//cout << "rhid " << rhid << endl;
//for(int rr = 0; rr < rhids.size(); rr++){
//	if(rhids[rr] == rhid) cout << "ID MATCH FOUND" << endl;
//
//
//}
                        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				//skip rhs that have already been looked at - avoids duplicates in SC
				auto jrhit = std::find(jrhids.begin(), jrhids.end(), rhid);
				if(jrhit != jrhids.end()) continue;
//cout << "rh passed duplicate check with eta " << _base->ECALRecHit_eta->at(rhidx) << endl;
				//if rh is in endcap, skip
				if(fabs(_base->ECALRecHit_eta->at(rhidx)) > 1.479) continue;
//cout << "rh passed eta req" << endl;
				//remove timing reco (ratio) failed fits
				if(_base->ECALRecHit_time->at(rhidx) == 0.) continue;
//cout << "rh passed timing reco check" << endl;
				//energy cut
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;				
//cout << "rh passed min energy" << endl;


				//spike rejection? - only studied for rhE > 4 GeV
				if(_spikes){
					if(1 - _base->ECALRecHit_swCross->at(rhidx) < 0.02*log10(_base->ECALRecHit_energy->at(rhidx))+0.02)
						continue;
				
				}
				//TOF from 0 to rh location
				drh = _base->ECALRecHit_0TOF->at(rhidx);
				//TOF from PV to rh location - use this to improve time covariance
				dpv = _base->ECALRecHit_pvTOF->at(rhidx); 
				if(_spatial_corr)
					timecorr = drh - dpv;
				else
					timecorr = 0;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh;
				float time = _timecalibTool->getCorrectedTime(_base->ECALRecHit_time->at(rhidx), _base->ECALRecHit_ampres->at(rhidx), _base->ECALRecHit_ID->at(rhidx), _base->Evt_run, _timecalibTag, _mctype);
				_rhIdToRes[_base->ECALRecHit_ID->at(rhidx)] = (double)_timecalibTool->getTimeResoltuion( _base->ECALRecHit_ampres->at(rhidx), _base->ECALRecHit_ID->at(rhidx), _base->Evt_run, _timecalibTag, _mctype);	
				/*	
				double time = _base->ECALRecHit_time->at(rhidx);
				if(_calib){
                          		calibfactor = GetTimeCalibrationFactor(_base->ECALRecHit_ID->at(rhidx), (int)_base->Evt_run);
				}
				else{
					calibfactor = 0;	
				}
				if(_timesmear){
					time = SmearRecHitTime(_base->ECALRecHit_ampres->at(rhidx), time);
				}
				*/
				time = time + timecorr;
				rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                _base->ECALRecHit_rhz->at(rhidx), (double)time);
                              //cout << "rh time " << rh.t() << endl; 
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
//cout << "rh passed in time enough req" << endl;
//	cout << "adding rh with x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) <<  " etajetpoint " << rh.eta() << " phi " << _base->ECALRecHit_phi->at(rhidx) << " phijp " << rh.phi() << " timecorr " << timecorr << " calib " << calibfactor << endl;			
				
				rhe = _base->ECALRecHit_energy->at(rhidx);
				//multiply energy by hitsAndFractions fraction
				//indexed within supercluster
				if(_applyFrac){
					rhe *= fracs[r];
					if(rhe < 0.1) continue;
				}					
				//cout << std::setprecision(10) << "rh e og: " << _base->ECALRecHit_energy->at(rhidx) << " frac: " << fracs[r] << " rhE*frac " << _base->ECALRecHit_energy->at(rhidx)*fracs[r] << " new rh e: " << rhe << endl;
				rh.SetEnergy(rhe);
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				//specify if this rechit as tripped a gain switch - if so, the time will be assigned a large uncertainty
				if(_base->ECALRecHit_hasGS1->at(rhidx) || _base->ECALRecHit_hasGS6->at(rhidx)){
					//cout << "RH HAS INVALID TIME with weight " << rh.GetWeight() << endl;
					rh.SetInvalidTime();
				}
	//cout << "adding rh with rhidx " << rhidx << " x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) << " phi " << _base->ECALRecHit_phi->at(rhidx) << " nrhs so far " << nrhs << " r " << r << " rhid " << rhid << " counts in SC " << count(rhs.begin(), rhs.end(), rhid) << " counts in ECAL " << count(rhids.begin(), rhids.end(), rhid) << endl;
				nrhs++; 
                                sc_rhs.push_back(rh);
                		jrhids.push_back(_base->ECALRecHit_ID->at(rhidx));
	//cout << "adding rh with rhidx " << rhidx << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) << " phi " << _base->ECALRecHit_phi->at(rhidx) << " energy " << rhe << " res " << _rhIdToRes[rh.rhId()] <<  endl;
		        }
			//else{
			//	cout << "no match found for rh #" << r << " with id " << rhid << " in ECAL rechit list" << endl;
			//}

                }

		Jet supercl(sc_rhs, vtx);
		supercl.SetUserIdx(sc);
		if(supercl.GetNRecHits() < 2) continue;
//cout << "sc passed min # rhs cut " << endl;
	//	cout << jrhids.size() << " nrhs in pho" << endl;
	//	for(auto rh : jrhids) cout << "rh id  " << rh << " count " << count(jrhids.begin(), jrhids.end(), rh) << endl;
		supercls.push_back(supercl);
		jrhids.clear();
		sc_rhs.clear();
        }

	return 0;
}
