#include "BaseProducer.hh"
#include <iomanip>

void BaseProducer::GetTrueJets(vector<Jet>& jets, int evt, double gev){
        if(gev == -1) gev = _gev;
        double px, py, pz, pt, phi, eta;
        jets.clear();
	//true = skip
	//false = keep (ok)
	bool hemVeto = false;	
	double minrhE = 0.5;
	if(evt > _nEvts) return;

        _base->GetEntry(evt);
        int nJets = (int)_base->Jet_energy->size();
        int nrhs, rhidx;
	bool jetid;
	double dr, deta, dphi, eme;
	double timecorr, drh, dpv, calibfactor;

	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;

	Point vtx(3);
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

		eme = _base->Jet_energy->at(j)*(_base->Jet_neEmEF->at(j)+_base->Jet_chEmEF->at(j));
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
				if(_base->ECALRecHit_energy->at(rhidx) < minrhE) continue;				
				//spike rejection? - only for rhE > 4 GeV
				if(_spikes && _data && _base->ECALRecHit_energy->at(rhidx) > 4){
					cout << "rejecting spikes" << endl;
					if( _base->ECALRecHit_swCross->at(rhidx) < 0.02*log10(_base->ECALRecHit_energy->at(rhidx))+0.02)
						continue;
				}
				
				//TOF from 0 to rh location
				drh = _base->ECALRecHit_0TOF->at(rhidx);
				//TOF from PV to rh location
				dpv = _base->ECALRecHit_pvTOF->at(rhidx); 
				timecorr = drh - dpv;
				//redo dr matching tighter - dr = 0.5
				dr = sqrt(deltaR2(_base->Jet_eta->at(j), _base->Jet_phi->at(j), _base->ECALRecHit_eta->at(rhidx), _base->ECALRecHit_phi->at(rhidx)));
				if(dr > 0.5) continue;				


				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh;
				if(_calibmap){
                          		calibfactor = GetTimeCalibrationFactor(_base->ECALRecHit_ID->at(rhidx));
					rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx) + timecorr - calibfactor);
				}
				else{
					calibfactor = 0;	
					rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx) + timecorr);
				}
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
	//cout << "adding rh with x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << rh.t() << " calib " << calibfactor << endl;			
                                
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				jet.AddRecHit(rh);
                        }

                }
//cout << "BaseProducer::GetTrueJets - jet #" << jets.size() << " of " << nJets << " j - " << j << " energy: " << _base->Jet_energy->at(j) << " phi: " << _base->Jet_phi->at(j) << " eta: " << _base->Jet_eta->at(j) << " has " << rhs.size() << " rhs - event #" << _base->Evt_event << " charged EMF: " << _base->Jet_chEmEF->at(j) << " neutral EMF: " << _base->Jet_neEmEF->at(j) << " charged HF: " << _base->Jet_chHEF->at(j) << " neutral HF: " << _base->Jet_neHEF->at(j) <<  endl;
		//jet.Print();
		if(jet.GetNRecHits() < _minnrhs) continue;
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
	Point vtx(3);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);
	for(int p = 0; p < nPhos; p++){
		//if selected photons # is already 2, return (only want 2 highest pt photons that pass selection)
		if(selPhoCount == 2) return;
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
                if(_base->Photon_pt->at(p) < 30.) continue;
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
				//TODO: removed when ntuples are fixed!
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
				timecorr = drh - dpv;
	//cout << "pho #" << p << " rh # " << r << " time " << _base->ECALRecHit_time->at(rhidx) << " drh " << drh << " dpv " << dpv << " saved time (with other factors) " << _base->ECALRecHit_time->at(rhidx) + timecorr - calibfactor << endl;

				

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh;
				if(_calibmap){
                          		calibfactor = GetTimeCalibrationFactor(_base->ECALRecHit_ID->at(rhidx));
					rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx) + timecorr - calibfactor);
				}
				else{	
					rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx) + timecorr);
				}
                               
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
	////cout << "adding rh with x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) <<  " etajetpoint " << rh.eta() << " phi " << _base->ECALRecHit_phi->at(rhidx) << " phijp " << rh.phi() << " timecorr " << timecorr << " calib " << calibfactor << endl;			
				
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
	//cout << "adding rh with rhidx " << rhidx << " x " << _base->ECALRecHit_rhx->at(rhidx) << " y " << _base->ECALRecHit_rhy->at(rhidx) << " z " << _base->ECALRecHit_rhz->at(rhidx) << " t " << _base->ECALRecHit_time->at(rhidx) << " eta " << _base->ECALRecHit_eta->at(rhidx) << " phi " << _base->ECALRecHit_phi->at(rhidx) << " nrhs so far " << nrhs << " r " << r << " rhid " << rhid << " counts in SC " << count(rhs.begin(), rhs.end(), rhid) << " counts in ECAL " << count(rhids.begin(), rhids.end(), rhid) << endl;
				nrhs++; 
                                pho.AddRecHit(rh);
                		jrhids.push_back(_base->ECALRecHit_ID->at(rhidx));
		        }

                }
		if(pho.GetNRecHits() < 2) continue;
		selPhoCount++;
	//	cout << jrhids.size() << " nrhs in pho" << endl;
	//	for(auto rh : jrhids) cout << "rh id  " << rh << " count " << count(jrhids.begin(), jrhids.end(), rh) << endl;

		phos.push_back(pho);
		jrhids.clear();
        }


}
