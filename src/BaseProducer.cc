#include "BaseProducer.hh"


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
	double dr, deta, dphi, eme;
	double timecorr, drh;

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
                if(fabs(_base->Jet_eta->at(j)) > 1.5) continue;
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
				//TOF from 0 to rh location
				drh = hypo(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), _base->ECALRecHit_rhz->at(rhidx))/_c;
				//TOF from PV to rh location
				timecorr = drh - hypo((_base->ECALRecHit_rhx->at(rhidx) - _base->PV_x), (_base->ECALRecHit_rhy->at(rhidx) - _base->PV_y), (_base->ECALRecHit_rhz->at(rhidx) - _base->PV_z))/_c;

                           
				//redo dr matching tighter - dr = 0.5
				dr = sqrt(deltaR2(_base->Jet_eta->at(j), _base->Jet_phi->at(j), _base->ECALRecHit_eta->at(rhidx), _base->ECALRecHit_phi->at(rhidx)));
				if(dr > 0.5) continue;				

				//remove timing reco (ratio) failed fits
				if(_base->ECALRecHit_time->at(rhidx) == 0.) continue;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
                               	//undo adjustment currently made in ntuples (traw - drh/c)
				//TODO: UNDO THIS LATER WHEN NTUPLES ARE UPDATED
				JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx) + timecorr + drh);
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
                                
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
        int nPhos = 2;//only take leading + subleading photons //(int)_base->Photon_energy->size();
        int nrhs, rhidx;
	double timecorr = 0; //to get photon time in "PV frame"
	double drh;

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
                pt = _base->Photon_pt->at(p);
                phi = _base->Photon_phi->at(p);
                eta = _base->Photon_eta->at(p);

                px = pt*cos(phi);
                py = pt*sin(phi);
                pz = pt*sinh(eta);

                nrhs = _base->Photon_rhIds->size();
	
		//hem veto?
		if(_year == 2018 && _data){
			hemVeto = ( (_base->Evt_run >= 319077) && (eta > -1.58) && (eta < -1.34) && (phi > 4.8) && (phi < 5.4) );
			//skip whole event
			if(hemVeto) return;
		}

		//Photon selection
                if(_base->Photon_pt->at(p) < 30.) continue;
		if(fabs(_base->Photon_eta->at(p)) > 1.5) continue;
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
		vector<unsigned int> rhs = _base->Photon_rhIds->at(p);
                for(int r = 0; r < rhs.size(); r++){
                        unsigned int rhid = rhs[r];
                        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				//TOF from 0 to rh location
				drh = hypo(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), _base->ECALRecHit_rhz->at(rhidx))/_c;
				//dont need (would bring photon into PV frame) - TOF from PV to rh location
				timecorr = drh;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
                               	//undo adjustment currently made in ntuples (traw - drh/c)
				//TODO: UNDO THIS LATER WHEN NTUPLES ARE UPDATED
				JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx) + timecorr + drh);
                               
				//rec hit selection
				if(fabs(rh.t()) > 20) continue;
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
                                pho.AddRecHit(rh);
                        }

                }
		if(pho.GetNRecHits() < 2) continue;
		phos.push_back(pho);
        }


}
