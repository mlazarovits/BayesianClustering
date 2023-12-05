#include "BaseProducer.hh"



void BaseProducer::GetTrueJets(vector<Jet>& jets, int evt){
        double px, py, pz, pt, phi, eta;
        jets.clear();

        if(evt > _nEvts) return;

        _base->GetEntry(evt);
        int nJets = (int)_base->Jet_energy->size();
        int nrhs, rhidx;
	bool jetid;

	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;

	Point vtx(3);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);
	for(int j = 0; j < nJets; j++){
                pt = _base->Jet_pt->at(j);
                phi = _base->Jet_phi->at(j);
                eta = _base->Jet_eta->at(j);

                px = pt*cos(phi);
                py = pt*sin(phi);
                pz = pt*sinh(eta);

                nrhs = _base->Jet_drRhIds->size();

		//Jet selection
		if(nrhs < 2) continue;
                if(_base->Jet_pt->at(j) < 30.) continue;
                if(fabs(_base->Jet_eta->at(j)) > 1.5) continue;

		//create jet id (based on tight 2017 Run II recommendations)
		//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
		
		jetid = (_base->Jet_neEmEF->at(j) < 0.9) && (_base->Jet_neHEF->at(j) < 0.9) && (_base->Jet_nConstituents->at(j) > 1)
                        && (_base->Jet_chHEF->at(j) > 0.0) && (_base->Jet_chHM->at(j) > 0);
                if(!jetid) continue;



                Jet jet(px, py, pz, _base->Jet_energy->at(j));
                jet.SetVertex(vtx);
		//set rec hits in jet
		vector<unsigned int> rhs = _base->Jet_drRhIds->at(j);
                vector<unsigned int> egidxs = _base->Jet_egIndxs->at(j);
                for(int r = 0; r < rhs.size(); r++){
                        unsigned int rhid = rhs[r];
                        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
                                JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx),
                                        _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx));
                                rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*_gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
                                jet.AddRecHit(rh);
                        }

                }
                jets.push_back(jet);
        }


}
