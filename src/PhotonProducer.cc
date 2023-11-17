#include "PhotonProducer.hh"
#include "Matrix.hh"

#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>
PhotonProducer::PhotonProducer(){ };



PhotonProducer::~PhotonProducer(){ 
//	delete _base;
//	delete _file;
}


PhotonProducer::PhotonProducer(TFile* file) : BaseProducer(file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	//TTree* tree = (TTree*)file->Get("tree/llpgtree");
	//_base = new ReducedBase(tree);
	//_nEvts = _base->fChain->GetEntries();

}


//get rechits for all photons in an event
void PhotonProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs, nphotons, nRHs_evt;
	rhs.clear();
	double etaMax = 0.5;
	double etaMin = -etaMax;  
	double phiMax = 2.;
	double phiMin = -2.8;
	int cnt = 0;
	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	nphotons = (int)_base->Photon_rhIds->size();
	nRHs_evt = (int)_base->ECALRecHit_ID->size();
	
	for(int p = 0; p < nphotons; p++){
		nRHs = (int)_base->Photon_rhIds->at(p).size();
		unsigned long long id;
		for(int r = 0; r < nRHs; r++){
			//add tof = d_pv to time to get correct RH time
			//t = rh_time - d_rh/c + d_pv/c
			id = _base->Photon_rhIds->at(p).at(r);
			for(int j = 0; j < nRHs_evt; j++){
				if(_base->ECALRecHit_ID->at(j) == id){
					//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
					JetPoint rh(_base->ECALRecHit_rhx->at(j), _base->ECALRecHit_rhy->at(j), _base->ECALRecHit_rhz->at(j), _base->ECALRecHit_time->at(j)+_base->ECALRecHit_TOF->at(j));
					rh.SetEnergy(_base->ECALRecHit_energy->at(j));
					rh.SetEta(_base->ECALRecHit_eta->at(j));
					rh.SetPhi(_base->ECALRecHit_phi->at(j));
					rh.SetRecHitId(id);
					rh.SetWeight(_base->ECALRecHit_energy->at(j)*_gev);	
					
					//cleaning cuts
					if(!cleanRH(rh)) break;
				
					rhs.push_back(rh);
					break;
				}
				else continue;
			}
		}
	
	}
}

//get rec hits for a particular photon for an event
void PhotonProducer::GetRecHits(vector<JetPoint>& rhs, int evt, int pho){
	double x, y, z, t, E, eta, phi;
	int nRHs, nRHs_evt;
	rhs.clear();
	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	//make sure photon number is in vector
	if(pho >= (int)_base->Photon_rhIds->size()) return;
	nRHs = (int)_base->Photon_rhIds->at(pho).size();
	nRHs_evt = (int)_base->ECALRecHit_ID->size();
	unsigned int id;
	
	
	for(int r = 0; r < nRHs; r++){
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		id = _base->Photon_rhIds->at(pho).at(r);
		for(int j = 0; j < nRHs_evt; j++){
			if(_base->ECALRecHit_ID->at(j) == id){
				//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
				JetPoint rh(_base->ECALRecHit_rhx->at(j), _base->ECALRecHit_rhy->at(j), _base->ECALRecHit_rhz->at(j), _base->ECALRecHit_time->at(j)+_base->ECALRecHit_TOF->at(j));
				rh.SetEnergy(_base->ECALRecHit_energy->at(j));
				rh.SetEta(_base->ECALRecHit_eta->at(j));
				rh.SetPhi(_base->ECALRecHit_phi->at(j));
				rh.SetRecHitId(id);
				rh.SetWeight(_base->ECALRecHit_energy->at(j)*_gev);	
				
				//cleaning cuts
				if(!cleanRH(rh)) break;
				
				rhs.push_back(rh);
				break;
			}
		}
	
	}
}

//get rec hits for a particular photon for an event
void PhotonProducer::GetRecHits(vector<Jet>& rhs, int evt, int pho){
	double x, y, z, t, E, eta, phi;
	rhs.clear();
	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	//make sure photon number is in vector
	if(pho >= (int)_base->Photon_rhIds->size()) return;
	int nRHs = (int)_base->Photon_rhIds->at(pho).size();
	int nRHs_evt = (int)_base->ECALRecHit_ID->size();
	unsigned int id;

	//set vertex info
	Point vtx = Point(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);

	for(int r = 0; r < nRHs; r++){
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		id = _base->Photon_rhIds->at(pho).at(r);
		for(int j = 0; j < nRHs_evt; j++){
			if(_base->ECALRecHit_ID->at(j) == id){
				//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
				JetPoint rh(_base->ECALRecHit_rhx->at(j), _base->ECALRecHit_rhy->at(j), _base->ECALRecHit_rhz->at(j), _base->ECALRecHit_time->at(j)+_base->ECALRecHit_TOF->at(j));
				rh.SetEnergy(_base->ECALRecHit_energy->at(j));
				rh.SetEta(_base->ECALRecHit_eta->at(j));
				rh.SetPhi(_base->ECALRecHit_phi->at(j));
				rh.SetRecHitId(id);
				rh.SetWeight(_base->ECALRecHit_energy->at(j)*_gev);	
				//cleaning cuts
				if(!cleanRH(rh)) break;
			
				Jet jet(rh);
				jet.SetVertex(vtx);	
				rhs.push_back(jet);
				break;
			}
		}
	
	}
}

//ctor from rec hit collection - integrating into ntuplizer - in CMSSW

void PhotonProducer::GetPrimaryVertex(Point& vtx, int evt){
	//reset to empty 3-dim point	
	vtx = Point(3);
	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);

}


