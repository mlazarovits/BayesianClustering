#include "PhotonProducer.hh"
#include "Clusterizer.hh"
#include "Matrix.hh"

#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>
PhotonProducer::PhotonProducer(){ };



PhotonProducer::~PhotonProducer(){ 
	_file->Close();
	delete _base;
	delete _file;
}


PhotonProducer::PhotonProducer(TFile* file) : BaseProducer(file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	//_file = file;
	//TTree* tree = (TTree*)file->Get("tree/llpgtree");
	//_base = new ReducedBase(tree);
	//_nEvts = _base->fChain->GetEntries();

}

void PhotonProducer::GetRecHits(vector<vector<JetPoint>>& rhs){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs, nRHs_evt, nphotons;
	rhs.clear();
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		nphotons = (int)_base->Photon_rhIds->size();
		nRHs_evt = (int)_base->ECALRecHit_ID->size();
		rhs.push_back({});
		for(int p = 0; p < nphotons; p++){
			nRHs = (int)_base->Photon_rhIds->at(p).size();
			unsigned long long id;
			for(int r = 0; r < nRHs; r++){
				//add tof = d_pv to time to get correct RH time
				//t = rh_time - d_rh/c + d_pv/c
				id = _base->Photon_rhIds->at(p).at(r);
				rh.SetRecHitId(id);
				for(int j = 0; j < nRHs_evt; j++){
					if(_base->ECALRecHit_ID->at(j) == id){
						//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
						rh = JetPoint(_base->ECALRecHit_rhx->at(j), _base->ECALRecHit_rhy->at(j), _base->ECALRecHit_rhz->at(j), _base->ECALRecHit_time->at(j)+_base->ECALRecHit_TOF->at(j));
						rh.SetEnergy(_base->ECALRecHit_energy->at(j));
						rh.SetEta(_base->ECALRecHit_eta->at(j));
						rh.SetPhi(_base->ECALRecHit_phi->at(j));
						
						//cleaning cuts
						if(!cleanRH(rh)) break;
						
						rhs[i].push_back(rh);
						break;
					}
					else continue;
				}
			}
		}
	}
}

//get rechits for all photons in an event
void PhotonProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs, nphotons, nRHs_evt;
	rhs.clear();
	double etaMax = 0.5;
	double etaMin = -etaMax;  
	double phiMax = 2.;
	double phiMin = -2.8;
	int cnt = 0;
	for(int i = 0; i < _nEvts; i++){
		if(i == evt){
			_base->GetEntry(i);
			nphotons = (int)_base->Photon_rhIds->size();
			nRHs_evt = (int)_base->ECALRecHit_ID->size();
			
			for(int p = 0; p < nphotons; p++){
				nRHs = (int)_base->Photon_rhIds->at(p).size();
				unsigned long long id;
				for(int r = 0; r < nRHs; r++){
					//add tof = d_pv to time to get correct RH time
					//t = rh_time - d_rh/c + d_pv/c
					id = _base->Photon_rhIds->at(p).at(r);
					rh.SetRecHitId(id);
					for(int j = 0; j < nRHs_evt; j++){
						if(_base->ECALRecHit_ID->at(j) == id){
							//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
							rh = JetPoint(_base->ECALRecHit_rhx->at(j), _base->ECALRecHit_rhy->at(j), _base->ECALRecHit_rhz->at(j), _base->ECALRecHit_time->at(j)+_base->ECALRecHit_TOF->at(j));
							rh.SetEnergy(_base->ECALRecHit_energy->at(j));
							rh.SetEta(_base->ECALRecHit_eta->at(j));
							rh.SetPhi(_base->ECALRecHit_phi->at(j));
							
							//cleaning cuts
							if(!cleanRH(rh)) break;
						
							rhs.push_back(rh);
							break;
						}
						else continue;
					}
				}
			}
			return;
		}
		else continue;
	}
}

//get rec hits for a particular photon for an event
void PhotonProducer::GetRecHits(vector<JetPoint>& rhs, int evt, int pho){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs, nRHs_evt;
	rhs.clear();
	for(int i = 0; i < _nEvts; i++){
		if(i == evt){
			_base->GetEntry(i);
	//		cout << (int)_base->Photon_energy->size() << " nphotons in evt " << evt << endl;
			//make sure photon number is in vector
			if(pho >= (int)_base->Photon_rhIds->size()) return;
			nRHs = (int)_base->Photon_rhIds->at(pho).size();
			nRHs_evt = (int)_base->ECALRecHit_ID->size();
			unsigned long long id;
			for(int r = 0; r < nRHs; r++){
				//add tof = d_pv to time to get correct RH time
				//t = rh_time - d_rh/c + d_pv/c
				id = _base->Photon_rhIds->at(pho).at(r);
				rh.SetRecHitId(id);
				for(int j = 0; j < nRHs_evt; j++){
					if(_base->ECALRecHit_ID->at(j) == id){
						//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
						rh = JetPoint(_base->ECALRecHit_rhx->at(j), _base->ECALRecHit_rhy->at(j), _base->ECALRecHit_rhz->at(j), _base->ECALRecHit_time->at(j)+_base->ECALRecHit_TOF->at(j));
						rh.SetEnergy(_base->ECALRecHit_energy->at(j));
						rh.SetEta(_base->ECALRecHit_eta->at(j));
						rh.SetPhi(_base->ECALRecHit_phi->at(j));
						
						//cleaning cuts
						if(!cleanRH(rh)) break;
						
						rhs.push_back(rh);
						break;
					}
					else continue;
				}
		
			}
			return;
		}
		else continue;
	}
}

//ctor from rec hit collection - integrating into ntuplizer - in CMSSW

void PhotonProducer::GetPrimaryVertex(Point& vtx, int evt){
	//reset to empty 3-dim point	
	vtx = Point(3);

	for(int i = 0; i < _nEvts; i++){
		if(i == evt){
			_base->GetEntry(i);
			vtx.SetValue(_base->PV_x, 0);
			vtx.SetValue(_base->PV_y, 1);
			vtx.SetValue(_base->PV_z, 2);
			return;
		}
		else continue;	
	}

}


