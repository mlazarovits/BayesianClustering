#include "JetProducer.hh"

#include "Clusterizer.hh"
#include "Matrix.hh"
#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>




JetProducer::JetProducer(){ };



JetProducer::~JetProducer(){ 
	_file->Close();
	delete _base;
	delete _file;
}


JetProducer::JetProducer(TFile* file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	_file = file;
	TTree* tree = (TTree*)file->Get("tree/llpgtree");
	_base = new ReducedBase(tree);
	_nEvts = _base->fChain->GetEntries();
}

void JetProducer::GetRecHits(vector<vector<JetPoint>>& rhs){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs;
	rhs.clear();
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		//TODO: switch to _base->nRHs when that's in the ntuples
		nRHs = (int)_base->ECALRecHit_ID->size();
		rhs.push_back({});
		for(int r = 0; r < nRHs; r++){
			//add tof = d_pv to time to get correct RH time
			//t = rh_time - d_rh/c + d_pv/c
			rh = JetPoint(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
			
			rh.SetEnergy(_base->ECALRecHit_energy->at(r));
			rh.SetEta(_base->ECALRecHit_eta->at(r));
			rh.SetPhi(_base->ECALRecHit_phi->at(r));
			rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
	
			rhs[i].push_back(rh);
		}
	}
}


void JetProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs;
	rhs.clear();
	double etaMax = 0.5;
	double etaMin = -etaMax;  
	double phiMax = 2.;
	double phiMin = -2.8;
	int cnt = 0;
	int nPho, nPhoRHs, rhId;
	vector<unsigned int> phoIds;
	

	if(evt > _nEvts) return;

	_base->GetEntry(evt);

	nRHs = (int)_base->ECALRecHit_ID->size();
	//clean for RHs in photon superclusters
	nPho = (int)_base->Photon_energy->size();
	for(int p = 0; p < nPho; p++){
		nPhoRHs = (int)_base->Photon_rhIds->at(p).size();
		for(int r = 0; r < nPhoRHs; r++) phoIds.push_back(_base->Photon_rhIds->at(p).at(r));
	}
	
	for(int r = 0; r < nRHs; r++){
		//clean out photon ids - if rh id is in phoIds continue
		rhId = _base->ECALRecHit_ID->at(r);	
		if(std::any_of(phoIds.begin(), phoIds.end(), [&](int iid){return iid == rhId;})) continue; 
	
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		rh = JetPoint(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
	
		rhs.push_back(rh);
	}	
}

//ctor from rec hit collection - integrating into ntuplizer - in CMSSW

void JetProducer::GetPrimaryVertex(Point& vtx, int evt){
	//reset to empty 3-dim point	
	vtx = Point(3);

	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);

}


