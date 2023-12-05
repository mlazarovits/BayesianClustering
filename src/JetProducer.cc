#include "JetProducer.hh"

#include "Clusterizer.hh"
#include "Matrix.hh"
#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>



JetProducer::JetProducer(){ 
}



JetProducer::~JetProducer(){ 
//	_file->Close();
//	delete _base;
//	delete _file;
}


JetProducer::JetProducer(TFile* file) : BaseProducer(file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	//_file = file;
	//cout << "b" << endl;
	//cout << "file name: " << file->GetName() << endl;
	//TTree* tree = (TTree*)file->Get("tree/llpgtree");
	//cout << "c" << endl;
	//_base = new ReducedBase(tree);
	//cout << "d" << endl;
	//_nEvts = _base->fChain->GetEntries();

}



void JetProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned int rhId;
	int nRHs;
	rhs.clear();
	

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nJets = (int)_base->Jet_energy->size();

	//actually get rhs for clustering
	nRHs = (int)_base->ECALRecHit_energy->size();
	vector<double> ws;
	//need to transfer from GeV (energy) -> unitless (number of points
	
	for(int r = 0; r < nRHs; r++){
		//clean out photon ids - if rh id is in phoIds continue
		rhId = _base->ECALRecHit_ID->at(r);
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		rh = JetPoint(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
		rh.SetWeight(_base->ECALRecHit_energy->at(r)*_gev);
	
		rhs.push_back(rh);
	}	
}




void JetProducer::GetRecHits(vector<Jet>& jets, int evt){
	double t, E, eta, phi;
	unsigned int rhId;
	jets.clear();

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nRHs = (int)_base->ECALRecHit_ID->size();

	Point vtx = Point(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	//make weights - E/e_avg
	vector<double> ws;

	for(int r = 0; r < nRHs; r++){
		//clean out photon ids - if rh id is in phoIds continue
		rhId = _base->ECALRecHit_ID->at(r);
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		JetPoint rh(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
		rh.SetWeight(_base->ECALRecHit_energy->at(r)*_gev);

		Jet j(rh);
		j.SetVertex(vtx);
		jets.push_back(j);
	}	
}


void JetProducer::GetPrimaryVertex(Point& vtx, int evt){
	//reset to empty 3-dim point	
	vtx = Point(3);

	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);

}




