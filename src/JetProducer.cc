#include "JetProducer.hh"
#include "ReducedBase.hh"
#include "Jet.hh"


JetProducer::JetProducer(){ };



JetProducer::~JetProducer(){ };



JetProducer::JetProducer(TFile* file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	TTree* tree = (TTree*)file->Get("tree/llpgtree");
	ReducedBase* base = new ReducedBase(tree);
	
	int nEvents = base->fChain->GetEntries();


	Jet jt;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs;
	for(int i = 0; i < nEvents; i++){
		base->GetEntry(i);
		//TODO: switch to base->nRHs when that's in the ntuples
		nRHs = (int)base->ERH_time->size();
		for(int r = 0; r < nRHs; r++){
			//jt = Jet(base->ERH_x->at(r), base->ERH_y->at(r), base->ERH_z->at(r), base->ERH_time->at(r));
			jt.SetEnergy(base->ERH_energy->at(r));
			jt.SetEta(base->ERH_eta->at(r));
			jt.SetPhi(base->ERH_phi->at(r));
			jt.SetRecHitId(base->ERH_ID->at(r));
		
			m_rechits.push_back(jt);
		}
	}
	file->Close();
}


//ctor from rec hit collection - integrating into ntuplizer

//make ctor that simulates rechits




