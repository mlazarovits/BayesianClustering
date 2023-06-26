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
		m_rechits.push_back({});
		for(int r = 0; r < nRHs; r++){
	
			//add tof = d_pv to time to get correct RH time
			//jt = Jet(base->ERH_x, base->ERH_y, base->ERH_z, base->ERH_time+base->ERH_TOF);
			
			jt.SetEnergy(base->ERH_energy->at(r));
			jt.SetEta(base->ERH_eta->at(r));
			jt.SetPhi(base->ERH_phi->at(r));
			jt.SetRecHitId(base->ERH_ID->at(r));
	
			m_rechits[i].push_back(jt);
		}
	}
	file->Close();
}


//ctor from rec hit collection - integrating into ntuplizer - in CMSSW


//make ctor that simulates rechits - see src/varGMM.C



//ctor for alice's pixels


