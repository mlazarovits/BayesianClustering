#include "PhotonProducer.hh"

PhotonProducer::PhotonProducer(){ };



PhotonProducer::~PhotonProducer(){ 
	m_file->Close();
	delete m_base;
	delete m_file;
}


PhotonProducer::PhotonProducer(TFile* file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	m_file = file;
	TTree* tree = (TTree*)file->Get("tree/llpgtree");
	m_base = new ReducedBase(tree);
	m_nEvts = m_base->fChain->GetEntries();
}

void PhotonProducer::GetRecHits(vector<vector<JetPoint>>& rhs){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs;
	rhs.clear();
	for(int i = 0; i < m_nEvts; i++){
		m_base->GetEntry(i);
		//TODO: switch to m_base->nRHs when that's in the ntuples
		nRHs = (int)m_base->ERH_time->size();
		rhs.push_back({});
		for(int r = 0; r < nRHs; r++){
			//add tof = d_pv to time to get correct RH time
			//t = rh_time - d_rh/c + d_pv/c
			rh = JetPoint(m_base->ERH_x->at(r), m_base->ERH_y->at(r), m_base->ERH_z->at(r), m_base->ERH_time->at(r)+m_base->ERH_TOF->at(r));
			
			rh.SetEnergy(m_base->ERH_energy->at(r));
			rh.SetEta(m_base->ERH_eta->at(r));
			rh.SetPhi(m_base->ERH_phi->at(r));
			rh.SetRecHitId(m_base->ERH_ID->at(r));
	
			rhs[i].push_back(rh);
		}
	}
}


void PhotonProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
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
	for(int i = 0; i < m_nEvts; i++){
		if(i == evt){
			m_base->GetEntry(i);
			nRHs = (int)m_base->Photon_rhIds->size();
			nRHs_evt = (int)m_base->ERH_ID->size();
			unsigned long long id;
			for(int r = 0; r < nRHs; r++){
				//add tof = d_pv to time to get correct RH time
				//t = rh_time - d_rh/c + d_pv/c
				id = m_base->Photon_rhIds->at(r);
				rh.SetRecHitId(id);
				for(int j = 0; j < nRHs_evt; j++){
					if(m_base->ERH_ID->at(j) == id){
						rh = JetPoint(m_base->ERH_x->at(j), m_base->ERH_y->at(j), m_base->ERH_z->at(j), m_base->ERH_time->at(j)+m_base->ERH_TOF->at(j));
						rh.SetEnergy(m_base->ERH_energy->at(j));
						rh.SetEta(m_base->ERH_eta->at(j));
						rh.SetPhi(m_base->ERH_phi->at(j));
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

	for(int i = 0; i < m_nEvts; i++){
		if(i == evt){
			m_base->GetEntry(i);
			vtx.SetValue(m_base->PV_x, 0);
			vtx.SetValue(m_base->PV_y, 1);
			vtx.SetValue(m_base->PV_z, 2);
			return;
		}
		else continue;	
	}

}


//make ctor that simulates rechits - see src/varGMM.C





