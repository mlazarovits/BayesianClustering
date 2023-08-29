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
	int nRHs, nRHs_evt, nphotons;
	rhs.clear();
	for(int i = 0; i < m_nEvts; i++){
		m_base->GetEntry(i);
		nphotons = (int)m_base->Photon_rhIds->size();
		nRHs_evt = (int)m_base->ECALRecHit_ID->size();
		rhs.push_back({});
		for(int p = 0; p < nphotons; p++){
			nRHs = (int)m_base->Photon_rhIds->at(p).size();
			unsigned long long id;
			for(int r = 0; r < nRHs; r++){
				//add tof = d_pv to time to get correct RH time
				//t = rh_time - d_rh/c + d_pv/c
				id = m_base->Photon_rhIds->at(p).at(r);
				rh.SetRecHitId(id);
				for(int j = 0; j < nRHs_evt; j++){
					if(m_base->ECALRecHit_ID->at(j) == id){
						rh = JetPoint(m_base->ECALRecHit_rhx->at(j), m_base->ECALRecHit_rhy->at(j), m_base->ECALRecHit_rhz->at(j), m_base->ECALRecHit_time->at(j)+m_base->ECALRecHit_TOF->at(j));
						rh.SetEnergy(m_base->ECALRecHit_energy->at(j));
						rh.SetEta(m_base->ECALRecHit_eta->at(j));
						rh.SetPhi(m_base->ECALRecHit_phi->at(j));
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
	for(int i = 0; i < m_nEvts; i++){
		if(i == evt){
			m_base->GetEntry(i);
			nphotons = (int)m_base->Photon_rhIds->size();
			nRHs_evt = (int)m_base->ECALRecHit_ID->size();
			
			for(int p = 0; p < nphotons; p++){
				nRHs = (int)m_base->Photon_rhIds->at(p).size();
				unsigned long long id;
				for(int r = 0; r < nRHs; r++){
					//add tof = d_pv to time to get correct RH time
					//t = rh_time - d_rh/c + d_pv/c
					id = m_base->Photon_rhIds->at(p).at(r);
					rh.SetRecHitId(id);
					for(int j = 0; j < nRHs_evt; j++){
						if(m_base->ECALRecHit_ID->at(j) == id){
							rh = JetPoint(m_base->ECALRecHit_rhx->at(j), m_base->ECALRecHit_rhy->at(j), m_base->ECALRecHit_rhz->at(j), m_base->ECALRecHit_time->at(j)+m_base->ECALRecHit_TOF->at(j));
							rh.SetEnergy(m_base->ECALRecHit_energy->at(j));
							rh.SetEta(m_base->ECALRecHit_eta->at(j));
							rh.SetPhi(m_base->ECALRecHit_phi->at(j));
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
	for(int i = 0; i < m_nEvts; i++){
		if(i == evt){
			m_base->GetEntry(i);
			//make sure photon number is in vector
			if(pho >= (int)m_base->Photon_rhIds->size()) return;
			nRHs = (int)m_base->Photon_rhIds->at(pho).size();
			nRHs_evt = (int)m_base->ECALRecHit_ID->size();
			unsigned long long id;
			for(int r = 0; r < nRHs; r++){
				//add tof = d_pv to time to get correct RH time
				//t = rh_time - d_rh/c + d_pv/c
				id = m_base->Photon_rhIds->at(pho).at(r);
				rh.SetRecHitId(id);
				for(int j = 0; j < nRHs_evt; j++){
					if(m_base->ECALRecHit_ID->at(j) == id){
						rh = JetPoint(m_base->ECALRecHit_rhx->at(j), m_base->ECALRecHit_rhy->at(j), m_base->ECALRecHit_rhz->at(j), m_base->ECALRecHit_time->at(j)+m_base->ECALRecHit_TOF->at(j));
						rh.SetEnergy(m_base->ECALRecHit_energy->at(j));
						rh.SetEta(m_base->ECALRecHit_eta->at(j));
						rh.SetPhi(m_base->ECALRecHit_phi->at(j));
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



//make cluster param histograms
void PhotonProducer::Skim(){
	TFile* ofile = new TFile("plots/photon_skims.root","RECREATE");

	int nPho;
	for(int i = 0; i < m_nEvts; i++){
		m_base->GetEntry(i);
		//nPho = m_base->

	}





}
