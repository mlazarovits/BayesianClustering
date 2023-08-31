#include "PhotonProducer.hh"
#include "Clusterizer.hh"
#include "Matrix.hh"

#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>
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
						//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
						rh = JetPoint(m_base->ECALRecHit_rhx->at(j), m_base->ECALRecHit_rhy->at(j), m_base->ECALRecHit_rhz->at(j), m_base->ECALRecHit_time->at(j)+m_base->ECALRecHit_TOF->at(j));
						rh.SetEnergy(m_base->ECALRecHit_energy->at(j));
						rh.SetEta(m_base->ECALRecHit_eta->at(j));
						rh.SetPhi(m_base->ECALRecHit_phi->at(j));
						
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
							//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
							rh = JetPoint(m_base->ECALRecHit_rhx->at(j), m_base->ECALRecHit_rhy->at(j), m_base->ECALRecHit_rhz->at(j), m_base->ECALRecHit_time->at(j)+m_base->ECALRecHit_TOF->at(j));
							rh.SetEnergy(m_base->ECALRecHit_energy->at(j));
							rh.SetEta(m_base->ECALRecHit_eta->at(j));
							rh.SetPhi(m_base->ECALRecHit_phi->at(j));
							
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
						//time = ECALRecHit_time + TOF = (rh_time - d_rh/c) + TOF
						rh = JetPoint(m_base->ECALRecHit_rhx->at(j), m_base->ECALRecHit_rhy->at(j), m_base->ECALRecHit_rhz->at(j), m_base->ECALRecHit_time->at(j)+m_base->ECALRecHit_TOF->at(j));
						rh.SetEnergy(m_base->ECALRecHit_energy->at(j));
						rh.SetEta(m_base->ECALRecHit_eta->at(j));
						rh.SetPhi(m_base->ECALRecHit_phi->at(j));
						
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


void PhotonProducer::CleaningSkim(){
	TFile* ofile = new TFile("plots/photon_cleaningSkims.root","RECREATE");
	//rh time
	TH1D* t_rh = new TH1D("t_rh","t_rh",50,-150, 150);
	//rh time vs. rh e
	TH2D* TvErh = new TH2D("TvErh","TvErh",50, -150, 150., 1000, 0., 1000);
	TH2D* TvErh_lowE = new TH2D("TvErh_lowE","TvErh_lowE",50, -150, 150., 1000, 0., 10.);
	TH2D* TvErh_bx = new TH2D("TvErh_bx","TvErh_bx",50, -150, 150., 100, 0., 1000);
	TH2D* TvErh_cut = new TH2D("TvErh_cut","TvErh_cut",50, -100, 100., 1000, 0., 1000);

	TH1D* phoE = new TH1D("phoE","phoE",1000, 0, 1000);

	int nPho;
	vector<JetPoint> rhs;
	for(int i = 0; i < m_nEvts; i++){
		m_base->GetEntry(i);
		nPho = (int)m_base->Photon_energy->size();
		for(int p = 0; p < nPho; p++){
			phoE->Fill(m_base->Photon_energy->at(p));
			//find subclusters for each photon
			GetRecHits(rhs, i, p);
			//cout << "evt #" << i << " photon #" << p << " nrhs: " << rhs.size() << endl;
			for(int r = 0; r < (int)rhs.size(); r++){
				t_rh->Fill(rhs[r].t());
				TvErh->Fill(rhs[r].t(), rhs[r].E());
				if(rhs[r].E() <= 10.0) TvErh_lowE->Fill(rhs[r].t(), rhs[r].E());
				if(fabs(rhs[r].t()) < 25) TvErh_bx->Fill(rhs[r].t(), rhs[r].E());
				if(rhs[r].E() < 3.0 && fabs(rhs[r].t()) > 50) continue;
				TvErh_cut->Fill(rhs[r].t(), rhs[r].E());

			}
		}
	}
	
	ofile->cd();
	ofile->WriteTObject(t_rh);
	ofile->WriteTObject(TvErh);
	ofile->WriteTObject(TvErh_lowE);
	ofile->WriteTObject(TvErh_bx);
	ofile->WriteTObject(TvErh_cut);
	ofile->WriteTObject(phoE);
	ofile->Close();

}



//make cluster param histograms
void PhotonProducer::Skim(){
	TFile* ofile = new TFile("plots/photon_skims.root","RECREATE");
	

	vector<TH1D*> TH1D_hists;
	//subcluster energy - average
	TH1D* e_avg = new TH1D("e_avg","e_avg",100,0.,50.);
	TH1D_hists.push_back(e_avg);
	//space slope
	TH1D* slope_space = new TH1D("slope_space","slope_space",50,-30,30);
	TH1D_hists.push_back(slope_space);
	//eta-time slope
	TH1D* slope_etaT = new TH1D("slope_etaT","slope_etaT",50,-2,2);
	TH1D_hists.push_back(slope_etaT);
	//phi-time slop
	TH1D* slope_phiT = new TH1D("slope_phiT","slope_phiT",50,-4,4);
	TH1D_hists.push_back(slope_phiT);
	//mean time - center in t
	TH1D* time_center = new TH1D("time_center","time_center",50,-30,30);
	TH1D_hists.push_back(time_center);
	//mean eta - center in eta
	TH1D* eta_center = new TH1D("eta_center","eta_center",50,-3.5,3.5);
	TH1D_hists.push_back(eta_center);
	//mean phi - center in phi
	TH1D* phi_center = new TH1D("phi_center","phi_center",50,-3.5,3.5);
	TH1D_hists.push_back(phi_center);
	
	//# of subclusters
	TH1I* nSubClusters = new TH1I("nSubClusters","nSubClusters",7,0,7.);
	//# of subclusters vs. photon reco energy
	TH2D* e_nSubClusters = new TH2D("e_nSubClusters","e_nSubClusters",50,0.,1000.,7,0.,7.);	

	

	int nPho;
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = -log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta,0,0);
	smear.SetEntry(dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	
	
	Clusterizer* algo = new Clusterizer();
	algo->SetAlpha(0.1);
	algo->SetThresh(1.);
	algo->SetMaxNClusters(5);
	algo->SetWeighted(true);
	algo->SetVerbosity(0);
	//algo->SetDataSmear(smear);


	GaussianMixture* gmm = new GaussianMixture();
	
	map<string, Matrix> params;

	vector<JetPoint> rhs;
	int nclusters;
	vector<double> eigenvals, avg_Es;
	vector<Matrix> eigenvecs;
	for(int i = 0; i < m_nEvts; i++){
		m_base->GetEntry(i);
		nPho = (int)m_base->Photon_energy->size();
		for(int p = 0; p < nPho; p++){
			//find subclusters for each photon
			GetRecHits(rhs, i, p);
			cout << "evt: " << i << " of " << m_nEvts << "  pho: " << p << " nrhs: " << rhs.size() << "\r" << flush;
	
			gmm = algo->FindSubjets(Jet(rhs));
			
			nclusters = gmm->GetNClusters();
			nSubClusters->Fill(nclusters);
			e_nSubClusters->Fill(m_base->Photon_energy->at(p), nclusters);
			
			gmm->GetAvgWeights(avg_Es);		


			for(int k = 0; k < nclusters; k++){
				params = gmm->GetParameters(k);
				eta_center->Fill(params["mean"].at(0,0));
				phi_center->Fill(params["mean"].at(1,0));
				time_center->Fill(params["mean"].at(2,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//largest eigenvalue is last
				//phi/eta
				slope_space->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(0,0));
        			//eta/time
				slope_etaT->Fill(eigenvecs[2].at(0,0)/eigenvecs[2].at(2,0));
				//phi/time
				slope_phiT->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(2,0));

				//average cluster energy
				e_avg->Fill(avg_Es[k]);
			}
	
			rhs.clear();
			eigenvals.clear();
			eigenvecs.clear();
		}

	}
	ofile->cd();
	for(int i = 0; i < (int)TH1D_hists.size(); i++) TH1D_hists[i]->Write();
	nSubClusters->Write();
	e_nSubClusters->Write();

}
