#ifndef JETSIMPRODUCER_HH
#define JETSIMPRODUCER_HH
#include "ReducedBaseSim.hh"
#include "Jet.hh"
#include "JetPoint.hh"
#include "TH2D.h"
#include "TFile.h"

class JetSimProducer{
	public:
		JetSimProducer();
		virtual ~JetSimProducer();

		JetSimProducer(TFile* file);

		void GetRecHits(vector<JetPoint>& rhs, int evt);
		void GetRecHits(vector<Jet>& rhs, int evt);
		void GetGenJets(vector<Jet>& genAK4jets, vector<Jet>& genAK8jets, vector<Jet>& genAK15jets, int evt);
		void GetGenParticles(vector<Jet>& genparts, int evt);
		void GetRecoJets(vector<Jet>& recoAK4jets, vector<Jet>& recoAK8jets, vector<Jet>& recoAK15jets, int evt);
		void GetPrimaryVertex(BayesPoint& vtx, int evt);
		ReducedBaseSim* GetBase(){ return _base; }
		void SetTransferFactor(double g){ _gev = g; }

		void SetRecoMinPt(double pt){_ptmin_reco = pt; }
		void SetRecoMinE(double e){ _Emin_reco = e;}
		void SetGenMinPt(double pt){_ptmin_gen = pt; }
		void SetGenMinE(double e){ _Emin_gen = e;}
		void SetMinRhE(double r){ _minrhE = r; }
		void SetMinNrhs(int r){ _minNrhs = r; }
		void SetMinNGenConsts(int c){ _nConstsmin = c; }
		void SortJets(vector<Jet>& jets);
		
		void PrintPreselection(){
			cout << "Default energy transfer factor: " << _gev << endl;
			cout << "Minimum reco pt: " << _ptmin_reco << endl;
			cout << "Minimum reco energy: " << _Emin_reco << endl;
			cout << "Minimum gen pt: " << _ptmin_gen << endl;
			cout << "Minimum gen energy: " << _Emin_gen << endl;
			cout << "Minimum rh energy: " << _minrhE << endl;
			cout << "Minimum # rh: " << _minNrhs << endl;
			cout << "Minimum # constituents (gen jets only): " << _nConstsmin << endl;
		}


		void EtaPhiMap(string fname, vector<Jet>& rhs){
			TFile* f = new TFile(fname.c_str(),"RECREATE");
			f->cd();

			double maxEta = -1.5;
			double minEta = 1.5;
	
			//double maxPhi = -4*atan(1); //pi
			//double minPhi = 4*atan(1);
			double maxPhi = 0; //pi
			double minPhi = 8*atan(1);
			double cell = 4*atan(1)/180;
			for(int r = 0; r < rhs.size(); r++){
				if(rhs[r].eta() > maxEta)
					maxEta = rhs[r].eta();
				if(rhs[r].eta() < minEta)
					minEta = rhs[r].eta();
				
				//if(rhs[r].phi_std() > maxPhi)
				//	maxPhi = rhs[r].phi_std();
				//if(rhs[r].phi_std() < minPhi)
				//	minPhi = rhs[r].phi_std();
				if(rhs[r].phi() > maxPhi)
					maxPhi = rhs[r].phi();
				if(rhs[r].phi() < minPhi)
					minPhi = rhs[r].phi();
			}
			cout << "minEta " << minEta << " maxEta " << maxEta << " minPhi " << minPhi << " maxPhi " << maxPhi << endl;	
			TH2D* etaPhiMap = new TH2D("etaPhiMap","etaPhiMap;eta;phi",(int)rhs.size(),minEta-cell,maxEta+cell,(int)rhs.size(),minPhi-cell,maxPhi+cell);
			double totE = 0;		
			for(int r = 0; r < rhs.size(); r++){
				cout << "r " << r << " eta " << rhs[r].eta() << " phi " << rhs[r].phi() << " E " << rhs[r].E() << endl;
				etaPhiMap->Fill(rhs[r].eta(), rhs[r].phi(), rhs[r].E());
				totE += rhs[r].E();
			}
			cout << "total E " << totE << " scaled total E " << totE*_gev << endl;
			etaPhiMap->Write();
			f->Close();
		}

		double dR(double eta1, double phi1, double eta2, double phi2){
			//phi wraparound
			double dphi = (phi1-phi2);
			dphi = acos(cos(dphi));
			return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
		}

	private:
		double _gev;
		ReducedBaseSim* _base = nullptr;
		int _nEvts, _minNrhs;
		double _ptmin_reco, _ptmin_gen, _minrhE, _Emin_reco, _Emin_gen;
		int _nConstsmin;
		double _c = 29.9792458; // speed of light in cm/ns
		vector<int> lepIds = {11, 12, 13, 14, 15, 16};

};
#endif
