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
		void GetGenJets(vector<Jet>& genjets, int evt);
		void GetRecoJets(vector<Jet>& recojets, int evt);
		void GetPrimaryVertex(BayesPoint& vtx, int evt);
		ReducedBaseSim* GetBase(){ return _base; }
		void SetTransferFactor(double g){ _gev = g; }

		void SetRecoPtMin(double pt){_recoptmin = pt; }
		void SetMinRhE(double r){ _minrhE = r; }

		void SortJets(vector<Jet>& jets);
		
		void PrintPreselection(){
			cout << "Default energy transfer factor: " << _gev << endl;
			cout << "Minimum reco pt: " << _recoptmin << endl;
			cout << "Minimum rh (barrel only) energy: " << _minrhE << endl;
		}


		void EtaPhiMap(string fname, vector<Jet>& rhs){
			TFile* f = new TFile(fname.c_str(),"RECREATE");
			f->cd();

			//TODO: change bin size to be 1 cell in det eta, phi
			//TH2D* etaPhiMap = new TH2D("etaPhiMap","etaPhiMap;eta;phi",40,-1.4775,1.7,40,-3.575,3.25);
			TH2D* etaPhiMap = new TH2D("etaPhiMap","etaPhiMap;eta;phi",40,-1.4775,1.7,40,-3.575,3.25);
			for(int r = 0; r < rhs.size(); r++){
				cout << "r " << r << " eta " << rhs[r].eta() << " phi_std " << rhs[r].phi_std() << " E " << rhs[r].E() << endl;
				etaPhiMap->Fill(rhs[r].eta(), rhs[r].phi_std(), rhs[r].E());
			}
			etaPhiMap->Write();
			f->Close();
		}

	private:
		double _gev;
		ReducedBaseSim* _base = nullptr;
		int _nEvts;
		double _recoptmin, _minrhE;
};
#endif
