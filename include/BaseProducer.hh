#ifndef BaseProducer_HH
#define BaseProducer_HH

#include "ReducedBase.hh"
#include "Jet.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TCanvas.h"

class BaseProducer{
	public:
		BaseProducer(){ 
			_gev = 1;
			_isocut = false;
			_minpt = 0;
			_mineme = 0;
			_minnrhs = 0;
		};
		BaseProducer(TFile* file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			if(gSystem->AccessPathName(file->GetName())){ cout << "Error: file " << file->GetName() << " doesn't exist." << endl; return; }
			TTree* tree = (TTree*)file->Get("tree/llpgtree");
			_base = new ReducedBase(tree);
			_nEvts = _base->fChain->GetEntries();
			//default to 1 GeV = 1 entry -> gev = 1
			_gev = 1;
			_isocut = false;
			_minpt = 30;
			_mineme = 20;
			_minnrhs = 15;
		}
		virtual ~BaseProducer(){ 
			delete _base;
			
		};

		//TODO remove jet points - use only jets
		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt) = 0;
		virtual void GetRecHits(vector<Jet>& rhs, int evt){};
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt, int obj) = 0;
		virtual void GetRecHits(vector<Jet>& rhs, int evt, int obj){};
		virtual void GetPrimaryVertex(Point& vtx, int evt) = 0;


		void GetTrueJets(vector<Jet>& jets, int evt);
		void GetTruePhotons(vector<Jet>& phos, int evt);

		bool _isocut;
		void SetIsoCut(){ _isocut = true; }		

		ReducedBase* GetBase(){ return _base; }


		ReducedBase* _base = nullptr;
		int _nEvts;

		//energy weight transfer factor
		//w = g*E -> g = N/GeV
		void SetTransferFactor(double g){ _gev = g; }
		double _gev;


		void SetMinPt(double p){ _minpt = p; }
		double _minpt;
		void SetMinNrhs(double p){ _minnrhs = p; }
		double _minnrhs;
		void SetMinEmE(double p){ _mineme = p; }
		double _mineme;

		double deltaR2(double e1, double p1, double e2, double p2){
			double de = e1 - e2;
			double dp = fabs(p1 - p2);
			double pi = acos(-1);
			if(dp > pi)
				dp -= 2*pi;
			return de*de + dp*dp;
		}



};
#endif
