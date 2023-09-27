#ifndef BaseProducer_HH
#define BaseProducer_HH

#include "ReducedBase.hh"
#include "JetPoint.hh"
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
		BaseProducer(){ };
		BaseProducer(TFile* file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			if(gSystem->AccessPathName(file->GetName())){ cout << "Error: file " << file->GetName() << " doesn't exist." << endl; return; }
			_file = file;
			TTree* tree = (TTree*)file->Get("tree/llpgtree");
			_base = new ReducedBase(tree);
		//	_base->GetEntry(0);
		//	cout << "base prod init - " << _base->Photon_energy->size() << endl;
			_nEvts = _base->fChain->GetEntries();
			_inclpho = false;
		}
		virtual ~BaseProducer(){ };

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		virtual void GetRecHits(vector<vector<JetPoint>>& rhs) = 0;
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt) = 0;
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt, int obj) = 0;
		virtual void GetPrimaryVertex(Point& vtx, int evt) = 0;

		ReducedBase* GetBase(){ return _base; }

		void SetInclPho(bool ph){ _inclpho = ph; }	

		TFile* _file;
		ReducedBase* _base = nullptr;
		int _nEvts;

		bool _inclpho;




};
#endif
