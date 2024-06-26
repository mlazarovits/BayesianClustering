#ifndef JETSIMPRODUCER_HH
#define JETSIMPRODUCER_HH
#include "ReducedBaseSim.hh"
#include "Jet.hh"
#include "JetPoint.hh"
#include "TFile.h"

class JetSimProducer{
	public:
		JetSimProducer();
		virtual ~JetSimProducer();

		JetSimProducer(TFile* file);

		void GetRecHits(vector<JetPoint>& rhs, int evt);
		void GetRecHits(vector<Jet>& rhs, int evt);
		void GetGenJets(vector<Jet>& genjets, int evt);
		void GetPrimaryVertex(Point& vtx, int evt);
		ReducedBaseSim* GetBase(){ return _base; }
		void SetTransferFactor(double g){ _gev = g; }

	private:
		double _gev;
		ReducedBaseSim* _base = nullptr;
		int _nEvts;
};
#endif
