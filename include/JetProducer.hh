#ifndef JETPRODUCER_HH
#define JETPRODUCER_HH

#include "JetPoint.hh"
#include "Jet.hh"
#include "TFile.h"
#include "BaseProducer.hh"

class JetProducer : public BaseProducer{
	public:
		JetProducer();
		virtual ~JetProducer();

		//get rechits from file to cluster
		JetProducer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		//make ctor that simulates rechits
		

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		void GetRecHits(vector<JetPoint>& rhs, int evt);
		void GetRecHits(vector<Jet>& rhs, int evt);
		void GetSimRecHits(vector<Jet>& rhs, int evt);
		void GetRecHits(vector<JetPoint>& rhs, int evt, int jet){ };
		void GetGenJets(vector<Jet>& genjets, int evt);
		void GetPrimaryVertex(Point& vtx, int evt);






};
#endif
