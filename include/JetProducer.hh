#ifndef JETPRODUCER_HH
#define JETPRODUCER_HH

#include "Jet.hh"
#include "TFile.h"
#include "TTree.h"

class JetProducer{
	public:
		JetProducer();
		virtual ~JetProducer();

		//get rechits from file to cluster
		JetProducer(TFile* file);
		//make ctor that simulates rechits

		//returns rec hits (jets)
		vector<Jet> GetRecHits(){ return m_rechits; }

	private:
		//individual rec hits (jets)
		vector<Jet> m_rechits;









};
#endif
