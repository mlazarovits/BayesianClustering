#ifndef JETPRODUCER_HH
#define JETPRODUCER_HH

#include "JetPoint.hh"
#include "ReducedBase.hh"
#include "TFile.h"

class JetProducer{
	public:
		JetProducer();
		virtual ~JetProducer();

		//get rechits from file to cluster
		JetProducer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		//make ctor that simulates rechits
		

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		void GetRecHits(vector<vector<JetPoint>>& rhs);
		void GetRecHits(vector<JetPoint>& rhs, int evt)
;
		void GetPrimaryVertex(Point& vtx, int evt);
	private:
		//individual rec hits (jets)
		//vector<vector<JetPoint>> m_rechits;
		ReducedBase* m_base = nullptr;
		int m_nEvts;
		TFile* m_file;







};
#endif
