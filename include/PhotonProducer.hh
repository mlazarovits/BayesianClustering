#ifndef PHOTONPRODUCER_HH
#define PHOTONPRODUCER_HH

#include "JetPoint.hh"
#include "TFile.h"
#include "TTree.h"
#include "ReducedBase.hh"

class PhotonProducer{
	public:
		PhotonProducer();
		virtual ~PhotonProducer();

		//get rechits from file to cluster
		PhotonProducer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		//make ctor that simulates rechits
		

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		void GetRecHits(vector<vector<JetPoint>>& rhs);
		void GetRecHits(vector<JetPoint>& rhs, int evt);
		void GetRecHits(vector<JetPoint>& rhs, int evt, int pho);
		void GetPrimaryVertex(Point& vtx, int evt);

		void Skim();
	private:
		//individual rec hits (jets)
		//vector<vector<JetPoint>> m_rechits;
		ReducedBase* m_base = nullptr;
		int m_nEvts;
		TFile* m_file;







};
#endif
