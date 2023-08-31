#ifndef PHOTONPRODUCER_HH
#define PHOTONPRODUCER_HH

#include "JetPoint.hh"
#include "ReducedBase.hh"
#include <TFile.h>

class PhotonProducer{
	public:
		PhotonProducer();
		virtual ~PhotonProducer();

		//get rechits from file to cluster
		PhotonProducer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		//true = keep
		//false = drop
		bool cleanRH(JetPoint rh){
			if(fabs(rh.t()) > 25 || rh.e() < 1.) return false;
			else return true;
		}	
	

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		void GetRecHits(vector<vector<JetPoint>>& rhs);
		void GetRecHits(vector<JetPoint>& rhs, int evt);
		void GetRecHits(vector<JetPoint>& rhs, int evt, int pho);
		void GetPrimaryVertex(Point& vtx, int evt);

		void CleaningSkim();
		void Skim();
	private:
		//individual rec hits (jets)
		//vector<vector<JetPoint>> m_rechits;
		ReducedBase* m_base = nullptr;
		int m_nEvts;
		TFile* m_file;







};
#endif
