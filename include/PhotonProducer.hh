#ifndef PHOTONPRODUCER_HH
#define PHOTONPRODUCER_HH

#include "JetPoint.hh"
#include <TFile.h>
#include "BaseProducer.hh"
class PhotonProducer : public BaseProducer{
	public:
		PhotonProducer();
		virtual ~PhotonProducer();

		//get rechits from file to cluster
		PhotonProducer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		//true = keep
		//false = drop
		bool cleanRH(JetPoint rh){
			if(fabs(rh.t()) > 25 && rh.e() < 1.) return false;
			else return true;
		}	
	

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		void GetRecHits(vector<vector<JetPoint>>& rhs);
		void GetRecHits(vector<JetPoint>& rhs, int evt);
		void GetRecHits(vector<JetPoint>& rhs, int evt, int pho);
		void GetPrimaryVertex(Point& vtx, int evt);

		void CleaningSkim();
		void Skim();




		//stacked photon LLID hists
		//list of photon ids
		vector<int> ids = {22, 32, 25, 35, 23, 24, 21, 31, 20, 30, 26, 27};
		vector<string> id_names = {"#Chi \rightarrow #gamma", "ISR", "Not SUSY"};
		vector<TH1D*> pho_llp_polang;
		vector<TH1D*> pho_llp_azimang;
		//add photon id to hist name/title


		void MakeIDHists(){
			for(int i = 0; i < (int)id_names.size(); i++) pho_llp_polang.push_back(new TH1D(("polar_ang"+std::to_string(i)).c_str(),("polar_ang"+std::to_string(ids[i])).c_str(),50,-3.5,3.5));	

			for(int i = 0; i < (int)id_names.size(); i++) pho_llp_azimang.push_back(new TH1D(("azimuth_ang"+std::to_string(i)).c_str(),("azimuth_ang"+std::to_string(ids[i])).c_str(),50,-3.5,3.5));	

		}



};
#endif
