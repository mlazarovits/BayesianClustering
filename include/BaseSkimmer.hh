#ifndef BaseSkimmer_HH
#define BaseSkimmer_HH

//#include "ReducedBase.hh"
#include "JetPoint.hh"
#include "BaseProducer.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>

using std::vector;
using std::string;
class BaseSkimmer{
	public:
		BaseSkimmer(){ 
			_gev = 1/10.;
			_data = false;
			_debug = false;
			_smear = true;
		};
		BaseSkimmer(TFile* file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			//getting the stuff below from producer in derived class
			//TTree* tree = (TTree*)file->Get("tree/llpgtree");
			//_base = new ReducedBase(tree);
			//_nEvts = _base->fChain->GetEntries();
			//_base->GetEntry(0);
			//cout << "base skim init - " << _base->Photon_energy->size() << endl;
		
			_gev = 1;
			_data = false;
			_debug = false;
			_smear = true;
			
			_hists1D.push_back(nSubClusters);
			_hists1D.push_back(time_center);
			_hists1D.push_back(eta_center);
			_hists1D.push_back(phi_center);
			_hists1D.push_back(objE);
			_hists1D.push_back(clusterE);


		}
		virtual ~BaseSkimmer(){ 
			delete _base;
			_hists1D.clear();
			_hists2D.clear();
		}

		virtual void Skim() = 0;

		ReducedBase* _base = nullptr;
		int _nEvts;
		BaseProducer* _prod;
		bool _data;
		bool _debug;
		int _evti, _evtj;
		string _cms_label, _oname;
		double _gev;

		
		void SetData(bool d){ _data = d; }
		void SetDebug(bool d){ _debug = d; }
		void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj; }
		void SetOutfile(string fname){ _oname = fname; }
		void SetTransferFactor(double gev){
			_gev = gev;
			_prod->SetTransferFactor(_gev);
		}

		void SetMinPt(double p){ _prod->SetMinPt(p); }
		void SetMinNrhs(double p){ _prod->SetMinNrhs(p); }
		void SetMinEmE(double p){ _prod->SetMinEmE(p); }

		void Profile2DHist(TH2D* inhist, TH1D* outhist, vector<TH1D*>& profs);

		vector<TH1D*> _hists1D;
		//0 - # of subclusters
		TH1D* nSubClusters = new TH1D("nSubClusters","nSubClusters",10,0,10.);
		//1 - mean time - center in t
		TH1D* time_center = new TH1D("time_center","time_center",50,-20,20);
		//2 - mean eta - center in eta
		TH1D* eta_center = new TH1D("eta_center","eta_center",50,-3.5,3.5);
		//3 - mean phi - center in phi
		TH1D* phi_center = new TH1D("phi_center","phi_center",50,-0.1,6.3);
		//4 - object energy
		TH1D* objE = new TH1D("objE","objE",50,0,1000);
		//5 - cluster energy
		TH1D* clusterE = new TH1D("clusterE","clusterE",10,0,1000);

	
		//two dimensional histograms
		vector<TH2D*> _hists2D;

		//reco object histograms
		//NOT in hist vectors
		TH2D* objE_clusterE = new TH2D("objE_clusterE","objE_clusterE;objE;clusterE",50,0,1050,50,0,1050);



		//struct for different types of plots (ie signal, ISR, fakes, etc.)
		struct plotCat{
			string legName;
			string plotName;
			
			vector<string> histcatnames = {"","lead","notlead"};
		
			vector<TH1D*> hists1D_nom;
			//for lead subcluster
			vector<TH1D*> hists1D_lead;
			//for !lead subcluster
			vector<TH1D*> hists1D_notlead;
			vector<vector<TH1D*>> hists1D;

			vector<TH2D*> hists2D_nom;
			//for lead subcluster
			vector<TH2D*> hists2D_lead;
			//for !lead subcluster
			vector<TH2D*> hists2D_notlead;
			vector<vector<TH2D*>> hists2D;
			vector<double> ids;
		
			plotCat(const vector<TH1D*>& in1dhists, const vector<TH2D*>& in2dhists, string plotname = "", string legname = ""){
				hists1D.push_back(hists1D_nom);
				hists1D.push_back(hists1D_lead);
				hists1D.push_back(hists1D_notlead);
				
				hists2D.push_back(hists2D_nom);
				hists2D.push_back(hists2D_lead);
				hists2D.push_back(hists2D_notlead);

				plotName = plotname;
				legName = legname;
				
				string name;
				//for each histogram (variable or correlation)
				for(int i = 0; i < (int)in1dhists.size(); i++){
					//create a clone for each type
					for(int j = 0; j < hists1D.size(); j++){
						//make sure they have the right add-on name (ie leading, !lead, etc)
						TH1D* hist = (TH1D*)in1dhists[i]->Clone();
						hists1D[j].push_back(hist);
						name = hist->GetName();
						if(!plotName.empty()) name += "_"+plotName;
						if(!histcatnames[j].empty()) name += "_"+histcatnames[j];
						hists1D[j][i]->SetName(name.c_str());
						if(!plotName.empty()) hists1D[j][i]->SetTitle("");
					}

				}
				//for each histogram
				for(int i = 0; i < (int)in2dhists.size(); i++){
					//create a clone for each type
					for(int j = 0; j < hists2D.size(); j++){
						TH2D* hist = (TH2D*)in2dhists[i]->Clone();
						hists2D[j].push_back(hist);
						name = hist->GetName();
						if(!plotName.empty()) name += "_"+plotName;
						name += "_"+histcatnames[j];
						hists2D[j][i]->SetName(name.c_str());
						if(!plotName.empty()) hists2D[j][i]->SetTitle("");
					}

				}
			}
			

		};
		vector<plotCat> plotCats;
		void MakeIDHists(string sample){
			//total
			plotCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			plotCats.push_back(tot);	
			
			if(sample.find("GMSB") != string::npos){
				//notSunm
				plotCat notSunm(_hists1D, _hists2D, "notSunm","notSunm");
				//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
				notSunm.ids = {29, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
				plotCats.push_back(notSunm);
				
				//signal
				plotCat sig(_hists1D, _hists2D, "chiGam","#Chi^{0} #rightarrow #gamma");
				sig.ids = {22};
				plotCats.push_back(sig);
			}
			else if(sample.find("JetHT") != string::npos){
				//data
				plotCat jetht(_hists1D, _hists2D, "JetHT", "JetHT");
				jetht.ids = {-999};
				plotCats.push_back(jetht);
			}
			else if(sample.find("GJets") != string::npos){
				//data
				plotCat gjets(_hists1D, _hists2D, "GJets", "GJets");
				gjets.ids = {-999};
				plotCats.push_back(gjets);
			}
			else return;

		}



		void SetCMSLabel(string lab){ _cms_label = lab; }

		void SetSmear(bool t){ _smear = t; }
		void SetTimeSmear(bool t){ _timesmear = t; }
		bool _smear, _timesmear;


};
#endif
