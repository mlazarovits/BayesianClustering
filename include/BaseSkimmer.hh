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
#include "SampleWeight.hh"
using weights = SampleWeight::weights;

using std::vector;
using std::string;
class BaseSkimmer{
	public:
		BaseSkimmer(){ 
			_gev = 1;
			_data = false;
			_debug = false;
			_smear = true;
			_timesmear = false;
			_skip = 1;
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
			_debug = false;
			_smear = true;
			_timesmear = false;
			_skip = 1;
			
			string filename = file->GetName();	
			if(filename.find("SIM") != string::npos)
				_data = false;
			else
				_data = true;
			
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
		double _c = 29.9792458; // speed of light in cm/ns
		
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
		void SetMinRhE(double r){ _prod->SetMinRhE(r); }

		void Profile2DHist(TH2D* inhist, TH1D* outhist, vector<TH1D*>& profs);

		vector<TH1D*> _hists1D;
		//0 - # of subclusters
		TH1D* nSubClusters = new TH1D("nSubClusters","nSubClusters",10,0,10.);
		//1 - mean time - center in t
		TH1D* time_center = new TH1D("timeCenter","timeCenter",50,-20,20);
		//2 - mean eta - center in eta
		TH1D* eta_center = new TH1D("etaCenter","etaCenter",50,-1.6,1.6);
		//3 - mean phi - center in phi
		TH1D* phi_center = new TH1D("phiCenter","phiCenter",50,-0.2,6.4);
		//4 - object energy
		TH1D* objE = new TH1D("objE","objE",50,0,1000);
		//5 - cluster energy
		TH1D* clusterE = new TH1D("clusterE","clusterE",10,0,1000);

	
		//two dimensional histograms
		vector<TH2D*> _hists2D;

		//reco object histograms
		//NOT in hist vectors
		TH2D* objE_clusterE = new TH2D("objE_clusterE","objE_clusterE;objE;clusterE",50,0,1050,50,0,1050);



		//for sample weights
		SampleWeight _swts;
		//weight to apply to all histograms
		double _weight;
		//skip for event loop
		int _skip;
		void SetSkip(int i){ _skip = i; _weight *= _skip; }
		

		//struct for different types of plots (ie signal, ISR, fakes, etc.)
		struct procCat{
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
		
			procCat(const vector<TH1D*>& in1dhists, const vector<TH2D*>& in2dhists, string plotname = "", string legname = "", bool leadsep = true){
				hists1D.push_back(hists1D_nom);
				hists2D.push_back(hists2D_nom);
				if(leadsep){
					hists1D.push_back(hists1D_lead);
					hists1D.push_back(hists1D_notlead);
					hists2D.push_back(hists2D_lead);
					hists2D.push_back(hists2D_notlead);
				}
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
						if(!histcatnames[j].empty()) name += "_"+histcatnames[j];
						hists2D[j][i]->SetName(name.c_str());
						if(!plotName.empty()) hists2D[j][i]->SetTitle("");
					}

				}
			}
				
			//reset proc cat hists to different hists
			void SetHists(const vector<TH1D*>& in1dhists, const vector<TH2D*>& in2dhists, bool leadSep = true){
				string name;
				int nhists;
			//cout << "SetHists for "  << plotName << endl;	
//cout << "pre clear " << hists1D.size() << " " << hists1D[0].size() << endl;
				//how many hist categories to loop over
				if(leadSep) nhists = hists1D.size();
				else nhists = 1; //just nominal
				hists1D.clear(); hists2D.clear();
//cout << "post clear " << hists1D.size() << " " << hists1D[0].size() << endl;
				for(int i = 0; i < nhists; i++){ hists1D.push_back({}); hists2D.push_back({}); }

				//for each histogram (variable or correlation)
				for(int i = 0; i < (int)in1dhists.size(); i++){
					//create a clone for each type
					for(int j = 0; j < nhists; j++){
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
//				cout << "SetHists - n hist: " << hists1D[0].size() << " " << in1dhists.size() << endl;
				//for each histogram
				for(int i = 0; i < (int)in2dhists.size(); i++){
					//create a clone for each type
					for(int j = 0; j < nhists; j++){
						TH2D* hist = (TH2D*)in2dhists[i]->Clone();
						hists2D[j].push_back(hist);
						name = hist->GetName();
						if(!plotName.empty()) name += "_"+plotName;
						if(!histcatnames[j].empty()) name += "_"+histcatnames[j];
						hists2D[j][i]->SetName(name.c_str());
						if(!plotName.empty()) hists2D[j][i]->SetTitle("");
					}

				}
//				cout << "2d SetHists - n hist: " << hists2D[0].size() << " " << in2dhists.size() << endl;
			}
			void AddHist(TH1D* inhist){
				//create a clone for each type
				int n1dhist;
				string name;
				for(int j = 0; j < hists1D.size(); j++){
					//make sure they have the right add-on name (ie leading, !lead, etc)
					TH1D* hist = (TH1D*)inhist->Clone();
					hists1D[j].push_back(hist);
					n1dhist = hists1D.size()-1;
					name = hist->GetName();
					if(!plotName.empty()) name += "_"+plotName;
					name += "_"+histcatnames[j];
					hists1D[j][n1dhist]->SetName(name.c_str());
					if(!plotName.empty()) hists1D[j][n1dhist]->SetTitle("");
				}
			}	
			void AddHist(TH2D* inhist){
				int n2dhist;
				string name;
				for(int j = 0; j < hists1D.size(); j++){
					TH2D* hist = (TH2D*)inhist->Clone();
					hists2D[j].push_back(hist);
					n2dhist = hists2D.size()-1;
					name = hist->GetName();
					if(!plotName.empty()) name += "_"+plotName;
					name += "_"+histcatnames[j];
					hists2D[j][n2dhist]->SetName(name.c_str());
					if(!plotName.empty()) hists2D[j][n2dhist]->SetTitle("");
				}
			}	

			

		};
		vector<procCat> _procCats;
		void MakeProcCats(string sample, bool leadsep = true){
			//total
			procCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			_procCats.push_back(tot);	
			//cout << "sample " << sample << endl;	
			if(sample.find("GMSB") != string::npos){
				//notSunm
				procCat notSunm(_hists1D, _hists2D, "notSunm","notSunm", leadsep);
				//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
				notSunm.ids = {97, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
				_procCats.push_back(notSunm);
				
				//signal
				//do string matching to find specific grid point
				string lambda, ctau;
				string lmatch = "L-";
				string sample_l = sample.substr(sample.find(lmatch));
				lambda = sample_l.substr(0,sample_l.find("_"));
				
				string ctmatch = "Ctau";
				string sample_ctau = sample.substr(sample.find(ctmatch));
				ctau = sample_ctau.substr(0,sample_ctau.find("_"));
				
				//cout << "lambda " << lambda << " ctau " << ctau << endl;

				string lfancy = lambda.substr(lambda.find("-")+1);
				//lfancy.insert(lfancy.find("TeV")-2," ");
				lfancy += " TeV";
				string ctfancy = ctau.substr(ctau.find("-")+1);
				//ctfancy.insert(ctfancy.find("cm")," ");
				ctfancy += " cm";

	
				string plotName = "chiGam_"+lambda+"_"+ctau;
				while(plotName.find("-") != string::npos)
					plotName.replace(plotName.find("-"),1,"");

				string legName = "#Chi^{0} #rightarrow #gamma, L = "+lfancy+" c#tau = "+ctfancy;
				procCat sig(_hists1D, _hists2D, plotName, legName, leadsep);
				sig.ids = {22};
				_procCats.push_back(sig);
			}
			else if(sample.find("JetHT") != string::npos){
				//data
				procCat jetht(_hists1D, _hists2D, "JetHTPD", "JetHTPD", leadsep);
				jetht.ids = {-999};
				_procCats.push_back(jetht);
			}
			else if(sample.find("MET") != string::npos){
				//data
				procCat met(_hists1D, _hists2D, "METPD", "METPD", leadsep);
				met.ids = {-999};
				_procCats.push_back(met);
			}
			else if(sample.find("EGamma") != string::npos || sample.find("DoubleEG") != string::npos){
				//data
				procCat egam(_hists1D, _hists2D, "DoubleEGPD", "DoubleEGPD", leadsep);
				egam.ids = {-999};
				_procCats.push_back(egam);
			}
			else if(sample.find("GJets") != string::npos){
				//data
				procCat gjets(_hists1D, _hists2D, "GJets", "GJets", leadsep);
				gjets.ids = {-999};
				_procCats.push_back(gjets);
			}
			else if(sample.find("ttbar") != string::npos){
				procCat ttbar(_hists1D, _hists2D, "ttbar", "t#bar{t}",leadsep);
				ttbar.ids = {-999};
				_procCats.push_back(ttbar);
			}
			else if(sample.find("QCD") != string::npos){
				procCat qcd(_hists1D, _hists2D, "QCD", "QCD multi-jets",leadsep);
				qcd.ids = {-999};
				_procCats.push_back(qcd);
			}
			else return;

		}



		void SetCMSLabel(string lab){ _cms_label = lab; }

		void SetSmear(bool t){ _smear = t; }
		void SetTimeSmear(bool t){ _timesmear = t; }
		void SetTimeCalibrationMap(string f){
                        if(gSystem->AccessPathName(f.c_str())){
                                cout << "Error: file " << f << " does not exist. Time calibration file could not be set." << endl;
                                return;
                        }
                        TFile* file = TFile::Open(f.c_str());
                        _prod->SetTimeCalibrationMap(file);
                }
		bool _smear, _timesmear;

		void SetSpikeRejection(bool s){ _prod->RejectSpikes(s); }

};
#endif
