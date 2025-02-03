#ifndef BaseSkimmer_HH
#define BaseSkimmer_HH

//#include "ReducedBase.hh"
#include "JetPoint.hh"
#include "BasePDFMixture.hh"
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
//frugally deep for loading NN model
#include "fdeep/fdeep.hpp"

using weights = SampleWeight::weights;
static bool Esort(JetPoint j1, JetPoint j2){ return (j1.E() > j2.E()); }
static bool ptsort(Jet j1, Jet j2){ return (j1.pt() > j2.pt()); }
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
			_weight = 1;
			_ngrid = 7;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
                	_BHFilter = beamHaloFilter(0); //default not applied

			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);

			
			cout << "Default NN model: small8CNN_EMultr.json" << endl;


		};
		BaseSkimmer(TFile* file){
			_gev = 1;
			_debug = false;
			_smear = true;
			_timesmear = false;
			_skip = 1;
			_weight = 1;
			_ngrid = 7;
                	_BHFilter = beamHaloFilter(0); //default not applied
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);
			
			
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

			cout << "Default NN model: small8CNN_EMultr.json" << endl;


		}
		
		BaseSkimmer(string flist){
			_gev = 1;
			_debug = false;
			_smear = true;
			_timesmear = false;
			_skip = 1;
			_weight = 1;
			_ngrid = 7;
                	_BHFilter = beamHaloFilter(0); //default not applied
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);
			
			
			if(flist.find("SIM") != string::npos)
				_data = false;
			else
				_data = true;
			
			_hists1D.push_back(nSubClusters);
			_hists1D.push_back(time_center);
			_hists1D.push_back(eta_center);
			_hists1D.push_back(phi_center);
			_hists1D.push_back(objE);
			_hists1D.push_back(clusterE);

			cout << "Default NN model: small8CNN_EMultr.json" << endl;


		}
		virtual ~BaseSkimmer(){ 
			delete _base;
			_hists1D.clear();
			_hists2D.clear();
		}

		virtual void Skim() = 0;

		//make tchain from filelist
		//debug after jet res for 1st exo talk are done
		TChain* MakeTChain(string flist){
			if(gSystem->AccessPathName(flist.c_str())){ 
				cout << "Error: file " << flist << " doesn't exist." << endl; 
				return nullptr; 
			}	
			std::ifstream infile(flist);
			TChain* ch = new TChain("tree/llpgtree");
			ch->SetTitle(flist.c_str());
			string file;
			cout << "TChaining files in " << flist << "..." << endl;
			while(std::getline(infile,file)){
				if( file[0] == '#' ) continue;
				if(gSystem->AccessPathName(file.c_str())){
					cout << "Skipping file " << file << " - not found." << endl;
					continue;
				}
				//std::cout << "--  adding file: " << file << std::endl;
				ch->Add(file.c_str());
			}
			cout << "Done TChaining" << endl;
			return ch;
		}

		ReducedBase* _base = nullptr;
		int _nEvts;
		BaseProducer* _prod;
		bool _data;
		bool _debug;
		int _evti, _evtj;
		string _cms_label, _oname;
		double _gev, _thresh, _alpha, _emAlpha;
		map<string, Matrix> _prior_params;
		void SetThresh(double t){ _thresh = t; }
		void SetBHCAlpha(double a){ _alpha = a; }
		void SetEMAlpha(double a){ _emAlpha = a; }
		void SetPriorParameters(map<string, Matrix> params){_prior_params = params;} 		

		double _c = 29.9792458; // speed of light in cm/ns
		struct DetIDStruct {
		        DetIDStruct() {}
		        DetIDStruct(const int ni1, const int ni2, const Int_t nTT, const Int_t & necal, const double eta, const double phi) : i1(ni1), i2(ni2), TT(nTT), ecal(necal), deteta(eta), detphi(phi){}
		        //Int_t i1; // EB: iphi, EE: ix
		        int i1;
		//      Int_t i2; // EB: ieta, EE: iy
		        int i2;
		        double deteta, detphi;
		        Int_t TT; // trigger tower
		        Int_t ecal; // EB, EM, EP
		};//<<>>struct DetIDStruct
	
		std::map<UInt_t,DetIDStruct> _detIDmap;		
		std::map<pair<int, int>, UInt_t> _ietaiphiID;	
		
		//this function and the corresponding DetIDStruct (above) are courtesy of Jack King 
		//https://github.com/jking79/LLPgammaAnalyzer/blob/master/macros/KUCMS_Skimmer/KUCMSHelperFunctions.hh  
		void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap, std::map<pair<int,int>, UInt_t> &iEtaiPhiToDetID ){
		        const std::string detIDConfigEB("info/fullinfo_v2_detids_EB.txt");
		        std::ifstream infile( detIDConfigEB, std::ios::in);
		
		        UInt_t cmsswId, dbID;
		        pair<int, int> ietaiphi;
		        int hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
		        double deteta, detphi;
		        std::string pos;
		
		        while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM >> detphi >> deteta){
		            //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << 0 << std::endl;
		            DetIDMap[cmsswId] = {iphi,ieta,TT25,0,deteta,detphi};
		            ietaiphi = make_pair(ieta, iphi);
		            iEtaiPhiToDetID[ietaiphi] = cmsswId;
		            //auto idinfo = DetIDMap[cmsswId];
		            //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
		        }//while (infile >>
		
		}//<<>>void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )


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
		void SetMaxRhE(double r){ _prod->SetMaxRhE(r); }
		void SetSpatialCorr(bool r){ _prod->SetSpatialCorr(r); }
		void SetCNNGrid(double n){ _ngrid = n; }

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
					n1dhist = hists1D[j].size()-1;
					name = hist->GetName();
					if(!plotName.empty()) name += "_"+plotName;
					if(!histcatnames[j].empty()) name += "_"+histcatnames[j];
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
					n2dhist = hists2D[j].size()-1;
					name = hist->GetName();
					if(!plotName.empty()) name += "_"+plotName;
					if(!histcatnames[j].empty()) name += "_"+histcatnames[j];
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
				
				//notSunm
				procCat notSunm(_hists1D, _hists2D, "notSunm","notSunm", leadsep);
				//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
				notSunm.ids = {97, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
				_procCats.push_back(notSunm);
			}
			else if(sample.find("SMS-GlGl") != string::npos){
				//data
				procCat glgl(_hists1D, _hists2D, "GluGluN2", "GluGluN2", leadsep);
				glgl.ids = {-999};
				_procCats.push_back(glgl);
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



		//create function to write photon subcluster variables to CSV file for MVA training
		//include column for process?
		string _csvname;
		ofstream _csvfile;
	
		vector<string> _inputs;
		int _ngrid;
		void SetObs(){
			//sample
			_inputs.push_back("sample");
			//event
			_inputs.push_back("event");
			//supercl
			_inputs.push_back("object");
			//subcl
			_inputs.push_back("subcl");
			//etacenter
			_inputs.push_back("eta_center");
			//phicenter
			_inputs.push_back("phi_center");
			//timecenter
			_inputs.push_back("time_center");
			//etasig
			_inputs.push_back("eta_sig");
			//phisig
			_inputs.push_back("phi_sig");
			//etaphicov
			_inputs.push_back("etaphi_cov");
			//timeetacov
			_inputs.push_back("timeeta_cov");
			//energy
			_inputs.push_back("energy");
			//sw+
			_inputs.push_back("sw+");
			//subcl major length
			_inputs.push_back("major_length");
			//subcl minor length
			_inputs.push_back("minor_length");
			//subcl max pt / total E
			_inputs.push_back("maxOvtotE");
			//R9
			_inputs.push_back("R9");
			//Sietaieta
			_inputs.push_back("Sietaieta");
			//Siphiiphi
			_inputs.push_back("Siphiiphi");
			//Smajor
			_inputs.push_back("Smajor");
			//Sminor
			_inputs.push_back("Sminor");
			//ecalPFClusterIsoOvPt
			_inputs.push_back("ecalPFClusterIsoOvPt");
			//hcalPFClusterIsoOvPt
			_inputs.push_back("hcalPFClusterIsoOvPt");
			//trkSumPtHollowConeDR03OvPt
			_inputs.push_back("trkSumPtHollowConeDR03OvPt");
			//trkSumPtSolidConeDR04
			_inputs.push_back("trkSumPtSolidConeDR04");
                	//ecalRHSumEtConeDR04
                	_inputs.push_back("ecalRHSumEtConeDR04");
                	//hadTowOverEM
                	_inputs.push_back("hadTowOverEM");
			//lead photon? 1 = yes, 0 = no
			_inputs.push_back("lead");
			//CNN inputs
			for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++)
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++)
					_inputs.push_back("CNNgrid_E_cell"+to_string(i)+"_"+to_string(j));
			for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++)
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++)
					_inputs.push_back("CNNgrid_t_cell"+to_string(i)+"_"+to_string(j));
			for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++)
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++)
					_inputs.push_back("CNNgrid_r_cell"+to_string(i)+"_"+to_string(j));
			//label
			_inputs.push_back("label");
			for(auto s : _inputs){
				if(s != "label") _csvfile << s << ","; 
				else _csvfile << s << endl;
			}

		}
		void WriteObs(map<string,double> inputs, string object){
			string samp = "";
			if(_oname.find("condor") == string::npos){
				if(_oname.find("MET") != string::npos)
					samp = "METPD";
				else if(_oname.find("JetHT") != string::npos)
					samp = "JetHTPD";
				else if(_oname.find("GMSB") != string::npos){
					samp = _oname.substr(_oname.find("GMSB"),_oname.find("cm") - _oname.find("GMSB") + 2);
				}
				else if(_oname.find("GJets") != string::npos)
					samp = _oname.substr(_oname.find("GJets"),_oname.find("_AODSIM") - _oname.find("GJets"));
				else if(_oname.find("QCD") != string::npos)
					samp = _oname.substr(_oname.find("QCD"),_oname.find("_AODSIM") - _oname.find("QCD"));
				else if(_oname.find("DEG") != string::npos)
					samp = _oname.substr(_oname.find("DEG"),_oname.find("_AODSIM") - _oname.find("DEG"));
				else samp = "notFound";
			}
			else{
				if(_oname.find("MET") != string::npos)
					samp = "METPD";
				else if(_oname.find("JetHT") != string::npos)
					samp = "JetHTPD";
				else if(_oname.find("GMSB") != string::npos){
					samp = _oname.substr(_oname.find("GMSB"),_oname.find("_"+object) - _oname.find("GMSB"));
				}
				else if(_oname.find("GJets") != string::npos)
					samp = _oname.substr(_oname.find("GJets"),_oname.find("_"+object) - _oname.find("GJets"));
				else if(_oname.find("QCD") != string::npos)
					samp = _oname.substr(_oname.find("QCD"),_oname.find("_"+object) - _oname.find("QCD"));
				else if(_oname.find("EGamma") != string::npos)
					samp = _oname.substr(_oname.find("EGamma"),_oname.find("_"+object) - _oname.find("EGamma"));
				else samp = "notFound";
			}
			_csvfile << samp;// << evt << "," << obj << "," << ncl;
			for(int d = 1; d < _inputs.size(); d++){
				//cout << _inputs[d] << ": " << inputs[_inputs[d]] << endl; 
				_csvfile << "," << inputs[_inputs[d]];
			}
			_csvfile << endl;
		}
		void GetNeighborE(vector<JetPoint>& rhs, int r, vector<pair<int, int>>& icoords, vector<double>& Es, bool skipCenter = false, int ngrid = 5){
			int ieta, iphi;
			JetPoint rh;
			int ngrid_thresh = ngrid/2;
			if(r == -1){
				//choose highest E rh as center
				sort(rhs.begin(), rhs.end(), Esort);
				rh = rhs[0];
			}
			else{
				rh = rhs[r];
			}
			int rh_ieta = _detIDmap[rh.rhId()].i2;
			int rh_iphi = _detIDmap[rh.rhId()].i1;
			icoords.clear(); Es.clear();
			int deta, dphi;
			for(int j = 0; j < (int)rhs.size(); j++){
				ieta = _detIDmap[rhs[j].rhId()].i2;
				iphi = _detIDmap[rhs[j].rhId()].i1;
				deta = ieta - rh_ieta;
				dphi = iphi - rh_iphi;
				//do wraparound
				if(dphi > 180)
					dphi = 360 - dphi;
				else if(dphi < -180)
					dphi = -(360 + dphi);
				if(fabs(ieta - rh_ieta) <= ngrid_thresh && fabs(iphi - rh_iphi) <= ngrid_thresh){
					if(skipCenter && (ieta - rh_ieta == 0) && (iphi - rh_iphi == 0)) continue;
					else{
						Es.push_back(rhs[j].E());
					}
					if(rh_ieta < 0)
						icoords.push_back(make_pair(-(ieta - rh_ieta), iphi - rh_iphi));
					else
						icoords.push_back(make_pair(ieta - rh_ieta, iphi - rh_iphi));
				}
				//else
				//	cout << "r " << j << " deta " << deta << " dphi " << dphi << " E " << rhs[j].E() << " ngrid thresh " << ngrid_thresh << endl;
			}

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
		enum beamHaloFilter{
                        notApplied = 0,
                        applied = 1,
                        invApplied = 2
                };
                beamHaloFilter _BHFilter;
                void SetBeamHaloFilter(int bh){
                        _BHFilter = beamHaloFilter(bh);
                        if(_BHFilter == notApplied) cout << "Not applying beam halo filter." << endl;
                        else if(_BHFilter == applied) cout << "Applying beam halo filter." << endl;
                        else cout << "Applying inverse beam halo filter." << endl;
                }

		void GetCenterXtal(vector<JetPoint>& rhs, JetPoint& center){
			//get center of pts in ieta, iphi -> max E point
			double maxE = 0;
			for(int r = 0; r < rhs.size(); r++){
				if(rhs[r].E() > maxE){
					maxE = rhs[r].E();
					center = rhs[r];
				}
			}
		}


		//write nxn grid of E, t, r_nk here
		//void MakeCNNInputGrid(BasePDFMixture* model, int k, vector<JetPoint>& rhs, JetPoint& center, vector<map<string,double>>& mapobs){
		void MakeCNNInputGrid(BasePDFMixture* model, int k, vector<JetPoint>& rhs, JetPoint& center, map<string,double>& mapobs){
			map<pair<int,int>, vector<double>> grid;
			//make sure ngrid is odd to include center crystal
			if(_ngrid % 2 == 0)
				_ngrid++;

			int ngrid_boundary = (_ngrid-1)/2;
			//set default channel values to 0
			//{E, t, r, Er}
			for(int i = -ngrid_boundary; i < ngrid_boundary+1; i++)
				for(int j = -ngrid_boundary; j < ngrid_boundary+1; j++)
					grid[make_pair(i,j)] = {0., 0., 0.};


			//get ngrid x ngrid around center point 
			int ieta, iphi;
			int rh_ieta = _detIDmap[center.rhId()].i2;
			int rh_iphi = _detIDmap[center.rhId()].i1;
			int deta, dphi;
			Matrix post = model->GetPosterior();
			for(int j = 0; j < (int)rhs.size(); j++){
				ieta = _detIDmap[rhs[j].rhId()].i2;
				iphi = _detIDmap[rhs[j].rhId()].i1;
				//do eta flip
				if(rh_ieta < 0)
					deta = -(ieta - rh_ieta);
				else
					deta = ieta - rh_ieta;
				dphi = iphi - rh_iphi; 
				//needs wraparound
				if(dphi > 180)
					dphi = 360 - dphi;
				else if(dphi < -180)
					dphi = -(360 + dphi);
				if(fabs(deta) <= ngrid_boundary && fabs(dphi) <= ngrid_boundary){
					//posterior is weighted s.t. sum_k post_nk = w_n = E*_gev, need to just have unweighted probs since E is already here
					//r_nk = post_nk/w_n s.t. sum_k (post_nk/w_n) = w_n/w_n = 1
					grid[make_pair(deta, dphi)] = {rhs[j].E(), rhs[j].t() - center.t(), post.at(j,k)/model->GetData()->at(j).w()};	
					//mapobs[k]["CNNgrid_E_cell"+to_string(deta)+"_"+to_string(dphi)] = rhs[j].E();
					//mapobs[k]["CNNgrid_t_cell"+to_string(deta)+"_"+to_string(dphi)] = rhs[j].t() - center.t();
					//mapobs[k]["CNNgrid_r_cell"+to_string(deta)+"_"+to_string(dphi)] = post.at(j,k)/model->GetData()->at(j).w();
					//mapobs[k]["CNNgrid_Er_cell"+to_string(deta)+"_"+to_string(dphi)] = rhs[j].E()*(post.at(j,k)/model->GetData()->at(j).w());
					//if(deta == 0 && dphi == -1) cout << "cell (" << deta << ", " << dphi << ") weights E = " << rhs[j].E() << ", t = " << rhs[j].t() - center.t() << ", r = " << mapobs[k]["CNNgrid_r_cell"+to_string(deta)+"_"+to_string(dphi)] << endl;

				}
			}
			pair<int, int> icoords_grid;
			for(int i = -(_ngrid-1)/2.; i < (_ngrid-1)/2+1; i++){
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
					icoords_grid = make_pair(i,j);
					mapobs["CNNgrid_E_cell"+to_string(i)+"_"+to_string(j)] = grid[icoords_grid][0];
					mapobs["CNNgrid_t_cell"+to_string(i)+"_"+to_string(j)] = grid[icoords_grid][1];
					mapobs["CNNgrid_r_cell"+to_string(i)+"_"+to_string(j)] = grid[icoords_grid][2];
					mapobs["CNNgrid_Er_cell"+to_string(i)+"_"+to_string(j)] = grid[icoords_grid][0]*grid[icoords_grid][2];
					//if(i == 0 && j == -1) cout << "cell (" << i << ", " << j << ") weights E = " << grid[icoords_grid][0] << ", t = " << grid[icoords_grid][1] << ", r = " << grid[icoords_grid][2] << endl; 
				}
			}
		}

		//pass name of json from frugally-deep-master/keras_export/convert_model.py
		//for DNN input
		fdeep::model _nnmodel = fdeep::load_model("json/small8CNN_EMultr.json");
		void SetNNModel(string fname){
			_nnmodel = fdeep::load_model(fname);
		}
		//sets which features to use in NN prediction - make sure this matches with the given model (json)
		vector<string> _nnfeatures;
		void SetNNFeatures(vector<string>& features){
			_nnfeatures = features;
		}
		//returns predicted class, by reference the discriminator value - DNN
		int DNNPredict(map<string,double>& obs, double& ovalue){
			vector<float> input_sample;
			//transform obs to input_sample
			for(int f = 0; f < _nnfeatures.size(); f++)
				input_sample.push_back(obs[_nnfeatures[f]]);		

			int size = input_sample.size();
			fdeep::tensor input_tensor = fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(size)), input_sample);
			
			//predict_class returns predicted class number and value of max output neuron
			pair<size_t, double> result = _nnmodel.predict_class_with_confidence({input_tensor});
			ovalue = result.second;
			return (int)result.first;
		}
	
	
		//for CNN input - (ngrid, ngrid, nchannels) grid - CNN
		//grid maps int pair to vector of channels
		int CNNPredict(map<string, double>& obs, double& ovalue){
			int nchan = _nnfeatures.size();
			fdeep::tensor_shape tensor_shape(_ngrid, _ngrid, nchan);
			fdeep::tensor input_tensor(tensor_shape, 0.0f);
			//transform grid to input_sample
			for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++){
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
					for(int c = 0; c < _nnfeatures.size(); c++){
						double val;
						if (_nnfeatures[c].find("Mult") != string::npos){
							string ch1(_nnfeatures[c]);
							ch1 = ch1[0];
							string ch2(_nnfeatures[c]);
							ch2 = ch2[ch2.size()-1];
							val = obs["CNNgrid_"+ch1+"_cell"+to_string(i)+"_"+to_string(j)]*obs["CNNgrid_"+ch2+"_cell"+to_string(i)+"_"+to_string(j)];
						}
						else
							val = obs["CNNgrid_"+_nnfeatures[c]+"_cell"+to_string(i)+"_"+to_string(j)];
						input_tensor.set(fdeep::tensor_pos(i+(_ngrid-1)/2, j+(_ngrid-1)/2, c), val);
					}
					
				}
			}	
	
			//predict_class returns predicted class number and value of max output neuron
			pair<size_t, double> result = _nnmodel.predict_class_with_confidence({input_tensor});
			ovalue = result.second;
			return (int)result.first;
		}

	//returns vector sum of given objects
	Jet VectorSum(vector<Jet>& jets){
		Jet ret;
		for(auto j : jets){
			ret.add(j);
		}
		return ret;
	}
};
#endif
