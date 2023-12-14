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
			_timesmear = true;
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
		
			_gev = 1/10.;
			_data = false;
			_debug = false;
			_timesmear = true;
			_hists1D.push_back(nSubClusters);
			_hists1D.push_back(time_center);
			_hists1D.push_back(eta_center);
			_hists1D.push_back(phi_center);
			_hists1D.push_back(slope_space);
			_hists1D.push_back(slope_etaT);
			_hists1D.push_back(slope_phiT);
			_hists1D.push_back(polar_ang);
			_hists1D.push_back(azimuth_ang);
			_hists1D.push_back(e_subcl);
			_hists1D.push_back(rotundity_3D);
			_hists1D.push_back(rotundity_2D);
			_hists1D.push_back(velocity);
			_hists1D.push_back(eigen2D_ratio);	
			_hists1D.push_back(clusterE);
                	_hists1D.push_back(etaSig); 
                	_hists1D.push_back(phiSig); 
                	_hists1D.push_back(timeSig);
			_hists1D.push_back(fracE);
			_hists1D.push_back(azimuth_ang_2D);
                	_hists1D.push_back(etaSig_pos); 
                	_hists1D.push_back(etaSig_neg); 
			_hists1D.push_back(etaphi_cov); 
			_hists1D.push_back(timeeta_cov);
			_hists1D.push_back(timephi_cov);
			_hists1D.push_back(objE);


			_hists2D.push_back(time_E);
			_hists2D.push_back(az_E);
			_hists2D.push_back(rot2D_E);
			_hists2D.push_back(eta_phi);
			_hists2D.push_back(t_eta);
			_hists2D.push_back(t_phi);
			_hists2D.push_back(t_mixcoeff);
			_hists2D.push_back(E_mixcoeff);
			_hists2D.push_back(rot3D_E);
			_hists2D.push_back(npts_E);
			_hists2D.push_back(nsubcl_nrhs);
			_hists2D.push_back(nsubcl_mmcoeff);
			_hists2D.push_back(etaSig_phiSig); 
                	_hists2D.push_back(timeSig_etaSig);
                	_hists2D.push_back(timeSig_phiSig);
			_hists2D.push_back(fracE_mmcoeff);
			_hists2D.push_back(nsubcl_fracE);
			_hists2D.push_back(timeSig_timeCenter);
			_hists2D.push_back(rot2D_az2D);
			_hists2D.push_back(timeSig_fracE);
			_hists2D.push_back(timeSig_totE);
			_hists2D.push_back(timeEtaCov_totE);
			_hists2D.push_back(timePhiCov_totE);
	

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
		void SetTransferFactor(double gev){ _gev = gev; _prod->SetTransferFactor(_gev); }

		vector<TH1D*> _hists1D;
		//0 - # of subclusters
		TH1D* nSubClusters = new TH1D("nSubClusters","nSubClusters",10,0,10.);
		//1 - mean time - center in t
		TH1D* time_center = new TH1D("time_center","time_center",50,-20,20);
		//2 - mean eta - center in eta
		TH1D* eta_center = new TH1D("eta_center","eta_center",50,-3.5,3.5);
		//3 - mean phi - center in phi
		TH1D* phi_center = new TH1D("phi_center","phi_center",50,-3.5,3.5);
		//4 - space slope
		TH1D* slope_space = new TH1D("slope_space","slope_space",50,-30,30);
		//5 - eta-time slope
		TH1D* slope_etaT = new TH1D("slope_etaT","slope_etaT",50,-2,2);
		//6 - phi-time slope
		TH1D* slope_phiT = new TH1D("slope_phiT","slope_phiT",50,-4,4);
		//7 - polar angle
		TH1D* polar_ang = new TH1D("polar_ang","polar_ang",50,-0.5,3.5);		
		//8 - azimuth angle
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",25,-3.5,3.5);		
		//9 - subcluster energy - total
		TH1D* e_subcl = new TH1D("e_subcl","e_subcl",50,0.,2000.);
		//10 - ellipsoid rotundity
		TH1D* rotundity_3D = new TH1D("rotundity_3D","rotundity_3D",10,0.74,1.01);
		//11 - spatial rotundity
		TH1D* rotundity_2D = new TH1D("rotundity_2D","rotundity_2D",20,0.4,1.1);
		//12 - velocity = z/r*k for k transfer factor to velocity units
		TH1D* velocity = new TH1D("velocity","velocity",31,-1.,30.);
		//13 - ratio of 2D eigenvals
		TH1D* eigen2D_ratio = new TH1D("eigen2D_ratio","eigen2D_ratio",50,0.,1.);
		//14 - cluster energy
		TH1D* clusterE = new TH1D("clusterE","clusterE",10,0,1000);
		//15 - eta sigma	
                TH1D* etaSig = new TH1D("etaSig","etaSig",25,0.01, 0.09);
		//16 - phi sigma	
                TH1D* phiSig = new TH1D("phiSig","phiSig",25,0.01,0.09);
		//17 - time sigma	
                TH1D* timeSig = new TH1D("timeSig","timeSig",25,0,25.);
		//18 - fraction of energy in cluster
		TH1D* fracE = new TH1D("fracE","fracE",50,0.,1.1);
		//19 - azimuth angle in 2D
		TH1D* azimuth_ang_2D = new TH1D("azimuth_ang_2D","azimuth_ang_2D",50,-3.5,3.5);		
		//20 - eta sigma for positive eta clusters
                TH1D* etaSig_pos = new TH1D("etaSig_pos","etaSig_pos",25,0.01, 0.06);
		//21 - eta sigma for negative eta clusters
                TH1D* etaSig_neg = new TH1D("etaSig_neg","etaSig_neg",25,0.01, 0.06);
		//22 - normalized covariance - eta/phi
		TH1D* etaphi_cov = new TH1D("etaphi_cov","etaphi_cov",10,0,0.01);
		//23 - normalized covariance - time/eta
		TH1D* timeeta_cov = new TH1D("timeeta_cov","timeeta_cov",10,0,0.01);
		//24 - normalized covariance - time/phi
		TH1D* timephi_cov = new TH1D("timephi_cov","timephi_cov",10,0,0.01);
		//25 - object energy
		TH1D* objE = new TH1D("objE","objE",50,0,100);


	
		//two dimensional histograms
		vector<TH2D*> _hists2D;
		//0 - time v subcl subcluster energy
		TH2D* time_E = new TH2D("time_subclE","time_subclE;time_center;E;a.u.", 50,-30,30,10,0,1000);
		//1 - azimuthal angle v subcl energy
		TH2D* az_E = new TH2D("az_subclE","az_subclE;azimuthal_angle;E;a.u.",50,-3.5,3.5,10,0,1000);
		//2 - rotundity (2D) v subcl energy
		TH2D* rot2D_E = new TH2D("rot2D_subclE","rot2D_subclE;rotundity2D;E;a.u.",50,0.4,1.1,10,0,1000);
		//3 - eta v phi
		TH2D* eta_phi = new TH2D("eta_phi","eta_phi;eta_center;phi_center",50,-3.5,3.5,50,-3.5,3.5);
		//4 - t v eta
		TH2D* t_eta = new TH2D("t_eta","t_eta;time_center;eta_center",50,-30,30,50,-3.5,3.5);
		//5 - t v phi
		TH2D* t_phi = new TH2D("t_phi","t_phi;time_center;phi_center;a.u.",50,-30,30,50,-3.5,3.5);
		//6 - time to mixing coeff
		TH2D* t_mixcoeff = new TH2D("t_mixcoeff","t_mixcoeff;time_center;mixing_coeff;a.u.",50,-20,20,50,0,1.1);
		//7 - subcl E vs mixing coeff
		TH2D* E_mixcoeff = new TH2D("subclE_mixcoeff","subclE_mixcoeff;E;mixing_coeff;a.u.",20,0,2000,20,0,1.1);	
		//8 - rotundity (3D) v subcluster energy
		TH2D* rot3D_E = new TH2D("rot3D_subclE","rot3D_subclE;rotundity3D;E;a.u.",50,0.4,1.1,50,0,100);
		//9 - npts effective (ie weighted) v subcl energy
		TH2D* npts_E = new TH2D("npts_E","npts_subclE;npts;E;a.u.",50,0.,100,10,0,1000);
		//10 - nsubclusters vs nrhs in cluster
		TH2D* nsubcl_nrhs = new TH2D("nsubcl_nrhs","nsubcl_nrhs;nsubclusters;nrhs",10,0,10,10,0,50);
		//11 - nsubclusters vs mm coeff
		TH2D* nsubcl_mmcoeff = new TH2D("nsubcl_mmcoeff","nsubcl_mmcoeff;nsubclusters;mmcoeff",10,0,10,20,0.,1.1);
                //12 - eta sigma v phi sigma
		TH2D* etaSig_phiSig = new TH2D("etaSig_phiSig","etaSig_phiSig;etaSig;phiSig;a.u.",25,0.01,0.09,25,0.01,0.09);
                //13 - time sigma v phi sigma
                TH2D* timeSig_etaSig = new TH2D("timeSig_etaSig","timeSig_etaSig;timeSig;etaSig;a.u.",25,0,25.,25,0.01,0.09);
                //14 - time sigma v phi sigma
                TH2D* timeSig_phiSig = new TH2D("timeSig_phiSig","timeSig_phiSig;timeSig;phiSig;a.u.",25,0,25.,25,0.01,0.09);
		//15 - fraction of energy in subcluster vs mm coeff of subcluster
		TH2D* fracE_mmcoeff = new TH2D("fracE_mmcoeff","fracE_mmcoeff;fracE;mmcoeff",20,0.,1.,20,0,1.1);
		//16 - number of subclusters vs fraction of energy in particular subcluster (really only applicable to lead subcluster)
		TH2D* nsubcl_fracE = new TH2D("nsubcl_fracE","nsubcl_fracE;nSubClusters;fracE",10,0,10.,20,0,1.1);
		//17 - time sigma vs time center
		TH2D* timeSig_timeCenter = new TH2D("timeSig_timeCenter","timeSig_timeCenter",25,0,25,50,-20,20);
		//18 - 2d rotundity vs 2d az angle
		TH2D* rot2D_az2D = new TH2D("rot2D_az2D","rot2D_az2D",50,0.4,1.1,50,-3.5,3.5);
		//19 - time variation vs frac E
		TH2D* timeSig_fracE = new TH2D("timeSig_fracE","timeSig_fracE;timeSig;fracE",25,0,25,20,0.,1.1);
		//20 - time sigma vs total E
		TH2D* timeSig_totE = new TH2D("timeSig_totE","timeSig_totE;timeSig;totE",25,0,25,20,0,2000);
		//21 - time-eta covariance vs total E
		TH2D* timeEtaCov_totE = new TH2D("timeEtaCov_totE","timeEtaCov_totE;timeEtaCov;totE",10,0,25,20,0,2000);
		//22 - time-phi covariance vs total E
		TH2D* timePhiCov_totE = new TH2D("timePhiCov_totE","timePhiCov_totE;timePhiCov;totE",10,0,25,20,0,2000);



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


		void SetCMSLabel(string lab){ _cms_label = lab; }

		void SetTimeSmear(bool t){ _timesmear = t; }
		bool _timesmear;


};
#endif
