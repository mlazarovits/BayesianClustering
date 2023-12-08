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
                	_hists1D.push_back(etaVar); 
                	_hists1D.push_back(phiVar); 
                	_hists1D.push_back(timeVar);
			_hists1D.push_back(fracE);
			_hists1D.push_back(azimuth_ang_2D);


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
			_hists2D.push_back(etaVar_phiVar); 
                	_hists2D.push_back(timeVar_etaVar);
                	_hists2D.push_back(timeVar_phiVar);
			_hists2D.push_back(fracE_mmcoeff);
			_hists2D.push_back(nsubcl_fracE);
	

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
		void SetTransferFactor(double gev){ _gev = gev; }

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
		//15 - eta variance	
                TH1D* etaVar = new TH1D("etaVar","etaVar",25,0, 1e-1);
		//16 - phi variance	
                TH1D* phiVar = new TH1D("phiVar","phiVar",25,0,1e-1);
		//17 - time variance	
                TH1D* timeVar = new TH1D("timeVar","timeVar",25,0,50.);
		//18 - fraction of energy in cluster
		TH1D* fracE = new TH1D("fracE","fracE",50,0.,1.1);
		//19 - azimuth angle in 2D
		TH1D* azimuth_ang_2D = new TH1D("azimuth_ang_2D","azimuth_ang_2D",50,-3.5,3.5);		

	
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
                //12 - eta variance v phi variance
		TH2D* etaVar_phiVar = new TH2D("etaVar_phiVar","etaVar_phiVar;etaVar;phiVar;a.u.",25,0,1e-1,25,0,1e-1);
                //13 - time variance v phi variance
                TH2D* timeVar_etaVar = new TH2D("timeVar_etaVar","timeVar_etaVar;timeVar;etaVar;a.u.",25,0,50.,25,0,1e-1);
                //14 - time variance v phi variance
                TH2D* timeVar_phiVar = new TH2D("timeVar_phiVar","timeVar_phiVar;timeVar;phiVar;a.u.",25,0,50.,25,0,1e-1);
		//15 - fraction of energy in subcluster vs mm coeff of subcluster
		TH2D* fracE_mmcoeff = new TH2D("fracE_mmcoeff","fracE_mmcoeff;fracE;mmcoeff",20,0.,1.,20,0,1.1);
		//16 - number of subclusters vs fraction of energy in particular subcluster (really only applicable to lead subcluster)
		TH2D* nsubcl_fracE = new TH2D("nsubcl_fracE","nsubcl_fracE;nSubClusters;fracE",10,0,10.,20,0,1.1);


		//reco object histograms
		//NOT in hist vectors
		TH1D* objE = new TH1D("objE","objE",50,0,100);
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


		void TDRMultiHist(vector<TH1D*> &hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, vector<string> label){
			if(hist.size() != label.size()){ cout << "Error: number of histograms and labels don't match." << endl; return;}
			if(hist.size() == 0 || label.size() == 0) return;
			can->cd();
			can->SetGridx(1);
			can->SetGridy(1);
			can->SetTitle("");
			TLegend* myleg = new TLegend(0.7, 0.7, 0.8, 0.8 );
			myleg->SetFillColor(0);
			myleg->SetBorderSize(0);
			myleg->SetTextFont(42);
			myleg->SetTextSize(0.04);
		
			//offset for log scale
			if(miny == 0) miny += 1e-6;
			vector<int> k = {kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack};	


			for( int i = 0 ; i < int(hist.size()); i++){
				hist[i]->UseCurrentStyle();
				hist[i]->SetStats(false);
				hist[i]->GetXaxis()->CenterTitle(true);
				hist[i]->GetXaxis()->SetTitle(xtit.c_str());
				hist[i]->GetYaxis()->CenterTitle(true);
				hist[i]->GetYaxis()->SetTitle(ytit.c_str());
				hist[i]->GetYaxis()->SetRangeUser(miny, maxy + maxy/10.);
				hist[i]->SetLineColor(i+2);
				hist[i]->SetMarkerStyle(i+25);
				hist[i]->SetMarkerColor(i+2);
				if( i == 0 ){
					hist[i]->Draw("ep");
				}else{
					hist[i]->Draw("epsame");
				}
				myleg->AddEntry( hist[i], (label[i]).c_str(), "p" );
				gPad->Update();
			}
			myleg->Draw("same"); 
			gPad->Update();
			string lat_cms = "#bf{CMS} #it{WIP} "+_cms_label;
			TLatex lat;
			lat.SetNDC();
			lat.SetTextSize(0.04);
			lat.SetTextFont(42);
			lat.DrawLatex(0.02,0.92,lat_cms.c_str());
			//lat.DrawLatex(0.4,0.92,plot_title.c_str());

			return;
		}
		void TDRHist(TH1D* hist, TCanvas* &can, string plot_title, string xtit, string ytit){
			can->cd();
			can->SetGridx(1);
			can->SetGridy(1);
			hist->SetTitle("");
			hist->UseCurrentStyle();
			hist->SetStats(false);
			hist->GetXaxis()->CenterTitle(true);
			hist->GetXaxis()->SetTitle(xtit.c_str());
			hist->GetYaxis()->CenterTitle(true);
			hist->GetYaxis()->SetTitle(ytit.c_str());
			hist->Draw();
			
			string lat_cms = "#bf{CMS} #it{WIP} "+_cms_label;
			TLatex lat;
			lat.SetNDC();
			lat.SetTextSize(0.04);
			lat.SetTextFont(42);
			lat.DrawLatex(0.02,0.92,lat_cms.c_str());
			//lat.DrawLatex(0.4,0.92,plot_title.c_str());

			return;
		}
		


		void TDR2DHist(TH2D* hist, TCanvas* &can, string xtit, string ytit, string title = ""){
			can->cd();
			//can->SetGridx(1);
			//can->SetGridy(1);
			hist->SetTitle("");
			hist->UseCurrentStyle();
			hist->SetStats(false);
			hist->GetXaxis()->CenterTitle(true);
			hist->GetXaxis()->SetTitle(xtit.c_str());
			hist->GetYaxis()->CenterTitle(true);
			hist->GetYaxis()->SetTitle(ytit.c_str());
			hist->Draw("colz");
			
			string lat_cms = "#bf{CMS} #it{WIP} "+_cms_label+" "+title;
			TLatex lat;
			lat.SetNDC();
			lat.SetTextSize(0.04);
			lat.SetTextFont(42);
			lat.DrawLatex(0.02,0.92,lat_cms.c_str());
			//lat.DrawLatex(0.4,0.92,plot_title.c_str());

			return;
		}


		void FindListHistBounds(vector<TH1D*>& hists, double& ymin, double& ymax){
			//insert to find max, min
			int N = (int)hists.size();
			if(N < 1){ ymax = 0; ymin = 0; return; }
			ymax = -999;
			ymin = 999;
			for(int i = 0; i < N; i++){
				if(hists[i]->GetMaximum() > ymax) ymax = hists[i]->GetMaximum();
				if(hists[i]->GetMinimum() < ymin) ymin = hists[i]->GetMinimum();

			}		

		}




};
#endif
