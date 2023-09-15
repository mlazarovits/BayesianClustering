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

using std::string;
class BaseSkimmer{
	public:
		BaseSkimmer(){ };
		BaseSkimmer(TFile* file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			_file = file;
			//getting the stuff below from producer in derived class
			//TTree* tree = (TTree*)file->Get("tree/llpgtree");
			//_base = new ReducedBase(tree);
			//_nEvts = _base->fChain->GetEntries();
			//_base->GetEntry(0);
			//cout << "base skim init - " << _base->Photon_energy->size() << endl;
		
			_data = false;
			_debug = false;	
			hists1D.push_back(nSubClusters);
			hists1D.push_back(time_center);
			hists1D.push_back(eta_center);
			hists1D.push_back(phi_center);
			hists1D.push_back(slope_space);
			hists1D.push_back(slope_etaT);
			hists1D.push_back(slope_phiT);
			hists1D.push_back(polar_ang);
			hists1D.push_back(azimuth_ang);
			hists1D.push_back(e_avg);
			hists1D.push_back(rotundity_3D);
			hists1D.push_back(rotundity_2D);
			hists1D.push_back(e_avg_lead);
			hists1D.push_back(e_avg_sublead);
			hists1D.push_back(velocity);
			hists1D.push_back(npts_lead);
			hists1D.push_back(fracpts_lead);
			hists1D.push_back(fracE_lead);
			hists1D.push_back(eigen2D_ratio);	
			hists1D.push_back(objE);
			hists1D.push_back(clusterE);
			hists1D.push_back(rotundity_2D_lead);
			hists1D.push_back(rotundity_2D_notlead);
			hists1D.push_back(time_center_lead);
			hists1D.push_back(time_center_sublead);

			hists2D.push_back(time_totE);
			hists2D.push_back(time_totE_lead);
			hists2D.push_back(time_totE_sublead);
			hists2D.push_back(az_totE);
			hists2D.push_back(rot2D_totE);
			hists2D.push_back(eta_phi);
			hists2D.push_back(t_eta);
			hists2D.push_back(t_phi);
			hists2D.push_back(nsubcl_fracElead);
			hists2D.push_back(objE_clusterE);
			hists2D.push_back(t_mixcoeff);



		}
		virtual ~BaseSkimmer(){ 
			_file->Close();
			//delete _base;
			hists1D.clear();
			hists2D.clear();
		}

		virtual void CleaningSkim() = 0;
		virtual void Skim() = 0;

		TFile* _file;
		ReducedBase* _base = nullptr;
		int _nEvts;
		BaseProducer* _prod;
		bool _data;
		bool _debug;


		void SetData(bool d){ _data = d; }
		void SetDebug(bool d){ _debug = d; }

		string _cms_label;

		vector<TH1D*> hists1D;
		//0 - # of subclusters
		TH1D* nSubClusters = new TH1D("nSubClusters","nSubClusters",7,0,7.);
		//1 - mean time - center in t
		TH1D* time_center = new TH1D("time_center","time_center",50,-30,30);
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
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",50,-3.5,3.5);		
		//9 - subcluster energy - average
		TH1D* e_avg = new TH1D("e_avg","e_avg",50,0.,50.);
		//10 -ellipsoid rotundity
		TH1D* rotundity_3D = new TH1D("rotundity_3D","rotundity_3D",20,0,1.1);
		//11 - spatial rotundity
		TH1D* rotundity_2D = new TH1D("rotundity_2D","rotundity_2D",20,0,1.1);
		//12 - leading subcluster energy - average
		TH1D* e_avg_lead = new TH1D("e_avg_lead","e_avg_lead",50,0.,50.);
		//13 - subleading subcluster energy - average
		TH1D* e_avg_sublead = new TH1D("e_avg_sublead","e_avg_sublead",50,0.,50.);
		//14 - velocity = z/r*k for k transfer factor to velocity units
		TH1D* velocity = new TH1D("velocity","velocity",50,-200,200);
		//15 - number of pts in lead subcluster
		TH1D* npts_lead = new TH1D("npts_lead","npts_lead",50,0,100);
		//16 - fraction of pts in lead subcluster
		TH1D* fracpts_lead = new TH1D("fracpts_lead","fracpts_lead",50,0,1);
		//17 - fraction of energy in lead subcluster
		TH1D* fracE_lead = new TH1D("fracE_lead","fracE_lead",50,0.,1.);
		//18 - ratio of 2D eigenvals
		TH1D* eigen2D_ratio = new TH1D("eigen2D_ratio","eigen2D_ratio",50,0.,1.);
		//19 - total object energy
		TH1D* objE = new TH1D("objE","objE",50,0,100);
		//20 - cluster energy
		TH1D* clusterE = new TH1D("clusterE","clusterE",50,0,100);
		//21 - spatial rotundity for lead subcluster
		TH1D* rotundity_2D_lead = new TH1D("rotundity_2D_lead","rotundity_2D_lead",20,0,1.1);
		//22 - spatial rotundity for !lead subcluster
		TH1D* rotundity_2D_notlead = new TH1D("rotundity_2D_notlead","rotundity_2D_lead",20,0,1.1);
		//23 - time center for lead subcluster
		TH1D* time_center_lead = new TH1D("time_center_lead","time_center_lead",50,-30,30);
		//24 - time center for sublead subcluster
		TH1D* time_center_sublead = new TH1D("time_center_sublead","time_center_sublead",50,-30,30);


		//two dimensional histograms
		vector<TH2D*> hists2D;
		//0 - time v tot subcluster energy
		TH2D* time_totE = new TH2D("time_totE","time_totE;time_center;totalE;a.u.", 50,-30,30,50,0,50);
		//1 - time v tot lead subcluster energy
		TH2D* time_totE_lead = new TH2D("time_totE_lead","time_totE_lead;time_center_lead;totalE_lead;a.u.", 50,-30,30,50,0,50);
		//2 - time v tot sublead subcluster energy
		TH2D* time_totE_sublead = new TH2D("time_totE_sublead","time_totE_sublead:time_center_sublead;totalE_sublead;a.u.", 50,-30,30,50,0,50);
		//3 - azimuthal angle v energy
		TH2D* az_totE = new TH2D("az_totE","az_totE:azimuthal_angl;totalE;a.u.",50,-3.5,3.5,50,0,50);
		//4 - rotundity (2D) v energy
		TH2D* rot2D_totE = new TH2D("rot2D_totE","rot2D_totE;rotundity2D;totalE;a.u.",50,0.,1.1,50,0,50);
		//5 - eta v phi
		TH2D* eta_phi = new TH2D("eta_phi","eta_phi;eta_center;phi_center;a.u.",50,-3.5,3.5,50,-3.5,3.5);
		//6 - t v eta
		TH2D* t_eta = new TH2D("t_eta","t_eta;time_center;eta_center;a.u.",50,-30,30,50,-3.5,3.5);
		//7 - t v phi
		TH2D* t_phi = new TH2D("t_phi","t_phi;time_center;phi_center;a.u.",50,-30,30,50,-3.5,3.5);
		//8 - nsubClusters vs fraction of lead subcl E
		TH2D* nsubcl_fracElead = new TH2D("nsubcl_fracElead","nsubcl_fracElead;nSubClusters;fracE_lead;a.u.",7,0,7,50,0,1);
		//9 - objE (true) to cluster E (algo)
		TH2D* objE_clusterE = new TH2D("objE_clusterE","objE_clusterE;objE;clusterE;a.u.",50,0,100,50,0,100);
		//10 - time to mixing coeff
		TH2D* t_mixcoeff = new TH2D("t_mixcoeff","t_mixcoeff;time_center;mixing_coeff;a.u.",50,-30,30,50,0,1.);
		
	
		//struct for different types of plots (ie signal, ISR, fakes, etc.)
		struct plotCat{
			string legName;
			string plotName;
			vector<TH1D*> hists1D;
			vector<TH2D*> hists2D;
			vector<double> ids;
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
		


		void TDR2DHist(TH2D* hist, TCanvas* &can, string xtit, string ytit){
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
			
			string lat_cms = "#bf{CMS} #it{WIP} "+_cms_label;
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
