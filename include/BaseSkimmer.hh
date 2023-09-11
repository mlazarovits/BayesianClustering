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
		}
		virtual ~BaseSkimmer(){ 
			_file->Close();
			delete _base;
		}

		virtual void CleaningSkim() = 0;
		virtual void Skim() = 0;

		TFile* _file;
		ReducedBase* _base = nullptr;
		int _nEvts;
		BaseProducer* _prod;

		string _cms_label;

		vector<TH1D*> hists1D;
		//# of subclusters
		TH1D* nSubClusters = new TH1D("nSubClusters","nSubClusters",7,0,7.);
		//mean time - center in t
		TH1D* time_center = new TH1D("time_center","time_center",50,-30,30);
		//mean eta - center in eta
		TH1D* eta_center = new TH1D("eta_center","eta_center",50,-3.5,3.5);
		//mean phi - center in phi
		TH1D* phi_center = new TH1D("phi_center","phi_center",50,-3.5,3.5);
		//space slope
		TH1D* slope_space = new TH1D("slope_space","slope_space",50,-30,30);
		//eta-time slope
		TH1D* slope_etaT = new TH1D("slope_etaT","slope_etaT",50,-2,2);
		//phi-time slop
		TH1D* slope_phiT = new TH1D("slope_phiT","slope_phiT",50,-4,4);
		//polar angle
		TH1D* polar_ang = new TH1D("polar_ang","polar_ang",50,-0.5,3.5);		
		//azimuth angle
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",50,-3.5,3.5);		
		//subcluster energy - average
		TH1D* e_avg = new TH1D("e_avg","e_avg",50,0.,50.);
		//ellipsoid rotundity
		TH1D* rotundity_3D = new TH1D("rotundity_3D","rotundity_3D",50,0,10);
		//spatial rotundity
		TH1D* rotundity_2D = new TH1D("rotundity_2D","rotundity_2D",50,0,10);
		//leading subcluster energy - average
		TH1D* e_avg_lead = new TH1D("e_avg_lead","e_avg_lead",50,0.,50.);
		//subleading subcluster energy - average
		TH1D* e_avg_sublead = new TH1D("e_avg_sublead","e_avg_sublead",50,0.,50.);


		//struct for different types of plots (ie signal, ISR, fakes, etc.)
		struct plotCat{
			string legName;
			string plotName;
			vector<TH1D*> hists1D;
			vector<double> ids;
		};



		void SetCMSLabel(string lab){ _cms_label = lab; }


		void TDRMultiHist(vector<TH1D*> &hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, vector<string> label){
			if(hist.size() != label.size()){ cout << "Error: number of histograms and labels don't match." << endl; return;}
			if(hist.size() == 0 || label.size() == 0) return;
			can->cd();
			can->SetGridx(1);
			can->SetGridy(1);
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
			lat.DrawLatex(0.4,0.92,plot_title.c_str());

			return;
		}
		void TDRHist(TH1D* hist, TCanvas* &can, string plot_title, string xtit, string ytit){
			can->cd();
			can->SetGridx(1);
			can->SetGridy(1);
			hist->UseCurrentStyle();
			hist->SetStats(false);
			hist->GetXaxis()->CenterTitle(true);
			hist->GetXaxis()->SetTitle(xtit.c_str());
			hist->GetYaxis()->CenterTitle(true);
			hist->GetYaxis()->SetTitle(ytit.c_str());
			hist->Draw();
			
			string lat_cms = "#bf{CMS} #it{WIP} 2017 GMSB";
			TLatex lat;
			lat.SetNDC();
			lat.SetTextSize(0.04);
			lat.SetTextFont(42);
			lat.DrawLatex(0.02,0.92,lat_cms.c_str());
			lat.DrawLatex(0.4,0.92,plot_title.c_str());

			return;
		}


		void FindListHistBounds(vector<TH1D*>& hists, double& ymin, double& ymax){
			//shellsort to find max, min
			int N = (int)hists.size()-1;
			if(N < 1){ ymax = 0; ymin = 0; return; }
			int i, j, h;
			TH1D* v = nullptr;
			for(h = 1; h <= N/9; h = 3*h+1) ;
			for( ; h > 0; h /= 3)
				for(i = h+1; i <= N; i += 1){
					v = hists[i]; j = i;
					while(j > h && hists[j - h]->GetMaximum() > v->GetMaximum()){ hists[j] = hists[j - h]; j -= h; }
					hists[j] = v;
				}
		
			ymax = hists[0]->GetMaximum();
			ymin = hists[N-1]->GetMinimum();

		}


};
#endif
