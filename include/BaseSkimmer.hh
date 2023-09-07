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
		TH1D* polar_ang = new TH1D("polar_ang","polar_ang",50,-3.5,3.5);		
		//azimuth angle
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",50,-3.5,3.5);		
		//subcluster energy - average
		TH1D* e_avg = new TH1D("e_avg","e_avg",100,0.,50.);



		//struct for different types of plots (ie signal, ISR, fakes, etc.)
		struct plotCat{
			string legName;
			string plotName;
			vector<TH1D*> hists1D;
			vector<double> ids;
		};




		void TDRMultiHist(vector<TH1D*> &hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, vector<string> label){
			cout << "TDRHist - name: " << xtit << endl;
			if(hist.size() != label.size()){ cout << "Error: number of histograms and labels don't match." << endl; return;}
			if(hist.size() == 0 || label.size() == 0) return;
			can->cd();
			TLegend* myleg = new TLegend( 0.6, 0.6, 0.8, 0.8 );
			myleg->SetFillColor(0);
			myleg->SetBorderSize(0);
			myleg->SetTextSize(0.025);
			myleg->SetHeader( (label[0]).c_str() );
			for( int i = 0 ; i < int(hist.size()); i++){
				cout << "i - " << i << " label - " << label[i] << endl; 
				hist[i]->UseCurrentStyle();
				hist[i]->SetStats(false);
				hist[i]->GetXaxis()->CenterTitle(true);
				hist[i]->GetXaxis()->SetTitle(xtit.c_str());
				hist[i]->GetYaxis()->CenterTitle(true);
				hist[i]->GetYaxis()->SetTitle(ytit.c_str());
				hist[i]->GetYaxis()->SetRangeUser(miny, maxy);
				hist[i]->SetLineColor(i+2);
				hist[i]->SetMarkerStyle(i+25);
				if( i == 0 ){
					hist[i]->Draw();
				}else{
					hist[i]->Draw("same");
				}
				myleg->AddEntry( hist[i], (label[i]).c_str(), "L" );
				gPad->Update();
			}
			myleg->Draw("same"); 
			gPad->Update();
			string lat_cms = "#bf{CMS} #it{Preliminary} (13 TeV)";
			TLatex lat;
			lat.SetNDC();
			lat.SetTextSize(0.05);
			lat.SetTextFont(42);
			lat.DrawLatex(0.15,0.9325,lat_cms.c_str());
			lat.DrawLatex((0.82),0.9325,plot_title.c_str());

			return;


		}




};
#endif
