#ifndef BaseProducer_HH
#define BaseProducer_HH

#include "ReducedBase.hh"
#include "JetPoint.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TCanvas.h"

class BaseProducer{
	public:
		BaseProducer(){ };
		BaseProducer(TFile* file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			_file = file;
			TTree* tree = (TTree*)file->Get("tree/llpgtree");
			_base = new ReducedBase(tree);
			_nEvts = _base->fChain->GetEntries();
		}
		virtual ~BaseProducer(){ };

		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		virtual void GetRecHits(vector<vector<JetPoint>>& rhs) = 0;
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt) = 0;
		virtual void GetPrimaryVertex(Point& vtx, int evt) = 0;

		virtual void CleaningSkim() = 0;
		virtual void Skim() = 0;

		TFile* _file;
		ReducedBase* _base = nullptr;
		int _nEvts;

		vector<TH1D*> TH1D_hists;
		//subcluster energy - average
		TH1D* e_avg = new TH1D("e_avg","e_avg",100,0.,50.);
		//space slope
		TH1D* slope_space = new TH1D("slope_space","slope_space",50,-30,30);
		//eta-time slope
		TH1D* slope_etaT = new TH1D("slope_etaT","slope_etaT",50,-2,2);
		//phi-time slop
		TH1D* slope_phiT = new TH1D("slope_phiT","slope_phiT",50,-4,4);
		//mean time - center in t
		TH1D* time_center = new TH1D("time_center","time_center",50,-30,30);
		//mean eta - center in eta
		TH1D* eta_center = new TH1D("eta_center","eta_center",50,-3.5,3.5);
		//mean phi - center in phi
		TH1D* phi_center = new TH1D("phi_center","phi_center",50,-3.5,3.5);
		//polar angle
		TH1D* polar_ang = new TH1D("polar_ang","polar_ang",50,-3.5,3.5);		
		//azimuth angle
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",50,-3.5,3.5);		


		//# of subclusters
		TH1I* nSubClusters = new TH1I("nSubClusters","nSubClusters",7,0,7.);
		//# of subclusters vs. photon reco energy
		TH2D* e_nSubClusters = new TH2D("e_nSubClusters","e_nSubClusters",50,0.,1000.,7,0.,7.);	



			


		void TDRMultiHist(vector<TH1D*> &hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, vector<string> label){
			can->cd();
			TLegend* myleg = new TLegend( 0.6, 0.6, 0.8, 0.8 );
			myleg->SetFillColor(0);
			myleg->SetBorderSize(0);
			myleg->SetTextSize(0.025);
			myleg->SetHeader( (label[0]).c_str() );
			for( int i = 0 ; i < int(hist.size()); i++){
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
					myleg->AddEntry( hist[i], (label[i+1]).c_str(), "L" );
					gPad->Update();
				}else{
					hist[i]->Draw("same");
					myleg->AddEntry( hist[i], (label[i+1]).c_str(), "L" );
					gPad->Update();
				}
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
