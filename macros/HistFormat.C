#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>

using std::string;


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


void TDRMultiHist(vector<TH1D*> &hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, vector<string> label, string cms_label){
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
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());
	//lat.DrawLatex(0.4,0.92,plot_title.c_str());

	return;
}

void TDR2DHist(TH2D* hist, TCanvas* &can, string xtit, string ytit, string cms_label, string title = ""){
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
	
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label+" "+title;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());
	//lat.DrawLatex(0.4,0.92,plot_title.c_str());

	return;
}



void TDRHist(TH1D* hist, TCanvas* &can, string plot_title, string xtit, string ytit, string cms_label){
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
	
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());
	//lat.DrawLatex(0.4,0.92,plot_title.c_str());

	return;
};


void GetHists(TDirectory* dir, string type, vector<TH1D*>& hists, vector<string>& labels){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th1d("TH1D");
	labels.clear();
	hists.clear();

	vector<string> procs = {"chiGam","notSunm"};
	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d){
			//get 1D histograms
			TH1D* hist = (TH1D*)key->ReadObj();
			if(hist){
				name = hist->GetName();
				if(type.empty()){
					if(name.find("lead") != string::npos) continue;
				}
				else{
					if(name.find("_"+type) == string::npos) continue;
				}
				hist->Scale(1./hist->Integral());
				hists.push_back(hist);
				for(int i = 0; i < procs.size(); i++){
					if(name.find(procs[i]) != string::npos){
						labels.push_back(procs[i]);
						break;
					}
				}
			}
		}
	}
}


string GetExtraLabel(string in_file){
	int idx = in_file.find("Skim_")+4;
	int idx2 = in_file.find("_emAlpha");
	//no extra string in between
	if(idx == idx2) return "";
	else return in_file.substr(idx+1,idx2-idx-1);

}

string GetCMSLabel(string in_file){
	//get version
	string cmslab;
	int idx = in_file.find("NperGeV0p");
	//based on NperGeV0pXXX being the last label before cms label
	//and 3 digits following "GeV0p"
	cmslab = in_file.substr(idx+13);
	cmslab = cmslab.substr(0,cmslab.find(".root"));
	return cmslab;
}

void HistFormat(string file){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str());
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name, xtitle, ytitle;
	string oname = f->GetName();
	oname = oname.substr(0,oname.find(".root"));
	oname = oname+"_formatted.root";
	TFile* ofile = new TFile(oname.c_str(),"RECREATE");

	vector<string> types = {"","lead","notlead"};

	string cmslab = GetCMSLabel(file);
	string extra = GetExtraLabel(file);
	if(!extra.empty()) cmslab += " "+extra;	

	TString th1d("TH1D");
	TString th2d("TH2D");
	TString tdir("TDirectoryFile");
	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d){
			//get 1D histograms
			TH1D* hist = (TH1D*)key->ReadObj();
			if(hist){
				name = hist->GetName();
				TCanvas *cv = new TCanvas(name.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				TDRHist(hist, cv, name, name, "a.u",cmslab);
				cv->Write(); 
			}
		}
		if(key->GetClassName() == th2d){
			//get 2D histograms
			TH2D* hist = (TH2D*)key->ReadObj();
			if(hist){
				name = hist->GetName();
				xtitle = hist->GetXaxis()->GetTitle();
				ytitle = hist->GetYaxis()->GetTitle();
				TCanvas *cv = new TCanvas(name.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				TDR2DHist(hist, cv, xtitle, ytitle, cmslab);
				cv->Write(); 
			}
		}
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = (TDirectory*)key->ReadObj();
			double ymin, ymax;
			vector<string> labels;
			if(dir){
				dir->cd();
				vector<TH1D*> hists;
				//loop throught types
				for(int i = 0; i < types.size(); i++){
					//loop through hists in dir
					name = dir->GetName();
					if(!types[i].empty()) name += "_"+types[i];
					GetHists(dir, types[i], hists, labels);
					if(hists.size() < 1) continue;
					FindListHistBounds(hists, ymin, ymax);
					if(ymin == 0 && ymax == 0) continue;
					TCanvas *cv = new TCanvas(name.c_str(), "");
					ofile->cd();
					//draw as tcanvases
					TDRMultiHist(hists, cv, name, name, "a.u", ymin, ymax, labels, cmslab);
					cv->Write(); 
				}
			}
		}
	}
	cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;


	ofile->Close();

};
