#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>
#include <sstream>
#include <iostream>
using std::string;


void GetRatioHistRange(TH1D* ratio_hist, double& ratio_axis_min, double& ratio_axis_max){
	// make y-axis in ratio plot be symmetric and cover all data and mc points
	// ratio_hist doesn't return the correct value for GetYaxis()->GetXmax/min() so we brute force it
	// start with getting the min and max deviation from the first bin (0 bin is underflow)
	ratio_axis_min = ratio_hist->GetBinContent(1)-ratio_hist->GetBinError(1);
	ratio_axis_max = ratio_hist->GetBinContent(1)+ratio_hist->GetBinError(1);
	
	// loop over remaining bins to find max deviation up and down for the data (b = nbins is last bin)
	for (int b = 2; b <= ratio_hist->GetNbinsX(); b++) {
	    if(ratio_hist->GetBinContent(b)-ratio_hist->GetBinError(b) < ratio_axis_min)
	      ratio_axis_min = ratio_hist->GetBinContent(b)-ratio_hist->GetBinError(b);
	    if(ratio_hist->GetBinContent(b)+ratio_hist->GetBinError(b) > ratio_axis_max)
	      ratio_axis_max = ratio_hist->GetBinContent(b)+ratio_hist->GetBinError(b);
	  }
	
	// round to nearest three decimal places
	ratio_axis_min = std::floor(ratio_axis_min*1000.)/1000.;
	ratio_axis_max = std::ceil(ratio_axis_max*1000.)/1000.;
	
	// make it so we have the same range from 1 for min and max
	//   example: 0.9, 1.15 becomes 0.85, 1.15
	if(fabs(1.-ratio_axis_min) > fabs(ratio_axis_max-1.)){
	  ratio_axis_max = 1.-ratio_axis_min+1.;
	}
	else if(fabs(ratio_axis_max-1.) > fabs(1.-ratio_axis_min)){
	  ratio_axis_min = 1.+1.-ratio_axis_max;
	}
	
	// add offset so there is a bit of white space between histograms and axis
	if(ratio_axis_max > 1.5){
	  ratio_axis_min -= 0.11;
	  ratio_axis_max += 0.11;
	}
	else if(ratio_axis_max > 1.25){
	  ratio_axis_min -= 0.06;
	  ratio_axis_max += 0.06;
	}
	else{
	  ratio_axis_min -= 0.02;
	  ratio_axis_max += 0.02;
	}
	
	// minimum deviation from 1 is 7% so we can clearly see 5% on the axis label (5% is the default minimum deviation to show)
	if(ratio_axis_min > 0.93 && ratio_axis_max < 1.07){
	  ratio_axis_min = 0.93;
	  ratio_axis_max = 1.07;
	}
}


void stringReplaceAll(string& input, string& oldstring, string& newstring){
	while(input.find(oldstring) != string::npos){
		input.replace(input.find(oldstring), oldstring.size(), newstring);
	}

}

string SignalLegEntry(string label){
	//cout << "label " << label << endl;
	string lambda, ctau;
	string lmatch = "_L";
	string sample_l = label.substr(label.find(lmatch)+2);
	lambda = sample_l.substr(0,sample_l.find("_"));
	
	string ctmatch = "Ctau";
	string sample_ctau = label.substr(label.find(ctmatch)+ctmatch.size());
	ctau = sample_ctau.substr(0,sample_ctau.find("_"));
	
	//cout << "sample_l "<< sample_l << " lambda " << lambda << " ctau " << ctau << endl;

	string lfancy = lambda.substr(lambda.find("_")+1);
	//lfancy.insert(lfancy.find("TeV")-2," ");
	lfancy += " TeV";
	string ctfancy = ctau.substr(ctau.find("_")+1);
	if(ctfancy.find("p") != string::npos)
		ctfancy.replace(ctfancy.find("p"),1,".");
	//ctfancy.insert(ctfancy.find("cm")," ");
	ctfancy += " cm";

	//string legName = "#Chi^{0} #rightarrow #gamma, L = "+lfancy+" c#tau = "+ctfancy;
	string legName = "GMSB, L = "+lfancy+" c#tau = "+ctfancy;
	return legName;
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


void TDRMultiHist(vector<TH1D*> hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, string cms_label, bool b_ratio = false){
	double hlo = 0.12;
	double hhi = 0.2;
	double hbo = 0.05;
	double hto = 0.07;
	double ratio_h = 0.25;
	if(can == nullptr) return;
	if(hist.size() == 0)return;
	string title, name, legentry;
	string canname = can->GetName();

	can->cd();
	can->SetGridx(1);
	can->SetGridy(1);
	can->SetTitle(xtit.c_str());
	//can->SetLeftMargin(hlo);
	//can->SetRightMargin(hhi);
	//can->SetBottomMargin(hbo);
	//can->SetTopMargin(hto);
	can->SetGridy();
	//can->SetLogy();
	can->Draw();
	can->cd();

	double pad_xlo, pad_xhi, pad_ylo, pad_yhi;
	if(b_ratio){
		pad_xlo = 0.;
		pad_ylo = 0.;
		pad_xhi = 1.;
		pad_yhi = 1.;
	}
	else{
		pad_xlo = 0.;
		pad_ylo = 0.;
		pad_xhi = 1.;
		pad_yhi = 1.;
	}	
	TPad* pad = new TPad(("pad_"+canname).c_str(),
			      ("pad_"+canname).c_str(),
			      pad_xlo, pad_ylo, pad_xhi, pad_yhi);
	
	if(b_ratio)
		pad->SetBottomMargin(hbo+ratio_h);

	if(b_ratio){
		pad->SetLeftMargin(hlo);
		//pad->SetRightMargin(hhi);
		pad->SetTopMargin(hto);
		pad->SetGridy();
		pad->SetGridx();
		//pad->SetLogy();
	}
	pad->Draw();
	pad->cd();



	TLegend* myleg = nullptr;
	if(cms_label.find("photon") != string::npos) myleg = new TLegend(0.6, 0.7, 0.9, 0.9);
	else{
		if(b_ratio)
			myleg = new TLegend(0.563,0.701,0.78,0.890);
		else
			myleg = new TLegend(0.73,0.675,0.89,0.882);
		//myleg->SetMargin(0.1);
	}
	myleg->SetFillColor(0);
	myleg->SetBorderSize(0);
	myleg->SetTextFont(42);
	if(hist.size() > 4) myleg->SetTextSize(0.025);
	else myleg->SetTextSize(0.04);
	
	//offset for log scale
	if(miny == 0) miny += 1e-6;

	//sort hists alphabetically
	map<string, TH1D*> nameToHist;
	for( int i = 0 ; i < int(hist.size()); i++){
		name = hist[i]->GetName();
		nameToHist[name] = hist[i];
		hist[i] = nullptr;
	}
	int h = 0;
	for(map<string, TH1D*>::iterator it = nameToHist.begin(); it != nameToHist.end(); it++){
		hist[h] = it->second;
		h++;
	}
	//make color map
	TColor::SetColorThreshold(0.09);
	map<string, int> labelToColor;
	map<string, int> labelToMark;
	int offset = 1;
	labelToColor["chiGam"] =  TColor::GetColor("#86bbd8");
	labelToColor["GMSB"] =  TColor::GetColor("#86bbd8");

	labelToColor["notSunm"] =  TColor::GetColor("#9e0059");
	labelToColor["GJets"] 	=  TColor::GetColor("#f6ae2d");
	labelToColor["JetHT"] 	=  TColor::GetColor("#3d348b");
	labelToColor["MET"]   	=  TColor::GetColor("#671E76");
	labelToColor["ttbar"] 	=  TColor::GetColor("#CA5743");
	labelToColor["QCD"]   	=  TColor::GetColor("#9E0059");
	labelToColor["DoubleEG"] = TColor::GetColor("#9e0059");
	//later colors to use
	//"#CA5743" - jasper (burnt orange)
	//"#BEB583" - sage

	labelToColor["!median"] = TColor::GetColor("#f7a278");
	labelToColor["!eAvg"] = TColor::GetColor("#6859f1");
	labelToColor["!mmAvg"] = TColor::GetColor("#52b788");
	labelToColor["!eMax"] = TColor::GetColor("#E2C2FF");

	//MC symbols - primary shapes
	labelToMark["!chiGam"] =  20;
	labelToMark["!GMSB"] =  20;
	labelToMark["!notSunm"] = 72;
	labelToMark["!GJets"] =   73;
	labelToMark["ttbar"] = 24;
	labelToMark["QCD"] = 25;
	//data symbols - some form of open cross
	labelToMark["!JetHT"] =   75;
	labelToMark["!MET"] =   85;
	labelToMark["!DoubleEG"] =   83;
	
	//signal point additions
	if(cms_label.find("photons") != string::npos){
		labelToMark["L150"] = -1;
		labelToMark["L350"] = 1;
		labelToMark["Ctau0p1"] = 1;
		labelToMark["Ctau200"] = 2;
	}

	labelToMark["median"] = 71;
	labelToMark["eAvg"] =   72; 
	labelToMark["mmAvg"] =  73;
	labelToMark["eMax"] =   74; 

	int col, mark;	
	for( int i = 0 ; i < int(hist.size()); i++){
		hist[i]->UseCurrentStyle();
		hist[i]->SetStats(false);
		hist[i]->GetXaxis()->CenterTitle(true);
		if(!b_ratio)
			hist[i]->GetXaxis()->SetTitle(xtit.c_str());
		else
			hist[i]->GetXaxis()->SetTitle("");
		hist[i]->GetYaxis()->CenterTitle(true);
		hist[i]->GetYaxis()->SetTitle(ytit.c_str());
		hist[i]->GetYaxis()->SetRangeUser(miny, maxy + maxy/10.);
		
		if(b_ratio)
			hist[i]->GetXaxis()->SetLabelOffset(0.5);

		title = hist[i]->GetTitle();
		legentry = title.empty() ? hist[i]->GetName() : title;
		//if a key from labeltocolor is in legentry, set that color
		for(map<string, int>::iterator it = labelToColor.begin(); it != labelToColor.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(legentry.find(match) != string::npos){
				col = it->second;
				if(title.empty()) legentry = match;
				break;
			}
			else col = 1;
		}
		for(map<string, int>::iterator it = labelToMark.begin(); it != labelToMark.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(legentry.find(match) != string::npos){
				mark = it->second;
				break;
			}
			else mark = 1;
		}
		if(cms_label.find("photon") != string::npos){
			//do different signal points
			if(legentry.find("chiGam") != string::npos){
				for(map<string, int>::iterator it = labelToMark.begin(); it != labelToMark.end(); it++){
					string match = it->first;
					if(match.find("!") != string::npos) continue; //already looped over above
					if(legentry.find(match) != string::npos){
						mark += it->second; //need to loop over all additions bc signal point has two (L, ctau)
					}
				}
			}
		}
		
		hist[i]->SetLineColor(col);
		//hist[i]->SetLineWidth(2);
		hist[i]->SetMarkerStyle(mark);
		hist[i]->SetMarkerColor(col);
		//hist[i]->SetMarkerSize(1);
		if(i == 0)
			title = hist[i]->GetTitle();
		hist[i]->SetTitle("");
		if( i == 0 ){
			hist[i]->Draw("ep");
		}else{
			hist[i]->Draw("epsame");
		}
		myleg->AddEntry( hist[i], legentry.c_str(), "p" );
		pad->Update();
	}


	title = title.substr(title.find("_")+1);
	if(title.find("chiGam") != string::npos){
		title = SignalLegEntry(title);
	}
	pad->cd();
	myleg->Draw("same"); 
	pad->Update();
	//do ratio if specified
	if(b_ratio){
		if(hist.size() > 2){
			cout << "taking ratio of first two hists: " << hist[0]->GetName() << " to " << hist[1]->GetName() << endl;
		}
		TH1D* ratio_hist = (TH1D*)hist[0]->Clone(("data_ratio_"+canname).c_str());
		ratio_hist->SetLineColor(kBlack);
		ratio_hist->SetMarkerColor(kBlack);
		ratio_hist->SetFillColor(kWhite);
		ratio_hist->SetMarkerStyle(8);
		ratio_hist->SetMarkerSize(0.7);
		//hist order should be data, MC so ratio is data/MC
		ratio_hist->Divide(hist[1]);
		double ratio_min, ratio_max;
		GetRatioHistRange(ratio_hist, ratio_min, ratio_max);

		ratio_hist->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
		can->cd();
		double offset = 0.02;
		TPad* pad_ratio = new TPad(("pad_"+canname).c_str(),
		    		 ("pad_"+canname).c_str(),
		    		 0.0, hbo, 1., hbo+ratio_h - offset);
		pad_ratio->SetBottomMargin(0.15);
		pad_ratio->SetLeftMargin(hlo);
		//pad_ratio->SetRightMargin(0.);
		pad_ratio->SetTopMargin(0.);
		pad_ratio->SetGridy();
		pad_ratio->SetGridx();
		pad_ratio->Draw();
		pad_ratio->cd();

		ratio_hist->GetXaxis()->SetLabelOffset(0.005);
		ratio_hist->GetXaxis()->SetLabelSize(0.16);
		ratio_hist->GetYaxis()->CenterTitle();
		ratio_hist->GetYaxis()->SetTitleFont(42);
		ratio_hist->GetYaxis()->SetTitleSize(0.16);
		ratio_hist->GetYaxis()->SetTitleOffset(0.29);
		ratio_hist->GetYaxis()->SetLabelFont(42);
		ratio_hist->GetYaxis()->SetLabelSize(0.16);
		ratio_hist->GetXaxis()->SetTickLength(0.08);
		ratio_hist->GetYaxis()->SetTitle("#frac{data}{MC}");

		ratio_hist->Draw("ep");
		pad_ratio->Update();
		can->cd();
		pad->cd();
	}
	string lat_cms = "#bf{CMS} #it{Work in Progress} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.03,0.94,lat_cms.c_str());
	TLatex lat1;
	lat1.SetNDC();
	lat1.SetTextSize(0.04);
	lat1.SetTextFont(42);
	if(b_ratio)
		lat1.DrawLatex(0.50,0.94,plot_title.c_str());
	
}

void TDR2DHist(TH2D* hist, TCanvas* &can, string xtit, string ytit, string cms_label, string title = ""){
	can->SetGridx();
	can->SetGridy();
	can->cd();
	hist->SetTitle("");
	hist->UseCurrentStyle();
	hist->SetStats(false);
	hist->GetXaxis()->CenterTitle(true);
	hist->GetXaxis()->SetTitle(xtit.c_str());
	hist->GetYaxis()->CenterTitle(true);
	hist->GetYaxis()->SetTitle(ytit.c_str());
	string histname = hist->GetName();
	if((hist->GetNbinsX() == 2 && hist->GetNbinsY() == 2) || histname.find("geoEavg_genDeltaTime_meanRecoGenDeltaT") != string::npos ){
		if(histname.find("geoEavg_genDeltaTime_meanRecoGenDeltaT") != string::npos)
			hist->SetMarkerSize(1.3);
		//count histograms
		else{
			hist->Scale(1./hist->Integral());
			hist->SetMarkerSize(3.);
		}
		hist->Draw("colz1text");
	}
	else hist->Draw("colz1");

	string name = hist->GetName();
	//if(name.find("genDeltaTime_meanRecoGenDeltaT") != string::npos) cout << "n entries: " << hist->GetEntries() << endl;
	
	string lat_cms = "#it{Work In Progress} "+cms_label+" "+title;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.025);
	lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());

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
	string name = hist->GetName();
	if(hist->Integral() != 0 && name.find("sigma") == string::npos) hist->Scale(1./hist->Integral());
	hist->Draw();
	
	string lat_cms = "#it{Work In Progress} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.025);
	lat.SetTextFont(42);
	lat.DrawLatex(0.1,0.92,lat_cms.c_str());

	return;
};


void GetHistsProc(TDirectory* dir, string& proc, vector<TH1D*>& hists){
	if(!dir) return;
	dir->cd();
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString tdir("TDirectoryFile");
	TString th1d("TH1D");
	hists.clear();
	while((key = (TKey*)iter())){
		if(key->GetClassName() == tdir){
			TDirectory* ddir = dynamic_cast<TDirectory*>(key->ReadObj());
			if(!ddir) continue;
			ddir->cd();
			TList* llist = ddir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name;
			while((kkey = (TKey*)iiter())){
				if(kkey->GetClassName() == th1d){
					TH1D* hist = dynamic_cast<TH1D*>(kkey->ReadObj());
					if(!hist) continue;
					name = hist->GetName();
					if(name.find(proc) == string::npos) continue;
					if(hist->Integral() != 0 && name.find("sigma") == string::npos && name.find("mean") == string::npos) hist->Scale(1./hist->Integral());
					hists.push_back(hist);
				}
			}
		}
	}	
}

void GetHistsProcMethods(TDirectory* dir, string& proc, vector<TH1D*>& hists, vector<string>& methods){
	if(!dir) return;
	dir->cd();
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString tdir("TDirectoryFile");
	TString th1d("TH1D");
	hists.clear();
cout << "GETHISTSPROCMETHOD with dir " << dir->GetName() << " for proc " << proc << endl;
	while((key = (TKey*)iter())){
		cout << "key " << key->GetName() << " class " << key->GetClassName() << endl;
		if(key->GetClassName() == th1d){
			for(int m = 0; m < methods.size(); m++){
				cout << "m " << m << endl;
				//TDirectory* ddir = dynamic_cast<TDirectory*>(key->ReadObj());
				TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
				if(!hist) continue;
				cout << "hist not null" << endl;
				string name = hist->GetName();
				cout << "name " << name << " method " << methods[m] << endl;
				if(name.find(methods[m]) == string::npos) continue; 
				cout << "methods[m] " << methods[m] << endl;
				if(name.find(proc) == string::npos) continue;
				hists.push_back(hist);
				/*
				ddir->cd();
				cout << "--in dir " << ddir->GetName() << endl;
				TList* llist = ddir->GetListOfKeys();
				TIter iiter(llist);
				TKey* kkey;
				while((kkey = (TKey*)iiter())){
					if(kkey->GetClassName() == th1d){
						TH1D* hist = dynamic_cast<TH1D*>(kkey->ReadObj());
						if(!hist) continue;
						name = hist->GetName();
						if(name.find(proc) == string::npos) continue;
						if(hist->Integral() != 0 && name.find("sigma") == string::npos && name.find("mean") == string::npos) hist->Scale(1./hist->Integral());
						cout << "getting hist " << name << " for method " << methods[m] << endl;
						hists.push_back(hist);
					}
				}
				*/
			}
		}
	}	
}


void GetHistsProcs(TDirectory* dir, vector<string>& procs, vector<TH1D*>& hists, bool scale = true){
	if(!dir) return;
	dir->cd();
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString tdir("TDirectoryFile");
	TString th1d("TH1D");
	hists.clear();
	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d){
			TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
			if(!hist) continue;
			name = hist->GetName();
			for(int s = 0; s < procs.size(); s++){
				if(name.find(procs[s]) != string::npos){
					if(hist->Integral() == 0) continue;
					if(name.find("sigma") == string::npos && name.find("mean") == string::npos && scale) hist->Scale(1./hist->Integral());
					if(!scale) hist->GetYaxis()->SetTitle("events");
					string title = hist->GetTitle();
					if(title.find(procs[s]) == string::npos) hist->SetTitle((title+"_"+procs[s]).c_str());
					hists.push_back(hist);

				}
				else continue;
			}
		}
	}
}


void GetHistsType(TDirectory* dir, string type, vector<TH1D*>& hists){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th1d("TH1D");
	hists.clear();

	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d){
			//get 1D histograms
			TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
			if(hist){
				name = hist->GetName();
				if(type.empty()){
					if(name.find("lead") != string::npos) continue;
				}
				else{
					if(name.find("_"+type) == string::npos) continue;
				}
				if(isnan(hist->Integral())) continue;
				if(hist->GetEntries() == 0) continue;
				if(hist->Integral() != 0 && name.find("sigma") == string::npos) hist->Scale(1./hist->Integral());
				hists.push_back(hist);
			}
		}
	}
}

void GetHists(TDirectory* dir, string proc, vector<TH1D*>& hists){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th1d("TH1D");
	hists.clear();

	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d){
			//get 2D histograms
			TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
			if(hist){
				name = hist->GetName();
				if(name.find(proc) == string::npos) continue;
				if(isnan(hist->Integral())) continue;
				if(hist->Integral() != 0) hist->Scale(1./hist->Integral());
				else continue;
				hists.push_back(hist);
			}
		}
	}
}
void GetHists(TDirectory* dir, string proc, vector<TH2D*>& hists){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th2d("TH2D");
	hists.clear();

	while((key = (TKey*)iter())){
		if(key->GetClassName() == th2d){
			//get 2D histograms
			TH2D* hist = dynamic_cast<TH2D*>(key->ReadObj());
			if(hist){
				name = hist->GetName();
				if(name.find(proc) == string::npos) continue;
				if(isnan(hist->Integral())) continue;
				hists.push_back(hist);
			}
		}
	}
}
void GetProcs(TDirectory* dir, vector<string>& procs){
	vector<TH1D*> hists;
	procs.clear();
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString tdir("TDirectoryFile");
	while((key = (TKey*)iter())){
		if(key->GetClassName() == tdir){
			TDirectory* ddir = dynamic_cast<TDirectory*>(key->ReadObj());
			GetHists(ddir, "", hists);
			for(int i = 0; i < hists.size(); i++){
				name = hists[i]->GetTitle();
				name = name.substr(name.find("_")+1);
				//check if name in vector already
				if(find(procs.begin(), procs.end(), name) == procs.end()) 
					procs.push_back(name);
			}
		}
	}

}


string GetExtraLabel(string in_file){
	int idx = in_file.find("Skim_")+4;
	int idx2 = in_file.find("_emAlpha");
	if(idx == string::npos || idx2 == string::npos) return "";
	//no extra string in between
	if(idx == idx2) return "";
	else return in_file.substr(idx+1,idx2-idx-1);

}

string GetCMSLabel(string in_file){
	//get version
	string cmslab;
	string year = "";
	if(in_file.find("condor") != string::npos){
		//get from v[0-9] to cm
		//get version
		std::smatch m;
		std::regex re("_v[0-9]+_");
		string version = "";
		std::regex_search(in_file,m,re);
		for(auto x : m) version += x;
		int vidx = in_file.find(version)+version.size();
		int cmidx = in_file.rfind("_AOD");
		cmslab = in_file.substr(vidx,cmidx-vidx);
		m = smatch();
		re = std::regex("_R[0-9]+");
		std::regex_search(in_file,m,re);
		for(auto x : m) year += x;
	}
	else{
		int idx = in_file.find("NperGeV0p");
		//based on NperGeV0pXXX being the last label before cms label
		//and 3 digits following "GeV0p"
		cmslab = in_file.substr(idx+13);
		cmslab = cmslab.substr(0,cmslab.find("_output"));
	}
	cmslab = cmslab.substr(0,cmslab.find(".root"));
	cmslab += year;	
	//remove directory prefixes
	int cnt = std::count(cmslab.begin(), cmslab.end(), '/');
	for(int i = 0; i < cnt; i++)
		cmslab = cmslab.substr(cmslab.find("/")+1);
	if(cmslab.find("condor_") != string::npos)
		cmslab = cmslab.substr(cmslab.find("condor_")+7);
	return cmslab;
	
}



bool HistCheck(TDirectory* dir){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	TString th1d("TH1D");
	TString th2d("TH2D");
	TString tdir("TDirectoryFile");
	
	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d)
			return true;
		else if(key->GetClassName() == th2d)
			return true;
		else if(key->GetClassName() == tdir)
			return false;
	}
	return false;


}



void ProcStackHists(string file, vector<string>& procs, string oname, string match="",bool b_ratio = false){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"READ");
	TFile* ofile = TFile::Open(oname.c_str(),"UPDATE");
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name, xtitle, ytitle;
	TString th1d("TH1D");
	TString th2d("TH2D");
	TString tdir("TDirectoryFile");

	string procsname = "";
	for(auto s : procs) procsname += "_"+s;

	double ymin, ymax;
	string ylab, xlab;
	while((key = (TKey*)iter())){
		name = key->GetName();
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = dynamic_cast<TDirectory*>(key->ReadObj());
			if(!dir) continue;
			name = dir->GetName();
			if(!match.empty() && name.find(match) == string::npos) continue;
			dir->cd();
			//each dir is full of histograms (stack by process)
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name = dir->GetName();
			vector<TH1D*> hists;
			cout << "\ndir name: " << name << endl;
			GetHistsProcs(dir, procs, hists,false);
			if(hists.size() > 1){
				cout << hists.size() << " hists got" << endl;
				for(auto h : hists){
					cout << h->GetName() << endl;
				}
				FindListHistBounds(hists, ymin, ymax);
				if(ymin == 0 && ymax == 0) continue;
				name = dir->GetName();
				string cvname = name+procsname+"_procStack";
				TCanvas *cv = new TCanvas(cvname.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				if(name.find("sigma") != string::npos || name.find("mean") != string::npos){
					continue;
				}
				else{
					xlab = name;
					ylab = "events"; 
				}
				TDRMultiHist(hists, cv, name, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", b_ratio);
				cout << "writing canvas (1D) " << cv->GetName() << endl;
				cv->Write(); 
			}
		
		}
	}
	ofile->Close();
	f->Close();

};


void Hist2D(string file, string proc, string method, string oname, string match){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"READ");
	TFile* ofile = TFile::Open(oname.c_str(),"UPDATE");
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name, xtitle, ytitle;
	//string oname = f->GetName();
	//oname = oname.substr(0,oname.find(".root"));
	//oname = oname+"_formatted.root";
	//TFile* ofile = new TFile(oname.c_str(),"RECREATE");

	string cmslab = "";
	if(proc == "QCD"){
		cmslab = "QCD Multijets, 2017";
	}
	else if(proc == "ttbar"){
		cmslab = "t#bar{t}";
	}
	else cmslab = "process";
	cmslab += ", "+method;
	//string cmslab = GetCMSLabel(file);
	//string extra = "";
	//if(file.find("Skim") != string::npos) extra = GetExtraLabel(file);
	//if(!extra.empty()) cmslab += " "+extra;	


	TString th2d("TH2D");
	TString tdir("TDirectoryFile");

	string dirname, ddirname, histname, xlab, ylab;

	while((key = (TKey*)iter())){
		name = key->GetName();
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = dynamic_cast<TDirectory*>(key->ReadObj());
			if(!dir) continue;
			name = dir->GetName();
			if(!match.empty() && name.find(match) == string::npos) continue;
			dir->cd();
			//if directory is full of directories (stack by method)
			//each dir is full of histograms (stack by process)
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name = dir->GetName();
			vector<TH2D*> hists;
			TH2D* hist;
			cout << "\ndir name: " << dir->GetName() << endl;

			//method in this subdir needs to match what's given
			if(name.find(method) == string::npos) continue;
			cout << "method " << method << endl;
			//we're in the directory with hists of one method split by procs
			GetHists(dir, proc, hists);
			if(hists.size() > 0){
				cout << "hists got" << endl;
				for(auto h : hists) cout << h->GetName() << endl;
				hist = hists[0];
				name = hist->GetName();
				string cvname = name;
				TCanvas *cv = new TCanvas(cvname.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				xlab = hist->GetXaxis()->GetTitle();
				ylab = hist->GetYaxis()->GetTitle();
				if(xlab.empty() && ylab.empty()){
					xlab = name.substr(0,name.find("_"));
					ylab = name.substr(xlab.size()+1,name.find("_")-1);
				}
				cout << "x " << xlab << " y " << ylab << " name " << name << endl;
				TDR2DHist(hist, cv, xlab, ylab, cmslab, "");
				cout << "writing canvas (2D) " << cv->GetName() << endl;
				cv->Write(); 
			}
		
		}
	}

	ofile->Write();
	ofile->Close();
	f->Close();



}



void HistFormat_DataMCComp(string file){
	string oname = file;
	oname = oname.substr(0,oname.find(".root"));

	oname = oname+"_formatted.root";
	TFile* ofile = new TFile(oname.c_str(),"RECREATE");
	ofile->Close();	

	vector<string> procs = {"DoubleEG","GJets"};
	//decay types

	ProcStackHists(file, procs, oname, "IsoBkgSel",true);
	

	cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;

}



