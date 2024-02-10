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


void TDRMultiHist(vector<TH1D*> hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, string cms_label){
	if(hist.size() == 0)return;
	can->cd();
	can->SetGridx(1);
	can->SetGridy(1);
	can->SetTitle("");
	TLegend* myleg = new TLegend(0.7, 0.7, 0.9, 0.9);
	myleg->SetFillColor(0);
	myleg->SetBorderSize(0);
	myleg->SetTextFont(42);
	myleg->SetTextSize(0.04);

	//offset for log scale
	if(miny == 0) miny += 1e-6;

	string name;
	string canname = can->GetName();

	string legentry;
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
	labelToColor["!chiGam"] =  TColor::GetColor("#86bbd8");
	labelToColor["!GMSB"] =  TColor::GetColor("#86bbd8");
	labelToColor["!notSunm"] = TColor::GetColor("#9e0059");
	labelToColor["!GJets"] =   TColor::GetColor("#f6ae2d");
	labelToColor["!JetHT"] =   TColor::GetColor("#3d348b");

	labelToColor["median"] = TColor::GetColor("#f7a278");
	labelToColor["eAvg"] = TColor::GetColor("#6859f1");
	labelToColor["mmAvg"] = TColor::GetColor("#52b788");

	labelToMark["chiGam"] =  71;
	labelToMark["GMSB"] =  71;
	labelToMark["notSunm"] = 72;
	labelToMark["GJets"] =   73;
	labelToMark["JetHT"] =   74;

	labelToMark["!median"] = 71;
	labelToMark["!eAvg"] =   72; 
	labelToMark["!mmAvg"] =  73;

	int col, mark;	
	for( int i = 0 ; i < int(hist.size()); i++){
		hist[i]->UseCurrentStyle();
		hist[i]->SetStats(false);
		hist[i]->GetXaxis()->CenterTitle(true);
		hist[i]->GetXaxis()->SetTitle(xtit.c_str());
		hist[i]->GetYaxis()->CenterTitle(true);
		hist[i]->GetYaxis()->SetTitle(ytit.c_str());
		hist[i]->GetYaxis()->SetRangeUser(miny, maxy + maxy/10.);
		

		legentry = hist[i]->GetTitle();
		if(canname.find("jet") != string::npos){
			if(legentry.find("notSunm") != string::npos) continue;
			if(legentry.find("chiGam") != string::npos){
				legentry.replace(legentry.find("chiGam"),6,"GMSB");
			}
		}		
		

		//if a key from labeltocolor is in legentry, set that color
		for(map<string, int>::iterator it = labelToColor.begin(); it != labelToColor.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(legentry.find(match) != string::npos){
				col = it->second;
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
		hist[i]->SetLineColor(col);
		//hist[i]->SetLineWidth(2);
		hist[i]->SetMarkerStyle(mark);
		hist[i]->SetMarkerColor(col);
		//hist[i]->SetMarkerSize(1);
		hist[i]->SetTitle("");
		if( i == 0 ){
			hist[i]->Draw("ep");
		}else{
			hist[i]->Draw("epsame");
		}
		myleg->AddEntry( hist[i], legentry.c_str(), "p" );
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
	return;
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
	if(hist->GetNbinsX() == 2 && hist->GetNbinsY() == 2){
		hist->SetMarkerSize(3.);
		hist->Draw("colztext");
	}
	else hist->Draw("colz");

	string name = hist->GetName();
	//if(name.find("genDeltaTime_meanRecoGenDeltaT") != string::npos) cout << "n entries: " << hist->GetEntries() << endl;
	
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label+" "+title;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
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
	
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
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
			TDirectory* ddir = (TDirectory*)key->ReadObj();
			if(!ddir) continue;
			ddir->cd();
			TList* llist = ddir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name;
			while((kkey = (TKey*)iiter())){
				if(kkey->GetClassName() == th1d){
					TH1D* hist = (TH1D*)kkey->ReadObj();
					if(!hist) continue;
					name = hist->GetName();
					if(name.find(proc) == string::npos) continue;
					hists.push_back(hist);
				}
			}
		}
	}	



}


void GetHists(TDirectory* dir, string type, vector<TH1D*>& hists){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th1d("TH1D");
	hists.clear();

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
				if(isnan(hist->Integral())) continue;
				if(hist->GetEntries() == 0) continue;
				if(hist->Integral() != 0 && name.find("sigma") == string::npos) hist->Scale(1./hist->Integral());
				hists.push_back(hist);
			}
		}
	}
}

void GetHists(TDirectory* dir, string type, vector<TH2D*>& hists){
	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th2d("TH2D");
	hists.clear();

	while((key = (TKey*)iter())){
		if(key->GetClassName() == th2d){
			//get 1D histograms
			TH2D* hist = (TH2D*)key->ReadObj();
			if(hist){
				name = hist->GetName();
				if(type.empty()){
					if(name.find("lead") != string::npos) continue;
				}
				else{
					if(name.find("_"+type) == string::npos) continue;
				}
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
			TDirectory* ddir = (TDirectory*)key->ReadObj();
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
	if(in_file.find("output") == string::npos){
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
	}
	else{
		int idx = in_file.find("NperGeV0p");
		//based on NperGeV0pXXX being the last label before cms label
		//and 3 digits following "GeV0p"
		cmslab = in_file.substr(idx+13);
		cmslab = cmslab.substr(0,cmslab.find("_output"));
	}
	cmslab = cmslab.substr(0,cmslab.find(".root"));
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

void HistFormat(string file){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"UPDATE");
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
	string extra = "";
	if(file.find("Skim") != string::npos) extra = GetExtraLabel(file);
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
			string keyname = key->GetName();
			TH2D* hist = (TH2D*)key->ReadObj();
			if(hist){
				name = hist->GetName();
				xtitle = hist->GetXaxis()->GetTitle();
				ytitle = hist->GetYaxis()->GetTitle();
				TCanvas *cv = new TCanvas(name.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				TDR2DHist(hist, cv, xtitle, ytitle, cmslab, hist->GetTitle());
				cv->Write(); 
			}
		}
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = (TDirectory*)key->ReadObj();
			double ymin, ymax;
			string ylab, xlab;
			if(!dir) continue;
			dir->cd();
			//cout << "dir name: " << dir->GetName() << endl;
			//get histograms (stack these)
			vector<TH1D*> hists;
			//loop through types for 1D hists
			for(int i = 0; i < types.size(); i++){
				//loop through hists in dir
				name = dir->GetName();
				if(!types[i].empty()) name += "_"+types[i];
				GetHists(dir, types[i], hists);
				if(hists.size() < 1) continue;
				FindListHistBounds(hists, ymin, ymax);
				if(ymin == 0 && ymax == 0) continue;
				TCanvas *cv = new TCanvas(name.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				if(name.find("sigma") != string::npos || name.find("mean") != string::npos){
					xlab = hists[0]->GetXaxis()->GetTitle();
					ylab = hists[0]->GetYaxis()->GetTitle();
				}
				else{
					xlab = name;
					ylab = "a.u."; 
				}
				TDRMultiHist(hists, cv, name, xlab, ylab, ymin, ymax, cmslab);
				cv->Write(); 
			}
			//for 2D hists
			vector<TH2D*> hists2D;
			//loop through hists in dir
			//if(!types[i].empty()) name += "_"+types[i];
			GetHists(dir, "", hists2D);
			string title;
			if(hists2D.size() > 0){
				for(int h = 0; h < hists2D.size(); h++){
					name = hists2D[h]->GetName();
					title = hists2D[h]->GetTitle();
					TCanvas *cv = new TCanvas(name.c_str(), "");
					ofile->cd();
					//draw as tcanvases
					xlab = hists2D[h]->GetXaxis()->GetTitle();
					ylab = hists2D[h]->GetYaxis()->GetTitle();
					TDR2DHist(hists2D[h], cv, xlab, ylab, cmslab, title);
					cv->Write(); 
				}
			}
			//if directory is full of directories (stack by method)
			//each dir is full of histograms (stack by process)
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name;
			while((kkey = (TKey*)iiter())){
				if(kkey->GetClassName() == th2d){
					TH2D* hist = (TH2D*)kkey->ReadObj();
					if(!hist) continue;	
					//cout << "dir: " << dir->GetName() << " getting kkey: " << kkey->GetName() << " hist " << hist->GetName() << " with entries " << hist->GetEntries() << endl;
				}
				//writing stack hist - same method different procs
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = (TDirectory*)kkey->ReadObj();
					if(!ddir) continue;
					name = ddir->GetName();
					//cout << " dir name: " << name << endl;
					ddir->cd();
					//we're in the directory with hists of one method split by procs
					GetHists(ddir, "", hists);
					if(hists.size() > 0){
						FindListHistBounds(hists, ymin, ymax);
						if(ymin == 0 && ymax == 0) continue;
						TCanvas *cv = new TCanvas(name.c_str(), "");
						ofile->cd();
						//draw as tcanvases
						if(name.find("sigma") != string::npos || name.find("mean") != string::npos){
							xlab = hists[0]->GetXaxis()->GetTitle();
							ylab = hists[0]->GetYaxis()->GetTitle();
						}
						else{
							xlab = name;
							ylab = "a.u."; 
						}
						TDRMultiHist(hists, cv, name, xlab, ylab, ymin, ymax, cmslab);
						cv->Write(); 
					}
					vector<TH2D*> hhists2D;
					//loop through hists in dir
					//if(!types[i].empty()) name += "_"+types[i];
					GetHists(ddir, "", hhists2D);
					if(hhists2D.size() > 0){
						for(int h = 0; h < hhists2D.size(); h++){
							if(hhists2D[h]->GetEntries() == 0) continue;
							name = hhists2D[h]->GetName();
							title = hhists2D[h]->GetTitle();
							TCanvas *cv = new TCanvas(name.c_str(), "");
							ofile->cd();
							//draw as tcanvases
							xlab = hhists2D[h]->GetXaxis()->GetTitle();
							ylab = hhists2D[h]->GetYaxis()->GetTitle();
							TDR2DHist(hhists2D[h], cv, xlab, ylab, cmslab, title);
							cv->Write(); 
						}
					}
				}

			}
			//writing methodStack hist - same proc different methods
			string dirname = dir->GetName();
			//only do for sigma plots for now - can remove this later to change
			if(dirname.find("sigma") == string::npos && dirname.find("mean") == string::npos) continue;	
			vector<string> procs;
			GetProcs(dir, procs);
			for(int p = 0; p < procs.size(); p++){
				GetHistsProc(dir, procs[p], hists);
				name = dirname+"_"+procs[p]+"_methodStack";
				if(name.find("jet") != string::npos && name.find("notSunm") != string::npos) continue;
				
				if(hists.size() < 1) continue;
				FindListHistBounds(hists, ymin, ymax);
				if(ymin == 0 && ymax == 0) continue;
				TCanvas *cv = new TCanvas(name.c_str(), "");
				ofile->cd();
				//draw as tcanvases
				if(name.find("sigma") != string::npos || name.find("mean") != string::npos){
					xlab = hists[0]->GetXaxis()->GetTitle();
					ylab = hists[0]->GetYaxis()->GetTitle();
				}
				else{
					xlab = name;
					ylab = "a.u."; 
				}
				TDRMultiHist(hists, cv, name, xlab, ylab, ymin, ymax, cmslab);
				cv->Write(); 
				
			}
		
		}
	}
	cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;

	ofile->Close();
	f->Close();

};
