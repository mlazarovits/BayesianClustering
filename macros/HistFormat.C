#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>

using std::string;


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


void TDRMultiHist(vector<TH1D*> hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, string cms_label){
	if(can == nullptr) return;
	if(hist.size() == 0)return;
	can->cd();
	can->SetGridx(1);
	can->SetGridy(1);
	can->SetTitle("");
	TLegend* myleg = nullptr;
	if(cms_label.find("photon") != string::npos) myleg = new TLegend(0.6, 0.7, 0.9, 0.9);
	else{
		myleg = new TLegend(0.737,0.675,0.878,0.882);
		//myleg->SetMargin(0.1);
	}
	myleg->SetFillColor(0);
	myleg->SetBorderSize(0);
	myleg->SetTextFont(42);
	if(hist.size() > 4) myleg->SetTextSize(0.025);
	else myleg->SetTextSize(0.04);
	
	//offset for log scale
	if(miny == 0) miny += 1e-6;

	string title;
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
	labelToColor["chiGam"] =  TColor::GetColor("#86bbd8");
	labelToColor["GMSB"] =  TColor::GetColor("#86bbd8");
	//TODO (maybe): set different signal grid points to different shades of above color	

	labelToColor["notSunm"] =  TColor::GetColor("#9e0059");
	labelToColor["GJets"] 	=  TColor::GetColor("#f6ae2d");
	labelToColor["JetHT"] 	=  TColor::GetColor("#3d348b");
	labelToColor["MET"]   	=  TColor::GetColor("#671E76");
	labelToColor["ttbar"] 	=  TColor::GetColor("#CA5743");
	labelToColor["QCD"]   	=  TColor::GetColor("#9E0059");
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
		hist[i]->GetXaxis()->SetTitle(xtit.c_str());
		hist[i]->GetYaxis()->CenterTitle(true);
		hist[i]->GetYaxis()->SetTitle(ytit.c_str());
		hist[i]->GetYaxis()->SetRangeUser(miny, maxy + maxy/10.);
		


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
		gPad->Update();
	}
	title = title.substr(title.find("_")+1);
	if(title.find("chiGam") != string::npos){
		title = SignalLegEntry(title);
	}
	myleg->Draw("same"); 
	gPad->Update();
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.025);
	lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());
	TLatex lat1;
	lat1.SetNDC();
	lat1.SetTextSize(0.04);
	lat1.SetTextFont(42);
	lat1.DrawLatex(0.60,0.92,title.c_str());
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
	
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label+" "+title;
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
	
	string lat_cms = "#bf{CMS} #it{WIP} "+cms_label;
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
					if(hist->Integral() != 0 && name.find("sigma") == string::npos && name.find("mean") == string::npos) hist->Scale(1./hist->Integral());
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
			//get 2D histograms
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
	TString tgraph("TGraph");
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
				//cout << "writing hist " << cv->GetName() << " in file " << endl;
				cv->Write(); 
			}
		}
		if(key->GetClassName() == tgraph){
			//get TGraphs
			string keyname = key->GetName();
			TGraph* gr = (TGraph*)key->ReadObj();
			if(gr){
				name = gr->GetName();
				xtitle = gr->GetXaxis()->GetTitle();
				ytitle = gr->GetYaxis()->GetTitle();
				TCanvas *cv = new TCanvas(name.c_str(), "");
				ofile->cd();
				gr->SetMarkerStyle(20);
				//draw as tcanvases
				gr->Draw("AP");
				//cout << "writing gr " << cv->GetName() << " in file " << endl;
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
			//cout << "---in dir " << dir->GetName() << endl;
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
			//cout << "getting 2D hists" << endl;
			//for 2D hists
			vector<TH2D*> hists2D;
			//loop through hists in dir
			//if(!types[i].empty()) name += "_"+types[i];
			GetHists(dir, "", hists2D);
			string title;
			//cout << dir->GetName() << " has " << hists2D.size() << " hists" << endl;
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
					//cout << "writing canvas " << cv->GetName() << " in dir " << dir->GetName() << endl;
					cv->Write("",TObject::kOverwrite); 
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
					//cout << " getting kkey: " << kkey->GetName() << " hist " << hist->GetName() << " with entries " << hist->GetEntries() << endl;
				}
				//writing stack hist - same method different procs
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = (TDirectory*)kkey->ReadObj();
					if(!ddir) continue;
					name = ddir->GetName();
					//cout << " dir name: " << name << endl;
					ddir->cd();
					//cout << " ---in ddir " << ddir->GetName() << endl;
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
						//cout << "writing canvas (1D) " << cv->GetName() << endl;
						cv->Write(); 
					}
					vector<TH2D*> hhists2D;
					//loop through hists in dir
					//if(!types[i].empty()) name += "_"+types[i];
					GetHists(ddir, "", hhists2D);
					//cout << "got hists from dir " << ddir->GetName() << endl;
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
							//cout << "writing canvas " << cv->GetName() << endl;
							cv->Write(); 
						}
					}
				}

			}
			//writing methodStack hist - same proc different methods
			string dirname = dir->GetName();
			//only do for sigma plots + profiles for now - can remove this later to change
			if(dirname.find("sigma") == string::npos && dirname.find("mean") == string::npos && dirname.find("profile") == string::npos) continue;	
			vector<string> procs;
			GetProcs(dir, procs);
			for(int p = 0; p < procs.size(); p++){
				if(procs[p].find("notSunm") != string::npos) continue;
				GetHistsProc(dir, procs[p], hists);
				name = dirname+"_"+procs[p]+"_methodStack";
				
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








