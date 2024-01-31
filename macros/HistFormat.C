#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>

using std::string;

void clearptrvec(vector<TH1D*>& vec){
	for(int i = 0; i < vec.size(); i++)
		delete vec[i];
	vec.clear();
}


void Profile2DHist(TH2D* inhist, TH1D* outhist, vector<TH1D*>& profs){
	int nbins = inhist->GetNbinsX();
	profs.clear();
	string profilename = "";
	string profiletitle = "";
	//skip overflow + underflow bins
	for(int i = 1; i < nbins; i++){
		TH1D* phist = (TH1D*)inhist->ProjectionY("tmp",i,i);
		profilename = inhist->GetName();
		profiletitle = inhist->GetTitle();
		profilename.insert(profilename.find("_"+profiletitle),"_bin"+std::to_string(i));
		phist->SetTitle(profiletitle.c_str());	
		phist->SetName(("profile_"+profilename).c_str());
		phist->GetXaxis()->SetTitle(inhist->GetXaxis()->GetTitle());	
		profs.push_back(phist);

		outhist->GetXaxis()->SetTitle(inhist->GetXaxis()->GetTitle());
		string ytitle = inhist->GetYaxis()->GetTitle();
		ytitle = "#sigma "+ytitle;
		outhist->GetYaxis()->SetTitle(ytitle.c_str());

		//get values for param init
		double mean = phist->GetMean();
		double stddev = phist->GetStdDev();
		double norm = phist->Integral();
		double low = phist->GetBinLowEdge(0);
		double high = -low;
		//check that initial parameter values are ok
		if( stddev >= 0.0 && norm > 0.){
			TF1* fit = new TF1("fit","gaus",low,high);
			//fit->SetParameter(0,norm);
			//fit->SetParameter(1,mean);
			//fit->SetParameter(2,stddev);
			//phist->Fit(fit->GetName(),"RBQ0");
			phist->Fit(fit->GetName(),"RQ0");
			
			double fit_stddev = fit->GetParameter(2);
			double fit_stddev_err = fit->GetParError(2);
			//set new contents
			outhist->SetBinContent(i, fit_stddev);
			outhist->SetBinError(i, fit_stddev_err);
			delete fit;
		}
	}
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
	if(hist->GetNbinsX() == 2 && hist->GetNbinsY() == 2)
		hist->Draw("colztext");
	else hist->Draw("colz");
	
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

void Profile2DHists(TFile* f){
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;

	TString th2d("TH2D");
	TString tdir("TDirectoryFile");
	string name, newname, histname, newhistname, histtitle, newhisttitle;
	string dirname, ddirname, dirname_new, ddirname_new;

	string filename = f->GetName();
	bool print = (filename.find("GMSB") != string::npos);

	if(filename.find("photons") != string::npos) return;
	f->cd();
	//loop over all dirs - top level is variable	
	while((key = (TKey*)iter())){
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = (TDirectory*)key->ReadObj();
			if(!dir) continue;
			dirname = dir->GetName();
			dir->cd();
			if(dirname.find("diffDeltaTime") != string::npos){
				dirname_new = dirname.replace(dirname.find("diff"), 4, "sigma");
			}
			else{
				continue;
			}
			//get total hist in this directory
			vector<TH2D*> hists;
			vector<TH1D*> profs;
			GetHists(dir, "", hists);
			if(hists.size() < 1) continue;
			//if(print) cout << "dir name: " << dirname << " has " << hists.size() << " hists" << endl;
			//write only total 1D profile from total 2D (doesn't end in process) to newname
			for(int h = 0; h < hists.size(); h++){
				histname = hists[h]->GetName();
				histname.replace(histname.find("diff"), 4, "sigma");
				//cout << " getting hist " << dirname_new << "/" << histname << endl;
				TH1D* outhist = (TH1D*)f->Get((dirname_new+"/"+histname).c_str()); 
				//if(print) cout << outhist->GetEntries() << " entries in " << outhist->GetName() << endl;
				if(outhist->GetEntries() > 1) continue;
				//cout << " writing sigma hist to " << dirname_new+"/"+outhist->GetName() << endl;
				//profile 1D histograms
				Profile2DHist(hists[h], outhist, profs);
				//if(print) cout << " inhist " << hists[h]->GetName() << " has " << hists[h]->GetEntries() << " entries and outhist " << histname << " has entries " << outhist->GetEntries() << endl;
				for(int b = 0; b < profs.size(); b++){
					f->cd();
					//write total profiles to 
					string profilename = profs[b]->GetName();
					string profilepath = profilename;
					string profiletitle = profs[b]->GetTitle();
					profilepath = profilepath.substr(0,profilepath.find("_"+profiletitle));
					profilepath = profilepath+"_stack";
					//if method in profile name - add next dir to profile path
					TH1D* prof = (TH1D*)f->Get((profilepath+"/"+profilename).c_str()); 
					TDirectory* profdir = (TDirectory*)f->Get(profilepath.c_str());
					profdir->cd();
					if(prof == nullptr){ cout << "profile null" << endl; continue; }
					//cout << "og entries " << prof->GetEntries() << endl;
					//if(prof->GetEntries() > 1) continue;
					prof = (TH1D*)profs[b]->Clone();
					//cout << "new entries " << prof->GetEntries() << endl;
				}
				for(int b = 0; b < profs.size(); b++)
					delete profs[b];
			}
			f->cd();
			dir->cd();
			//if directory is full of directories (stack by method)
			//each dir is full of histograms (stack by process)
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			//loop over process directories
			while((kkey = (TKey*)iiter())){
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = (TDirectory*)kkey->ReadObj();
					if(!ddir) continue;
					ddirname = ddir->GetName();
					ddirname_new = ddirname;
					ddir->cd();
					//we're in the directory with hists of one method split by procs
					//get 2D proc split hists
					vector<TH2D*> hhists;
					vector<TH1D*> pprofs;
					GetHists(ddir, "", hhists);
					if(hhists.size() < 1) continue;
					//cout << " dir name: " << ddirname << " has " << hhists.size() << " hists" << endl;
					ddirname_new.replace(ddirname_new.find("diff"), 4, "sigma");
					ddirname_new = dirname_new+"/"+ddirname_new;
					//want to profile the rest of the histograms
					for(int h = 0; h < hhists.size(); h++){
						//cout << " i: " << i << " hist: " << hists[i]->GetName() << " newname: " << newname << " newhistname: " << newhistname << endl;
						histname = hhists[h]->GetName();
						newhistname = histname;
						newhistname.replace(newhistname.find("diff"), 4, "sigma");
						//cout << "  getting hist " << ddirname_new+"/"+newhistname << endl;
						TH1D* outhist2 = (TH1D*)f->Get((ddirname_new+"/"+newhistname).c_str());
						if(outhist2 == nullptr){ continue; }
						//only fill if empty
						if(outhist2->GetEntries() > 1){ continue; }
						//cout << "  i: " << h << " inhist " << hists[h]->GetName() << " has " << hists[h]->GetEntries() << " entries and already has profiles, outhist " << newhistname << " has entries " << outhist2->GetEntries() << endl; continue; 	
						//profile 1D histograms
						Profile2DHist(hhists[h], outhist2, pprofs);
	//cout << "  i: " << h << " inhist has " << hists[h]->GetEntries() << " entries and " << profs.size() << " profiles, outhist " << newhistname << " has entries " << outhist2->GetEntries() << endl;
						for(int b = 0; b < pprofs.size(); b++){
							f->cd();
							//write by process profiles to 
							string profilename = pprofs[b]->GetName();
							string profilepath = profilename;
							string profiletitle = pprofs[b]->GetTitle();
							profilepath = profilepath.substr(0,profilepath.find("_"+profiletitle));
							profilepath = profilepath+"_stack/"+profilepath;
							profilepath = profilepath+"_"+profiletitle.substr(0,profiletitle.rfind("_")+1)+"procs";
							//if method in profile name - add next dir to profile path
							//cout << "Getting: " << profilepath+"/"+profilename << endl;
							TH1D* prof = (TH1D*)f->Get((profilepath+"/"+profilename).c_str());
							TDirectory* profdir1 = (TDirectory*)f->Get((profilepath.substr(0,profilepath.rfind("/"))).c_str());
							TDirectory* profdir2 = (TDirectory*)f->Get((profilepath).c_str());
							//if(print) cout << "Setting profile " << prof->GetName() << " in dir " << profdir1->GetName() << "/" << profdir2->GetName() << endl;
							profdir1->cd();
							profdir2->cd();
							if(prof == nullptr){ cout << "profile null" << endl; continue; }
							//if(print) cout << "og entries " << prof->GetEntries() << " new prof: " << pprofs[b]->GetEntries() << endl;
							//if(prof->GetEntries() > 1) continue;
							prof = (TH1D*)pprofs[b]->Clone();
							//if(print) cout << "new entries " << prof->GetEntries() << endl;
						}
						for(int b = 0; b < pprofs.size(); b++)
							delete pprofs[b];
					}
				}

			}
		}
	}
	f->Write("",TObject::kOverwrite);
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

	//profile relevant 2D histograms first
	Profile2DHists(f);

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
			//loop through types
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
				if(name.find("sigma") != string::npos){
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
			//if directory is full of directories (stack by method)
			//each dir is full of histograms (stack by process)
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name;
			while((kkey = (TKey*)iiter())){
				//writing stack hist - same method different procs
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = (TDirectory*)kkey->ReadObj();
					if(!ddir) continue;
					name = ddir->GetName();
					//cout << " dir name: " << name << endl;
					ddir->cd();
					//we're in the directory with hists of one method split by procs
					GetHists(ddir, "", hists);
					if(hists.size() < 1) continue;
					FindListHistBounds(hists, ymin, ymax);
					if(ymin == 0 && ymax == 0) continue;
					TCanvas *cv = new TCanvas(name.c_str(), "");
					ofile->cd();
					//draw as tcanvases
					if(name.find("sigma") != string::npos){
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
			//writing stack hist - same proc different methods
			string dirname = dir->GetName();
			//only do for sigma plots for now - can remove this later to change
			if(dirname.find("sigma") == string::npos) continue;	
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
				if(name.find("sigma") != string::npos){
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
