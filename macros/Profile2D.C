#include <string>
#include <vector>

using std::string;
using std::vector;


void Profile2DHist(TH2D* inhist, TH1D* sig_outhist, TH1D* mean_outhist, vector<TH1D*>& profs){
	int nbins = inhist->GetNbinsX();
	profs.clear();
	string profilename = "";
	string profiletitle = "";
	//skip overflow + underflow bins
	cout << "\nPROFILING HIST " << inhist->GetName() << " has " << nbins << " bins" << " sighist bins " << sig_outhist->GetNbinsX() << endl;
	for(int i = 1; i < nbins+1; i++){
		TH1D* phist = (TH1D*)inhist->ProjectionY("tmp",i,i);
		if(!phist) continue;
		if(phist->GetEntries() == 0 || phist->Integral() == 0) continue;
		//cout << "bin #" << i << " nentries in profile " << phist->GetEntries() << " integral " << phist->Integral() << endl;
		profilename = inhist->GetName();
		profiletitle = inhist->GetTitle();
		profilename.insert(profilename.find("_"+profiletitle),"_bin"+std::to_string(i));
		phist->SetTitle(profiletitle.c_str());	
		phist->SetName(("profile_"+profilename).c_str());
		phist->GetXaxis()->SetTitle(inhist->GetXaxis()->GetTitle());
		profs.push_back(phist);

		sig_outhist->GetXaxis()->SetTitle(inhist->GetXaxis()->GetTitle());
		string ytitle = inhist->GetYaxis()->GetTitle();
		ytitle = "#sigma "+ytitle;
		sig_outhist->GetYaxis()->SetTitle(ytitle.c_str());
	
		if(mean_outhist){	
			mean_outhist->GetXaxis()->SetTitle(inhist->GetXaxis()->GetTitle());
			ytitle = inhist->GetYaxis()->GetTitle();
			ytitle = "#mu "+ytitle;
			mean_outhist->GetYaxis()->SetTitle(ytitle.c_str());
		}
		//get values for param init
		double mean = phist->GetMean();
		double stddev = phist->GetStdDev();
		double norm = phist->Integral();
		double low = phist->GetBinLowEdge(0);
		double high = -low;
		if(profilename.find("gamPV") != string::npos) high = phist->GetBinLowEdge(phist->GetNbinsX());
		cout << "norm " << norm << " mean " << mean << " stddev " << stddev << " low " << low << " high " << high << " " << phist->GetBinLowEdge(phist->GetNbinsX()-1)+phist->GetBinWidth(phist->GetNbinsX()-1) << endl;
		int ngoodbins = 0;
		for(int b = 0; b < phist->GetNbinsX(); b++){ if(phist->GetBinContent(b) > 0) ngoodbins++; }//cout << "bin #" << b << " error " << phist->GetBinError(b) << " content " << phist->GetBinContent(b) << endl;
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
			double fit_mean = fit->GetParameter(1);
			double fit_mean_err = fit->GetParError(1);
			//set new contents
			sig_outhist->SetBinContent(i, fit_stddev);
			sig_outhist->SetBinError(i, fit_stddev_err);
			cout << "bin #" << i << " sig " << fit_stddev << " err " << fit_stddev_err << " npts in fit " << fit->GetNumberFitPoints() << " nentries in profile " << phist->GetEntries() << " ngoodentries " << ngoodbins << endl;
			if(mean_outhist){
				//cout << "bin #" << i << " mean " << fit_mean << " nentries in profile " << phist->GetEntries() << endl;
				mean_outhist->SetBinContent(i, fit_mean);
				mean_outhist->SetBinError(i, fit_mean_err);
			}
			delete fit;
		}
	}
	//cout << "return" << endl;
}


void GetHists(TDirectory* dir, vector<TH1D*>& hists){
	dir->cd();

	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th1d("TH1D");
	hists.clear();
	TH1D* hist = nullptr;
	while((key = (TKey*)iter())){
		if(key->GetClassName() == th1d){
			name = key->GetName();
			hist = (TH1D*)dir->Get((name).c_str());
			if(!hist) continue;
			//cout << "getting hist " << dir->GetName() << "/" << hist->GetName() << " with entries " << hist->GetEntries() << endl;
			if(hist) hists.push_back(hist);
		}
	}
}

void GetHists(TDirectory* dir, vector<TH2D*>& hists){
	dir->cd();

	TList* list = dir->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name;
	TString th2d("TH2D");
	hists.clear();
	TH2D* hist = nullptr;
	while((key = (TKey*)iter())){
		if(key->GetClassName() == th2d){
			name = key->GetName();
			hist = (TH2D*)dir->Get((name).c_str());
			if(!hist) continue;
			if(hist) hists.push_back(hist);
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
	string sig_dirname, sig_ddirname, sig_dirname_new, sig_ddirname_new;
	string mean_dirname, mean_ddirname, mean_dirname_new, mean_ddirname_new;
	string dirname, ddirname;

	string filename = f->GetName();
	bool print = (filename.find("GMSB") != string::npos);

	//only do for jets
	if(filename.find("photons") != string::npos) return;
	f->cd();
	//loop over all dirs - top level is variable	
	while((key = (TKey*)iter())){
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = dynamic_cast<TDirectory*>(key->ReadObj());
			if(!dir) continue;
			dirname = dir->GetName();
			//cout << "In dir " << dirname << endl;
			dir->cd();
			if(dirname.find("diffDeltaTime") == string::npos){
				continue;
			}
			sig_dirname = dirname;
			mean_dirname = dirname;
			sig_dirname_new = sig_dirname.replace(sig_dirname.find("diff"), 4, "sigma");
			mean_dirname_new = mean_dirname.replace(mean_dirname.find("diff"), 4, "mean");
			//get total hist in this directory
			vector<TH2D*> hists;
			vector<TH1D*> profs;
			GetHists(dir, hists);
			//cout << "dir name: " << dirname << " has " << hists.size() << " hists" << endl;
			if(hists.size() < 1) continue;
			//write only total 1D profile from total 2D (doesn't end in process) to newname
			for(int h = 0; h < hists.size(); h++){
				histname = hists[h]->GetName();
				histname.replace(histname.find("diff"), 4, "sigma");
				TH1D* sig_outhist = (TH1D*)f->Get((sig_dirname_new+"/"+histname).c_str()); 
				if(sig_outhist->GetEntries() > 1) continue;
				histname = hists[h]->GetName();
				histname.replace(histname.find("diff"), 4, "mean");
				TH1D* mean_outhist = nullptr;
				if(dirname.find("recoGen") != string::npos) mean_outhist = (TH1D*)f->Get((mean_dirname_new+"/"+histname).c_str()); 
				if(mean_outhist && mean_outhist->GetEntries() > 1) continue;
				//if(print) cout << sig_outhist->GetEntries() << " entries in " << sig_outhist->GetName() << endl;
				//if(print && mean_outhist) cout << mean_outhist->GetEntries() << " entries in " << mean_outhist->GetName() << endl;
				//profile 1D histograms
				Profile2DHist(hists[h], sig_outhist, mean_outhist, profs);
				//cout << "1 - " << profs.size() << " profiles for " << hists[h]->GetNbinsX() << " bins" << endl;
				//if(print && mean_outhist) cout << " inhist " << hists[h]->GetName() << " has " << hists[h]->GetEntries() << " entries and mean_outhist " << histname << " has entries " <<  sig_outhist->GetEntries() << " mean: " << mean_outhist->GetEntries() << endl;
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
					if(prof == nullptr){ cout << "profile " << profilepath+"/"+profilename << " null" << endl; continue; }
					TDirectory* profdir = (TDirectory*)f->Get(profilepath.c_str());
					profdir->cd();
					//cout << "og entries " << prof->GetEntries() << " in " << prof->GetName() << endl;
					//if(prof->GetEntries() > 1) continue;
					prof = (TH1D*)profs[b]->Clone();
					//cout << "new entries " << prof->GetEntries() <<" in " << prof->GetName() <<  endl;
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
					sig_ddirname_new = ddirname;
					mean_ddirname_new = ddirname;
					ddir->cd();
					//we're in the directory with hists of one method split by procs
					//get 2D proc split hists
					vector<TH2D*> hhists;
					vector<TH1D*> pprofs;
					GetHists(ddir, hhists);
					if(hhists.size() < 1) continue;
				//	cout << " dir name: " << ddirname << " has " << hhists.size() << " hists" << endl;
					sig_ddirname_new.replace(sig_ddirname_new.find("diff"), 4, "sigma");
					sig_ddirname_new = sig_dirname_new+"/"+sig_ddirname_new;
					mean_ddirname_new.replace(mean_ddirname_new.find("diff"), 4, "mean");
					mean_ddirname_new = mean_dirname_new+"/"+mean_ddirname_new;
					//want to profile the rest of the histograms
					for(int h = 0; h < hhists.size(); h++){
						//cout << " h: " << h << " hist: " << hists[h]->GetName() << " newname: " << newname << " newhistname: " << newhistname << endl;
						histname = hhists[h]->GetName();
						newhistname = histname;
						newhistname.replace(newhistname.find("diff"), 4, "sigma");
						//cout << "  getting hist " << ddirname_new+"/"+newhistname << endl;
						TH1D* sig_outhist2 = (TH1D*)f->Get((sig_ddirname_new+"/"+newhistname).c_str());
						newhistname = histname;
						newhistname.replace(newhistname.find("diff"), 4, "mean");
						TH1D* mean_outhist2 = nullptr;
						if(ddirname.find("recoGen") != string::npos) mean_outhist2 = (TH1D*)f->Get((mean_ddirname_new+"/"+newhistname).c_str());
						if(sig_outhist2 == nullptr){ continue; }
						//only fill if empty
						if(sig_outhist2->GetEntries() > 1){ continue; }
						if(mean_outhist2 && mean_outhist2->GetEntries() > 1){ continue; }
						//cout << "  i: " << h << " inhist " << hists[h]->GetName() << " has " << hists[h]->GetEntries() << " entries and already has profiles, outhist " << newhistname << endl; 	
						//profile 1D histograms
						Profile2DHist(hhists[h], sig_outhist2, mean_outhist2, pprofs);
						//cout << "2 - " << pprofs.size() << " profiles for " << hhists[h]->GetNbinsX() << " bins" << endl;
						for(int b = 0; b < pprofs.size(); b++){
							//cout << "pprof #" << b << endl;
							f->cd();
							//write by process profiles to 
							string profilename = pprofs[b]->GetName();
							string profilepath = profilename;
							string profiletitle = pprofs[b]->GetTitle();
							//cout << "profiletitle " << profiletitle << endl;
							//cout << "profilepath og " << profilepath << endl;
							profilepath = profilepath.substr(0,profilepath.find("_"+profiletitle));
							//cout << "profilepath 1  " << profilepath << endl;
							profilepath = profilepath+"_stack/"+profilepath;
							//cout << "profilepath 2  " << profilepath << endl;
							//profilepath = profilepath+"_"+profiletitle.substr(0,profiletitle.rfind("_")+1)+"procStack";
							profilepath += ddirname.substr(ddirname.find_last_of("_",ddirname.find("_procStack")-1));
							//cout << "profilepath 3  " << profilepath << endl;
							//cout << "ddirname " << ddirname << endl;
							//if method in profile name - add next dir to profile path
							//cout << "Getting: " << profilepath << " / " << profilename << endl;
							TH1D* prof = (TH1D*)f->Get((profilepath+"/"+profilename).c_str());
							if(prof == nullptr){ cout << "profile " << profilepath+"/"+profilename << " null" << endl; continue; }
							TDirectory* profdir1 = (TDirectory*)f->Get((profilepath.substr(0,profilepath.rfind("/"))).c_str());
							TDirectory* profdir2 = (TDirectory*)f->Get((profilepath).c_str());
							//cout << "Setting profile " << prof->GetName() << " in dir " << profdir1->GetName() << "/" << profdir2->GetName() << endl;
							profdir1->cd();
							profdir2->cd();
							//cout << "og entries " << prof->GetEntries() << " new prof: " << pprofs[b]->GetEntries() << endl;
							//if(prof->GetEntries() > 1) continue;
							prof = (TH1D*)pprofs[b]->Clone();
							//cout << "new entries " << prof->GetEntries() << endl;
						}
						for(int b = 0; b < pprofs.size(); b++){
							//cout << " 1 - delete prof #" << b << endl;
							delete pprofs[b];
							//cout << " 2 - delete prof #" << b << endl;
						}
					}
				}

			}
		}
	}
	//cout << "writing" << endl;
	f->Write("",TObject::kOverwrite);
}





double GetGausMean(TH1D* inhist){
	//cout << "inhist " << inhist->GetName() << " entries " << inhist->GetEntries() << endl;
	//get values for param init
	double mean = inhist->GetMean();
	double stddev = inhist->GetStdDev();
	double norm = inhist->Integral();
	double low = inhist->GetBinLowEdge(0);
	double high = -low;
	//check that initial parameter values are ok
	if( stddev >= 0.0 && norm > 0.){
		TF1* fit = new TF1("fit","gaus",low,high);
		//fit->SetParameter(0,norm);
		//fit->SetParameter(1,mean);
		//fit->SetParameter(2,stddev);
		//inhist->Fit(fit->GetName(),"RBQ0");
		inhist->Fit(fit->GetName(),"RQ0");
		return fit->GetParameter(1);
	}
	else return -999;	
}


void Make2DHist(TFile* f, string histdirname){
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;

	TString th2d("TH2D");
	TString th1d("TH1D");
	TString tdir("TDirectoryFile");
	string name, newname, histname, newhistname, histtitle, newhisttitle;
	string dirname, ddirname;

	string filename = f->GetName();
	bool print = (filename.find("GMSB") != string::npos);


	//get histograms
	TDirectory* stackdir = (TDirectory*)f->Get((histdirname+"_stack").c_str());
	vector<TH2D*> outhists;
	//gDirectory->pwd();
	GetHists(stackdir,outhists);
	//cout << "got " << outhists.size() << " hists for " << stackdir->GetName() << endl;
	return;
	//outhist_procStack[i][j] - i: tr method, j: process
	vector<vector<TH2D*>> outhists_procStack;
	TList* liststack = stackdir->GetListOfKeys();
	TIter iterstack(liststack);
	TKey* keystack;
	TDirectory* procstackdir = nullptr;
	string keyname;
	vector<TDirectory*> procstackdirs;
	while((keystack = (TKey*)iterstack())){
		if(keystack->GetClassName() == tdir){
			keyname = keystack->GetName();
			procstackdir = (TDirectory*)f->Get((histdirname+"_stack/"+keyname).c_str());
			if(!procstackdir) continue;
			procstackdirs.push_back(procstackdir);
			outhists_procStack.push_back({});
			GetHists(procstackdirs[procstackdirs.size()-1],outhists_procStack[outhists_procStack.size()-1]);
		}

	}
	//cout << outhists.size() << " # outhists" << endl;
//cout << "A" << endl;
	int xbin, ybin;
	string ybinname;
	const char* match = "bin";
	TDirectory* dir = nullptr;
	while((key = (TKey*)iter())){
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			dir = (TDirectory*)key->ReadObj();
			if(!dir) continue;
			dirname = dir->GetName();
			//grab only profile hists
			if(dirname.find("profile") == string::npos) continue;
			//needs to be geoEavg profiles to match 2D hist
			if(dirname.find("geoEavg") == string::npos) continue;
			//needs to be binned profile - 2D = 2 "bin"s
			if(std::count(dirname.begin(), dirname.end(), *match) != 2) continue;
			dir->cd();
			xbin = std::stod(dirname.substr(dirname.find("_bin")+4,dirname.find("_stack")-dirname.find("_bin")-4));
			ybinname = dirname.substr(0,dirname.find("_bin"+to_string(xbin)));
			ybin = std::stod(ybinname.substr(ybinname.find("bin")+3));
			//should be in bin (alphabetical) order
			//do process integrated 2D plots
			vector<TH1D*> hists;
			//get profiles - for all time reco methods
			GetHists(dir, hists);
			//cout << "dirname: " << dirname << " xbin: " << xbin << " ybin: " << ybin << " dir name: " << dirname << " has " << hists.size() << " hists" << endl;
			if(hists.size() < 1) continue;
			//loop through profiles
			double mean;
			string trMethod;
			for(int h = 0; h < hists.size(); h++){
				mean = GetGausMean(hists[h]);
				if(mean == -999) continue;	
				//get TH2D for time reco h
				//cout << "h " << h << " hist 1D " << hists[h]->GetName() << " mean: " << mean << " hist 2D " << outhists[h]->GetName() << endl;	
				if(!outhists[h]){ cout << "1 - outhist " << h <<  " null" << endl; continue; }
				//fill TH2D for bin ebin, genbin with mean
				outhists[h]->SetBinContent(xbin, ybin, mean);
				//cout << "xbin: " << xbin << " ybin: " << ybin << " mean: " << mean << " outhist " << outhists[h]->GetName() << " current # entries " << outhists[h]->GetEntries() << endl;
			}
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			int itr = 0;
			//loop over process directories
			while((kkey = (TKey*)iiter())){
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = (TDirectory*)kkey->ReadObj();
					if(!ddir) continue;
					ddirname = ddir->GetName();
					//grab only profile hists
					if(ddirname.find("profile") == string::npos) continue;
					//needs to be geoEavg profiles to match 2D hist
					if(ddirname.find("geoEavg") == string::npos) continue;
					//needs to be binned profile - 2D = 2 "bin"s
					if(std::count(ddirname.begin(), ddirname.end(), *match) != 2) continue;
					ddir->cd();
					xbin = std::stod(dirname.substr(ddirname.find("_bin")+4,ddirname.find("_stack")-ddirname.find("_bin")-4));
					ybinname = ddirname.substr(0,ddirname.find("_bin"+to_string(xbin)));
					ybin = std::stod(ybinname.substr(ybinname.find("bin")+3));
					//these are the hists for processes so should be the same bin numbers from parent dir
					vector<TH1D*> hhists;
					//get profiles - for all time reco methods
					//cout << "in dir " << ddirname << endl;	
					//string procName;
					//procName = see root files for examples
					GetHists(ddir, hhists);
					if(hhists.size() < 1) continue;
					for(int h = 0; h < hhists.size(); h++){
						mean = GetGausMean(hhists[h]);
						if(mean == -999) continue;	
						//get TH2D for time reco h
						if(!outhists_procStack[itr][h]){ cout << "2 - outhist null" << endl; continue; }
						//fill TH2D for bin ebin, genbin with mean
						outhists_procStack[itr][h]->SetBinContent(xbin, ybin, mean);
//						cout << "xbin: " << xbin << " ybin: " << ybin << " mean: " << mean << " outhist: " << outhists_procStack[itr][h]->GetName() << " current # entries " << outhists_procStack[itr][h]->GetEntries() << endl;
					}
					itr++;
				}
			}
		}
	}
	list = f->GetListOfKeys();
	iter = TIter(list);
	key = nullptr;
	while((key = (TKey*)iter())){
		if(key->GetClassName() == tdir){
			dir = (TDirectory*)key->ReadObj();
			if(!dir) continue;
			dirname = dir->GetName();
			if(dirname.find(histdirname) == string::npos) continue;
			dir->cd();
			//gDirectory->pwd();
			for(auto hist : outhists) hist->Write("",TObject::kOverwrite); 
	
		}
	}

	for(int i = 0; i < outhists_procStack.size(); i++){
		procstackdirs[i]->cd();
		//cout << "writing process hists to dir ";
		//gDirectory->pwd();
		for(int j = 0; j < outhists_procStack[i].size(); j++){
			//cout << "Writing " << outhists_procStack[i][j]->GetName() << " with " << outhists_procStack[i][j]->GetEntries() << " entries" << endl;
			outhists_procStack[i][j]->Write("", TObject::kOverwrite);
		}
	}
}


void Profile2D(string file){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"UPDATE");
	string name, xtitle, ytitle;
	vector<string> types = {"","lead","notlead"};
//cout << "opened file" << endl;
	//profile relevant 2D histograms first
	Profile2DHists(f);
//cout << "profiled hists" << endl;
	Make2DHist(f,"geoEavg_genDeltaTime_meanRecoGenDeltaT");
//cout << "made 2d hist" << endl;
	f->Close();

	cout << "Wrote profiles for " << f->GetName() << endl;

}
