#include <string>
#include <TFile.h>
using std::string;


void MakePDFs(string file, string odir, string histname, bool logy = false, bool logz = false, bool lead = false){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	if(file.find("formatted") == string::npos && (histname.find("jetcan") == string::npos) && histname.find("evtdisp") == string::npos){
		cout << "Please give formatted file with TCanvases." << endl;
		return;
	}
	string odirname = file.substr(0,file.find_last_of("_"));
	odirname = odirname.substr(odirname.find("/")+1);
	odir = "plots_local/"+odirname+"/"+odir+"/";
	if(gSystem->AccessPathName(odir.c_str())){
		gSystem->mkdir(odir.c_str(),kTRUE);
	}
	TFile* f = TFile::Open(file.c_str());
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	TCanvas* cv;
	string name;
	TString cvstring("TCanvas");
	while((key = (TKey*)iter())){
		if(key->GetClassName() != cvstring) continue;
		cv = (TCanvas*)key->ReadObj();
		name = cv->GetName();
		if(name.find(histname) != string::npos){
			if(lead){
				if(name.find("_lead") == string::npos) continue;
			}
			if(logy){
				cv->SetLogy();
				name += "log";	
			}
			if(logz){
				cv->SetLogz();
				name += "log";	
			}
			cv->SaveAs((odir+name+".pdf").c_str());
		}
		
	}
	cout << "Writing plots to " << odir << endl;
};
