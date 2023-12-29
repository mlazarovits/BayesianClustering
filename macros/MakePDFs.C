#include <string>
#include <TFile.h>
using std::string;


void MakePDFs(string file, string odir, string histname){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	string odirname = file.substr(0,file.find_last_of("."));
	odirname = odirname.substr(odirname.find("/")+1);
	odir = "plots_local/"+odirname+"/"+odir+"/";
	cout << "Writing plots to " << odir << endl;
	if(gSystem->AccessPathName(odir.c_str())){
		gSystem->mkdir(odir.c_str(),kTRUE);
	}
	TFile* f = TFile::Open(file.c_str());
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	TCanvas* cv;
	string name;
	while((key = (TKey*)iter())){
		cv = (TCanvas*)key->ReadObj();
		name = cv->GetName();
		if(name.find(histname) != string::npos){
			cv->SaveAs((odir+name+".pdf").c_str());
		}
		
	}
};
