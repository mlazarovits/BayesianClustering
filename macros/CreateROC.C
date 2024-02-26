#include <string>
#include "TColor.h"
#include "TFile.h"

using std::string;


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


//void CreateROC(string file="skims/condor_photons_defaultv2_GJets_AOD_v14_GJets_HT-400To600.root", string hist="timeeta_cov_stack/timeeta_cov_GJets"){
void CreateROC(vector<string> filenames={}, string hist="", string odir="./"){
if(filenames.empty()){ cout << "Please provide file name" << endl; return; }
if(hist == ""){ cout << "Please provide hist name" << endl; return; }
vector<int> cols;
cols.push_back(TColor::GetColor("#86bbd8"));
cols.push_back(TColor::GetColor("#f6ae2d"));
	
string histname;
if(hist.find("/") != string::npos) histname = hist.substr(hist.find("/")+1);
else histname = hist;


TLegend* leg = new TLegend(0.11,0.76,0.26,0.89);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextFont(42);
leg->SetTextSize(0.02);

TCanvas* cv = new TCanvas((histname+"_cv").c_str(),(histname+"_cv").c_str(),800,600);
cv->cd();
cv->SetGridx();
cv->SetGridy();
cv->SetTitle((histname+" ROC").c_str());

for(int f = 0; f < filenames.size(); f++){
	if(gSystem->AccessPathName(filenames[f].c_str())){
		cout << "File " << filenames[f] << " does not exist" << endl;
		return;
	}
	
	
	TFile* file = TFile::Open(filenames[f].c_str(), "READ");
	string cms_label = GetCMSLabel(filenames[f]);
	
	vector<string> procs;
	procs.push_back("GJets");
	procs.push_back("chiGam");
	
	
	vector<TH1D*> hists;
	vector<double> tots;
	int nBins = 0;
	for(int p = 0; p < procs.size(); p++){
		hists.push_back((TH1D*)file->Get((hist+"_stack/"+hist+"_"+procs[p]).c_str()));
		if(!hists[p]){ cout << "Histogram " << hist+"_"+procs[p] << " not found in file " << file->GetName() << endl; return; }
		nBins = hists[p]->GetNbinsX()+1;
		tots.push_back(double(hists[p]->Integral()));
	}
	
	cout << "Getting ROC for " << histname << endl;
	if(hists[0]->GetNbinsX()+1 != nBins){ cout << "# bins do not match." << endl; return;}
	int nCuts = nBins;
	vector<double> xs, ys;
	int bin;
	for(int i = 1; i < nCuts+1; i++){
		bin = i*nBins/nCuts;
		xs.push_back(double(hists[0]->Integral(0,bin))/tots[0]);
		ys.push_back(double(hists[1]->Integral(0,bin))/tots[1]);
		cout << "i: "<< i << " integral to bin " << bin << " x: " << hists[0]->Integral(0,bin) << " total " << tots[0] << " y: " << hists[1]->Integral(0,bin) << " total: " << tots[1] << endl;
	
	}
	cout << "npts " << xs.size() << " " << ys.size() << " nbins " << nBins << endl;
	
	TGraph* roc = new TGraph(nCuts, &xs[0], &ys[0]);
	roc->SetTitle("");
	roc->GetXaxis()->SetTitle("Bkg Efficiency");
	roc->GetYaxis()->SetTitle("Signal Efficiency");
	roc->SetMarkerSize(1.2);
	roc->SetMarkerStyle(22+f);
	roc->SetMarkerColor(cols[f]);
	roc->SetLineColor(cols[f]);
	roc->GetXaxis()->CenterTitle(true);
	roc->GetYaxis()->CenterTitle(true);
	roc->GetXaxis()->SetRangeUser(0.,1.);
	roc->GetYaxis()->SetRangeUser(0.,1.);
	leg->AddEntry(roc, cms_label.c_str(), "l");
	
	if(f == 0) roc->Draw("AC"); 
	else roc->Draw("sameC");
}
leg->Draw("same");


string lat_cms = "#bf{CMS} #it{WIP} "+histname+"ROC";
TLatex lat;
lat.SetNDC();
lat.SetTextSize(0.025);
lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());
string odirname = "";
for(int f = 0; f < filenames.size(); f++){
	odirname += filenames[f].substr(0,filenames[f].find_last_of("."));
	if(f != filenames.size()-1) odirname += "_";
}
odirname = odirname.substr(odirname.find("/")+1);
odir = "plots_local/"+odirname+"/"+odir+"/";
if(gSystem->AccessPathName(odir.c_str())){
	gSystem->mkdir(odir.c_str(),kTRUE);
}
cout << "Saving as " << odir+histname+"_stackROC.pdf" << endl;
cv->SaveAs((odir+histname+"_stackROC.pdf").c_str());


};
