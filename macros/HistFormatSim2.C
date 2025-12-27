#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>
#include <sstream>
#include <iostream>
using std::string;

enum plotFormat{
	allStack = 0, //default original
	procStack = 1,
	methodStack = 2,
	decayStack = 3,
	diFileStack = 4,
	anyStack = 5
};

map<string, string> hist_titles;
map<string, string> proc_titles;
map<string, string> pt_threshes;
map<vector<string>, string> group_hist_titles;
vector<string> relvars = {"subclRelEtaVar","subclRelPhiVar","subclRelTimeVar"};
void BuildMaps(){
	hist_titles["Mass"] = "Jet Mass [GeV]";
	hist_titles["Jet_pt"] = "Jet p_{T} [GeV]";
	hist_titles["JetPt"] = "Jet p_{T} [GeV]";
	hist_titles["Energy"] = "Jet Energy [GeV]";
	hist_titles["TimeCenter"] = "Jet Time [ns]";
	hist_titles["EtaCenter"] = "Jet Pseudorapidity (#eta)";
	hist_titles["PhiCenter"] = "Jet Azimuthal Angle (#phi)";
	hist_titles["nJets"] = "Number of Jets";
	hist_titles["nSubclustersJet"] = "Number of Subclusters";
	hist_titles["nSubclusters"] = "Number of Subclusters";
	hist_titles["JetSize"] = "Jet Size";
	
	group_hist_titles[relvars] = "Relative Variances";

	proc_titles["ttbar"] = "t#bar{t}";
	proc_titles["QCD"] = "QCD dijets";
	proc_titles["singleW"] = "single W^{#pm}";
	proc_titles["Wgluon"] = "W^{#pm}+gluon";

	pt_threshes["ttbar"] = "250";
	pt_threshes["singleW"] = "175";
	pt_threshes["QCD"] = "175";
	pt_threshes["Wgluon"] = "175";
	
}

void stringReplaceAll(string& input, string& oldstring, string& newstring){
	while(input.find(oldstring) != string::npos){
		input.replace(input.find(oldstring), oldstring.size(), newstring);
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



void TDRMultiHist(vector<TH1D*> hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, string cms_label, plotFormat pf = allStack, vector<string> legentries = {}){
	if(can == nullptr) return;
	if(hist.size() == 0)return;
	can->cd();
	can->SetGridx(1);
	can->SetGridy(1);
	can->SetTitle("");
	TLegend* myleg = nullptr;
	//if(pf == 1) myleg = new TLegend(0.512,0.684,0.876,0.882);
	//if(pf == 1) myleg = new TLegend(0.512,0.684,0.876,0.882);
	//myleg = new TLegend(0.601,0.675,0.803,0.873);
	if(hist.size() != 6)
		myleg = new TLegend(0.691,0.675,0.873,0.873);
	else{
		myleg = new TLegend(0.471,0.665,0.873,0.873);
		myleg->SetNColumns(3);
	}
	myleg->SetFillColor(0);
	myleg->SetBorderSize(0);
	myleg->SetTextFont(42);
	if(hist.size() > 4) myleg->SetTextSize(0.03);
	//else if(pf == 1) myleg->SetTextSize(0.03);	
	else myleg->SetTextSize(0.04);
	
	//offset for log scale
	if(miny == 0) miny += 1e-6;

	string title, name, histtitle;
	string canname = can->GetName();
	if(canname[0] == '_'){
		canname = canname.substr(1);
		can->SetName(canname.c_str());
	}
	if(xtit[0] == '_'){
		xtit = xtit.substr(1);
	}

	string legentry;
	//sort hists alphabetically/numerically
	map<string, TH1D*> nameToHist;
	for( int i = 0 ; i < int(hist.size()); i++){
		name = hist[i]->GetName();
		if(name.find("AK") != string::npos){
			if(name[name.find("AK")+2] == '1') //push to end
				name = name.substr(0,name.find("AK")+2) + "9" + name.substr(name.find("AK")+2); 
		}
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
	//default
	if(pf == 0){
		labelToColor["singleW"] =   TColor::GetColor("#f6ae2d");
		labelToColor["QCD"]   	=  TColor::GetColor("#9E0059");

		labelToColor["!reco"] = TColor::GetColor("#f7a278");
		labelToColor["!AK4"] = TColor::GetColor("#f7a278");
		labelToColor["!BHC"] = TColor::GetColor("#6859f1");
		labelToColor["!BHCnoPU"] = TColor::GetColor("#9E0059");

		labelToMark["!singleW"] = 72;
		labelToMark["!QCD"] =   73;
		labelToMark["!Wgluon"] = 74;
		
		labelToMark["reco"] = 71;
		labelToMark["AK4"] = 71;
		labelToMark["BHC"] =   72; 
		labelToMark["BHCnoPU"] =   70; 


	}
	//procStack formatting
	if(pf == 1){
		labelToColor["singleW"] =   TColor::GetColor("#f6ae2d");
		labelToColor["QCD"]   	=  TColor::GetColor("#9E0059");
		
		labelToMark["singleW"] = 72;
		labelToMark["QCD"] =   73;
		//TODO: add for ttbar + Wgluon
		labelToMark["Wgluon"] = 74;
	}
	//methodStack formatting
	else if(pf == 2){
		labelToColor["genAK4"] = TColor::GetColor("#DB95D0");
		labelToColor["genAK8"] = TColor::GetColor("#A73993");
		labelToColor["genAK15"] = TColor::GetColor("#2E0F2B");
		labelToMark["genAK4"] = 74;
		labelToMark["genAK8"] = 75;
		labelToMark["genAK15"] = 76;

		labelToColor["recoAK4"] = TColor::GetColor("#F7A278");
		labelToColor["recoAK8"] = TColor::GetColor("#F36E2B");
		labelToColor["recoAK15"] = TColor::GetColor("#AE410A");
	
		labelToMark["recoAK4"] = 71;
		labelToMark["recoAK8"] = 72;
		labelToMark["recoAK15"] = 73;
		
		labelToColor["BHC"] = TColor::GetColor("#6859f1");
		labelToColor["BHCnoPU"] = TColor::GetColor("#9E0059");
		labelToMark["BHC"] =   78; 
		labelToMark["BHCnoPU"] =   77; 

	}
	//decay stack formatting
	else if(pf == 3){
		labelToColor["_b"] = TColor::GetColor("#F5B700");
		labelToColor["qg"] = TColor::GetColor("#306B34");
		labelToColor["lep"] = TColor::GetColor("#6859f1");
		
		labelToColor["fullHad"] = TColor::GetColor("#F5B700");
		labelToColor["semiLep"] = TColor::GetColor("#306B34");
		labelToColor["fullLep"] = TColor::GetColor("#6859f1");
		
		labelToMark["_b"] = 114;
		labelToMark["qg"] = 115;
		labelToMark["lep"] = 116;
		
		labelToMark["fullHad"] = 114;
		labelToMark["semiLep"] = 115;
		labelToMark["fullLep"] = 116;

	}
	else if(pf == 5){
		//use hist names as labels
		int nlabel = 70;
		int ncolor = 3;
		vector<string> ncolors = {"#3EB8F4","#0781C5","#25008B","#30B08E","#5CD6AB","#AEECCB"};
		//TODO - pull colors from consistent color palette
		for( int i = 0 ; i < int(hist.size()); i++){
			if(legentries.size() == hist.size())
				legentry = legentries[i];
			else{
				legentry = hist[i]->GetTitle();
				if(legentry.find(xtit) != string::npos){
					legentry = legentry.substr(0,legentry.find(xtit));
				}
			}
			int coloridx = i;
			if(i >= ncolors.size())
				coloridx = 0;
			labelToMark[legentry] = nlabel+i;
			labelToColor[legentry] = TColor::GetColor(ncolors[coloridx].c_str());
cout << "color " << ncolors[coloridx] << " i " << i << " coloridx " << coloridx << " legentry " << legentry << " color " << labelToColor[legentry] << endl;	
//cout << "legentry " << legentry << " title " << xtit << " title " << hist[i]->GetTitle() << " mark " << labelToMark[legentry] << " col " << labelToColor[legentry] << endl;
		}
		//for(auto leg : legentries) cout << "tdrmulti leg entries " << leg << endl;
	}


	int col, mark;
	bool highpt = false;
	bool lowpt = false;
	vector<string> legentries2;
	for( int i = 0 ; i < int(hist.size()); i++){
		//cout << "i " << i << " hists size " << hist.size() << endl;
		hist[i]->UseCurrentStyle();
		hist[i]->SetStats(false);
		hist[i]->GetXaxis()->CenterTitle(true);
		hist[i]->GetXaxis()->SetTitle(xtit.c_str());
//cout << "title " << xtit << " canname " << canname << " y title " << ytit << " histname " << hist[i]->GetName() << endl;
	if(pf == 5) cout << "name " << hist[i]->GetName() << " legentry " << legentries[i] << endl;
		hist[i]->GetYaxis()->CenterTitle(true);
		hist[i]->GetYaxis()->SetTitle(ytit.c_str());
		//else hist[i]->GetYaxis()->SetTitle("#sigma #Delta t (ns)");
		

//		cout << "miny " << miny << " max " << 3*maxy << endl;
		hist[i]->GetYaxis()->SetRangeUser(1e-4, 1.5*maxy);

		hist[i]->GetYaxis()->SetLabelFont(132);
		hist[i]->GetXaxis()->SetLabelFont(132);
		hist[i]->GetYaxis()->SetTitleFont(132);
		hist[i]->GetXaxis()->SetTitleFont(132);
		hist[i]->GetYaxis()->SetTitleSize(0.04);
		hist[i]->GetXaxis()->SetTitleOffset(1.05);
		hist[i]->GetXaxis()->SetTitleSize(0.04);


		if(pf == 0){
			legentry = hist[i]->GetTitle();
			title = legentry;	 
			histtitle = hist[i]->GetTitle();
			legentry = hist[i]->GetTitle(); 
			title = title.substr(title.find("_")+1);
		}
		else if(pf == 1){
			legentry = hist[i]->GetTitle();
			title = legentry;	 
			histtitle = hist[i]->GetTitle();
			legentry = legentry.substr(legentry.find("_")+1);
			title = title.substr(0,title.find("_"));
		}
		else if(pf == 2){ 
			legentry = hist[i]->GetTitle();
			title = legentry;	 
			histtitle = hist[i]->GetTitle();
			legentry = legentry.substr(0,legentry.find("_"));
			title = title.substr(title.find("_")+1);
		}
		else if(pf == 3){ 
			legentry = hist[i]->GetTitle();
			title = legentry;	 
			histtitle = hist[i]->GetTitle();
			string title = hist[i]->GetTitle();
			string name = hist[i]->GetName();
			legentry = name.substr(0,name.find(title));
			legentry = legentry.substr(legentry.rfind("_")+1);
			string sub1 = title.substr(0,title.find("_procStack"));
			legentry = name.substr(0,name.find(legentry)-1);
		}
		else if(pf == 4){ 
			legentry = hist[i]->GetTitle();
			title = legentry;	 
			histtitle = hist[i]->GetTitle();
			string title = hist[i]->GetTitle();
			legentry = title.substr(title.rfind("_")+1);
		}
		else if(pf == 5){
			legentry = legentries[i];//hist[i]->GetTitle();
			if(legentry.find(xtit) != string::npos){
				legentry = legentry.substr(0,legentry.find(xtit));
			}
		}
		else continue;
		//remove PD 
		if(legentry.find("PD") != string::npos)
			legentry = legentry.substr(0,legentry.find("PD"));		


		//if a key from labeltocolor is in legentry, set that color
		for(map<string, int>::iterator it = labelToColor.begin(); it != labelToColor.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(pf != 5) match += "Jet";
			if(legentry.find(match) != string::npos){
				col = it->second;
				break;
			}
			else col = 1;
		}
		for(map<string, int>::iterator it = labelToMark.begin(); it != labelToMark.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(pf != 3 && pf != 5){
				match += "Jet_";
				if(histtitle.find(match) != string::npos){
					mark = it->second;
					break;
				}
				else mark = 1;
			}
			else{
				if(legentry.find(match) != string::npos){
					mark = it->second;
					legentry = match;	
					if(strncmp(&legentry[0],"_",1) == 0) legentry = legentry.substr(1,legentry.size());
					break;
				}
				else mark = 1;

			}
		}
		string name = hist[i]->GetName();
		string type = "";
		if(name.find("_lead") != string::npos){
			type = "high p_{T}";
			mark -= 4;
			highpt = true;
			col += 2;	
		}
		else if(name.find("_notlead") != string::npos){
			type = "low p_{T}";
			mark += 1;
			lowpt = true;
			col -= 1;	
		}
		else{
			lowpt = false;
			highpt = false;
		}
		legentry = legentry.substr(0,legentry.find("Jet"));
		if(type != "") legentry += " "+type;

		legentries2.push_back(legentry);
cout << "col " << col << " mark " << mark << endl;
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

	}
	for( int i = 0 ; i < int(hist.size()); i++){
		string legentry = legentries2[i];
		if(highpt && !lowpt){
			legentry = legentry.substr(0,legentry.find(" high"));
		}	
		if(!highpt && lowpt){
			legentry = legentry.substr(0,legentry.find(" low"));
		}	
		if(legentry == "highMass") legentry = "high mass";
		if(legentry == "lowMass") legentry = "low mass";
		if(legentry == "Wmass") legentry = "W mass";
		if(legentry.find("reco") != string::npos) legentry = legentry.substr(legentry.find("reco")+4);
		myleg->AddEntry( hist[i], legentry.c_str(), "p" );

		gPad->Update();
	}	
	//cout << "plot title " << plot_title << endl;
	//draw mass lines
	if(plot_title.find("W^{#pm}") != string::npos && xtit == "Jet Mass [GeV]"){
		TLine* wmass = new TLine(80.4,1e-4, 80.4,1.5*maxy);
		wmass->SetLineStyle(6);
		wmass->Draw("same");

	}
	if(plot_title.find("t#bar{t}") != string::npos && xtit == "Jet Mass [GeV]"){
		TLine* wmass = new TLine(80.4,1e-4, 80.4,1.5*maxy);
		wmass->SetLineStyle(6);
		wmass->Draw("same");
		
		TLine* topmass = new TLine(172.,1e-4, 172.,1.5*maxy);
		topmass->SetLineStyle(6);
		topmass->Draw("same");

	}
	//cout << "highpt " << highpt << " lowpt " << lowpt << endl;
	myleg->SetTextFont(132);
	myleg->SetMargin(0.3);
	myleg->Draw("same"); 
	gPad->Update();

	string proc = canname.substr(canname.find("_")+1, canname.find("_",canname.find("_")+1) - canname.find("_")-1);
	string pt_thresh = pt_threshes[proc];
	if(highpt && !lowpt){
		string jetsel_str = "#font[132]{Jet p_{T} #geq "+pt_thresh+" GeV}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.69,0.62,jetsel_str.c_str());
	}
	if(highpt && lowpt){
		string jetsel_str = "#font[132]{#splitline{low p_{T}: Jet p_{T} < "+pt_thresh+" GeV}{high p_{T}: Jet p_{T} #geq "+pt_thresh+" GeV}}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.61,0.56,jetsel_str.c_str());
		//adjust legend bounds
		myleg->SetX1NDC(0.60);
		myleg->SetY1NDC(0.65);
		myleg->SetX2NDC(0.85);
		myleg->SetY2NDC(0.87);
		gPad->Modified();
	}
	if(canname.find("highMass") != string::npos){
		string jetsel_str = "#font[132]{high mass: Jet mass > 100 GeV}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.59,0.57,jetsel_str.c_str());

	}
	if(canname.find("Wmass") != string::npos){
		string jetsel_str = "#font[132]{W mass: 70 GeV < Jet mass < 90 GeV}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.51,0.52,jetsel_str.c_str());

	}
	if(canname.find("lowMass") != string::npos){
		string jetsel_str = "#font[132]{low mass: Jet mass < 50 GeV}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.614,0.47,jetsel_str.c_str());

	}
	

	
	//string lat_cms = "#bf{Pythia 8} event generation, #sqrt{s} = 13 TeV"+cms_label;
	string lat_cms = "#font[22]{Pythia 8} #font[132]{event generation, #sqrt{s} = 13 TeV"+cms_label+"}";
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.03,0.92,lat_cms.c_str());
	TLatex lat1;
	lat1.SetNDC();
	lat1.SetTextSize(0.04);
	lat1.SetTextFont(42);
	plot_title = "#font[132]{"+plot_title+"}";
	lat1.DrawLatex(0.75,0.92,plot_title.c_str());
	
	//draw sigma formula
	if(canname.find("sigma") != string::npos){
		TLatex sigFormula;
		sigFormula.SetNDC();
		sigFormula.SetTextSize(0.03);
		sigFormula.SetTextFont(42);
		string xstr = hist[0]->GetXaxis()->GetTitle();
		string sigstr = hist[0]->GetYaxis()->GetTitle();
		if(xstr.find("(GeV)") != string::npos) xstr.replace(xstr.find("(GeV)"),5,"");
		string paramsStr = "("+sigstr+")^{2} = #frac{N^{2}}{("+xstr+")^{2}} + #frac{S^{2}}{("+xstr+")} + 2C^{2}";
		//string paramsStr = "("+sigstr+")^{2} = #frac{N^{2}}{("+xstr+")^{2}} + C^{2}";
		cout << "formula Y: " << 0.4+(hist.size()+1)*0.05 << endl;
		sigFormula.DrawLatex(0.4, 0.4+(hist.size()+1)*0.05, paramsStr.c_str());
	}

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
	hist->GetZaxis()->SetTitle("a.u.");
	hist->GetZaxis()->CenterTitle(true);
	hist->GetYaxis()->SetLabelFont(132);
	hist->GetXaxis()->SetLabelFont(132);
	hist->GetZaxis()->SetLabelFont(132);
	hist->GetYaxis()->SetTitleFont(132);
	hist->GetXaxis()->SetTitleFont(132);
	hist->GetZaxis()->SetTitleFont(132);
	hist->GetYaxis()->SetTitleSize(0.04);
	hist->GetXaxis()->SetTitleOffset(1.05);
	hist->GetXaxis()->SetTitleSize(0.04);
	hist->GetZaxis()->SetTitleSize(0.04);
	string histname = hist->GetName();
	hist->Scale(1./hist->Integral());
	hist->Draw("colz1");
	string name = hist->GetName();
	//if(name.find("genDeltaTime_meanRecoGenDeltaT") != string::npos) cout << "n entries: " << hist->GetEntries() << endl;
	
	string lat_cms = "#font[22]{Pythia 8} #font[132]{event generation, #sqrt{s} = 13 TeV}";
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.03,0.92,lat_cms.c_str());

	string canname = can->GetName();
	string proc = canname.substr(canname.find("_")+1, canname.find("_",canname.find("_")));
	string pt_thresh = pt_threshes[proc];
	if(cms_label.find("high pt") != string::npos){
		cms_label = cms_label.substr(0,cms_label.find(" high pt"));
		string jetsel_str = "#font[132]{Jet p_{T} #geq "+pt_thresh+" GeV}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.67,0.85,jetsel_str.c_str());
	}
	if(cms_label.find("low pt") != string::npos){
		cms_label = cms_label.substr(0,cms_label.find(" low pt"));
		string jetsel_str = "#font[132]{Jet p_{T} < "+pt_thresh+" GeV}";
		TLatex jetsel;	
		jetsel.SetNDC();
		jetsel.SetTextSize(0.04);
		jetsel.SetTextFont(42);
		jetsel.DrawLatex(0.67,0.85,jetsel_str.c_str());
	}
	cms_label = "#font[132]{"+cms_label+"}";
	
	TLatex lat1;
	lat1.SetNDC();
	lat1.SetTextSize(0.04);
	lat1.SetTextFont(42);
	lat1.DrawLatex(0.69,0.92,cms_label.c_str());

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


void GetHistsProcs(TDirectory* dir, vector<string>& procs, vector<TH1D*>& hists){
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
					if(hist->Integral() != 0 && name.find("sigma") == string::npos && name.find("mean") == string::npos) hist->Scale(1./hist->Integral());
					string title = hist->GetTitle();
					if(title.find(procs[s]) == string::npos) hist->SetTitle((title+"_"+procs[s]).c_str());
					if(hist->Integral() != 0) hist->Scale(1./hist->Integral());
					else continue;
					hists.push_back(hist);

				}
				else continue;
			}
		}
	}
}


void GetHistsType(TDirectory* dir, string proc, string type, vector<TH1D*>& hists){
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
				if(name.find(proc) == string::npos) continue;
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

void ProcStackHists(string file, vector<string>& procs, string method, string oname, string match=""){
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
//	string oname = f->GetName();
//	oname = oname.substr(0,oname.find(".root"));
//	oname = oname+"_formatted.root";
	//TFile* ofile = new TFile(oname.c_str(),"RECREATE");

	//string cmslab = method;
	//string cmslab = GetCMSLabel(file);
	//string extra = "";
	//if(file.find("Skim") != string::npos) extra = GetExtraLabel(file);
	//if(!extra.empty()) cmslab += " "+extra;	

	TString th1d("TH1D");
	TString th2d("TH2D");
	TString tdir("TDirectoryFile");

	string procsname = "";
	for(auto s : procs) procsname += "_"+s;

	double ymin, ymax;
	string ylab, xlab;
	//cout << "method " << method << endl;
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
			vector<TH1D*> hists;
			//cout << "\ndir name: " << name << endl;
			if(name.find(method) == string::npos) continue;
			//we're in the directory with hists of one method split by procs
			GetHistsProcs(dir, procs, hists);
			if(hists.size() > 0){
				cout << "hists got" << endl;
				for(auto h : hists){
					cout << h->GetName() << endl;
				}
				FindListHistBounds(hists, ymin, ymax);
				if(ymin == 0 && ymax == 0) continue;
				name = dir->GetName();
				string cvname = name+procsname+"_procStack_"+method;
				TCanvas *cv = new TCanvas(cvname.c_str(), "");
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
				TDRMultiHist(hists, cv, method, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", procStack);
				cout << "writing canvas (1D) " << cv->GetName() << endl;
				cv->Write(); 
			}
		
		}
	}
	//cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;
	ofile->Write();
	ofile->Close();
	f->Close();

};


void Hist2D(string file, string proc, string method, string oname, string obs, string type = ""){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"READ");
	string histname = method+"Jet_"+obs;
	if(type != "")
		histname += "_"+type;
	TH2D* h = (TH2D*)f->Get(histname.c_str());
	if(h == nullptr) return;
	h->Scale(1./h->Integral());
	string name = h->GetName();
	string cvname = name;
	TCanvas *cv = new TCanvas(cvname.c_str(), "");
	//draw as tcanvases
	string xlab = hist_titles[obs.substr(obs.find("_")+1)];
	string ylab = hist_titles[obs.substr(0,obs.find("_"))];
	string methodname = method;
	if(methodname.find("reco") != string::npos)
		methodname = methodname.substr(4); //remove reco
	string cmslab = proc_titles[proc]+", "+methodname;
	//cout << "x " << xlab << " y " << ylab << " name " << name << " cmslab " << cmslab+type << endl;
	TDR2DHist(h, cv, xlab, ylab, cmslab+type, "");
	cv->SetRightMargin(0.15);
	cout << "writing canvas (2D) " << cv->GetName() << endl;
	
	TFile* ofile = TFile::Open(oname.c_str(), "UPDATE");
	ofile->cd();
	cv->Write(); 
	ofile->Write();
	ofile->Close();
	f->Close();


}


//method in this case would be BHC vs reco
void MethodStackHists(string file, string proc, vector<string> methods, string oname, string obs ="", vector<string> type = {}){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"READ");
	vector<TH1D*> hists;
	string histname;
	for(int m = 0; m < methods.size(); m++){
		string hname = methods[m]+"Jet_"+obs;
		if(type.size() > 0){
			for(int t = 0; t < type.size(); t++){
				histname = hname+"_"+type[t];
cout << "histname " << histname << endl;
				TH1D* h = (TH1D*)f->Get(histname.c_str());
				if(h == nullptr) continue;
				h->Scale(1./h->Integral());
				hists.push_back(h);
			}
		}
		else{
			histname = hname;
			TH1D* h = (TH1D*)f->Get(histname.c_str());
			if(h == nullptr) continue;
			h->Scale(1./h->Integral());
			hists.push_back(h);
		}
	}
	double ymin, ymax;
	string name = obs+"_"+proc+"_methodStack";
	if(find(methods.begin(), methods.end(), "BHCnoPU") != methods.end()) name += "_puComp";
	//if(type != "") name += "_"+type;
	if(type.size() > 0){
		for(int t = 0; t < type.size(); t++){
			name += "_"+type[t];
		}
	}
	if(hists.size() < 1) return;
	FindListHistBounds(hists, ymin, ymax);
	if(ymin == 0 && ymax == 0) return;
	TCanvas *cv = new TCanvas(name.c_str(), "");
	string xlab = hist_titles[obs];
	string ylab = "a.u."; 
	string cmslab = proc_titles[proc];
	
	//for(auto h : hists){ cout << "have hist " << h->GetName() << endl;  }
	//cout << "xlab " << xlab << " name " << name << " proc " << proc << endl;
	TDRMultiHist(hists, cv, cmslab, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", methodStack);
	TFile* ofile = TFile::Open(oname.c_str(), "UPDATE");
	ofile->cd();
	cv->Write(); 
		cout << "writing canvas (1D) " << cv->GetName() << endl;
				
	//cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;
	ofile->Write();
	ofile->Close();
	f->Close();

};

void AnyStackHists(string file, string proc, string method, vector<string> obses, string oname, vector<string> types = {}, vector<string> legentries = {}){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"READ");
	vector<TH1D*> hists;
	for(int o = 0; o < obses.size(); o++){
		string histname = method+"Jet_"+obses[o];
		if(types.size() > 0){
			for(int t = 0; t < types.size(); t++){
				string hname = histname+"_"+types[t];
				cout << "hname " << hname << endl;
				TH1D* h = (TH1D*)f->Get(hname.c_str());
				if(h == nullptr) continue;
				h->Scale(1./h->Integral());
				hists.push_back(h);
			}
		}
		else{
			TH1D* h = (TH1D*)f->Get(histname.c_str());
			if(h == nullptr) continue;
			h->Scale(1./h->Integral());
			hists.push_back(h);
		}
	}
	double ymin, ymax;
	if(hists.size() < 1) return;
	FindListHistBounds(hists, ymin, ymax);
	if(ymin == 0 && ymax == 0) return;
	string xlab = group_hist_titles[obses];
	string ylab = "a.u."; 
	string cmslab = proc_titles[proc];
	string name = xlab+"_"+proc+"_anyStack";
	while(name.find(" ") != string::npos)
		name = name.replace(name.find(" "),1,"");
	TCanvas *cv = new TCanvas(name.c_str(), "");
	
	//for(auto h : hists){ cout << "have hist " << h->GetName() << endl;  }
	//cout << "xlab " << xlab << " name " << name << " proc " << proc << endl;
	TDRMultiHist(hists, cv, cmslab, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", anyStack,legentries);
	TFile* ofile = TFile::Open(oname.c_str(), "UPDATE");
	ofile->cd();
	cv->Write(); 
		cout << "writing canvas (1D) " << cv->GetName() << endl;
				
	//cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;
	ofile->Write();
	ofile->Close();
	f->Close();
};








void FileStackHists(vector<string>& files, vector<string>& labels, string proc, string method, string oname, string match="",string plottitle=""){
	TFile* ofile = TFile::Open(oname.c_str(),"UPDATE");
	string cmslab = "";
	if(proc == "JetHTPD"){
		cmslab = "JetHT, Run F 2017";
	}
	else if(proc == "DEGPD"){
		cmslab = "DoubleEG, Run F 2017";
	}
	else cmslab = "";
	cmslab += " "+method;
	vector<TH1D> hists;
	//GetFileLabels(files,labels);
	for(int f = 0; f < files.size(); f++){
		if(gSystem->AccessPathName(files[f].c_str())){
			cout << "File " << files[f] << " does not exist." << endl;
			return;
		}
		TFile* file = TFile::Open(files[f].c_str(),"READ");
		TList* list = file->GetListOfKeys();
		TIter iter(list);
		TKey* key;
		string name, xtitle, ytitle;
//		string oname = f->GetName();
//		oname = oname.substr(0,oname.find(".root"));
//		oname = oname+"_formatted.root";
		//TFile* ofile = new TFile(oname.c_str(),"RECREATE");

		//string cmslab = method;
		//string cmslab = GetCMSLabel(file);
		//string extra = "";
		//if(file.find("Skim") != string::npos) extra = GetExtraLabel(file);
		//if(!extra.empty()) cmslab += " "+extra;	

		TString th1d("TH1D");
		TString th2d("TH2D");
		TString tdir("TDirectoryFile");

		string dirname, ddirname, histname, histtitle;

		while((key = (TKey*)iter())){
			name = key->GetName();
			//skip these dirs
			if(name.find("genDeltaTpvGambin") != string::npos) continue;
			if(key->GetClassName() == tdir){
				TDirectory* dir = dynamic_cast<TDirectory*>(key->ReadObj());
				double ymin, ymax;
				string ylab, xlab;
				if(!dir) continue;
				dirname = dir->GetName();
				if(dirname.find(match) == string::npos) continue;
				cout << "\ndir name: " << dir->GetName() << endl;
				dir->cd();
				TList* llist = dir->GetListOfKeys();
				TIter iiter(llist);
				TKey* kkey;
				//get subdir - per method all procs
				while((kkey = (TKey*)iiter())){
					if(kkey->GetClassName() == tdir){
						TDirectory* ddir = dynamic_cast<TDirectory*>(kkey->ReadObj());
						if(!ddir) continue;
						ddirname = ddir->GetName();
						//make sure this is for the specified method
						if(ddirname.find(method) == string::npos) continue;
						ddir->cd();
						cout << " ---in ddir " << ddir->GetName() << endl;
						//we're in the directory with hists of one method split by procs
						histname = ddirname;
						histname = histname.substr(0,histname.find("_procStack"));
						histname += "_"+proc;
						TH1D* hist = (TH1D*)file->Get((dirname+"/"+ddirname+"/"+histname).c_str());
						if(hist) cout << "got histogram " << hist->GetName() << " " << hist->GetEntries() << " " << hist->GetTitle() << endl;
						else cout << "hist null " << dirname+"/"+ddirname+"/"+histname << endl;
						histtitle = hist->GetTitle();
						histname += "_"+labels[f];
						hist->SetTitle((histtitle+"_"+labels[f]).c_str());	
						hist->SetName(histname.c_str());	
						//if(match.find("profile") != string::npos) hist->Scale(1./hist->Integral());
						if(hist) hists.push_back(*hist);
					}
				}	


			}
		}
		file->Close();
	}
		double ymin, ymax;
		string ylab, xlab;
		vector<TH1D*> histsp;
		if(hists.size() > 0){
			//make into pointers because...yeah...
			for(vector<TH1D>::iterator h = hists.begin(); h != hists.end(); h++) histsp.push_back(&(*h));
			FindListHistBounds(histsp, ymin, ymax);
			if(ymin == 0 && ymax == 0) return;
			string name = match+"_"+plottitle+"_"+proc+"_"+method;
			TCanvas *cv = new TCanvas(name.c_str(), "");
			ofile->cd();
			//draw as tcanvases
			if(match.find("profile") != string::npos){
				xlab = histsp[0]->GetName();
				ylab = histsp[0]->GetYaxis()->GetTitle();

			}
			else{
				xlab = histsp[0]->GetXaxis()->GetTitle();
				ylab = histsp[0]->GetYaxis()->GetTitle();
			}
			TDRMultiHist(histsp, cv, cmslab, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", diFileStack);
			cout << "writing canvas (1D) " << cv->GetName() << endl;
			cv->Write(); 
		}
		ofile->Write();
		ofile->Close();
};


void HistFormatSim2(string file){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	BuildMaps();
	string match1 = "condorSim_";
	string match2 = "_defaultv";
	string proc = file.substr(file.find(match1)+match1.size(), file.find(match2) - file.find(match1) - match2.size() - 1);
	string oname = file;
	oname = oname.substr(0,oname.find(".root"));

	oname = oname+"_formatted.root";
	TFile* ofile = new TFile(oname.c_str(),"RECREATE");
	ofile->Close();	

	vector<string> jettypes_recoBHC = {"recoAK4","recoAK8", "recoAK15", "BHC"};
	vector<string> jettypes_recoAK4BHC = {"recoAK4", "BHC"};
	vector<string> jettypes_recoAK8BHC = {"recoAK8", "BHC"};
	vector<string> jettypes_recoAK15BHC = {"recoAK15", "BHC"};
	vector<string> jettypes_BHCPU = {"BHC","BHCnoPU"};
	
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "nJets"); 
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "nJets",{"lead"}); 
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "JetSize"); 
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "JetSize",{"lead"}); 
	
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "EtaCenter"); 
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "PhiCenter"); 
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "TimeCenter"); 
	MethodStackHists(file, proc, jettypes_recoAK8BHC, oname, "TimeCenter",{"lead","notlead"});
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "Energy"); 
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "Mass",{"lead"});
	MethodStackHists(file, proc, jettypes_recoBHC, oname, "nSubclusters",{"lead"});

		
	MethodStackHists(file, proc, jettypes_BHCPU, oname, "Mass",{"lead"});
	MethodStackHists(file, proc, jettypes_BHCPU, oname, "nSubclusters",{"lead"});

	Hist2D(file, proc, "BHC", oname, "Mass_JetSize");
	Hist2D(file, proc, "recoAK4", oname, "Mass_JetSize");
	Hist2D(file, proc, "recoAK8", oname, "Mass_JetSize");
	Hist2D(file, proc, "recoAK15", oname, "Mass_JetSize");
	Hist2D(file, proc, "BHC", oname, "Mass_JetSize","lead");
	Hist2D(file, proc, "recoAK4", oname, "Mass_JetSize","lead");
	Hist2D(file, proc, "recoAK8", oname, "Mass_JetSize","lead");
	Hist2D(file, proc, "recoAK15", oname, "Mass_JetSize","lead");
	Hist2D(file, proc, "BHC", oname, "Mass_JetPt");
	Hist2D(file, proc, "recoAK4", oname, "Mass_JetPt");

	vector<string> legentries = {"#sigma^{2*}_{#eta}, E* #geq 0.5","#sigma^{2*}_{#eta}, E* < 0.5","#sigma^{2*}_{#phi}, E* #geq 0.5", "#sigma^{2*}_{#phi}, E* < 0.5", "#sigma^{2*}_{time}, E* #geq 0.5","#sigma^{2*}_{time}, E* < 0.5"};

	if(proc == "ttbar"){
		AnyStackHists(file, proc, "BHC", relvars, oname, {"lead_TopMatched_ge2Subcls_subclRelEge0p5","lead_TopMatched_ge2Subcls_subclRelElt0p5"}, legentries);
	}
	/*
	if(proc.find("W") != string::npos){
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_nSubclusters","lead");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_nSubclusters");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_subclusterEnergy");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_subclusterMass");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_subclusterLeadInvMass");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_subclusterEtaCenter");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_subclusterPhiCenter");
		MethodStackHists(file, proc, jettypes_recoAK15BHC, oname, "W_subclusterTimeCenter");
		
		AnyStackHists(file, proc, "BHC", {"W_highMass_nSubclustersJet","W_lowMass_nSubclustersJet","W_Wmass_nSubclustersJet"},oname,{"lead"},{"highMass","lowMass","Wmass"});


		Hist2D(file, proc, "BHC", oname, "W_dRGenPartons");
		Hist2D(file, proc, "recoAK4", oname, "W_dRGenPartons");
		Hist2D(file, proc, "recoAK8", oname, "W_dRGenPartons");
		Hist2D(file, proc, "recoAK15", oname, "W_dRGenPartons");
		//Hist2D(file, proc, "BHC", oname, "BHCJetW_subclEnergy_subclLeadIdx");
		
		AnyStackHists(file, proc, "BHC", {"W_highMass_partonMatchSubclPt_","W_highMass_partonNoMatchSubclPt_"}, oname, {"lead"},{"partonMatch","partonNoMatch"});
		AnyStackHists(file, proc, "BHC", {"W_highMass_partonMatchSubclSize_","W_highMass_partonNoMatchSubclSize_"}, oname, {"lead"},{"partonMatch","partonNoMatch"});
		AnyStackHists(file, proc, "BHC", {"W_highMass_partonMatchSubclPtOvJetPt","W_highMass_partonNoMatchSubclPtOvJetPt"}, oname, {"lead"},{"partonMatch","partonNoMatch"});
		AnyStackHists(file, proc, "BHC", {"W_highMass_partonMatchSubclTimeSigOvJetTimeSig","W_highMass_partonNoMatchSubclTimeSigOvJetTimeSig"}, oname, {"lead"},{"partonMatch","partonNoMatch"});
		Hist2D(file, proc, "BHC", oname, "W_highMass");

		//PU cleaning hists - labels (last arg) need to match (or be in) given hist names
		AnyStackHists(file, proc, "BHC", {"subclTimeCenter_PUlike","subclTimeCenter_PUcleaned"},oname,{"lead"},{"PUlike","PUcleaned"});
		AnyStackHists(file, proc, "BHC", {"subclusterTimeSig_PUlike","subclusterTimeSig_PUcleaned"},oname,{"lead"},{"PUlike","PUcleaned"});


		//causes segfault - need to rerun with _nom for the first hist
		AnyStackHists(file, proc, "BHC",{"BHCJetW_subclParton_dR_nom","BHCJetW_subclParton_dR_PUcleaned"},oname,{"lead"},{"","PUcleaned"});
		AnyStackHists(file, proc, "BHC",{"BHCJetW_subclParton_dR"},oname,{"lead"},{"","PUcleaned"});
		AnyStackHists(file, proc, "BHC",{"BHCJet_mass_","BHCJet_PUremoved_mass","BHCJet_PUdownweighted_mass","W_subclusterLeadInvMass"},oname,{"lead"},{"_mass","PUremoved","PUdownweighted","subclusterLeadInvMass"});
	}
	if(proc.find("gluon") != string::npos){
		MethodStackHists(file, proc, jettypes_recoBHC, oname, "genGluon_Eratio","lead");
		MethodStackHists(file, proc, jettypes_recoBHC, oname, "genGluon_dR","lead");
		MethodStackHists(file, proc, {"BHC"}, oname, "Gluon_nSubclusters","lead");
		MethodStackHists(file, proc, {"BHC"}, oname, "BHCJetGluon_subclParton_dR");
		MethodStackHists(file, proc, {"BHC"}, oname, "BHCJetGluon_subclParton_Eratio");
		MethodStackHists(file, proc, {"BHC"}, oname, "BHCJetGluon_subclParton_Eratio","lead");
		MethodStackHists(file, proc, {"BHC"}, oname, "BHCJetGluon_subclParton_dR","lead");
	}

	if(proc == "ttbar"){	
		MethodStackHists(file, proc, {"BHC"}, oname, "Top_nSubclusters","lead");
		MethodStackHists(file, proc, jettypes_recoAK4BHC, oname, "Top_subclusterMass");
		MethodStackHists(file, proc, jettypes_recoAK4BHC, oname, "Top_subclusterLeadInvMass");
	}
	*/
	cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;

}


