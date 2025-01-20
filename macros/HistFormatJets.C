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
	dijetRecoGenStack = 3,
	diFileStack = 4
};


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
	//string legName = "GMSB, L = "+lfancy+" c#tau = "+ctfancy;
	string legName = "SMS #tilde{g}#tilde{g}";//, L = "+lfancy+" c#tau = "+ctfancy;
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



void TDRMultiHist(vector<TH1D*> hist, TCanvas* &can, string plot_title, string xtit, string ytit, double miny, double maxy, string cms_label, plotFormat pf = allStack){
	if(can == nullptr) return;
	if(hist.size() == 0)return;
	can->cd();
	can->SetGridx(1);
	can->SetGridy(1);
	can->SetTitle("");
	TLegend* myleg = nullptr;
	//if(pf == 1) myleg = new TLegend(0.512,0.684,0.876,0.882);
	//if(pf == 1) myleg = new TLegend(0.512,0.684,0.876,0.882);
	myleg = new TLegend(0.143,0.686,0.345,0.884);
	//else myleg = new TLegend(0.696,0.696,0.875,0.882);
	myleg->SetFillColor(0);
	myleg->SetBorderSize(0);
	myleg->SetTextFont(42);
	if(hist.size() > 4) myleg->SetTextSize(0.025);
	//else if(pf == 1) myleg->SetTextSize(0.03);	
	else myleg->SetTextSize(0.04);
	
	//offset for log scale
	if(miny == 0) miny += 1e-6;

	string title, name, histtitle;
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
	//default
	if(pf == 0){
		//labelToColor["GMSB"] =  TColor::GetColor("#86bbd8");
		labelToColor["SMS"] =  TColor::GetColor("#86bbd8");
		labelToColor["c#tau = 200"] = TColor::GetColor("#55DBCB");
		labelToColor["c#tau = 0.1"] = TColor::GetColor("#86bbd8");
		labelToColor["c#tau = 1000"] = TColor::GetColor("#39A2AE");
		labelToColor["c#tau = 800"] = TColor::GetColor("#39A2AE");
		labelToColor["notSunm"] = TColor::GetColor("#9e0059");
		labelToColor["GJets"] =   TColor::GetColor("#f6ae2d");
		labelToColor["QCD"] =   TColor::GetColor("#f6ae2d");
		labelToColor["JetHT"] =   TColor::GetColor("#3d348b");
		labelToColor["MET"] =   TColor::GetColor("#671E76");
		labelToColor["DoubleEG"] = TColor::GetColor("#9e0059");

		labelToColor["!median"] = TColor::GetColor("#f7a278");
		labelToColor["!eAvg"] = TColor::GetColor("#6859f1");
		labelToColor["!mmAvg"] = TColor::GetColor("#52b788");
		labelToColor["!eMax"] = TColor::GetColor("#E2C2FF");

		//MC symbols - primary shapes
		labelToMark["!chiGam"] =  20;
		labelToMark["!GMSB"] =  20;
		labelToMark["!GluGlu"] =  20;
		labelToMark["!SMS"] =  20;
		labelToMark["!Ctau200"] = -1;
		labelToMark["!Ctau0p1"] = 1;
		labelToMark["!Ctau1000"] = 2;
		labelToMark["!Ctau800"] = 2;
		labelToMark["!notSunm"] = 72;
		labelToMark["!GJets"] =   73;
		labelToMark["!QCD"] =   73;
		//data symbols - some form of open cross
		labelToMark["!JetHT"] =   75;
		labelToMark["!MET"] =   85;
		labelToMark["!DoubleEG"] =   83;
		
		labelToMark["median"] = 71;
		labelToMark["eAvg"] =   72; 
		labelToMark["mmAvg"] =  73;
		labelToMark["eMax"] =   74; 


	}
	//procStack formatting
	if(pf == 1){
		labelToColor["chiGam"] =  TColor::GetColor("#86bbd8");
		labelToColor["GluGlu"] =  TColor::GetColor("#86bbd8");
		labelToColor["SMS"] =  TColor::GetColor("#86bbd8");
		//labelToColor["GMSB"] =  TColor::GetColor("#86bbd8");
		labelToColor["c#tau = 200"] = TColor::GetColor("#55DBCB");
		labelToColor["c#tau = 0.1"] = TColor::GetColor("#86bbd8");
		labelToColor["c#tau = 1000"] = TColor::GetColor("#39A2AE");
		labelToColor["c#tau = 800"] = TColor::GetColor("#39A2AE");
		labelToColor["GJets"] =   TColor::GetColor("#f6ae2d");
		labelToColor["QCD"] =   TColor::GetColor("#f6ae2d");
		labelToColor["JetHT"] =   TColor::GetColor("#3d348b");
		labelToColor["MET"] =   TColor::GetColor("#671E76");
		labelToColor["DoubleEG"] = TColor::GetColor("#9e0059");

		//MC symbols - primary shapes
		//labelToMark["chiGam"] =  20;
		//labelToMark["GMSB"] =  20;
		labelToMark["SMS"] =  20;
		labelToMark["GluGlu"] =  20;
		labelToMark["Ctau200"] = 2;
		labelToMark["Ctau0p1"] = 3;
		labelToMark["Ctau1000"] = 4;
		labelToMark["Ctau800"] = 4;
		labelToMark["notSunm"] = 72;
		labelToMark["GJets"] =   73;
		labelToMark["QCD"] =   73;
		//data symbols - some form of open cross
		labelToMark["JetHT"] =   75;
		labelToMark["MET"] =   85;
		labelToMark["DoubleEG"] =   83;



	}
	//methodStack formatting
	else if(pf == 2){
		labelToColor["median"] = TColor::GetColor("#f7a278");
		labelToColor["eAvg"] = TColor::GetColor("#6859f1");
		labelToColor["mmAvg"] = TColor::GetColor("#52b788");
		labelToColor["eMax"] = TColor::GetColor("#E2C2FF");
	
		labelToMark["median"] = 71;
		labelToMark["eAvg"] =   72; 
		labelToMark["mmAvg"] =  73;
		labelToMark["eMax"] =   74; 

	}
	//dijetRecoGenStack formatting
	else if(pf == 3){
		labelToColor["dijets"] = TColor::GetColor("#776DA7"); 
		labelToColor["recoGen"] = TColor::GetColor("#D05340");
		labelToColor["gamPV"] = TColor::GetColor("#48A9A6");
		
		labelToMark["dijets"] = 104;
		labelToMark["recoGen"] = 105;
		labelToMark["gamPV"] = 106;

	}
	//diFileStack formatting
	else if(pf == 4){
		labelToColor["preCalib"] = TColor::GetColor("#F5B700");
		labelToColor["postCalib"] = TColor::GetColor("#306B34");
		labelToColor["wSpikes"] = TColor::GetColor("#F5B700");
		labelToColor["woSpikes"] = TColor::GetColor("#306B34");
		
		labelToMark["preCalib"] = 114;
		labelToMark["postCalib"] = 115;
		labelToMark["wSpikes"] = 114;
		labelToMark["woSpikes"] = 115;

	}
	else{
		cout << "plotFormat option " << pf << " not available." << endl;
		return;
	}


	int col, mark;	
	for( int i = 0 ; i < int(hist.size()); i++){
		cout << "i " << i << " hists size " << hist.size() << endl;
		hist[i]->UseCurrentStyle();
		hist[i]->SetStats(false);
		hist[i]->GetXaxis()->CenterTitle(true);
		hist[i]->GetXaxis()->SetTitle(xtit.c_str());
cout << "title " << xtit << " canname " << canname << endl;
		hist[i]->GetYaxis()->CenterTitle(true);
		if(pf != 3) hist[i]->GetYaxis()->SetTitle(ytit.c_str());
		else hist[i]->GetYaxis()->SetTitle("#sigma #Delta t (ns)");
		cout << "miny " << miny << " max " << 3*maxy << endl;
		hist[i]->GetYaxis()->SetRangeUser(0, 1.5*maxy);
		if(canname.find("meanDeltaTime") == string::npos) hist[i]->GetYaxis()->SetRangeUser(0, 2.5*maxy);
		else hist[i]->GetYaxis()->SetRangeUser(miny, 1.5*maxy);
		

		legentry = hist[i]->GetTitle();
		title = legentry;	 
		histtitle = hist[i]->GetTitle();
		if(pf == 0){
			legentry = hist[i]->GetTitle(); 
			title = title.substr(title.find("_")+1);
			//if(title.find("chiGam") != string::npos){
			if(title.find("SMS") != string::npos){
				title = SignalLegEntry(title);
			}
		}
		else if(pf == 1){
			legentry = legentry.substr(legentry.find("_")+1);
			title = title.substr(0,title.find("_"));
			//if(legentry.find("chiGam") != string::npos){
			if(legentry.find("SMS") != string::npos){
				legentry = SignalLegEntry(legentry);
			}
		}
		else if(pf == 2){ 
			legentry = legentry.substr(0,legentry.find("_"));
			title = title.substr(title.find("_")+1);
			//if(title.find("chiGam") != string::npos){
			if(title.find("SMS") != string::npos){
				title = SignalLegEntry(title);
			}
		}
		else if(pf == 3){ 
			string title = hist[i]->GetTitle();
			title = "_"+title;
			string name = hist[i]->GetName();
			legentry = name.substr(0,name.find(title));
			legentry = legentry.substr(legentry.rfind("_")+1);
		}
		else if(pf == 4){ 
			string title = hist[i]->GetTitle();
			legentry = title.substr(title.rfind("_")+1);
		}
		else continue;
		//remove PD 
		if(legentry.find("PD") != string::npos)
			legentry = legentry.substr(0,legentry.find("PD"));		

		cout << " legentry " << legentry << " title " << title << " histtitle " << histtitle << endl;
		//if a key from labeltocolor is in legentry, set that color
		for(map<string, int>::iterator it = labelToColor.begin(); it != labelToColor.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(legentry.find(match) != string::npos){
				col = it->second;
				break;
			}
			else col = 1;
			//if(match.find("chiGam") != string::npos) cout << "legentry " << legentry << " match " << match << " col " << col << endl;
		}
		for(map<string, int>::iterator it = labelToMark.begin(); it != labelToMark.end(); it++){
			string match = it->first;
			if(match.find("!") != string::npos) match = match.substr(match.find("!")+1);
			if(pf != 3){
				if(histtitle.find(match) != string::npos){
					mark = it->second;
					break;
				}
				else mark = 1;
			}
			else{
				if(legentry.find(match) != string::npos){
					mark = it->second;
					break;
				}
				else mark = 1;

			}
		}
		//cout << "hist " << hist[i]->GetName() << " col " << col << " mark " << mark << endl;
		cout << "hist " << histtitle << " col " << col << " mark " << mark << endl;
		
		hist[i]->SetLineColor(col);
		//hist[i]->SetLineWidth(2);
		hist[i]->SetMarkerStyle(mark);
		hist[i]->SetMarkerColor(col);
		//hist[i]->SetMarkerSize(1);
		if(i == 0)
			title = hist[i]->GetTitle();
		hist[i]->SetTitle("");
		if( i == 0 ){
			if(canname.find("meanDeltaTime") != string::npos) hist[i]->GetYaxis()->SetTitleOffset(1.2);	
			hist[i]->Draw("ep");
		}else{
			hist[i]->Draw("epsame");
		}
		myleg->AddEntry( hist[i], legentry.c_str(), "p" );

		gPad->Update();

		if(canname.find("sigma") != string::npos && hist[i]->GetEntries() > 3){
			cout << "do fit for sigma" << endl;
			//string formula = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))";
			string formula = "sqrt((([0]*[0])/(x*x))+([1]*[1]/x)+(2*[2]*[2]))";
			TFormula* form = new TFormula("resFormula",formula.c_str());
	
			double xlo = hist[i]->GetXaxis()->GetBinLowEdge(hist[i]->GetXaxis()->GetFirst());
			double xhi = hist[i]->GetXaxis()->GetBinUpEdge (hist[i]->GetXaxis()->GetLast());
			TF1* fit = new TF1("fit",form->GetName(),xlo,xhi); 
			fit->SetLineColor(hist[i]->GetLineColor());
			hist[i]->Fit("fit","RQM0");
			fit->Draw("same");
			//myleg->AddEntry(fit,(legentry+" fit").c_str(),"l");
			//gPad->Update();
			//draw fit parameters on plot
			TLatex fitparams;
			fitparams.SetNDC();
			fitparams.SetTextSize(0.03);
			fitparams.SetTextFont(42);
			fitparams.SetTextColor(col);
			double val0 = fit->GetParameter(0);
			double val1 = fit->GetParameter(1);
			double val2 = fit->GetParameter(2);
			double err0 = fit->GetParError(0);
			double err1 = fit->GetParError(1);
			double err2 = fit->GetParError(2);
			std::ostringstream ss;
			ss << setprecision(2); 
			if(fabs(val0) < 9e-3)
				ss << std::scientific;
			ss << "N = " << val0;
			if(fabs(err0) < 9e-3)
				ss << std::scientific;
			else
				ss << std::fixed;	
			ss << " #pm " << err0 << " [GeV*ns],";
		
			if(fabs(val1) < 9e-3)
				ss << std::scientific;
			else
				ss << std::fixed;	
			ss << " S = " << val1;
			if(fabs(err1) < 9e-3)
				ss << std::scientific;
			else
				ss << std::fixed;	
			ss << " #pm " << err1 << " [#sqrt{GeV}*ns],";
			
			if(fabs(val2) < 9e-3)
				ss << std::scientific; 
			else
				ss << std::fixed;	
			ss << " C = " << fabs(val2);
			if(fabs(err2) < 9e-3)
				ss << std::scientific; 
			else
				ss << std::fixed;	
			ss << " #pm " << err2 << " [ns]";
			
			string teststr = ss.str();
			cout << "params " << i << " Y: " << 0.3+(hist.size()+1)*0.05-i*0.05 << endl;
			fitparams.DrawLatex(0.2,0.3+(hist.size()+1)*0.05-i*0.05,teststr.c_str());
		}
	}
	myleg->Draw("same"); 
	if(canname.find("meanDeltaTime") != string::npos){
		gPad->Update();
		myleg->SetX1NDC(0.13);
		myleg->SetY1NDC(0.18);
		myleg->SetX2NDC(0.33);
		myleg->SetY2NDC(0.38);
	} 
	gPad->Update();
	string lat_cms = "#bf{CMS} #it{Work in Progress} "+cms_label;
	TLatex lat;
	lat.SetNDC();
	lat.SetTextSize(0.04);
	lat.SetTextFont(42);
	lat.DrawLatex(0.02,0.92,lat_cms.c_str());
	TLatex lat1;
	lat1.SetNDC();
	lat1.SetTextSize(0.04);
	lat1.SetTextFont(42);
	lat1.DrawLatex(0.50,0.92,plot_title.c_str());

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
	string histname = hist->GetName();
	if(histname.find("Neighbors") != string::npos && histname.find("norm") == string::npos)
		hist->Scale(1./hist->Integral());

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
		if(key->GetClassName() == tdir){
			cout << "key " << key->GetName() << endl;
			for(int m = 0; m < methods.size(); m++){
				cout << "m " << m << endl;
				TDirectory* ddir = dynamic_cast<TDirectory*>(key->ReadObj());
				if(!ddir) continue;
				cout << "ddir not null" << endl;
				string name = ddir->GetName();
				cout << "name " << name << " method " << methods[m] << endl;
				if(name.find(methods[m]) == string::npos) continue; 
				cout << "methods[m] " << methods[m] << endl;
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
					hists.push_back(hist);

				}
				else continue;
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

//void HistFormat(string file){
void AllHists(string file){
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
			TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
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
			TH2D* hist = dynamic_cast<TH2D*>(key->ReadObj());
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
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = dynamic_cast<TDirectory*>(key->ReadObj());
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
					TH2D* hist = dynamic_cast<TH2D*>(kkey->ReadObj());
					if(!hist) continue;	
					//cout << " getting kkey: " << kkey->GetName() << " hist " << hist->GetName() << " with entries " << hist->GetEntries() << endl;
				}
				//writing stack hist - same method different procs
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = dynamic_cast<TDirectory*>(kkey->ReadObj());
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

	//ofile->Close();
	f->Close();

};









//void ProcStackHists(string file, vector<string>& procs, string method, TFile *ofile, string match=""){
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

	while((key = (TKey*)iter())){
		name = key->GetName();
		//skip these dirs
		if(name.find("genDeltaTpvGambin") != string::npos) continue;
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
			string name;
			vector<TH1D*> hists;
			//cout << "\ndir name: " << dir->GetName() << endl;

			while((kkey = (TKey*)iiter())){
				//writing stack hist - same method different procs
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = dynamic_cast<TDirectory*>(kkey->ReadObj());
					if(!ddir) continue;
					name = ddir->GetName();
					//method in this subdir needs to match what's given
					if(name.find(method) == string::npos) continue;
					ddir->cd();
					//cout << " ---in ddir " << ddir->GetName() << endl;
					//we're in the directory with hists of one method split by procs
					GetHistsProcs(ddir, procs, hists);
					if(hists.size() > 0){
						//cout << "hists got" << endl;
						//for(auto h : hists) cout << h->GetName() << endl;
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
		
		}
	}
	//cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;
	ofile->Write();
	ofile->Close();
	f->Close();

};



void MethodStackHists(string file, string proc, vector<string>& methods, string oname, string match ="",string year = ""){
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	TFile* f = TFile::Open(file.c_str(),"READ");
	TFile* ofile = TFile::Open(oname.c_str(), "UPDATE");
	TList* list = f->GetListOfKeys();
	TIter iter(list);
	TKey* key;
	string name, xtitle, ytitle;
	//string oname = f->GetName();
	//oname = oname.substr(0,oname.find(".root"));
	//oname = oname+"_formatted.root";
	//TFile* ofile = new TFile(oname.c_str(),"RECREATE");

	string cmslab = "";
	if(proc == "GJets"){
		cmslab = "GJets, HT 600 to Inf";
	}
	if(proc == "QCD"){
		cmslab = "QCD Multijets";
	}
	else if(proc == "JetHT"){
		cmslab = "JetHT, Run F";
	}
	else if(proc == "DoubleEG"){
		cmslab = "DoubleEG, Run F";
	}
	//string cmslab = GetCMSLabel(file);
	//string extra = "";
	//if(file.find("Skim") != string::npos) extra = GetExtraLabel(file);
	//if(!extra.empty()) cmslab += " "+extra;	

	cmslab += " "+year;
	TString th1d("TH1D");
	TString th2d("TH2D");
	TString tdir("TDirectoryFile");

	string methodsname = "";
	for(auto s : methods) methodsname += "_"+s;
	
	while((key = (TKey*)iter())){
		name = key->GetName();
		//skip these dirs
		if(name.find("genDeltaTpvGambin") != string::npos) continue;
		if(key->GetClassName() == tdir){
			//get stack histograms - in directory
			TDirectory* dir = dynamic_cast<TDirectory*>(key->ReadObj());
			double ymin, ymax;
			string ylab, xlab;
			if(!dir) continue;
			name = dir->GetName();
			if(!match.empty() && name.find(match) == string::npos) continue;
			//if(name.find("geoEavg_sigmaDeltaTime_recoGen") == string::npos) continue;
			cout << "\ndir name: " << dir->GetName() << endl;
			dir->cd();
			cout << "---in dir " << dir->GetName() << endl;
			cout << "dir name: " << dir->GetName() << endl;
			//get histograms (stack these)
			vector<TH1D*> hists;
			//if directory is full of directories (stack by method)
			//each dir is full of histograms (stack by process)
			TList* llist = dir->GetListOfKeys();
			TIter iiter(llist);
			TKey* kkey;
			string name;
			
			//writing methodStack hist - same proc different methods
			string dirname = dir->GetName();
			//only do for sigma plots + profiles for now - can remove this later to change
			if(dirname.find("sigma") == string::npos && dirname.find("mean") == string::npos && dirname.find("profile") == string::npos) continue;	
			cout << "about to gethistsprocmethods" << endl;
			GetHistsProcMethods(dir, proc, hists, methods);
			for(auto h : hists) cout << "got hist " << h->GetName() << endl;
			name = dirname+"_"+proc+"_methodStack"+methodsname;
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
			TDRMultiHist(hists, cv, cmslab, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", methodStack);
			cv->Write(); 
			cout << "writing canvas (1D) " << cv->GetName() << endl;
				
		}
		cout << "end for loop - name " << name << endl;
	}
	cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;
	ofile->Write();
	ofile->Close();
	f->Close();

};



void Hist2D(string file, string proc, string method, string oname, string match, string year){
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
	if(proc == "GJets"){
		cmslab = "GJets";
	}
	else if(proc == "QCD"){
		cmslab = "QCD Multijets,";
	}
	else if(proc == "JetHT"){
		cmslab = "JetHT, Run F";
	}
	else if(proc == "DoubleEG"){
		cmslab = "DoubleEG, Run F";
	}
	else cmslab = "process";
	cmslab += " "+year;
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
			string name;
			vector<TH2D*> hists;
			TH2D* hist;
			cout << "\ndir name: " << dir->GetName() << endl;

			while((kkey = (TKey*)iiter())){
				//writing stack hist - same method different procs
				if(kkey->GetClassName() == tdir){
					TDirectory* ddir = dynamic_cast<TDirectory*>(kkey->ReadObj());
					if(!ddir) continue;
					name = ddir->GetName();
					//method in this subdir needs to match what's given
					if(name.find(method) == string::npos) continue;
					ddir->cd();
					cout << " ---in ddir " << ddir->GetName() << endl;
					//we're in the directory with hists of one method split by procs
					GetHists(ddir, proc, hists);
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
		
		}
	}

	ofile->Write();
	ofile->Close();
	f->Close();



}


void ResolutionStackHists(string file, string proc, string method, string oname, string year){
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
	if(proc == "GJets"){
		cmslab = "GJets HT 600 to Inf";
	}
	else if(proc == "QCD"){
		cmslab = "QCD Multijets, HT 500 to 700";
	}
	else if(proc == "JetHTPD"){
		cmslab = "JetHT, Run F";
	}
	else if(proc == "DEGPD"){
		cmslab = "DoubleEG, Run F";
	}
	else cmslab = "process";
	cmslab += " "+year;
	cmslab += ", "+method;
	//string cmslab = GetCMSLabel(file);
	//string extra = "";
	//if(file.find("Skim") != string::npos) extra = GetExtraLabel(file);
	//if(!extra.empty()) cmslab += " "+extra;	


	TString th1d("TH1D");
	TString th2d("TH2D");
	TString tdir("TDirectoryFile");

	//dir names for recogen + dijet resolutions
	string recogenDir = "geoEavg_sigmaDeltaTime_recoGen_stack";
	string dijetDir = "geoAvgEecal_sigmaDeltaTime_dijets_stack";
	string gampvDir = "geoEavg_sigmaDeltaTime_gamPV";
	vector<TH1D*> hists;

	string dirname;
	string ddirname;
	string histname;

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
			if(dirname.find(recogenDir) == string::npos && dirname.find(dijetDir) == string::npos && dirname.find(gampvDir) == string::npos) continue;
			//cout << "\ndir name: " << dir->GetName() << endl;
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
					//cout << " ---in ddir " << ddir->GetName() << endl;
					//we're in the directory with hists of one method split by procs
					histname = ddirname;
					histname = histname.substr(0,histname.find("_procStack"));
					histname += "_"+proc;
					TH1D* hist = (TH1D*)f->Get((dirname+"/"+ddirname+"/"+histname).c_str());
					if(hist) cout << "got histogram " << hist->GetName() << " " << hist->GetEntries() << endl;
					else cout << "hist null " << dirname+"/"+ddirname+"/"+histname << endl;
					if(hist) hists.push_back(hist);
				}
			}	


		}
	}
	double ymin, ymax;
	string ylab, xlab;
	if(hists.size() > 0){
		FindListHistBounds(hists, ymin, ymax);
		if(ymin == 0 && ymax == 0) return;
		string name = "geoEavg_sigmaDeltaTime_dijetRecoGenGamPV_stack_"+proc+"_"+method;
		TCanvas *cv = new TCanvas(name.c_str(), "");
		ofile->cd();
		//draw as tcanvases
		xlab = hists[0]->GetXaxis()->GetTitle();
		ylab = hists[0]->GetYaxis()->GetTitle();
		TDRMultiHist(hists, cv, cmslab, xlab, ylab, ymin-fabs(ymin*0.5), ymax, "", dijetRecoGenStack);
		cout << "writing canvas (1D) " << cv->GetName() << endl;
		cv->Write(); 
	}
	ofile->Write();
	ofile->Close();
	f->Close();
}




void FileStackHists(vector<string>& files, vector<string>& labels, string proc, string method, string oname, string match="",string plottitle="", string year = ""){
	TFile* ofile = TFile::Open(oname.c_str(),"UPDATE");
	string cmslab = "";
	if(proc == "JetHTPD"){
		cmslab = "JetHT, Run F";
	}
	else if(proc == "DEGPD"){
		cmslab = "DoubleEG, Run F";
	}
	else cmslab = "";
	cmslab += " "+year;
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


void HistFormatJets(string file, string file2 = ""){
	string oname = file;
	oname = oname.substr(0,oname.find(".root"));
	if(!file2.empty()){
		string oname2 = file2;
		oname2 = oname2.substr(0,oname2.find(".root"));
		oname2 = oname2.substr(oname2.find("/")+1);
		oname += "_"+oname2;
	}
	oname = oname+"_formatted.root";
	TFile* ofile = new TFile(oname.c_str(),"RECREATE");
	ofile->Close();	

	string year = "";
	if(oname.find("_R17") != string::npos) year = "2017";
	if(oname.find("_R18") != string::npos) year = "2018";

	if(file2.empty()){
		vector<string> med_eAvg = {"median","eAvg"};
		vector<string> chiGam_QCD = {"chiGam","QCD"};
		vector<string> GluGlu_QCD = {"GluGlu","QCD"};
		//PV dijets for median + eAvg for QCD
		MethodStackHists(file, "QCD", med_eAvg, oname, "geoAvgEecal", year);
		//PV dijets for median + eAvg for JetHT
		MethodStackHists(file, "JetHT", med_eAvg, oname, "geoAvgEecal", year);
		MethodStackHists(file, "DoubleEG", med_eAvg, oname, "geoAvgEecal", year);
		//recoGen and gamPV for med + eAvg for QCD
		MethodStackHists(file, "QCD", med_eAvg, oname, "geoEavg", year);

		//same method in legend, only 1 proc in plot label
		vector<string> jetHT_QCD = {"JetHT","QCD"};
		vector<string> jetHT_QCD_DEG = {"JetHT","QCD","DoubleEG"};
		vector<string> DEG_QCD = {"DoubleEG","QCD"};
		//PV dijets for data (JetHT, year) + MC for median
		ProcStackHists(file, jetHT_QCD, "median", oname,"geoAvgEecal");
		////PV dijets for data (DEG) + MC for median
		ProcStackHists(file, DEG_QCD, "median", oname,"geoAvgEecal");
		////PV dijets for data (JetHT) + MC for eAvg
		ProcStackHists(file, jetHT_QCD, "eAvg", oname,"geoAvgEecal");
		////PV dijets for data (DEG) + MC for eAvg
		ProcStackHists(file, DEG_QCD, "eAvg", oname,"geoAvgEecal");
		///PV dijets for data (JetHT) + MC for mmAvg
		ProcStackHists(file, jetHT_QCD, "mmAvg", oname,"geoAvgEecal");
		///PV dijets for data (DEG) + MC for mmAvg
		ProcStackHists(file, DEG_QCD, "mmAvg", oname,"geoAvgEecal");
		//recoGen and gamPV for QCD + GMSB for eAvg
		ProcStackHists(file, GluGlu_QCD, "eAvg", oname, "geoEavg");
		//gamPV for QCD + JetHT
		ProcStackHists(file, jetHT_QCD, "eAvg", oname, "sigmaDeltaTime_gamPV");
		ProcStackHists(file, jetHT_QCD, "median", oname, "sigmaDeltaTime_gamPV");
		ProcStackHists(file, jetHT_QCD, "eAvg", oname, "diffDeltaTime_gamPV");
		ProcStackHists(file, jetHT_QCD, "median", oname, "diffDeltaTime_gamPV");
		//gamPV for QCD + DEG
		ProcStackHists(file, DEG_QCD, "median", oname, "sigmaDeltaTime_gamPV");
		ProcStackHists(file, DEG_QCD, "eAvg", oname, "sigmaDeltaTime_gamPV");
		ProcStackHists(file, DEG_QCD, "eAvg", oname, "diffDeltaTime_gamPV");
		ProcStackHists(file, DEG_QCD, "median", oname, "diffDeltaTime_gamPV");
		
		
		//PV dijets + reocGen + gamPV for QCD eAvg
		ResolutionStackHists(file, "QCD", "eAvg", oname, year);

		//jet properties hists - eavg vs med in data
		Hist2D(file, "JetHT", "eAvg", oname, "jetTime_Energy", year);
		Hist2D(file, "JetHT", "med", oname, "jetTime_Energy", year);
		Hist2D(file, "DoubleEG", "eAvg", oname, "jetTime_Energy", year);
		Hist2D(file, "DoubleEG", "med", oname, "jetTime_Energy", year);
		Hist2D(file, "QCD", "eAvg", oname, "jetTime_Energy", year);
		Hist2D(file, "QCD", "med", oname, "jetTime_Energy", year);

		Hist2D(file, "JetHT", "med", oname, "rhTime_Energy", year);
		Hist2D(file, "DoubleEG", "med", oname, "rhTime_Energy", year);
		Hist2D(file, "QCD", "med", oname, "rhTime_Energy", year);
		
		Hist2D(file, "JetHT", "med", oname, "rhTime_eta", year);
		Hist2D(file, "DoubleEG", "med", oname, "rhTime_eta", year);
		Hist2D(file, "QCD", "med", oname, "rhTime_eta", year);
		
		Hist2D(file, "JetHT", "med", oname, "rhPhi_eta", year);
		Hist2D(file, "DoubleEG", "med", oname, "rhPhi_eta", year);
		Hist2D(file, "QCD", "med", oname, "rhPhi_eta", year);
	
		Hist2D(file, "JetHT", "med", oname, "swCross_rhTime", year);
		Hist2D(file, "DoubleEG", "med", oname, "swCross_rhTime", year);
		Hist2D(file, "QCD", "med", oname,   "swCross_rhTime", year);
		
		Hist2D(file, "JetHT", "med", oname, "swCross_rhEnergy", year);
		Hist2D(file, "DoubleEG", "med", oname, "swCross_rhEnergy", year);
		Hist2D(file, "QCD", "med", oname,   "swCross_rhEnergy", year);
		
		Hist2D(file, "JetHT", "med", oname, "kWeird", year);
		Hist2D(file, "DoubleEG", "med", oname, "kWeird", year);
		Hist2D(file, "QCD", "med", oname,   "kWeird", year);

		Hist2D(file, "JetHT", "med", oname, "rhEta_rhPhi", year);
		Hist2D(file, "DoubleEG", "med", oname, "rhEta_rhPhi", year);
		Hist2D(file, "QCD", "med", oname,   "rhEta_rhPhi", year);
		
		Hist2D(file, "JetHT", "med", oname, "rhEovP_dRtrack", year);
		Hist2D(file, "DoubleEG", "med", oname, "rhEovP_dRtrack", year);
		Hist2D(file, "QCD", "med", oname,   "rhEovP_dRtrack", year);
		
		Hist2D(file, "JetHT", "med", oname, "rhTime_rhEta", year);
		Hist2D(file, "DoubleEG", "med", oname, "rhTime_rhEta", year);
		Hist2D(file, "QCD", "med", oname,   "rhTime_rhEta", year);
	
		Hist2D(file, "JetHT",    "med", oname, "Neighbors", year);
		Hist2D(file, "DoubleEG", "med", oname, "Neighbors", year);
		
		ProcStackHists(file, jetHT_QCD_DEG, "median", oname, "swCross_rhTime");
		ProcStackHists(file, jetHT_QCD_DEG, "median", oname, "LHratio");
		//ProcStackHists(file, jetHT_QCD, "eMax", oname, "jetPt");
		//ProcStackHists(file, jetHT_QCD, "eMax", oname, "jetEta");
		//ProcStackHists(file, jetHT_QCD, "eMax", oname, "jetPhi");
		//ProcStackHists(file, jetHT_QCD, "eMax", oname, "jetNrhs");
		//ProcStackHists(file, jetHT_QCD, "mmAvg", oname, "jetNSubclusters");
		//ProcStackHists(file, jetHT_QCD, "eAvg", oname, "jetTime");
		//ProcStackHists(file, jetHT_QCD, "med", oname, "jetTime");


	}
	//pre + post calibration for JetHT
	//pre + post calibraiton for DoubleEG
	vector<string> files;
	files.push_back(file);
	if(!file2.empty()){
		files.push_back(file2);
		vector<string> labels = {"postCalib","preCalib"};
		//for(int i = 0; i < labels.size(); i++) cout << "file " << files[i] << " has label " << labels[i] << endl;
		//FileStackHists(files,labels,"JetHTPD","eAvg",oname,"profile_geoAvgEecal","PrePostCalibration");
		//FileStackHists(files,labels,"DoubleEGPD","eAvg",oname,"geoAvgEecal_sigmaDeltaTime_dijets","PrePostCalibration");
		labels = {"wSpikes","woSpikes"};
		for(int i = 0; i < labels.size(); i++) cout << "file " << files[i] << " has label " << labels[i] << endl;
		for(int i = 0; i < 6; i++)
			FileStackHists(files,labels,"JetHTPD","eAvg",oname,"profile_geoAvgEecal_diffDeltaTime_dijets_bin"+std::to_string(i+1),"wWoSpikes");
		FileStackHists(files,labels,"JetHTPD","eAvg",oname,"geoAvgEecal_sigmaDeltaTime_dijets","wWoSpikes");
		//FileStackHists(files,labels,"DoubleEGPD","eAvg",oname,"geoAvgEecal_sigmaDeltaTime_dijets","wWoSpikes");
	}
	cout << "Wrote formatted canvases to: " << ofile->GetName() << endl;


};
