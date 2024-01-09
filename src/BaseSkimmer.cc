#include "BaseSkimmer.hh"
#include "TF1.h"


//consider changing outplot to TGraph (hist is ok for now)
void BaseSkimmer::Profile2DHist(TH2D* inhist, TH1D* outhist, vector<TH1D*>& profs){
	int nbins = inhist->GetNbinsX();
	profs.clear();
	string profilename = "";
	//skip overflow + underflow bins
	for(int i = 1; i < nbins; i++){
		TH1D* phist = (TH1D*)inhist->ProjectionY("tmp",i,i);
		profilename = inhist->GetName();
		profilename = profilename.substr(0,profilename.rfind("_"));
		profilename = profilename.substr(0,profilename.rfind("_"));
		phist->SetTitle(("profile_"+profilename+"_bin"+std::to_string(i)).c_str());	
		phist->SetName(("profile_"+profilename+"_bin"+std::to_string(i)).c_str());
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



