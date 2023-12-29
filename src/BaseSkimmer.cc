#include "BaseSkimmer.hh"
#include "TF1.h"

void BaseSkimmer::Profile2DHist(TH2D* inhist, TH1D* outhist, vector<TH1D*>& profs, double range){
	int nbins = inhist->GetNbinsX();
	profs.clear();
	//skip overflow + underflow bins
	for(int i = 1; i < nbins; i++){
		TH1D* phist = (TH1D*)inhist->ProjectionY("tmp",i,i);
		phist->SetTitle(("profile_bin"+std::to_string(i)).c_str());	
		phist->SetName(("profile_bin"+std::to_string(i)).c_str());
		phist->GetXaxis()->SetTitle("diffDeltaT_recoGen");	
		profs.push_back(phist);
		//get values for param init
		double mean = phist->GetMean();
		double stddev = phist->GetStdDev();
		double norm = phist->GetBinContent(phist->GetMaximumBin());
		double high = mean + range*stddev;
		double low = mean - range*stddev;
		//check that initial parameter values are ok
		if( stddev > 0.0 && norm > 0.){
			TF1* fit = new TF1("fit","gaus",low,high);
			fit->SetParameter(0,norm);
			fit->SetParameter(1,mean);
			fit->SetParameter(2,stddev);
			phist->Fit(fit->GetName(),"RBQ0");
			
			double fit_stddev = fit->GetParameter(2);
			double fit_stddev_err = fit->GetParError(2);
			//set new contents
			outhist->SetBinContent(i, fit_stddev);
			outhist->SetBinError(i, fit_stddev_err);
		cout << "bin " << i << ": stddev " << stddev << " fit_stddev " << fit_stddev << " err: " << fit_stddev_err << " norm " << norm << " phist integral: " << phist->Integral() << endl; 	
			delete fit;
		}
	}
		
}
