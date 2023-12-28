#include "BaseSkimmer.hh"
#include "TF1.h"

void BaseSkimmer::Profile2DHist(TH2D* inhist, TH1D* prof, double range = 0.2){
	int nbins = inhist->GetNbinsX();
	//skip overflow + underflow bins
	for(int i = 1; i < nbins; i++){
		TH1D* phist = (TH1D*)inhist->ProjectionY("tmp",i,i);
		//get values for param init
		double mean = phist->GetMean();
		double stddev = phist->GetStdDev();
		double norm = phist->GetBinContent(phist->GetMaximumBin());
		double high = mean + range*stddev;
		double low = mean - range*stddev;
		//check that initial parameter values are ok
		if( stddev > 0.0 && norm > 1){
			TF1* fit = new TF1("fit","gaus",low,high);
			fit->SetParameter(0,norm);
			fit->SetParameter(1,mean);
			fit->SetParameter(2,stddev);
			phist->Fit(fit->GetName(),"RBQ0");
			
			double fit_stddev = fit->GetParameter(2);
			double fit_stddev_err = fit->GetParError(2);
			//set new contents
			prof->SetBinContent(i, fit_stddev);
			prof->SetBinError(i, fit_stddev_err);
		
			delete fit;
		}
	}
		
}
