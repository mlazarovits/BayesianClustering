#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include <TFile.h>
#include "JetProducer.hh"
#include "TSystem.h"

using plotCat = BaseSkimmer::plotCat;
class JetSkimmer : public BaseSkimmer{
	public:
		JetSkimmer();
		virtual ~JetSkimmer();

		//get rechits from file to cluster
		JetSkimmer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		
		void CleaningSkim(){ };
		void Skim();

		//jet specific quantities
		TH1D* nClusters = new TH1D("nClusters","nClusters",20,0,20);
		TH1D* nTrueJets = new TH1D("nTrueJets","nTrueJets",20,0,20);
		//tPV = tJet - dRH/c (clock offset) + dPV/c (TOF - time to travel offset)
		TH1D* tPV = new TH1D("tPV","tPV",100,-1.,1.);
		//difference in tPV between two back-to-back jets
		TH1D* tPV_res = new TH1D("tPV_res",100,-10.,10.);

		//all hists referenced here are in hists1D
		void FillTotalHists(BasePDFMixture* model){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es;
			vector<Matrix> eigenvecs;
			int nclusters = model->GetNClusters();
			nSubClusters->Fill(nclusters);
			//e_nSubClusters->Fill(_base->Photon_energy->at(p), nclusters);
			model->GetAvgVarWeights(avg_Es);		
			double theta, phi, r, id, npts;
			//for energy weights since w_n = E_n*(N/sum_n E_n) s.t. sum_n w_n = N
			npts = (double)gmm->GetData()->GetNPoints();	
			nClusters->Fill((double)nclusters);
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				params = model->GetParameters(k);
				eta_center->Fill(params["mean"].at(0,0));
				phi_center->Fill(params["mean"].at(1,0));
				time_center->Fill(params["mean"].at(2,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//largest eigenvalue is last
				//phi/eta
				slope_space->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(0,0));
        			//eta/time
				slope_etaT->Fill(eigenvecs[2].at(0,0)/eigenvecs[2].at(2,0));
				//phi/time
				slope_phiT->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(2,0));
				//polar angle
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				r = sqrt(eigenvecs[2].at(0,0)*eigenvecs[2].at(0,0) + eigenvecs[2].at(1,0)*eigenvecs[2].at(1,0) + eigenvecs[2].at(2,0)*eigenvecs[2].at(2,0));
				theta = acos( eigenvecs[2].at(2,0) / r );
				polar_ang->Fill(theta);
				//azimuthal angle
				//phi = arctan(y/x)
				phi = atan2(eigenvecs[2].at(1,0) , eigenvecs[2].at(0,0));
				azimuth_ang->Fill(phi);
				
				//average cluster energy
				e_avg->Fill(avg_Es[k]);
				
				

			}
		}
		void WriteHists(TFile* ofile){
			string name;

			ofile->cd();
			for(int i = 0; i < (int)hists1D.size(); i++){
				TCanvas* cv = new TCanvas(name.c_str(), "");
				TDRHist(hists1D[i], cv, name, name, "a.u.");
				//write cv to file			
			//	cv->SaveAs((fname+"/"+name+".pdf").c_str());
				cv->Write();

			}
			ofile->Close();

		}





};
#endif
