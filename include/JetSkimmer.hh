#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "TFile.h"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"

class JetSkimmer : public BaseSkimmer{
	public:
		JetSkimmer();
		virtual ~JetSkimmer();

		//get rechits from file to cluster
		JetSkimmer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		
		void CleaningSkim(){ };
		void Skim();


		//all hists referenced here are in hists1D
		void FillTotalHists(BasePDFMixture* model){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es;
			vector<Matrix> eigenvecs;
			int nclusters = model->GetNClusters();
			nSubClusters->Fill(nclusters);
			//e_nSubClusters->Fill(_base->Photon_energy->at(p), nclusters);
			//model->GetAvgWeights(avg_Es);		
			double theta, phi, r, id;
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
				phi = atan(eigenvecs[2].at(1,0) / eigenvecs[2].at(0,0));
				azimuth_ang->Fill(phi);
				
				//average cluster energy
				e_avg->Fill(avg_Es[k]);
			}
		}





};
#endif
