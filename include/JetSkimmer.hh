#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include <TFile.h>
#include <TGraph.h>
#include "TSystem.h"
#include "BaseTree.hh"

using node = BaseTree::node;
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
		TH1D* tPV_res_avg = new TH1D("tPV_res_avg","tPV_res_avg",100,-10.,10.);
		TH1D* tPV_res_lead = new TH1D("tPV_res_lead","tPV_res_lead",100,-10.,10.);
		TH2D* e_nRhs = new TH2D("e_nRhs","e_nRhs",100,0,500,100,0,100);
		TH1D* t_rhs = new TH1D("t_rhs","t_rhs",100,-30,30); 
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* comptime_nrhs = new TGraph();
		

		//this is for one jet
		//all hists referenced here are in hists1D
		void FillModelHists(BasePDFMixture* model, double transf = 1.){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es, npts_unwt;
			vector<Matrix> eigenvecs;
			double theta, phi, r, id, npts, E_k;

			int nclusters = model->GetNClusters();
			nSubClusters->Fill(nclusters);
			
			model->GetAvgVarWeights(avg_Es);
			model->GetNormsUnwt(npts_unwt);
		
			nClusters->Fill((double)nclusters);
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				E_k = avg_Es[k]*transf*npts_unwt[k];

				params = model->GetPriorParameters(k);
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
				e_tot->Fill(E_k);
			

			}
			//tPV->Fill(time/time_denom);
		}
		


		//find back to back jets
		void FillPVHists(vector<node*> tree){
			int njets = (int)tree.size(); 
			double pi = acos(-1);

			double dtime, dphi, dr, phi1, t1, phi2, t2;
			//find pairs of jets to calculate resolution	
			//need to be back to back
			//time of subclusters is measured as center
			for(int i = 0; i < njets; i++){
				for(int j = i; j < njets; j++){
					//averaged over subclusters, weighted by mixing coeff
					CalcAvg(tree[i],phi1,t1); CalcAvg(tree[j],phi2,t2);
					if(fabs(phi1-phi2) < pi-0.1 && fabs(phi1-phi2) > pi+0.1) tPV_res_avg->Fill(t1-t2);
					//lead subcluster
					CalcLead(tree[i],phi1,t1); CalcLead(tree[j],phi2,t2);
					if(fabs(phi1-phi2) < pi-0.1 && fabs(phi1-phi2) > pi+0.1) tPV_res_lead->Fill(t1-t2);
			
				}
			}



		}


		void CalcLead(node* nnode, double& phi, double& t){
			BasePDFMixture* model = nnode->model;
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			vector<int> idxs;
			model->SortIdxs(idxs);

			map<string, Matrix> params = model->GetPriorParameters(idxs[kmax-1]);
			phi = params["pi"].at(0,0)*params["mean"].at(1,0);
			t = params["pi"].at(0,0)*params["mean"].at(2,0);
		}
		



		void CalcAvg(node* nnode, double& phi, double& t){
			BasePDFMixture* model = nnode->model;
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double ws;
			double pi;
			map<string, Matrix> params;
			for(int k = 0; k < kmax; k++){
				params = model->GetPriorParameters(k);
				phi += params["pi"].at(0,0)*params["mean"].at(1,0);
				t += params["pi"].at(0,0)*params["mean"].at(2,0);
				ws += params["pi"].at(0,0);
			}
			phi /= ws;
			t /= ws; 
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
			TCanvas* cv = new TCanvas("rh_time", "");
			TDRHist(t_rhs, cv, "rh_time", "rh_time", "a.u.");
			cv->Write();
			ofile->Close();

		}




};
#endif
