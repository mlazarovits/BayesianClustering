#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include "BayesCluster.hh"
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
		
		
		void Skim();
		void GMMOnly(){ _mmonly = true; }


		//true jet hists
		TH1D* nTrueJets = new TH1D("nTrueJets","nTrueJets",20,0,20);
		TH1D* rhTime = new TH1D("rhTime","rhTime",100,-30,30); 
		TH1D* pT_twoHardJets = new TH1D("pT_twoHardJets","pT Two Hardest Jets",100,0,3000);
		TH2D* e_nRhs = new TH2D("e_nRhs","e_nRhs",100,0,500,100,0,100);
		

		//tPV = tJet - dRH/c (clock offset) + dPV/c (TOF - time to travel offset)
		TH1D* PVtime_median = new TH1D("PVtime_median","PVtime_median",100,-10,10);	
		TH1D* PVtime_eAvg = new TH1D("PVtime_eAvg","PVtime_eAvg",100,-10,10);	
		//difference in tPV between two back-to-back jets
		//previous time definitions
		TH1D* PVtimeDiff_median = new TH1D("PVtimeDiff_median","PVtimeDiff_median",100,-20,20);	
		TH1D* PVtimeDiff_eAvg = new TH1D("PVtimeDiff_eAvg","PVtimeDiff_eAvg",100,-20,20);	

		//mm only jet plots
		//mm only time (from true jets)
		TH1D* PVtime_mmAvg = new TH1D("PVtime_mmAvg", "PVtime_mmAvg",100,-10,10);	
		//mm only time (from true jets)
		TH1D* PVtimeDiff_mmAvg = new TH1D("PVtimeDiff_mmAvg", "PVtimeDiff_mmAvg",100,-20,20);	
		TH1D* nSubClusters_mm = new TH1D("nSubClusters_mm","nSubClusters_mm",20,0,20);


		//comparing predicted jets + true jets
		//TH2D* nSubClusters_nConstituents = new TH2D("nSubClusters_nConstituents", "nSubClusters_nConstituents",50,0,20,50,0,20);

		//predicted jet plots
		TH1D* nClusters = new TH1D("nClusters","nClusters",20,0,20);
		TH1D* PVtime_median_pred = new TH1D("PVtime_median_pred","PVtime_median_pred",100,-10,10);	
		TH1D* PVtime_eAvg_pred = new TH1D("PVtime_eAvg_pred","PVtime_eAvg_pred",100,-10,10);	
		TH1D* PVtime_mmAvg_pred = new TH1D("PVtime_mmAvg_pred", "PVtime_mmAvg_pred",100,-10,10);	
		//difference in tPV between two back-to-back jets
		//previous time definitions
		TH1D* PVtimeDiff_median_pred = new TH1D("PVtimeDiff_median_pred","PVtimeDiff_median_pred",100,-10,10);	
		TH1D* PVtimeDiff_eAvg_pred = new TH1D("PVtimeDiff_eAvg_pred","PVtimeDiff_eAvg_pred",100,-10,10);	
		//mm only time (from true jets)
		TH1D* PVtimeDiff_mmAvg_pred = new TH1D("PVtimeDiff_mmAvg_pred", "PVtimeDiff_mmAvg_pred",100,-10,10);	
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();
	

		void FillTrueJetHists(const vector<Jet>& jets){
			int njets = _base->Jet_energy->size();	
			nTrueJets->Fill((double)njets);
		
			for(int j = 0; j < jets.size(); j++){
				e_nRhs->Fill(jets[j].E(), jets[j].GetNRecHits());
				
			}

			FillPVHists_TrueJets(jets);	
		}


		void FillMMOnlyJetHists(const vector<Jet>& jets, double alpha = 0.1, double emAlpha = 0.5){
			vector<GaussianMixture*> models;
			Matrix smear = Matrix(3,3);
			double dphi = acos(-1)/360.; //1 degree in radians
			double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
			//diagonal matrix
			smear.SetEntry(deta*deta,0,0);
			smear.SetEntry(dphi*dphi,1,1);
			smear.SetEntry(1.,2,2); //energy dependent smear in time, set in BayesCluster	
			for(int i = 0; i < jets.size(); i++){
				vector<JetPoint> rhs = jets[i].GetJetPoints();
				vector<Jet> rhs_jet;
				for(int r = 0; r < rhs.size(); r++){
					rhs_jet.push_back( Jet(rhs[r]) );
					rhs_jet[r].SetVertex(jets[i].GetVertex());
				}
				BayesCluster* algo = new BayesCluster(rhs_jet);
				algo->SetDataSmear(smear);
				//set time resolution smearing
				algo->SetTimeResSmear(0.2, 0.3*_gev);
				algo->SetThresh(1.);
				algo->SetAlpha(alpha);
				algo->SetSubclusterAlpha(emAlpha);
				algo->SetVerbosity(0);
				GaussianMixture* gmm = algo->SubCluster();
				nSubClusters_mm->Fill(gmm->GetNClusters());
				models.push_back(gmm);
			}		
			//need vector of all jets/models
			//to do back-to-back matching for resolution
			FillPVHists_MMOnlyJets(jets, models);	
			
	
		}	


		void FillPredJetHists(const vector<node*>& trees){
			int njets = 0;
			vector<node*> cleaned_trees;
			for(int i = 0; i < trees.size(); i++){
				if(trees[i] == nullptr) continue;
				//check for mirrored point - would be double counted
				if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0) continue;	
				FillModelHists(trees[i]->model);
				njets++;
				cleaned_trees.push_back(trees[i]);
			}
			nClusters->Fill(njets);
			FillPVHists_PredJets(cleaned_trees);
		}
	

		//this is for one jet
		//all hists referenced here are in hists1D
		void FillModelHists(BasePDFMixture* model){
			map<string, Matrix> params;
			vector<double> eigenvals, norms;
			vector<Matrix> eigenvecs;
			double theta, phi, r, id, npts, E_k;

			int nclusters = model->GetNClusters();
			nSubClusters->Fill(nclusters);
			model->GetNorms(norms);
		
			nClusters->Fill((double)nclusters);
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				E_k = norms[k]/_gev;

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
				clusterE->Fill(E_k);
			

			}
		}
		


		//find back to back jets
		void FillPVHists_PredJets(const vector<node*>& trees){
			int njets = (int)trees.size(); 
			double pi = acos(-1);

			double dtime, dphi, dr, phi1, t1, phi2, t2;
			//find pairs of jets to calculate resolution	
			//need to be back to back
			//time of subclusters is measured as center
			for(int i = 0; i < njets; i++){
				if(trees[i] == nullptr) continue;
				CalcMMAvgPhiTime(trees[i]->model, phi1, t1);
				PVtimeDiff_mmAvg_pred->Fill(t1);	
				for(int j = i+1; j < njets; j++){
					if(trees[i] == nullptr) continue;
					CalcMMAvgPhiTime(trees[j]->model, phi2, t2);
					//dphi within [pi-0.1,pi+0.1]
					dphi = fabs(phi1 - phi2);
					if(dphi < pi-0.1 || dphi > pi+0.1) continue;
					//median time
					//energy-weighted average
					//mm average over subclusters
					PVtimeDiff_mmAvg_pred->Fill(t1 - t2);	
				}
			}

		}
		//find back to back jets
		void FillPVHists_TrueJets(const vector<Jet>& jets){
			int njets = (int)jets.size(); 
			double pi = acos(-1);
			double dphi, phi1, phi2;
			//find pairs of jets to calculate resolution	
			//need to be back to back
			//time of subclusters is measured as center
			for(int i = 0; i < njets; i++){
				phi1 = jets[i].phi_02pi();
				//fill pv time histograms
				//median time
				PVtime_median->Fill( CalcMedianTime(jets[i]));
				//energy-weighted average
				PVtime_eAvg->Fill( CalcEAvgTime(jets[i]));	
				for(int j = i+1; j < njets; j++){
					phi2 = jets[j].phi_02pi();
					//dphi within [pi-0.1,pi+0.1]
					dphi = fabs(phi1 - phi2);
					if(dphi < pi-0.1 || dphi > pi+0.1) continue;
					//fill pv time difference histograms
					//median time
					PVtimeDiff_median->Fill( CalcMedianTime(jets[i]) - CalcMedianTime(jets[j]) );
					//energy-weighted average
					PVtimeDiff_eAvg->Fill( CalcEAvgTime(jets[i]) - CalcEAvgTime(jets[j]) );	
				}
			}

		}
		//find back to back jets
		void FillPVHists_MMOnlyJets(const vector<Jet>& jets, const vector<GaussianMixture*>& models){
			int njets = (int)jets.size();
			if(njets != models.size()){ cout << njets << " jets and " << models.size() << " models" << endl; return;} 
			double pi = acos(-1);
			double time;
			double dphi, phi1, phi2;
			Jet jet_i, jet_j;
			GaussianMixture model_i, model_j;
			map<double,double> pt_time;
			//find hardest two jets to calculate resolution	
			//need to be back to back
			//time of subclusters is measured as center
			for(int i = 0; i < njets; i++){
				//fill pv time histograms
				//mm average over subclusters
				time = CalcMMAvgTime(models[i]);
				pt_time[jets[i].pt()] = time;
				PVtime_mmAvg->Fill( time );
			}
			map<double,double>::reverse_iterator it = pt_time.rbegin();	
			double dtime = it->second;
			pT_twoHardJets->Fill(it->first);
			it++;
			pT_twoHardJets->Fill(it->first);
			dtime -= it->second;	
			//not needed for GMSB
			////dphi within [pi-0.1,pi+0.1]
			//phi1 = jets[jet1].phi_02pi();
			//phi2 = jets[jet2].phi_02pi();
			//dphi = fabs(phi1 - phi2);
			//if(dphi < pi-0.1 || dphi > pi+0.1) return;
			//fill pv time difference histograms
			//mm average over subclusters
			PVtimeDiff_mmAvg->Fill( dtime );

		}


		double CalcMedianTime(const Jet& j){
			vector<double> times;
			vector<JetPoint> rhs = j.GetJetPoints();
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++)
				times.push_back(rhs[i].t());
			//even - return average of two median times
			if(nrhs % 2 == 0)
				return (times[int(double(nrhs)/2.)] + times[int(double(nrhs)/2.)-1]/2.);
			//odd - return median
			else
				return times[int(double(nrhs)/2.)];

		}		


		double CalcEAvgTime(const Jet& j){
			double t, norm;
			vector<JetPoint> rhs = j.GetJetPoints();
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++){
				norm += rhs[i].E();
				t += rhs[i].E()*rhs[i].t();
			}
			return t/norm;
		}




		double CalcMMAvgTime(BasePDFMixture* model){
			int kmax = model->GetNClusters();
			double t = 0;
			double pi, norm;
			map<string, Matrix> params;
			for(int k = 0; k < kmax; k++){
				params = model->GetPriorParameters(k);
				t += params["pi"].at(0,0)*params["mean"].at(2,0);
				norm += params["pi"].at(0,0);
			}
			return t/norm; 
		}

		void CalcMMAvgPhiTime(BasePDFMixture* model, double& phi, double& t){
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double pi, ws;
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
		
		void CalcEAvgPhiTime(BasePDFMixture* model, double& phi, double& t){
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double pi, ws;
			map<string, Matrix> params;
			for(int i = 0; i < model->GetData()->GetNPoints(); i++){
				phi += model->GetData()->at(i).w()/_gev*model->GetData()->at(i).Value(1);
				t += model->GetData()->at(i).w()/_gev*model->GetData()->at(i).Value(2);
				ws += model->GetData()->at(i).w()/_gev;
			}
			phi /= ws;
			t /= ws; 
		}



		void WriteHists(TFile* ofile){
			string name;

			ofile->cd();
			//for condor skims, histograms need to be written because these
			//can be hadded together (TCanvases can't)
			for(int i = 0; i < (int)_hists1D.size(); i++){
				//name = hists1D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRHist(hists1D[i], cv, name, name, "a.u.");
				//write cv to file			
			//	cv->SaveAs((fname+"/"+name+".pdf").c_str());
				//cv->Write();
				_hists1D[i]->Write();
			}
			for(int i = 0; i < (int)_hists2D.size(); i++){
				//name = hists2D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDR2DHist(hists2D[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				_hists2D[i]->Write();
			}
			for(int i = 0; i < (int)graphs.size(); i++){
				//name = graphs[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRGraph(graphs[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				graphs[i]->Write();
			}
			ofile->Close();

		}


		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else return; 
		}


		private:
			bool _mmonly;
			enum Strategy{
				//Delauney strategy - NlnN time - for 2pi cylinder
				NlnN = 0,
				//traditional strategy - N^2 time
				N2 = 1
			};
		
			//clustering strategy - N^2 or NlnN
			Strategy _strategy;

			vector<TGraph*> graphs;
	
		void TDRGraph(TGraph* gr, TCanvas* &can, string plot_title, string xtit, string ytit){
			can->cd();
			can->SetGridx(1);
			can->SetGridy(1);
			gr->SetTitle("");
			gr->UseCurrentStyle();
			gr->GetXaxis()->CenterTitle(true);
			gr->GetXaxis()->SetTitle(xtit.c_str());
			gr->GetYaxis()->CenterTitle(true);
			gr->GetYaxis()->SetTitle(ytit.c_str());
			gr->Draw();
			
			string lat_cms = "#bf{CMS} #it{WIP} "+_cms_label;
			TLatex lat;
			lat.SetNDC();
			lat.SetTextSize(0.04);
			lat.SetTextFont(42);
			lat.DrawLatex(0.02,0.92,lat_cms.c_str());
			//lat.DrawLatex(0.4,0.92,plot_title.c_str());

			return;
		}
};
#endif
