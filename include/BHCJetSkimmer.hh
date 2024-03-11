#ifndef BHCJETSKIMMER_HH
#define BHCJETSKIMMER_HH
#include "JetSkimmer.hh"

class BHCJetSkimmer : public JetSkimmer{
	public:
		BHCJetSkimmer(){
			_evti = 0;
			_evtj = 0;
			_gev = 1./10.;
		}

		virtual ~BHCJetSkimmer(){ }

		BHCJetSkimmer(TFile* file) : JetSkimmer(file){
			/*
			_prod = new JetProducer(file);
			_prod->SetIsoCut();
		

			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/bhcjet_skims.root";
			_gev = 1./10.;
				
			objE_clusterE->SetTitle("jetE_clusterE");
			objE_clusterE->SetName("jetE_clusterE");
			*/
	
			nrhs_comptime->SetName("nrhs_comptime");
			nrhs_comptime->SetTitle("nrhs_comptime");
			graphs.push_back(nrhs_comptime);

		}
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else return; 
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
				//PVdeltaT_jet_mmAvg_pred->Fill(t1);	
				for(int j = i+1; j < njets; j++){
					if(trees[i] == nullptr) continue;
					CalcMMAvgPhiTime(trees[j]->model, phi2, t2);
					//dphi within [pi-0.1,pi+0.1]
					dphi = fabs(phi1 - phi2);
					if(dphi < pi-0.1 || dphi > pi+0.1) continue;
					//median time
					//energy-weighted average
					//mm average over subclusters
					//PVdeltaT_jet_mmAvg_pred->Fill(t1 - t2);	
				}
			}

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
		void WriteOutput(TFile* ofile){
			WriteHists(ofile);
			string name;
			ofile->cd();
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

		//predicted jet plots
		TH1D* nClusters = new TH1D("nClusters","nClusters",20,0,20);
		TH1D* PVtime_pred = new TH1D("PVtime_pred","PVtime_pred",100,-10,10);	
		//difference in tPV between two back-to-back jets
		//previous time definitions
		TH1D* deltaT_jet_pred = new TH1D("PVdeltaT_jet_pred","PVdeltaT_jet_pred",100,-10,10);	
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();

	private:
		vector<TGraph*> graphs;
		vector<Jet> _phos; //photons for event
		enum Strategy{
			//Delauney strategy - NlnN time - for 2pi cylinder
			NlnN = 0,
			//traditional strategy - N^2 time
			N2 = 1,
		};
		//clustering strategy - N^2 or NlnN
		Strategy _strategy;
};
#endif
