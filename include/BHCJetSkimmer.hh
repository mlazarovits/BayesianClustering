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
