#ifndef BHCJETSKIMMER_HH
#define BHCJETSKIMMER_HH
#include "JetSimProducer.hh"
#include "BaseSkimmer.hh"
#include "BaseTree.hh"
#include "TGraph.h"

using node = BaseTree::node;
using procCat = BaseSkimmer::procCat;

class BHCJetSkimmer{
	public:
		BHCJetSkimmer(){
			_evti = 0;
			_evtj = 0;
			_gev = 1./10.;
			_oname = "";
		}

		virtual ~BHCJetSkimmer(){ }

		BHCJetSkimmer(TFile* file){
			_prod = new JetSimProducer(file);

			_oname = "";
			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_gev = 1./10.;
				
	
			nrhs_comptime->SetName("nrhs_comptime");
			nrhs_comptime->SetTitle("nrhs_comptime");
			graphs.push_back(nrhs_comptime);

			_hists1D.push_back(nClusters);
			_hists1D.push_back(nSubClusters);	
			_hists1D.push_back(predJet_Energy);
			_hists1D.push_back(predJet_EtaCenter);
			_hists1D.push_back(predJet_PhiCenter);
			_hists1D.push_back(predJet_TimeCenter);

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
			for(int p = 0; p < _procCats.size(); p++){
				cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				for(int i = 0; i < trees.size(); i++){
					if(trees[i] == nullptr) continue;
					//check for mirrored point - would be double counted
					if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0) continue;	
					FillModelHists(trees[i]->model, p);
					njets++;
					cleaned_trees.push_back(trees[i]);
				}
				_procCats[p].hists1D[0][0]->Fill(njets);
			}
		}
		
		//all hists referenced here are in hists1D
		void FillModelHists(BasePDFMixture* model, int p){
			map<string, Matrix> params;
			vector<double> eigenvals, norms;
			vector<Matrix> eigenvecs;
			double theta, phi, r, id, npts, E_k;

			int nclusters = model->GetNClusters();
			
			nSubClusters->Fill(nclusters);
			model->GetNorms(norms);
		
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				E_k = norms[k]/_gev;

				params = model->GetPriorParameters(k);
				_procCats[p].hists1D[0][3]->Fill(params["mean"].at(0,0));
				_procCats[p].hists1D[0][4]->Fill(params["mean"].at(1,0));
				_procCats[p].hists1D[0][5]->Fill(params["mean"].at(2,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//total cluster energy
			}
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
		void MakeProcCats(string sample, bool leadsep = true){
			//total
			procCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			_procCats.push_back(tot);	
			
			if(sample.find("GMSB") != string::npos){
				//notSunm
				procCat notSunm(_hists1D, _hists2D, "notSunm","notSunm", leadsep);
				//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
				notSunm.ids = {97, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
				_procCats.push_back(notSunm);
				
				//signal
				procCat sig(_hists1D, _hists2D, "chiGam","#Chi^{0} #rightarrow #gamma", leadsep);
				sig.ids = {22};
				_procCats.push_back(sig);
			}
			else if(sample.find("JetHT") != string::npos){
				//data
				procCat jetht(_hists1D, _hists2D, "JetHT", "JetHT", leadsep);
				jetht.ids = {-999};
				_procCats.push_back(jetht);
			}
			else if(sample.find("GJets") != string::npos){
				//data
				procCat gjets(_hists1D, _hists2D, "GJets", "GJets", leadsep);
				gjets.ids = {-999};
				_procCats.push_back(gjets);
			}
			else if(sample.find("ttbar") != string::npos){
				procCat ttbar(_hists1D, _hists2D, "ttbar", "t#bar{t}",leadsep);
				ttbar.ids = {-999};
				_procCats.push_back(ttbar);
			}
			else if(sample.find("QCD") != string::npos){
				procCat qcd(_hists1D, _hists2D, "QCD", "QCD multi-jets",leadsep);
				qcd.ids = {-999};
				_procCats.push_back(qcd);
			}
			else return;

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

			ofile->Close();

		}
		
		void WriteTimeRecoCatStack(TFile* ofile){
			ofile->cd();
			string name, dirname, histname;
			//write 1D hists
			//variables
			int nhists = _procCats[0].hists1D[0].size();
			for(int i = 0; i < nhists; i++){
				name = _procCats[0].hists1D[0][i]->GetName();
				//dirname = name.substr(0,name.rfind("_"+trs[0].methodName));
			////cout << "i: " << i << " name " << name << " making dir " << dirname+"_stack" << endl;
				TDirectory* dir = ofile->mkdir((name+"_stack").c_str());
				dir->cd();
				//make process breakdown directory
				TDirectory *dir2 = dir->mkdir((dirname+"_procStack").c_str());
				//cout << "  making dir " << dir2->GetName() << endl;
				dir2->cd();
				for(int p = 0; p < _procCats.size(); p++){
				//loop over processes
			//		cout << "    proc " << trs[j].procCats[p].plotName << " hist " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " entries " << trs[j].procCats[p].hists1D[0][i]->GetEntries() << endl;			
					//histname = procCats[p].hists1D[0][i]->GetName();
					if(_procCats[p].hists1D[0][i] == nullptr) continue;
					if(_procCats[p].hists1D[0][i]->GetEntries() == 0){ continue; }
					//cout << "  n hists " << trs[j].procCats[0].hists1D[0].size() << endl;
					//cout << "writing " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;;
					_procCats[p].hists1D[0][i]->Write();

				} 
			}
			//write 2D hists
			
		}


		
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();
		
		//predicted jet plots
		TH1D* nClusters = new TH1D("nPredJets","nPredJets",20,0,20);
		TH1D* nSubClusters = new TH1D("nPredSubClusters","nPredSubClusters",50,0,50);
		TH1D* predJet_Energy = new TH1D("predJet_Energy","predJet_Energy",20,0,1000);
		TH1D* predJet_TimeCenter = new TH1D("predJet_timeCenter","predJet_timeCenter",50,-20,20);
		TH1D* predJet_EtaCenter = new TH1D("predJet_etaCenter","predJet_etaCenter",50,-1.6,1.6);
		TH1D* predJet_PhiCenter = new TH1D("predJet_phiCenter","predJet_phiCenter",50,-3.2,3.2);

		void SetSmear(bool t){ _smear = t; }
		void SetTimeSmear(bool t){ _timesmear = t; }
		void SetOutfile(string fname){ _oname = fname; }
		void SetTransferFactor(double gev){
			_gev = gev;
			_prod->SetTransferFactor(_gev);
		}
		void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj; }

	private:
		string _oname;
		vector<TH1D*> _hists1D;
		vector<TH2D*> _hists2D;
		vector<TGraph*> graphs;
		vector<Jet> _phos; //photons for event
		vector<procCat> _procCats;
		bool _smear, _timesmear;
		enum Strategy{
			//Delauney strategy - NlnN time - for 2pi cylinder
			NlnN = 0,
			//traditional strategy - N^2 time
			N2 = 1,
		};
		//clustering strategy - N^2 or NlnN
		Strategy _strategy;
		ReducedBaseSim* _base = nullptr;
		int _nEvts;
		JetSimProducer* _prod = nullptr;
		bool _data;
		bool _debug;
		int _evti, _evtj;
		double _gev;
		double _c = 29.9792458; // speed of light in cm/ns
		
};
#endif
