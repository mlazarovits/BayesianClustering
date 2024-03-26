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
				
	
			graphs.push_back(nrhs_comptime);

			_hists1D.push_back(nClusters);
			_hists1D.push_back(nSubClusters);	
			_hists1D.push_back(predJet_Energy);
			_hists1D.push_back(predJet_EtaCenter);
			_hists1D.push_back(predJet_PhiCenter);
			_hists1D.push_back(predJet_TimeCenter);
			_hists1D.push_back(predJet_dR);

		}
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else return; 
		}
		
		void FillPredJetHists(const vector<node*>& trees){
			int njets;
			vector<node*> cleaned_trees;
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = 0;
				for(int i = 0; i < trees.size(); i++){
					if(trees[i] == nullptr) continue;
					//check for mirrored point - would be double counted
					if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0) continue;	
					cout << "jet #" << njets << " has " << trees[i]->model->GetNClusters() << " subclusters" << endl;
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
			
			_procCats[p].hists1D[0][1]->Fill(nclusters);
			model->GetNorms(norms);
	
			//center is mm weighted avg of subclusters
			double ceta = 0;
			double cphi = 0;
			double ctime = 0;
			double norm = 0;
			double Etot = 0;
			PointCollection* pts = model->GetData();

	
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				E_k = norms[k]/_gev;
				Etot += E_k;

				params = model->GetPriorParameters(k);
				ceta += params["pi"].at(0,0)*params["mean"].at(0,0);
				cphi += params["pi"].at(0,0)*params["mean"].at(1,0);
				//cout << "pi " << params["pi"].at(0,0) << " ctime " << params["mean"].at(0,2) << " cphi " << params["mean"].at(0,1) << " ceta " << params["mean"].at(0,0) << endl;
				//cout << "switched indices - pi " << params["pi"].at(0,0) << " ctime " << params["mean"].at(2,0) << " cphi " << params["mean"].at(1,0) << " ceta " << params["mean"].at(0,0) << endl;
				ctime += params["pi"].at(0,0)*params["mean"].at(2,0);
				norm += params["pi"].at(0,0);
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//total cluster energy
			}
			_procCats[p].hists1D[0][2]->Fill(Etot);
			_procCats[p].hists1D[0][3]->Fill(ceta/norm);
			_procCats[p].hists1D[0][4]->Fill(cphi/norm);
			_procCats[p].hists1D[0][5]->Fill(ctime/norm);
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
		void WriteOutput(TFile* ofile){
			WriteStackHists(ofile);
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
				if(_hists1D[i]->GetEntries() == 0) continue;
				_hists1D[i]->Write();
			}
			for(int i = 0; i < (int)_hists2D.size(); i++){
				//name = hists2D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDR2DHist(hists2D[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				if(_hists2D[i]->GetEntries() == 0) continue;
				_hists2D[i]->Write();
			}


		}
		
		void WriteStackHists(TFile* ofile){
			ofile->cd();
			string name, dirname, histname;
			//write 1D hists
			//variables
			int nhists = _procCats[0].hists1D[0].size();
			for(int i = 0; i < nhists; i++){
				name = _procCats[0].hists1D[0][i]->GetName();
				_procCats[0].hists1D[0][i]->Write();
				//dirname = name.substr(0,name.rfind("_"+trs[0].methodName));
			////cout << "i: " << i << " name " << name << " making dir " << dirname+"_stack" << endl;
				//TDirectory* dir = ofile->mkdir((name+"_stack").c_str());
				//dir->cd();
				//make process breakdown directory
				TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
				//cout << "  making dir " << dir2->GetName() << endl;
				dir2->cd();
				for(int p = 1; p < _procCats.size(); p++){
				//loop over processes
			//		cout << "    proc " << trs[j].procCats[p].plotName << " hist " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " entries " << trs[j].procCats[p].hists1D[0][i]->GetEntries() << endl;			
					//histname = procCats[p].hists1D[0][i]->GetName();
					if(_procCats[p].hists1D[0][i] == nullptr) continue;
					if(_procCats[p].hists1D[0][i]->GetEntries() == 0){ continue; }
					//cout << "  n hists " << trs[j].procCats[0].hists1D[0].size() << endl;
					//cout << "writing " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;;
					_procCats[p].hists1D[0][i]->Write();

				}
				ofile->cd(); 
			}
			//write 2D hists
			
		}


		
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();
		
		//predicted jet plots
		//0
		TH1D* nClusters = new TH1D("nPredJets","nPredJets",20,0,20);
		//1
		TH1D* nSubClusters = new TH1D("nPredSubClusters","nPredSubClusters",50,0,50);
		//2
		TH1D* predJet_subClusterEnergy = new TH1D("predJet_subClusterEnergy","predJet_subClusterEnergy",20,0,1000);
		//3
		TH1D* predJet_subClusterTimeCenter = new TH1D("predJet_subClustertimeCenter","predJet_subClustertimeCenter",50,-20,20);
		//4
		TH1D* predJet_subClusterEtaCenter = new TH1D("predJet_subClusteretaCenter","predJet_subClusteretaCenter",50,-1.6,1.6);
		//5
		TH1D* predJet_subClusterPhiCenter = new TH1D("predJet_subClusterphiCenter","predJet_subClusterphiCenter",50,0.,6.3);
		//6
		TH1D* predJet_dR = new TH1D("predJet_dR","predJet_dR",50,0,5);

		void SetSmear(bool t){ _smear = t; }
		void SetTimeSmear(bool t){ _timesmear = t; }
		void SetOutfile(string fname){ _oname = fname; }
		void SetTransferFactor(double gev){
			_gev = gev;
			_prod->SetTransferFactor(_gev);
		}
		void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj; }
		void SetVerbosity(int verb){_verb = verb;}
		int _verb;

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
