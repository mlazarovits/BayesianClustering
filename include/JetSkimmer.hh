#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include "BayesCluster.hh"
#include "JetProducer.hh"
#include "PhotonProducer.hh"
#include <TFile.h>
#include <TGraph.h>
#include "TSystem.h"
#include "BaseTree.hh"

using node = BaseTree::node;
using plotCat = BaseSkimmer::plotCat;
class JetSkimmer : public BaseSkimmer{
	public:
		JetSkimmer(){
			_evti = 0;
			_evtj = 0;
		};
		virtual ~JetSkimmer(){ };

		//get rechits from file to cluster
		JetSkimmer(TFile* file) : BaseSkimmer(file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time		

			_prod = new JetProducer(file);
			_prod->SetTransferFactor(_gev);
			_prod->SetIsoCut();
		

			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/jet_skims_"+_cms_label+".root";
				
			objE_clusterE->SetTitle("jetE_clusterE");
			objE_clusterE->SetName("jetE_clusterE");
			//true jet hists
			_hists1D.push_back(nTrueJets);
			_hists1D.push_back(TrueJet_pT); 
			_hists1D.push_back(TrueJet_nRhs); 
			_hists1D.push_back(TrueJet_EmE); 
			_hists1D.push_back(TrueJet_nConstituents);
			_hists1D.push_back(TrueJet_twoHardestpT);	
			_hists1D.push_back(nSubClusters_mm);
			
			_timeHists.push_back(PVtime);
			_timeHists.push_back(deltaT_jet);
			_timeHists.push_back(deltaT_pvGam);	
			//key figures of merit	
			_timeHists.push_back(diffDeltaT_recoGen);
			_timeHists.push_back(ptAvg_resDiffDeltaT);
			
			//_hists2D.push_back(erhs_trhs);		
			
			//predicted jets - from BHC
			//_hists1D.push_back(nClusters);
			//_hists1D.push_back(rhTime);
			//_hists1D.push_back(comptime);
			//_hists1D.push_back(PVtime_median_pred);
			//_hists1D.push_back(PVtime_eAvg_pred);
			//_hists1D.push_back(PVtime_mmAvg_pred);
			//_hists1D.push_back(PVdeltaT_jet_median_pred);
			//_hists1D.push_back(PVdeltaT_jet_eAvg_pred);
			//_hists1D.push_back(PVdeltaT_jet_mmAvg_pred);
			//_hists2D.push_back(e_nRhs);

			nrhs_comptime->SetName("nrhs_comptime");
			nrhs_comptime->SetTitle("nrhs_comptime");
			graphs.push_back(nrhs_comptime);
			MakeTimeRecoCatHists();


		};
		//ctor from rec hit collection - integrating into ntuplizer

		enum TimeStrategy{
			med = 0, 
			eavg = 1, 
			mmavg = 2
		};		
		//struct for different types of time reco (ie median, eAvg, mmAvg)
		struct timeRecoCat{
			vector<TH1D*> hists1D;
			vector<TH2D*> hists2D;
			
			string methodName;

			timeRecoCat(const vector<TH1D*>& in1dhists, const vector<TH2D*>& in2dhists, const TimeStrategy& ts){
				if(ts == med)
					methodName = "median";
				else if(ts == eavg)
					methodName = "eAvg";
				else if(ts == mmavg)
					methodName = "mmAvg";
				else
					cout << "Error: time strategy " << ts << " not supported." << endl;		
				string name;
				//for each histogram (variable or correlation)
				for(int i = 0; i < (int)in1dhists.size(); i++){
						//make sure they have the right add-on name (ie leading, !lead, etc)
						TH1D* hist = (TH1D*)in1dhists[i]->Clone();
						name = hist->GetName();
						if(!methodName.empty()) name += "_"+methodName;
						hist->SetName(name.c_str());
						if(!methodName.empty()) hist->SetTitle("");
						hists1D.push_back(hist);
					}

				//for each histogram
				for(int i = 0; i < (int)in2dhists.size(); i++){
					TH2D* hist = (TH2D*)in2dhists[i]->Clone();
					name = hist->GetName();
					if(!methodName.empty()) name += "_"+methodName;
					hist->SetName(name.c_str());
					if(!methodName.empty()) hist->SetTitle("");
					hists2D.push_back(hist);

				}
			}
			


		};

		void Skim();
		vector<TH1D*> _timeHists;	

		//true jet hists
		TH1D* nTrueJets = new TH1D("nTrueJets","nTrueJets",20,0,20);
		TH1D* rhTime = new TH1D("rhTime","rhTime",100,-30,30); 
		TH1D* TrueJet_pT = new TH1D("TrueJet_pT","TrueJet_pT",100,0,1000);
		TH1D* TrueJet_nRhs = new TH1D("TrueJet_nRhs","TrueJet_nRhs",25,0,100);
		TH1D* TrueJet_EmE = new TH1D("TrueJet_EmE","TrueJet_EmE",50,0,600);
		TH1D* TrueJet_nConstituents = new TH1D("TrueJet_nConstituents","TrueJet_nConstituents",20,0,50);
		TH1D* TrueJet_twoHardestpT =  new TH1D("TrueJet_twoHardestpT","TrueJet_twoHardestpT",100,0,1000);	

		TH2D* e_nRhs = new TH2D("e_nRhs","e_nRhs",100,0,500,100,0,100);
		TH2D* erhs_trhs = new TH2D("erhs_trhs","erhs_trhs",100,0,4,100,-100,100);
		
		TH1D* nSubClusters_mm = new TH1D("nSubClusters_mm","nSubClusters_mm",20,0,20);
		

		//0 - pv time
		TH1D* PVtime = new TH1D("PVtime", "PVtime",100,-10,10);	
		//1 - delta t between jets (pv time frame)
		TH1D* deltaT_jet = new TH1D("PVdeltaT_jet", "PVdeltaT_jet",50,-3,3);	
		//2 - delta t between pv and photon 
		TH1D* deltaT_pvGam = new TH1D("deltaT_pvGam","deltaT_pvGam",25,-3,3);	
		//3 - difference in deltaT_pvGam between gen and reco
		TH1D* diffDeltaT_recoGen = new TH1D("diffDeltaT_recoGen","diffDeltaT_recoGen",50,-3,3);
		//4 - resolution of above difference as a function of pT avg of jets that go into PV time calculation
		//may need to adjust binning
		TH1D* ptAvg_resDiffDeltaT = new TH1D("ptAvg_resDiffDeltaT","ptAvg_resDiff",50,0,100);

		//comparing predicted jets + true jets
		//TH2D* nSubClusters_nConstituents = new TH2D("nSubClusters_nConstituents", "nSubClusters_nConstituents",50,0,20,50,0,20);

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


		vector<timeRecoCat> trCats;
		void MakeTimeRecoCatHists(){
			timeRecoCat trmed(_timeHists, _hists2D, med);
			timeRecoCat treavg(_timeHists, _hists2D, eavg);
			timeRecoCat trmmavg(_timeHists, _hists2D, mmavg);
		
			trCats.push_back(trmed);
			trCats.push_back(treavg);
			trCats.push_back(trmmavg);	
		}

	

		void FillTrueJetHists(const vector<Jet>& jets){
			int njets = jets.size();	
			nTrueJets->Fill((double)njets);
		
			double eECAL = 0;
			int ijet = 0;
			vector<JetPoint> rhs;
			for(int j = 0; j < njets; j++){
				e_nRhs->Fill(jets[j].E(), jets[j].GetNRecHits());
				objE->Fill(jets[j].E());
				TrueJet_pT->Fill(jets[j].pt());
				TrueJet_nRhs->Fill(jets[j].GetNRecHits());
	
				ijet = jets[j].GetUserIdx();
				eECAL = _base->Jet_neEmEF->at(ijet) + _base->Jet_chEmEF->at(ijet);
				eECAL *= jets[j].E();
				TrueJet_EmE->Fill(eECAL);
				TrueJet_nConstituents->Fill(_base->Jet_nConstituents->at(ijet));			
			
				
				rhs.clear();
				rhs = jets[j].GetJetPoints();
				for(int r = 0; r < rhs.size(); r++){
					erhs_trhs->Fill(rhs[r].E(), rhs[r].t());
				}		
				if(njets < 2) continue;
				if(j == 0 || j == 1){	
					TrueJet_twoHardestpT->Fill(jets[j].pt());
					TrueJet_twoHardestpT->Fill(jets[j].pt());
				}
			}
		}

		void FillPVTimeHists(vector<Jet>& jets, int tr_idx, const Matrix& smear = Matrix(), double emAlpha = 0.5, double alpha = 0.1, double tres_c = 0.2, double tres_n = 0.3){
			double time = -999;
			int njets = jets.size();
			double gamtime = -999;
			for(int j = 0; j < njets; j++){
				if(tr_idx == med) time = CalcMedianTime(jets[j]);
				else if(tr_idx == eavg) time = CalcEAvgTime(jets[j]);
				else if(tr_idx == mmavg){
					GaussianMixture* gmm = _subcluster(jets[j], smear, emAlpha, alpha, tres_c, tres_n);
					nSubClusters_mm->Fill(gmm->GetNClusters());
					time = CalcMMAvgTime(gmm);
					jets[j].SetJetTime(time);
				}
				//fill pv time - 0
				trCats[tr_idx].hists1D[0]->Fill( time );
				
				//fill deltaT_pvGam - 2
				if(_phos.size() < 1) continue;
				//only fill for two leading photons + this jet
				gamtime = CalcPhotonTime(TimeStrategy(tr_idx), _phos[0]);

			}
			if(njets < 2) return;
			pair<Jet,Jet> hardjets;
			if(jets[0].pt() < jets[1].pt()){
				FindHardestJetPair(jets, hardjets);
			}
			//pt sorted
			else{
				hardjets = std::make_pair(jets[0],jets[1]);
			}
			double dtime = GetDTime(hardjets);
			//fill deltaT_jets - 1
			trCats[tr_idx].hists1D[1]->Fill(dtime);	
		
			//fill difference in deltaT_pvGam of reco and gen - 3

			//fill res (sigma from gaussian fit) for deltaT_recoGen as a function of ptAvg of jets that go into pv time calc - 4

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
				
				//total cluster energy
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
		void FindHardestJetPair(const vector<Jet>& injets, pair<Jet,Jet>& outjets){
			int njets = (int)injets.size(); 
			map<double,Jet> pt_jet;
			map<double, int> pt_idx;
			for(int i = 0; i < njets; i++){
				pt_jet[injets[i].pt()] = injets[i];
				pt_idx[injets[i].pt()] = i;
			}
			map<double,Jet>::reverse_iterator it = pt_jet.rbegin();	
			map<double,int>::reverse_iterator it_idx = pt_idx.rbegin();	
			outjets.first = it->second;
			outjets.first.SetUserIdx(it_idx->second);
			it++;
			it_idx++;
			outjets.second = it->second;
			outjets.second.SetUserIdx(it_idx->second);
			//if needing to find "true" pair
			//calculate dphi here
			//if dphi doesn't satisfy back-to-back pi requirement (see below) check next two jets (ie it->second and it++; it->second;)
			//continue until pair is found

		}
		double GetDTime(pair<Jet,Jet>& jets){
			double pi = acos(-1);
			//not needed for GMSB
			if(_data){
				//dphi within [pi-0.1,pi+0.1]
				double phi1 = jets.first.phi_02pi();
				double phi2 = jets.second.phi_02pi();
				double dphi = fabs(phi1 - phi2);
				if(dphi < pi-0.1 || dphi > pi+0.1) return -999;
			}
			double t1 = jets.first.time();
			double t2 = jets.second.time();
			return t1 - t2;
		}

	


		double CalcMedianTime(Jet& j){
			vector<double> times;
			double time = -999;
			vector<JetPoint> rhs = j.GetJetPoints();
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++)
				times.push_back(rhs[i].t());
			//even - return average of two median times
			if(nrhs % 2 == 0)
				time = (times[int(double(nrhs)/2.)] + times[int(double(nrhs)/2.)-1]/2.);
			//odd - return median
			else
				time = times[int(double(nrhs)/2.)];
			j.SetJetTime(time);
			return time;
		}		


		double CalcEAvgTime(Jet& j){
			double t = 0;
			double norm = 0;
			double time = -999;
			vector<JetPoint> rhs = j.GetJetPoints();
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++){
				norm += rhs[i].E();
				t += rhs[i].E()*rhs[i].t();
			}
			time = t/norm;
			j.SetJetTime(time);
			return time;
		}


		GaussianMixture* _subcluster(const Jet& jet, const Matrix& smear = Matrix(), double emAlpha = 0.5, double alpha = 0.1, double tres_c = 0.2, double tres_n = 0.3){
			vector<JetPoint> rhs = jet.GetJetPoints();
			vector<Jet> rhs_jet;
			for(int r = 0; r < rhs.size(); r++){
				rhs_jet.push_back( Jet(rhs[r]) );
				rhs_jet[r].SetVertex(jet.GetVertex());
			}
			BayesCluster* algo = new BayesCluster(rhs_jet);
			algo->SetDataSmear(smear);
			//set time resolution smearing
			if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n);
			algo->SetThresh(1.);
			algo->SetAlpha(alpha);
			algo->SetSubclusterAlpha(emAlpha);
			algo->SetVerbosity(0);
			GaussianMixture* gmm = algo->SubCluster();
			return gmm;
		}



		double CalcMMAvgTime(BasePDFMixture* model){
			int kmax = model->GetNClusters();
			double t = 0;
			double norm = 0;
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

		double CalcPhotonTime(const TimeStrategy& ts, const Jet& pho){
			

			return -999;
		}

		void WriteTimeRecoCat1D(TFile* ofile, const timeRecoCat& tr){
			ofile->cd();
			string name;
			vector<TH1D*> hists1D = tr.hists1D;
			//write 1D hists
			for(int i = 0; i < (int)hists1D.size(); i++){
				if(tr.hists1D[i] == nullptr) continue;
				name = tr.hists1D[i]->GetName();
				if(tr.hists1D[i]->GetEntries() == 0){ continue; }//cout << "Histogram: " << name << " not filled." << endl; continue; }
				tr.hists1D[i]->Write();
			}
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
			erhs_trhs->Write();		
			for(int i = 0; i < (int)graphs.size(); i++){
				//name = graphs[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRGraph(graphs[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				graphs[i]->Write();
			}
			for(int i = 0; i < trCats.size(); i++)
				WriteTimeRecoCat1D(ofile, trCats[i]);

			ofile->Close();

		}


		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else if(i == 2) _strategy = MM;
			else return; 
		}


		private:
			enum Strategy{
				//Delauney strategy - NlnN time - for 2pi cylinder
				NlnN = 0,
				//traditional strategy - N^2 time
				N2 = 1,
				//mm only
				MM = 2
			};
		
		
			//clustering strategy - N^2 or NlnN
			Strategy _strategy;
			vector<TGraph*> graphs;
			vector<Jet> _phos; //photons for event

};
#endif
