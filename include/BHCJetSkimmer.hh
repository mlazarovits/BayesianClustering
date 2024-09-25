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
			_radius = 0;
			_smear = true;
			_timesmear = false;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_thresh = 1.;
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
			_radius = 0;
			_smear = true;
			_timesmear = false;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_thresh = 1.;
				
	
			graphs.push_back(nrhs_comptime);

			_hists1D.push_back(nClusters);
			_hists1D.push_back(nSubClusters);	
			_hists1D.push_back(predJet_subClusterEnergy);
			_hists1D.push_back(predJet_subClusterEtaCenter);
			_hists1D.push_back(predJet_subClusterPhiCenter);
			_hists1D.push_back(predJet_subClusterTimeCenter);
			_hists1D.push_back(predJet_jetSize);
			_hists1D.push_back(predJet_energy);
			_hists1D.push_back(predJet_pt);
			_hists1D.push_back(predJet_mass);
			_hists1D.push_back(jetGenE_sigmaDeltaPt_predGen);
			_hists1D.push_back(predGen_nJets);
			_hists1D.push_back(predJet_subClusterEtaVar);
			_hists1D.push_back(predJet_subClusterPhiVar);
			_hists1D.push_back(predJet_subClusterTimeVar);
			_hists1D.push_back(predJet_subClusteretaPhiCov);
			_hists1D.push_back(predJet_subClustertimeEtaCov);
			_hists1D.push_back(predJet_subClustertimePhiCov);
			_hists1D.push_back(nRecoJets);
			_hists1D.push_back(recoJet_jetSize);
			_hists1D.push_back(recoJet_energy);
			_hists1D.push_back(recoJet_pt);
			_hists1D.push_back(recoJet_mass);
			_hists1D.push_back(jetGenE_sigmaDeltaPt_recoGen);
			_hists1D.push_back(recoGen_nJets);
			_hists1D.push_back(recoGen_jetPtRatio);
			_hists1D.push_back(recoGen_jetERatio);
			_hists1D.push_back(recoJet_Wmass);
			_hists1D.push_back(recoJet_topmass);
			_hists1D.push_back(predJet_Wmass);
			_hists1D.push_back(predJet_topmass);
			_hists1D.push_back(predJet_Wmass_pTjjl100);
			_hists1D.push_back(predJet_Wmass_pTjjge100);
			_hists1D.push_back(recoJet_Wmass_pTjjl100);
			_hists1D.push_back(recoJet_Wmass_pTjjge100);
			_hists1D.push_back(predJet_etaVar_Wjj);
			_hists1D.push_back(predJet_etaVar_W_pTl100);
			_hists1D.push_back(predJet_etaVar_W_pTge100);
			_hists1D.push_back(predJet_etaVar_notWjj);
			_hists1D.push_back(predJet_phiVar_Wjj);
			_hists1D.push_back(predJet_phiVar_W_pTl100);
			_hists1D.push_back(predJet_phiVar_W_pTge100);
			_hists1D.push_back(predJet_phiVar_notWjj);
			_hists1D.push_back(predJet_timeVar_Wjj);
			_hists1D.push_back(predJet_timeVar_W_pTl100);
			_hists1D.push_back(predJet_timeVar_W_pTge100);
			_hists1D.push_back(predJet_timeVar_notWjj);
			_hists1D.push_back(predJet_EtaVar);
			_hists1D.push_back(predJet_PhiVar);
			_hists1D.push_back(predJet_TimeVar);
			_hists1D.push_back(predJet_etaPhiCov);
			_hists1D.push_back(predJet_timeEtaCov);
			_hists1D.push_back(predJet_timePhiCov);
			_hists1D.push_back(predJet_WjjDr);
			_hists1D.push_back(recoJet_WjjDr);
			_hists1D.push_back(predJet_Wpt);
			_hists1D.push_back(recoJet_Wpt);

			_hists2D.push_back(jetGenE_diffDeltaPt_predGen);
			_hists2D.push_back(jetGenE_diffDeltaPt_recoGen);
			_hists2D.push_back(genPt_recoPt);
			_hists2D.push_back(genJetMass_recoGenPtRatio);
			_hists2D.push_back(recoJetMass_recoJetPt);
			_hists2D.push_back(recoJetMass_recoJetSize);
			_hists2D.push_back(recoJetInvMassW_recoJetPairjetSize); 
			_hists2D.push_back(predJetMass_predJetPt);
			_hists2D.push_back(predJetMass_predJetSize);
			_hists2D.push_back(predJetInvMassW_predJetPairjetSize); 
			_hists2D.push_back(predJetPt_predJetSize);

		}
		void SetMinRhE(double r){ _prod->SetMinRhE(r); }
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else return; 
		}
	

		void CleanTrees(const vector<node*>& trees){
			_trees.clear();
			for(int i = 0; i < trees.size(); i++){
				if(trees[i] == nullptr) continue;
				//check for mirrored point - would be double counted
				if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0) continue;	
				_trees.push_back(trees[i]);
			}
		}


		void TreesToJets(){
			_predJets.clear();
			vector<JetPoint> rhs;
			double x, y, z, eta, phi, t, theta, px, py, pz;
			BayesPoint vertex({_pvx, _pvy, _pvz});
			for(int i = 0; i < _trees.size(); i++){
				//get points from tree
				PointCollection* pc = _trees[i]->points;
				//at least 2 points (rhs)
				if(pc->GetNPoints() < 2) continue;
				rhs.clear();
				//loop over points
				cout << "TREE " << i << endl;
				for(int p = 0; p < pc->GetNPoints(); p++){
					//declare JetPoint with x, y, z, t
					eta = pc->at(p).at(0);
					phi = pc->at(p).at(1);
					t = pc->at(p).at(2);
					x = _radius*cos(phi);
					y = _radius*sin(phi);
					theta = 2*atan2(1,exp(eta));
					z = _radius/tan(theta); 

					JetPoint jp(x, y, z, t);
					//add JetPoint to list of rhs
					jp.SetEnergy(pc->at(p).w()/_gev);
					jp.SetWeight(pc->at(p).w());
					jp.SetEta(eta);
					jp.SetPhi(phi);
					//cout << "rh - eta " << jp.eta() << " jp phi " << jp.phi() << " energy " << jp.e() << endl;
					rhs.push_back(jp);
				}
				//create new Jet
				Jet predJet(rhs, BayesPoint({_pvx, _pvy, _pvz}));
				//cout << "additive p: (" << pxtot << ", " << pytot << ", " << pztot << ", " << e << ")" << " mass " << sqrt((e*e) - (pxtot*pxtot + pytot*pytot + pztot*pztot)) << endl;
				//cout << "jet p: (" << predJet.px() << ", " << predJet.py() << ", " << predJet.pz() << ", " << predJet.e() << ") mass " << predJet.mass() << endl;
				//set PV info
				//predJet.SetVertex(Point({_pvx,_pvy,_pvz}));
				//set constituents (subclusters) here with model from tree
				int nsubclusters = _trees[i]->model->GetNClusters();
				double Ek; 
				vector<double> norms;
				_trees[i]->model->GetNorms(norms);
				cout << "jet " << _predJets.size() << " norms size " << norms.size() << " nsubcluster " << nsubclusters << endl;
				map<string,Matrix> params;
				for(int k = 0; k < nsubclusters; k++){
					params = _trees[i]->model->GetPriorParameters(k);
					Ek = norms[k]/_gev;
					cout << "mu" << endl; params["mean"].Print();
					Jet jet(params["mean"], params["cov"], Ek, params["pi"].at(0,0));
					predJet.AddConstituent(jet);
				}
				//put pt cut in for predjets of 20 GeV
				if(predJet.pt() < 20) continue; 
				//add Jet to jets	
				_predJets.push_back(predJet);	
			}
			cout << _predJets.size() << " pred jets" << endl;
			for(auto j : _predJets) cout << "jet px " << j.px() << " py " << j.py() << " pz " << j.pz() << " E " << j.E() << " mass " << j.mass() << endl;
		}

		void FillPredJetHists(){
			int njets = _predJets.size();
			double wmass, topmass, dr, dr_pair;
			int tmass_idx;
			pair<int,int> wmass_idxs;
			Jet w, top;
			for(int p = 0; p < _procCats.size(); p++){
				//if(p != 0) cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				_procCats[p].hists1D[0][0]->Fill(njets);
				cout << "# pred jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][11]->Fill(njets - (int)_genjets.size());
				FindResonances(_predJets,wmass_idxs,tmass_idx);
				//W mass
				if(wmass_idxs.first != -999){
					w = _predJets[wmass_idxs.first];
					w.add(_predJets[wmass_idxs.second]);
					wmass = w.mass();
					dr_pair = dR(_predJets[wmass_idxs.first].eta(), _predJets[wmass_idxs.first].phi(), _predJets[wmass_idxs.second].eta(), _predJets[wmass_idxs.second].phi());
					_procCats[p].hists1D[0][29]->Fill(wmass);
					Matrix cov1 = _predJets[wmass_idxs.first].GetCovariance();
					Matrix cov2 = _predJets[wmass_idxs.second].GetCovariance();
					Matrix covw = w.GetCovariance();
					_procCats[p].hists1D[0][35]->Fill(cov1.at(0,0));
					_procCats[p].hists1D[0][35]->Fill(cov2.at(0,0));
					_procCats[p].hists1D[0][39]->Fill(cov1.at(1,1));
					_procCats[p].hists1D[0][39]->Fill(cov2.at(1,1));
					_procCats[p].hists1D[0][43]->Fill(cov1.at(2,2));
					_procCats[p].hists1D[0][43]->Fill(cov2.at(2,2));


					_procCats[p].hists1D[0][53]->Fill(dr_pair);
					_procCats[p].hists2D[0][9]->Fill(wmass, dr_pair);
					_procCats[p].hists1D[0][55]->Fill(w.pt());
					if(w.pt() < 100){
						_procCats[p].hists1D[0][31]->Fill(w.mass());
						_procCats[p].hists1D[0][36]->Fill(covw.at(0,0));
						_procCats[p].hists1D[0][40]->Fill(covw.at(1,1));
						_procCats[p].hists1D[0][44]->Fill(covw.at(2,2));
					}
					else{
						_procCats[p].hists1D[0][32]->Fill(w.mass());
						_procCats[p].hists1D[0][37]->Fill(covw.at(0,0));
						_procCats[p].hists1D[0][41]->Fill(covw.at(1,1));
						_procCats[p].hists1D[0][45]->Fill(covw.at(2,2));
					}
				}
				//top mass
				if(tmass_idx != -999){
					top = _predJets[tmass_idx];
					top.add(w);
					topmass = top.mass();
					_procCats[p].hists1D[0][30]->Fill(topmass);
				}
				njets = _predJets.size();
				for(int j = 0; j < _predJets.size(); j++){
					//if(p != 0) cout << "pred jet #" << j << " phi " << _predJets[j].phi() << " eta " << _predJets[j].eta() << " energy " << _predJets[j].E() <<  " mass " << _predJets[j].mass() << " nConstituents " << _predJets[j].GetNConstituents() << " nRhs " << _predJets[j].GetNRecHits() << " pt " << _predJets[j].pt() << endl;
					//cout << "calc jet size for jet " << j << endl;
					dr = CalcJetSize(_predJets[j]);
					//cout << "pred jet j " << j << " dr " << dr << " n constituents " << _predJets[j].GetNConstituents() << endl;
					Matrix cov = _predJets[j].GetCovariance();
				
					_procCats[p].hists1D[0][7]->Fill(_predJets[j].e());
					_procCats[p].hists1D[0][8]->Fill(_predJets[j].pt());
					_procCats[p].hists1D[0][9]->Fill(_predJets[j].mass());
					if(dr != -999) _procCats[p].hists1D[0][6]->Fill(dr);	
					_procCats[p].hists1D[0][47]->Fill(cov.at(0,0));	
					_procCats[p].hists1D[0][48]->Fill(cov.at(1,1));	
					_procCats[p].hists1D[0][49]->Fill(cov.at(2,2));	
					_procCats[p].hists1D[0][50]->Fill(cov.at(0,1));	
					_procCats[p].hists1D[0][51]->Fill(cov.at(0,2));	
					_procCats[p].hists1D[0][52]->Fill(cov.at(2,1));	
					
					_procCats[p].hists2D[0][7]->Fill(_predJets[j].mass(), _predJets[j].pt());
					_procCats[p].hists2D[0][8]->Fill(_predJets[j].mass(), dr);
					_procCats[p].hists2D[0][10]->Fill(_predJets[j].pt(), dr);


					//get subcluster information
					for(auto subcl : _predJets[j].GetConstituents()){
						Matrix subcl_cov = subcl.GetCovariance();
						_procCats[p].hists1D[0][12]->Fill(subcl_cov.at(0,0));
						_procCats[p].hists1D[0][13]->Fill(subcl_cov.at(1,1));
						_procCats[p].hists1D[0][14]->Fill(subcl_cov.at(2,2));
						_procCats[p].hists1D[0][15]->Fill(subcl_cov.at(0,1));
						_procCats[p].hists1D[0][16]->Fill(subcl_cov.at(0,2));
						_procCats[p].hists1D[0][17]->Fill(subcl_cov.at(1,2));
					}



					if(j != wmass_idxs.first && j != wmass_idxs.second){
						_procCats[p].hists1D[0][38]->Fill(cov.at(0,0));
						_procCats[p].hists1D[0][42]->Fill(cov.at(1,1));
						_procCats[p].hists1D[0][46]->Fill(cov.at(2,2));
					}
				}
			}
		}
	
		void FillRecoJetHists(){
			int njets;
			double wmass, topmass, dr, dr_pair;
			pair<int, int> wmass_idxs;
			int tmass_idx;
			Jet w, top;
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = _recojets.size();
				for(int j = 0; j < _recojets.size(); j++){
					dr = CalcJetSize(_recojets[j]);
					cout << "reco jet #" << j << " phi " << _recojets[j].phi() << " eta " << _recojets[j].eta() << " energy " << _recojets[j].E() <<  " mass " << _recojets[j].mass() << " nConstituents " << _recojets[j].GetNConstituents() << " nRhs " << _recojets[j].GetNRecHits() << " pt " << _recojets[j].pt() << " dr " << dr << endl;
					_procCats[p].hists1D[0][19]->Fill(_recojets[j].e());
					_procCats[p].hists1D[0][20]->Fill(_recojets[j].e());
					//dr hist is #19
					_procCats[p].hists1D[0][21]->Fill(_recojets[j].pt());
					_procCats[p].hists1D[0][22]->Fill(_recojets[j].mass());
					_procCats[p].hists2D[0][4]->Fill(_recojets[j].mass(), _recojets[j].pt());
					if(dr != -999) _procCats[p].hists2D[0][5]->Fill(_recojets[j].mass(), dr);
				}
				_procCats[p].hists1D[0][18]->Fill(njets);
				//cout << "hist name " << _procCats[p].hists1D[0][18]->GetName() << " nentries " << _procCats[p].hists1D[0][18]->GetEntries() << endl;
				cout << "# reco jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][24]->Fill(njets - (int)_genjets.size());
				FindResonances(_recojets,wmass_idxs,tmass_idx);
				if(wmass_idxs.first != -999){
					w = _recojets[wmass_idxs.first];
					w.add(_recojets[wmass_idxs.second]);
					wmass = w.mass();
					dr_pair = dR(_recojets[wmass_idxs.first].eta(), _recojets[wmass_idxs.first].phi(), _recojets[wmass_idxs.second].eta(), _recojets[wmass_idxs.second].phi());
					_procCats[p].hists1D[0][29]->Fill(wmass);
					_procCats[p].hists1D[0][54]->Fill(dr_pair);
					_procCats[p].hists1D[0][56]->Fill(w.pt());

					_procCats[p].hists2D[0][6]->Fill(wmass,dr_pair);
					if(w.pt() < 100) _procCats[p].hists1D[0][33]->Fill(w.mass());
					else _procCats[p].hists1D[0][34]->Fill(w.mass());
				}
				if(tmass_idx != -999){
					top = _recojets[tmass_idx];
					top.add(w);
					topmass = top.mass();
					_procCats[p].hists1D[0][28]->Fill(topmass);
				}
			}
		}

		void FillResolutionHists(){
			//need to gen match jets to find difference in pt bw reco - gen
			//dr match (maybe do dE match later)
			for(int p = 0; p < _procCats.size(); p++){
				double dr;
				int bestIdx = 0;
				//cout << "proc " << p << " # reco jets " << _recojets.size() << " # gen jets " << _genjets.size() << endl;
				//reco jets
				for(int j = 0; j < _recojets.size(); j++){
					dr = 999;
					for(int g = 0; g < _genjets.size(); g++){
						if(dR(_recojets[j].eta(), _recojets[j].phi(), _genjets[g].eta(), _genjets[g].phi()) < dr){
							dr = dR(_recojets[j].eta(), _recojets[j].phi(), _genjets[g].eta(), _genjets[g].phi());
							bestIdx = g;
						}
					}
					//cout << "jet #" << j << " has best match with gen jet #" << bestIdx << " with dr " << dr << " reco E " << _recojets[j].E() << " gen energy " << _genjets[bestIdx].E() << " reco pt " << _recojets[j].pt() << " gen pt " << _genjets[bestIdx].pt() << endl;
					_procCats[p].hists2D[0][1]->Fill(_genjets[bestIdx].E(), _recojets[j].pt() - _genjets[bestIdx].pt());
					//cout << "hist has " << _procCats[p].hists2D[0][1]->GetEntries() << " entries" << endl;
					_procCats[p].hists1D[0][25]->Fill(_recojets[j].pt()/_genjets[bestIdx].pt());
					_procCats[p].hists1D[0][26]->Fill(_recojets[j].e()/_genjets[bestIdx].e());
					_procCats[p].hists2D[0][2]->Fill(_genjets[bestIdx].pt(), _recojets[j].pt());
					_procCats[p].hists2D[0][3]->Fill(_genjets[bestIdx].mass(), _recojets[j].pt()/_genjets[bestIdx].pt());
				}
				//predicted jets
				for(int j = 0; j < _predJets.size(); j++){
					dr = 999;
					for(int g = 0; g < _genjets.size(); g++){
						if(dR(_predJets[j].eta(), _predJets[j].phi(), _genjets[g].eta(), _genjets[g].phi()) < dr){
							dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genjets[g].eta(), _genjets[g].phi());
							bestIdx = g;
						}
					}
					//cout << "jet #" << j << " has best match with gen jet #" << bestIdx << " with dr " << dr << " reco E " << _predJets[j].E() << " gen energy " << _genjets[bestIdx].E() << endl;
					_procCats[p].hists2D[0][0]->Fill(_genjets[bestIdx].E(), _predJets[j].pt() - _genjets[bestIdx].pt());
				}
			}


		}

		void FillModelHists(){
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				for(int i = 0; i < _trees.size(); i++){
					//keep checks in just in case
					if(_trees[i] == nullptr) continue;
					//check for mirrored point - would be double counted
					if(_trees[i]->points->mean().at(1) > 2*acos(-1) || _trees[i]->points->mean().at(1) < 0) continue;
					FillModelHists(_trees[i]->model, p);
				}
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
	
		//use subclusters
		//can change to rhs later
		double CalcJetSize(const Jet& jet){
			int nSCs = jet.GetNConstituents();
			//if no subclusters (ie reco jet), take dR to be max distance between rechits
			if(nSCs < 1){
				double maxsize = -999;
				double size;
				vector<JetPoint> rhs = jet.GetJetPoints();
				if(rhs.size() < 2) return maxsize;
				for(int i = 0; i < rhs.size(); i++){
					for(int j = i+1; j < rhs.size(); j++){
						size = dR(rhs[i].eta(), rhs[i].phi(), rhs[j].eta(), rhs[j].phi());
						//cout << "i " << i << " j " << j << " dr " << dr << " eta " << rhs[i].eta() << " " << rhs[j].eta() << " phi "<< rhs[i].phi() << " " << rhs[j].phi() << endl;
						if(size > maxsize) maxsize = size;
					}
				}	
				return maxsize;
			}
			//if 1 subcluster, take dR to be 1 sigma (for pred jets)
			if(nSCs < 2){
				Matrix cov, mu;
				jet.GetClusterParams(mu, cov);
				return sqrt(cov.at(0,0) + cov.at(1,1));

			}

			//else, weighted average of subcluster sizes
			double size = 0;
			double norm = 0;
			Jet subjet;
			for(int i = 0; i < nSCs; i++){
				//get center in eta, phi for subcluster i and j
				subjet = jet.GetConstituent(i);
				size += CalcJetSize(subjet)*subjet.GetCoefficient();
				norm += subjet.GetCoefficient();
			}	
			return size/norm;
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
			WriteEmptyProfiles(ofile);
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
				//make process breakdown directory - not making these profiles rn
				TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
				//cout << "  making dir " << dir2->GetName() << " name " << name << endl;
				dir2->cd();
				for(int p = 1; p < _procCats.size(); p++){
				//loop over processes
					//if(name.find("profile") != string::npos) continue;
					if(_procCats[p].hists1D[0][i] == nullptr) continue;
					histname = _procCats[p].hists1D[0][i]->GetName();
					//cout << "    proc " << _procCats[p].plotName << " hist " << _procCats[p].hists1D[0][i]->GetName() << " " << _procCats[p].hists1D[0][i]->GetTitle() << " entries " << _procCats[p].hists1D[0][i]->GetEntries() << endl;			
					if(_procCats[p].hists1D[0][i]->GetEntries() == 0 && ((histname.find("sigma") == string::npos && histname.find("mean") == string::npos) && histname.find("profile") == string::npos)){ continue; }
					if(histname.find("profile") != string::npos){
						histname = histname.substr(0,histname.rfind("_"));
						_procCats[p].hists1D[0][i]->SetName(histname.c_str());
						_procCats[p].hists1D[0][i]->SetTitle(_procCats[p].plotName.c_str());
					}
					//cout << "  n hists " << _procCats[0].hists1D[0].size() << endl;
					//cout << "i " << i << " p " << p << " writing " << _procCats[p].hists1D[0][i]->GetName() << " " << _procCats[p].hists1D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;;
					_procCats[p].hists1D[0][i]->Write();

				}
				ofile->cd(); 
			}
			//cout << "2D hists" << endl;
			//write 2D hists
			nhists = _procCats[0].hists2D[0].size();
			for(int i = 0; i < nhists; i++){
				name = _procCats[0].hists2D[0][i]->GetName();
				histname = _procCats[0].hists2D[0][i]->GetName();
				//cout << "writing for " << name << " i " << i << " entries " << _procCats[0].hists2D[0][i]->GetEntries() << endl;
				//write total method histogram outside process directory
				if(_procCats[0].hists2D[0][i] == nullptr) continue;
				//cout << "passed null" << endl;
				if(_procCats[0].hists2D[0][i]->GetEntries() == 0 && histname.find("sigma") == string::npos){ continue; }
				//cout << "writing hist " << _procCats[0].hists2D[0][i]->GetName() << endl;
				_procCats[0].hists2D[0][i]->Write();
				//write method as directory within directory
				TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
				//cout << "  making dir " << dir2->GetName() << endl;
				dir2->cd();
				for(int p = 1; p < _procCats.size(); p++){
					//loop over processes
					if(_procCats[p].hists2D[0][i] == nullptr) continue;
					if(_procCats[p].hists2D[0][i]->GetEntries() == 0 && dirname.find("meanRecoGenDeltaT") == string::npos){ continue;}// cout << "Histogram for proc " << _procCats[p].hists2D[0][i]->GetName() << " not filled." << endl; continue; }
					//check if data can be run
					//cout << "writing " << _procCats[p].hists2D[0][i]->GetName() <<  " " <<  _procCats[p].hists2D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;
					histname = _procCats[p].hists2D[0][i]->GetName();
					_procCats[p].hists2D[0][i]->SetTitle(_procCats[p].plotName.c_str());
					_procCats[p].hists2D[0][i]->Write();
				} 
				ofile->cd(); 
			}
		}


		void WriteEmptyProfiles(TFile* ofile){
			ofile->cd();
			string name, addname;
			int i, j, nbins, nprofs;
			string profname, histname, match;
			for(int p = 0; p < _procCats.size(); p++){
				for(int i = 0; i < _procCats[p].hists2D[0].size(); i++){
					histname = _procCats[p].hists2D[0][i]->GetName();
					//make sure taking info from 2D diff histogram
					if(histname.find("diffDelta") == string::npos) continue;
					//make sure profiles get written
					nprofs = _procCats[p].hists2D[0][i]->GetNbinsX();
					for(int k = 1; k < nprofs+1; k++){
						histname = _procCats[p].hists2D[0][i]->GetName();
						profname = "profile_"+histname;
						//profname.insert(profname.size(),"_bin"+std::to_string(k));
			//cout << "p " << p << " profname " << profname << " plotname " << _procCats[p].plotName << endl;
						if(!_procCats[p].plotName.empty()) profname.insert(profname.find("_"+_procCats[p].plotName),"_bin"+std::to_string(k));
						else profname = profname += "_bin"+std::to_string(k);
						nbins = _procCats[p].hists2D[0][i]->GetNbinsY();
						TH1D* prof = new TH1D(profname.c_str(), profname.c_str(), nbins, _procCats[0].hists2D[0][i]->GetYaxis()->GetBinLowEdge(1), _procCats[0].hists2D[0][i]->GetYaxis()->GetBinUpEdge(nbins));
						prof->SetTitle(_procCats[p].plotName.c_str());
						prof->GetXaxis()->SetTitle(_procCats[p].hists2D[0][i]->GetYaxis()->GetTitle());	
						//cout << "adding hist " << prof->GetName() <<  " " << prof->GetTitle() << endl;
						_procCats[p].AddHist(prof);	
						//cout << "current list of hists " << endl;
						//for(int i = 0; i < _procCats[p].hists1D[0].size(); i++) cout << "i " << i << " p " << p << " " << _procCats[p].hists1D[0][i]->GetName() << endl;
					}	
				}				
		}	
	}
		
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();
		
		//predicted jet plots
		vector<double> xbins_recoGenPt = {0, 20, 30, 50, 100};
		//0
		TH1D* nClusters = new TH1D("nPredJets","nPredJets",10,0,10);
		//1
		TH1D* nSubClusters = new TH1D("predJet_nSubClusters","predJet_nSubClusters",10,0,10);
		//2
		TH1D* predJet_subClusterEnergy = new TH1D("predJet_subClusterEnergy","predJet_subClusterEnergy",20,0,500);
		//3
		TH1D* predJet_subClusterTimeCenter = new TH1D("predJet_subClusterTimeCenter","predJet_subClusterTimeCenter",25,0,15);
		//4
		TH1D* predJet_subClusterEtaCenter = new TH1D("predJet_subClusterEtaCenter","predJet_subClusterEtaCenter",25,-1.8,1.8);
		//5
		TH1D* predJet_subClusterPhiCenter = new TH1D("predJet_subClusterPhiCenter","predJet_subClusterPhiCenter",25,-0.1,6.3);
		//6
		TH1D* predJet_jetSize = new TH1D("predJet_jetSize","predJet_jetSize",50,0,2);
		//7
		TH1D* predJet_energy = new TH1D("predJet_energy","predJet_energy",50,0,400);
		//8
		TH1D* predJet_pt = new TH1D("predJet_pt","predJet_pt",50,0,150);
		//9
		TH1D* predJet_mass = new TH1D("predJet_mass","predJet_mass",50,0,200);
		//10 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_predGen = new TH1D("jetGenE_sigmaDeltaPt_predGen","jetGenE_sigmaDeltaPt_predGen",5,0,100);
		//11 - # pred jets - # gen jets
		TH1D* predGen_nJets = new TH1D("predGen_diffNJets","predGen_diffNJets",20,-10,10);
		//for subclusters
		//12 - eta sigma
		TH1D* predJet_subClusterEtaVar = new TH1D("predJet_subClusterEtaVar","predJet_subClusterEtaVar",25,0.,3.);
		//13 - phi sigma
		TH1D* predJet_subClusterPhiVar = new TH1D("predJet_subClusterPhiVar","predJet_subClusterPhiVar",25,0.,3.);
		//14 - time sigma
		TH1D* predJet_subClusterTimeVar = new TH1D("predJet_subClusterTimeVar","predJet_subClusterTimeVar",25,0.,5.);
		//15 - eta-phi covariance
		TH1D* predJet_subClusteretaPhiCov = new TH1D("predJet_subClusteretaPhiCov","predJet_subClusteretaPhiCov",25,-1.5,1.5);
		//16 - time-eta covariance
		TH1D* predJet_subClustertimeEtaCov = new TH1D("predJet_subClustertimeEtaCov","predJet_subClustertimeEtaCov",25,-1.5,1.5);
		//17 - time-phi covariance
		TH1D* predJet_subClustertimePhiCov = new TH1D("predJet_subClustertimePhiCov","predJet_subClustertimePhiCov",25,-1.5,1.5);
		//18
		TH1D* nRecoJets = new TH1D("nRecoJets","nRecoJets",10,0,10);
		//19
		TH1D* recoJet_jetSize = new TH1D("recoJet_jetSize","recoJet_jetSize",50,0,1);
		//20
		TH1D* recoJet_energy = new TH1D("recoJet_energy","recoJet_energy",50,0,200);
		//21
		TH1D* recoJet_pt = new TH1D("recoJet_pt","recoJet_pt",50,0,200);
		//22
		TH1D* recoJet_mass = new TH1D("recoJet_mass","recoJet_mass",60,0,60);
		//23 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_recoGen = new TH1D("jetGenE_sigmaDeltaPt_recoGen","jetGenE_sigmaDeltaPtOvJetGenE_recoGen",4,&xbins_recoGenPt[0]);
		//24 - # reco jets - # gen jets
		TH1D* recoGen_nJets = new TH1D("recoGen_diffNJets","recoGen_diffNJets",20,-10,10);
		//25 - reco jet pt/gen jet pt
		TH1D* recoGen_jetPtRatio = new TH1D("recoGen_jetPtRatio","recoGen_jetPtRatio",20,0,1.5);
		//26 - reco jet e - gen jet e		
		TH1D* recoGen_jetERatio = new TH1D("recoGen_jetERatio","recoGen_jetERatio",20,-10,10);
		//27 - reco jet W invariant mass
		TH1D* recoJet_Wmass = new TH1D("recoJet_Wmass","recoJet_Wmass;Invariant mass for best W candidate",50,50,250);
		//28 - reco jet top invariant mass
		TH1D* recoJet_topmass = new TH1D("recoJet_topmass","recoJet_topmass;Invariant mass for best top candidate",25,100,250);
		//29 - pred jet W invariant mass
		TH1D* predJet_Wmass = new TH1D("predJet_Wmass","predJet_Wmass;Invariant mass for best W candidate",50,50,250);
		//30 - pred jet top invariant mass
		TH1D* predJet_topmass = new TH1D("predJet_topmass","predJet_topmass;Invariant mass for best top candidate",25,100,250);
		//31 - pred jet W invariant mass, pT_jj < 100
		TH1D* predJet_Wmass_pTjjl100 = new TH1D("predJet_Wmass_pTjjl100","predJet_Wmass_pTjjl100;Invariant mass_pTjjl100 for best W candidate",50,0,250);
		//32 - pred jet W invariant mass, pT_jj >= 100
		TH1D* predJet_Wmass_pTjjge100 = new TH1D("predJet_Wmass_pTjjge100","predJet_Wmass_pTjjge100;Invariant mass_pTjjge100 for best W candidate",50,0,250);
		//33 - reco jet W invariant mass, pT_jj < 100
		TH1D* recoJet_Wmass_pTjjl100 = new TH1D("recoJet_Wmass_pTjjl100","recoJet_Wmass_pTjjl100;Invariant mass_pTjjl100 for best W candidate",50,0,250);
		//34 - reco jet W invariant mass, pT_jj >= 100
		TH1D* recoJet_Wmass_pTjjge100 = new TH1D("recoJet_Wmass_pTjjge100","recoJet_Wmass_pTjjge100;Invariant mass_pTjjge100 for best W candidate",50,0,250);
		//35 - W jet pair eta variance
		TH1D* predJet_etaVar_Wjj = new TH1D("predJet_etaVar_Wjj","predJet_etaVar_Wjj",50,0,1);
		//36 - W jet eta variance, pT < 100 (resolved)
		TH1D* predJet_etaVar_W_pTl100 = new TH1D("predJet_etaVar_W_pTl100","predJet_etaVar_W_pTl100",50,0,1);
		//37 - W jet eta variance, pT >= 100 (boosted)
		TH1D* predJet_etaVar_W_pTge100 = new TH1D("predJet_etaVar_W_pTge100","predJet_etaVar_W_pTge100",50,0,1);
		//38 - !W jet pair eta variance
		TH1D* predJet_etaVar_notWjj = new TH1D("predJet_etaVar_notWjj","predJet_etaVar_notWjj",50,0,1);
		//39 - W jet pair phi variance
		TH1D* predJet_phiVar_Wjj = new TH1D("predJet_phiVar_Wjj","predJet_phiVar_Wjj",50,0,1);
		//40 - W jet phi variance, pT < 100 (resolved)
		TH1D* predJet_phiVar_W_pTl100 = new TH1D("predJet_phiVar_W_pTl100","predJet_phiVar_W_pTl100",50,0,1);
		//41 - W jet phi variance, pT >= 100 (boosted)
		TH1D* predJet_phiVar_W_pTge100 = new TH1D("predJet_phiVar_W_pTge100","predJet_phiVar_W_pTge100",50,0,1);
		//42 - !W jet pair phi variance
		TH1D* predJet_phiVar_notWjj = new TH1D("predJet_phiVar_notWjj","predJet_phiVar_notWjj",50,0,1);
		//43 - W jet pair time variance
		TH1D* predJet_timeVar_Wjj = new TH1D("predJet_timeVar_Wjj","predJet_timeVar_Wjj",50,0,1);
		//44 - W jet time variance, pT < 100 (resolved)
		TH1D* predJet_timeVar_W_pTl100 = new TH1D("predJet_timeVar_W_pTl100","predJet_timeVar_W_pTl100",50,0,1);
		//45 - W jet time variance, pT >= 100 (boosted)
		TH1D* predJet_timeVar_W_pTge100 = new TH1D("predJet_timeVar_W_pTge100","predJet_timeVar_W_pTge100",50,0,1);
		//46 - !W jet pair time variance
		TH1D* predJet_timeVar_notWjj = new TH1D("predJet_timeVar_notWjj","predJet_timeVar_notWjj",50,0,1);
		//47 - eta sigma for jet
		TH1D* predJet_EtaVar = new TH1D("predJet_EtaVar","predJet_EtaVar",25,0.,1.);
		//48 - phi sigma for jet
		TH1D* predJet_PhiVar = new TH1D("predJet_PhiVar","predJet_PhiVar",25,0., 1.);
		//49 - time sigma for jet
		TH1D* predJet_TimeVar = new TH1D("predJet_TimeVar","predJet_TimeVar",25,0.,1.);
		//50 - eta-phi covariance for jet
		TH1D* predJet_etaPhiCov = new TH1D("predJet_etaPhiCov","predJet_etaPhiCov",25,-1.,1.);
		//51 - time-eta covariance for jet
		TH1D* predJet_timeEtaCov = new TH1D("predJet_timeEtaCov","predJet_timeEtaCov",25,-1.,1.);
		//52 - time-phi covariance for jet
		TH1D* predJet_timePhiCov = new TH1D("predJet_timePhiCov","predJet_timePhiCov",25,-0.75,0.75);
		//53 - dr between predicted jets for W candidate
		TH1D* predJet_WjjDr = new TH1D("predJet_WjjDr","predJet dR between W candidate jets",25,0,3.5);
		//54 - dr between reco jets for W candidate
		TH1D* recoJet_WjjDr = new TH1D("recoJet_WjjDr","recoJet dR between W candidate jets",25,0,3.5);
		//55 - pt of W candidate from pred jets
		TH1D* predJet_Wpt = new TH1D("predJet_Wpt","predJet W candidate pt",25,0,250.);
		//56 - pt of W candidate from reco jets
		TH1D* recoJet_Wpt = new TH1D("recoJet_Wpt","recoJet W candidate pt",25,0,250.);


		//2D plots
		//0 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_predGen = new TH2D("jetGenE_diffDeltaPt_predGen","jetGenE_diffDeltaPt_predGen;jet_{gen} E (GeV);#Delta p_{T}_{pred, gen} (GeV)",5,0,100,50,-50,50);

		//1 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_recoGen = new TH2D("jetGenE_diffDeltaPt_recoGen","jetGenE_diffDeltaPt_recoGen;jet_{gen} E (GeV);#Delta p_{T}_{reco, gen} (GeV)",4,&xbins_recoGenPt[0],50,-50,50);
		//2 - 2D histogram of gen pT vs reco pT
		TH2D* genPt_recoPt = new TH2D("genPt_recoPt","genPt_recoPt;genpt;recopt",50,5,50,50,5,50);
		//3 - gen jet mass vs reco/gen pt
		TH2D* genJetMass_recoGenPtRatio = new TH2D("genJetMass_recoGenPtRatio","genJetMass_recoGenPtRatio",20,0,10,20,0.8,1.2);
		//4 - reco jet mass vs reco jet pt
		TH2D* recoJetMass_recoJetPt = new TH2D("recoJetMass_recoJetPt","recoJetMass_recoJetPt;recoJetMass;recoJetPt",50,0,250,50,0,250);
		//5 - reco jet mass vs reco jet jetSize
		TH2D* recoJetMass_recoJetSize = new TH2D("recoJetMass_recoJetSize","recoJetMass_recoJetSize;recoJetMass;recoJetSize",50,0,250,50,0,1);
		//6 - reco m_jj ~ W mass vs jet pair jetSize
		TH2D* recoJetInvMassW_recoJetPairjetSize = new TH2D("recoJetInvMassW_recoJetPairjetSize","recoJetInvMassW_recoJetPairjetSize;reco m_jj;jetSize_jj",50,0,250,50,0,1); 
		//7 - pred jet mass vs pred jet pt
		TH2D* predJetMass_predJetPt = new TH2D("predJetMass_predJetPt","predJetMass_predJetPt;predJetMass;predJetPt",50,0,250,50,0,250);
		//8 - pred jet mass vs pred jet jetSize
		TH2D* predJetMass_predJetSize = new TH2D("predJetMass_predJetSize","predJetMass_predJetSize;predJetMass;predJetSize",50,0,250,50,0,1);
		//9 - pred m_jj ~ W mass vs jet pair jetSize 
		TH2D* predJetInvMassW_predJetPairjetSize = new TH2D("predJetInvMassW_predJetPairjetSize","predJetInvMassW_predJetPairjetSize;pred m_jj;jetSize_jj",50,0,250,50,0,1.); 
		//10 - pred jet pt vs pred jet jetSize
		TH2D* predJetPt_predJetSize = new TH2D("predJetPt_predJetSize","predJetPt_predJetSize;predJetPt;predJetSize",50,0,250,50,0,1);
	
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

		//void FindResonances(vector<Jet>& jets, double& wmass, double& tmass){
		void FindResonances(vector<Jet>& jets, pair<int,int>& wmass_idxs, int& tmass_idx){
			//find W candidates
			double mW = 80.377;
			double diff = 999;
			double mij;
			int j1, j2;
			if(jets.size() < 2){
				wmass_idxs = make_pair(-999,-999);
				tmass_idx = -999;
				return;
			}
			for(int i = 0; i < jets.size(); i++){
				for(int j = i+1; j < jets.size(); j++){
					mij = jets[i].invMass(jets[j]);
					if(fabs(mij - mW) < diff){
						diff = fabs(mij - mW);
						j1 = i;
						j2 = j;
					}
				}
			}
			wmass_idxs = make_pair(j1, j2);
			if(jets.size() < 3){
				tmass_idx = -999;
				return;
			}
			//find top candidates
			double mTop = 172.69;
			mij = -999;
			diff = 999;
			cout << "j1 " << j1 << " make w jet from jet " << jets[j1].px()  << endl;
			Jet wJet = jets[j1];
			wJet.add(jets[j2]);
			for(int i = 0; i < jets.size(); i++){
				if(i != j1 && i != j2){
					mij = jets[i].invMass(wJet);
					if(fabs(mij - mTop) < diff){
						diff = fabs(mij - mTop);
						tmass_idx = i;
					}
				}
			}

		}

	void SetAlpha(double a){_alpha = a;}
	void SetSubclusterAlpha(double a){_emAlpha = a; }
	void SetThreshold(double t){ _thresh = t; }


	private:
		string _oname;
		vector<TH1D*> _hists1D;
		vector<TH2D*> _hists2D;
		vector<TGraph*> graphs;
		vector<Jet> _phos; //photons for event
		vector<procCat> _procCats;
		vector<node*> _trees;
		vector<Jet> _predJets;
		vector<Jet> _genjets;
		vector<Jet> _recojets;
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
		double _radius; //radius of detector set by rhs in event (used for constructing jets)
		double _pvx, _pvy, _pvz;	
		double _alpha, _emAlpha, _thresh;

		double dR(double eta1, double phi1, double eta2, double phi2){
			//phi wraparound
			double dphi = (phi1-phi2);
			dphi = acos(cos(dphi));
			return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
		}
		
};
#endif
