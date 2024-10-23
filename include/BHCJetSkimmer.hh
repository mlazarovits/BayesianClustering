#ifndef BHCJETSKIMMER_HH
#define BHCJETSKIMMER_HH
#include "JetSimProducer.hh"
#include "BaseSkimmer.hh"
#include "BaseTree.hh"
#include "TGraph.h"
#include <set>

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
			_hists1D.push_back(nSubclusters);	
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
			_hists1D.push_back(recoJet_dR_b);
			_hists1D.push_back(BHCJet_dR_b);
			_hists1D.push_back(recoJet_dR_qg);
			_hists1D.push_back(BHCJet_dR_qg);
			_hists1D.push_back(recoJet_dR_lep);
			_hists1D.push_back(BHCJet_dR_lep);
			_hists1D.push_back(recoJet_genOvRecoE_b);
			_hists1D.push_back(BHCJet_genOvRecoE_b);
			_hists1D.push_back(recoJet_genOvRecoE_qg);
			_hists1D.push_back(BHCJet_genOvRecoE_qg);
			_hists1D.push_back(recoJet_genOvRecoE_lep);
			_hists1D.push_back(BHCJet_genOvRecoE_lep);
			_hists1D.push_back(qType_recoJet);
			_hists1D.push_back(qType_BHCJet);
			_hists1D.push_back(reco_nJets_fullHad);
			_hists1D.push_back(reco_nJets_semiLep);
			_hists1D.push_back(reco_nJets_fullLep);
			_hists1D.push_back(BHC_nJets_fullHad);
			_hists1D.push_back(BHC_nJets_semiLep);
			_hists1D.push_back(BHC_nJets_fullLep);


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
			_hists2D.push_back(prednSubclusters_jetSize);
			_hists2D.push_back(recoJet_genOvRecoE_dR_b);
			_hists2D.push_back(BHCJet_genOvRecoE_dR_b);
			_hists2D.push_back(recoJet_genOvRecoE_dR_qg);
			_hists2D.push_back(BHCJet_genOvRecoE_dR_qg);
			_hists2D.push_back(recoJet_genOvRecoE_dR_lep);
			_hists2D.push_back(BHCJet_genOvRecoE_dR_lep);
			_hists2D.push_back(recoJet_dRquark_Wenergy);
			_hists2D.push_back(BHCJet_dRquark_Wenergy);

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
				/*
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
				}*/
				//create new Jet
				Jet predJet(_trees[i]->model, BayesPoint({_pvx, _pvy, _pvz}), _gev, _radius);
				/*
				Jet predJet(rhs, BayesPoint({_pvx, _pvy, _pvz}));
				//cout << "additive p: (" << pxtot << ", " << pytot << ", " << pztot << ", " << e << ")" << " mass " << sqrt((e*e) - (pxtot*pxtot + pytot*pytot + pztot*pztot)) << endl;
				//cout << "jet p: (" << predJet.px() << ", " << predJet.py() << ", " << predJet.pz() << ", " << predJet.e() << ") mass " << predJet.mass() << endl;
				//set PV info
				//predJet.SetVertex(Point({_pvx,_pvy,_pvz}));
				//set constituents (subclusters) here with model from tree
				int nsubclusters = _trees[i]->model->GetNClusters();
				Matrix post = _trees[i]->model->GetPosterior(); //responsibilities of each point to every subcluster
				cout << "posterior has dims " << post.GetDims()[0] << " " << post.GetDims()[1] << endl;
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
				*/
				//put pt cut in for predjets of 20 GeV
				if(predJet.pt() < 20) continue; 
				//add Jet to jets	
				_predJets.push_back(predJet);	
			}
			cout << _predJets.size() << " pred jets" << endl;
			for(auto j : _predJets) cout << "pred jet px " << j.px() << " py " << j.py() << " pz " << j.pz() << " E " << j.E() << " m2 " << j.m2() << " mass " << j.mass() << endl;
		}

		void FillPredJetHists(){
			int njets = _predJets.size();
			int tmass_idx, id, widx;
			double wmass, topmass, dr, dr_pair, jetsize;
			pair<int, int> wmass_idxs;
			Jet w, top, genpart, genW;
			for(int p = 0; p < _procCats.size(); p++){
				//if(p != 0) cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				_procCats[p].hists1D[0][0]->Fill(njets);
				cout << "# pred jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][11]->Fill(njets - (int)_genjets.size());
				FindResonances(_predJets,wmass_idxs,tmass_idx);
				//W mass
				cout << "wmass idxs " << wmass_idxs.first << " " << wmass_idxs.second << endl;
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
					cout << "adding w jet with pt " << w.pt() << endl;
					top = _predJets[tmass_idx];
					top.add(w);
					topmass = top.mass();
					_procCats[p].hists1D[0][30]->Fill(topmass);
				}
				njets = _predJets.size();
				int nsubs;
				for(int j = 0; j < _predJets.size(); j++){
					//if(p != 0) cout << "pred jet #" << j << " phi " << _predJets[j].phi() << " eta " << _predJets[j].eta() << " energy " << _predJets[j].E() <<  " mass " << _predJets[j].mass() << " nConstituents " << _predJets[j].GetNConstituents() << " nRhs " << _predJets[j].GetNRecHits() << " pt " << _predJets[j].pt() << endl;
					dr = CalcJetSize(_predJets[j]);
					cout << "calc jet size for BHC jet " << j << ": " << dr << endl;
					nsubs = _predJets[j].GetNConstituents();
					_procCats[p].hists2D[0][11]->Fill(nsubs,dr);
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
					
					//TODO: update to exclusive gen match
					//do gen matching
					/*
					GenMatchJet(_predJets[j],genpart,id);
					cout << "found gen match to pred jet " << j << " with id " << id << " and dr " << dR(genpart.eta(), genpart.phi(), _predJets[j].eta(), _predJets[j].phi()) << " and gen/pred E ratio " << genpart.E()/_predJets[j].E() << " user idx " << genpart.GetUserIdx() << endl;
					dr = dR(genpart.eta(), genpart.phi(), _predJets[j].eta(), _predJets[j].phi());
					widx = _base->genpart_momIdx->at(genpart.GetUserIdx()); //this branch only filled for products of Ws
					cout << "widx " << widx << endl;
					if(widx != -1) _procCats[p].hists2D[0][19]->Fill(dr,_base->genpart_energy->at(widx));
					//fill gen match hists
					_procCats[p].hists1D[0][70]->Fill(id);
					cout << "fill dr, E ratio plots" << endl;
					//b's
					if(fabs(id) == 5){
						_procCats[p].hists1D[0][58]->Fill(dr);
						_procCats[p].hists1D[0][64]->Fill(genpart.E()/_predJets[j].E());
					}
					//g's, q's (not b's)
					else if(fabs(id) == 1 || fabs(id) == 2 || fabs(id) == 3 || fabs(id) == 4){
						_procCats[p].hists1D[0][60]->Fill(dr);
						_procCats[p].hists1D[0][66]->Fill(genpart.E()/_predJets[j].E());

					}
					//leptons
					else if(fabs(id) == 13 || fabs(id) == 14){
						_procCats[p].hists1D[0][62]->Fill(dr);
						_procCats[p].hists1D[0][68]->Fill(genpart.E()/_predJets[j].E());

					}
					else{ }
					*/
				}
				cout << "fill bhc njet plots" << endl;
				//fill njets based on top decay type
				//0 -> fully had
				//1 -> fully lep
				//2 -> semi lep 
				if(_topDecayType == 0)
					_procCats[p].hists1D[0][74]->Fill(njets);
				else if(_topDecayType == 1)
					_procCats[p].hists1D[0][76]->Fill(njets);
				else if(_topDecayType == 2)
					_procCats[p].hists1D[0][75]->Fill(njets);
				else{ }
			}
		}
	
		void FillRecoJetHists(){
			int njets, tmass_idx, widx, genidx, id;
			double wmass, topmass, dr, dr_pair, jetsize;
			pair<int, int> wmass_idxs;
			Jet w, top, genW;
			//do gen matching
			vector<int> genMatchIdxs; //one per jet, follows same indexing as jets
			GenMatchJet(_recojets,genMatchIdxs);
			cout << "final best matches" << endl;
			for(int b = 0; b < genMatchIdxs.size(); b++){
				if(genMatchIdxs[b] != -1) cout << " jet " << b << " is exclusively matched to gen particle " << genMatchIdxs[b] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs[b]), _base->genpart_phi->at(genMatchIdxs[b]), _recojets[b].eta(), _recojets[b].phi()) << endl;
				 else cout << " jet " << b << " could not be gen matched" << endl;

			}
			
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = _recojets.size();
				for(int j = 0; j < _recojets.size(); j++){
					jetsize = CalcJetSize(_recojets[j]);
					if(p == 0) cout << "calc reco jet size " << j << ": " << jetsize << endl;
					if(p == 0) cout << "reco jet #" << j << " phi " << _recojets[j].phi() << " eta " << _recojets[j].eta() << " energy " << _recojets[j].E() <<  " mass " << _recojets[j].mass() << " nConstituents " << _recojets[j].GetNConstituents() << " nRhs " << _recojets[j].GetNRecHits() << " pt " << _recojets[j].pt() << " jetsize " << jetsize << endl;
					_procCats[p].hists1D[0][19]->Fill(jetsize);
					_procCats[p].hists1D[0][20]->Fill(_recojets[j].e());
					_procCats[p].hists1D[0][21]->Fill(_recojets[j].pt());
					_procCats[p].hists1D[0][22]->Fill(_recojets[j].mass());
					_procCats[p].hists2D[0][4]->Fill(_recojets[j].mass(), _recojets[j].pt());
					if(jetsize != -999) _procCats[p].hists2D[0][5]->Fill(_recojets[j].mass(), jetsize);
			
					//fill gen match histograms
					//if no gen match, skip
					if(genMatchIdxs[j] == -1) continue;
					genidx = genMatchIdxs[j];	
					if(p == 0) cout << "found gen match to reco jet " << j << " with id " << _base->genpart_id->at(genidx) << " and dr " << dR(_base->genpart_eta->at(genidx), _base->genpart_phi->at(genidx), _recojets[j].eta(), _recojets[j].phi()) << " and gen/reco E ratio " << _base->genpart_energy->at(genidx)/_recojets[j].E() << " user idx " << genidx << endl;
					//get matched gen particle here
					dr = dR(_base->genpart_eta->at(genidx), _base->genpart_phi->at(genidx), _recojets[j].eta(), _recojets[j].phi());
					id = _base->genpart_id->at(genidx);
					widx = _base->genpart_momIdx->at(genidx); //this branch only filled for products of Ws
					if(p == 0)cout << "widx " << widx << endl;
					if(widx != -1) _procCats[p].hists2D[0][18]->Fill(dr,_base->genpart_energy->at(widx));
					//fill gen match hists
					_procCats[p].hists1D[0][69]->Fill(id);
					//b's
					if(fabs(id) == 5){
						_procCats[p].hists1D[0][57]->Fill(dr);
						_procCats[p].hists1D[0][63]->Fill(_base->genpart_energy->at(genidx)/_recojets[j].E());
					}
					//g's, q's (not b's)
					else if(fabs(id) == 1 || fabs(id) == 2 || fabs(id) == 3 || fabs(id) == 4){
						_procCats[p].hists1D[0][59]->Fill(dr);
						_procCats[p].hists1D[0][65]->Fill(_base->genpart_energy->at(genidx)/_recojets[j].E());

					}
					//leptons
					else if(fabs(id) == 13 || fabs(id) == 14){
						_procCats[p].hists1D[0][61]->Fill(dr);
						_procCats[p].hists1D[0][67]->Fill(_base->genpart_energy->at(genidx)/_recojets[j].E());

					}
					else{ }
				}
				_procCats[p].hists1D[0][18]->Fill(njets);
				//fill njets based on top decay type
				//0 -> fully had
				//1 -> fully lep
				//2 -> semi lep 
				if(_topDecayType == 0)
					_procCats[p].hists1D[0][71]->Fill(njets);
				else if(_topDecayType == 1)
					_procCats[p].hists1D[0][73]->Fill(njets);
				else if(_topDecayType == 2)
					_procCats[p].hists1D[0][72]->Fill(njets);
				else{ }

				//cout << "hist name " << _procCats[p].hists1D[0][18]->GetName() << " nentries " << _procCats[p].hists1D[0][18]->GetEntries() << endl;
				if(p == 0)cout << "# reco jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][24]->Fill(njets - (int)_genjets.size());


				FindResonances(_recojets,wmass_idxs,tmass_idx);
				if(p == 0)cout << "wmass idxs " << wmass_idxs.first << " " << wmass_idxs.second << " top idx " << tmass_idx << endl;
				if(wmass_idxs.first != -999){
					w = _recojets[wmass_idxs.first];
					w.add(_recojets[wmass_idxs.second]);
					wmass = w.mass();
					dr_pair = dR(_recojets[wmass_idxs.first].eta(), _recojets[wmass_idxs.first].phi(), _recojets[wmass_idxs.second].eta(), _recojets[wmass_idxs.second].phi());
					_procCats[p].hists1D[0][27]->Fill(wmass);
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
			double jetsize;
			//if no subclusters (ie reco jet), take jet size to be sqrt(etavar + phivar) from discrete covariance 
			if(nSCs < 1){
				Matrix recocov = Matrix(3,3);
				vector<JetPoint> rhs = jet.GetJetPoints();
				int nrhs = rhs.size();
				
				vector<double> rhv, jetv;
				jetv = {jet.eta(), jet.phi_std(), jet.time()};
				double diffi, diffj;
				for(int r = 0; r < nrhs; r++){
					rhv = {rhs[r].eta(), rhs[r].phi(), rhs[r].time()};
					for(int i = 0; i < 3; i++){
						for(int j = 0; j < 3; j++){
							//phi wraparound
							diffi = rhv[i] - jetv[i];
							diffj = rhv[j] - jetv[j];
							if(i == 1){
								if(diffi > 4*atan(1)) diffi -= 8*atan(1);
							}
							if(j == 1){
								if(diffj > 4*atan(1)) diffj -= 8*atan(1);
							}
						
							recocov.SetEntry( recocov.at(i,j) + diffi*diffj, i, j );
						}
					}
				}
				recocov.mult(recocov,1./nrhs);
				jetsize = sqrt(recocov.at(0,0) + recocov.at(1,1));
				/* //old reco jet size definition
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
				*/
			}
			else{
				Matrix cov, mu;
				jet.GetClusterParams(mu, cov);
				cout << "pred eta var " << cov.at(0,0) << " phi var " << cov.at(1,1) << endl;
				jetsize = sqrt(cov.at(0,0) + cov.at(1,1));
			}	
			return jetsize;
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
		TH1D* nClusters = new TH1D("nBHCJets","nBHCJets",10,0,10);
		//1
		TH1D* nSubclusters = new TH1D("BHCJet_nSubclusters","BHCJet_nSubclusters",10,0,10);
		//2
		TH1D* predJet_subClusterEnergy = new TH1D("BHCJet_subClusterEnergy","BHCJet_subClusterEnergy",20,0,500);
		//3
		TH1D* predJet_subClusterTimeCenter = new TH1D("BHCJet_subClusterTimeCenter","BHCJet_subClusterTimeCenter",25,0,15);
		//4
		TH1D* predJet_subClusterEtaCenter = new TH1D("BHCJet_subClusterEtaCenter","BHCJet_subClusterEtaCenter",25,-1.8,1.8);
		//5
		TH1D* predJet_subClusterPhiCenter = new TH1D("BHCJet_subClusterPhiCenter","BHCJet_subClusterPhiCenter",25,-0.1,6.3);
		//6
		TH1D* predJet_jetSize = new TH1D("BHCJet_jetSize","BHCJet_jetSize",50,0,1);
		//7
		TH1D* predJet_energy = new TH1D("BHCJet_energy","BHCJet_energy",50,0,400);
		//8
		TH1D* predJet_pt = new TH1D("BHCJet_pt","BHCJet_pt",50,0,500);
		//9
		TH1D* predJet_mass = new TH1D("BHCJet_mass","BHCJet_mass",50,0,150);
		//10 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_predGen = new TH1D("jetGenE_sigmaDeltaPt_predGen","jetGenE_sigmaDeltaPt_predGen",5,0,100);
		//11 - # pred jets - # gen jets
		TH1D* predGen_nJets = new TH1D("BHCGen_diffNJets","BHCGen_diffNJets",20,-10,10);
		//for subclusters
		//12 - eta sigma
		TH1D* predJet_subClusterEtaVar = new TH1D("BHCJet_subClusterEtaVar","BHCJet_subClusterEtaVar",25,0.,3.);
		//13 - phi sigma
		TH1D* predJet_subClusterPhiVar = new TH1D("BHCJet_subClusterPhiVar","BHCJet_subClusterPhiVar",25,0.,3.);
		//14 - time sigma
		TH1D* predJet_subClusterTimeVar = new TH1D("BHCJet_subClusterTimeVar","BHCJet_subClusterTimeVar",25,0.,5.);
		//15 - eta-phi covariance
		TH1D* predJet_subClusteretaPhiCov = new TH1D("BHCJet_subClusteretaPhiCov","BHCJet_subClusteretaPhiCov",25,-1.5,1.5);
		//16 - time-eta covariance
		TH1D* predJet_subClustertimeEtaCov = new TH1D("BHCJet_subClustertimeEtaCov","BHCJet_subClustertimeEtaCov",25,-1.5,1.5);
		//17 - time-phi covariance
		TH1D* predJet_subClustertimePhiCov = new TH1D("BHCJet_subClustertimePhiCov","BHCJet_subClustertimePhiCov",25,-1.5,1.5);
		//18
		TH1D* nRecoJets = new TH1D("nRecoJets","nRecoJets",10,0,10);
		//19
		TH1D* recoJet_jetSize = new TH1D("recoJet_jetSize","recoJet_jetSize",50,0,1);
		//20
		TH1D* recoJet_energy = new TH1D("recoJet_energy","recoJet_energy",50,0,400);
		//21
		TH1D* recoJet_pt = new TH1D("recoJet_pt","recoJet_pt",50,0,500);
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
		TH1D* recoJet_Wmass = new TH1D("recoJet_Wmass","recoJet_Wmass;Invariant mass for best W candidate",50,0,500);
		//28 - reco jet top invariant mass
		TH1D* recoJet_topmass = new TH1D("recoJet_topmass","recoJet_topmass;Invariant mass for best top candidate",50,0,500);
		//29 - pred jet W invariant mass
		TH1D* predJet_Wmass = new TH1D("BHCJet_Wmass","BHCJet_Wmass;Invariant mass for best W candidate",50,0,500);
		//30 - pred jet top invariant mass
		TH1D* predJet_topmass = new TH1D("BHCJet_topmass","BHCJet_topmass;Invariant mass for best top candidate",50,0,500);
		//31 - pred jet W invariant mass, pT_jj < 100
		TH1D* predJet_Wmass_pTjjl100 = new TH1D("BHCJet_Wmass_pTjjl100","BHCJet_Wmass_pTjjl100;Invariant mass_pTjjl100 for best W candidate",50,0,500);
		//32 - pred jet W invariant mass, pT_jj >= 100
		TH1D* predJet_Wmass_pTjjge100 = new TH1D("BHCJet_Wmass_pTjjge100","BHCJet_Wmass_pTjjge100;Invariant mass_pTjjge100 for best W candidate",50,0,500);
		//33 - reco jet W invariant mass, pT_jj < 100
		TH1D* recoJet_Wmass_pTjjl100 = new TH1D("recoJet_Wmass_pTjjl100","recoJet_Wmass_pTjjl100;Invariant mass_pTjjl100 for best W candidate",50,0,500);
		//34 - reco jet W invariant mass, pT_jj >= 100
		TH1D* recoJet_Wmass_pTjjge100 = new TH1D("recoJet_Wmass_pTjjge100","recoJet_Wmass_pTjjge100;Invariant mass_pTjjge100 for best W candidate",50,0,500);
		//35 - W jet pair eta variance
		TH1D* predJet_etaVar_Wjj = new TH1D("BHCJet_etaVar_Wjj","BHCJet_etaVar_Wjj",50,0,1);
		//36 - W jet eta variance, pT < 100 (resolved)
		TH1D* predJet_etaVar_W_pTl100 = new TH1D("BHCJet_etaVar_W_pTl100","BHCJet_etaVar_W_pTl100",50,0,1);
		//37 - W jet eta variance, pT >= 100 (boosted)
		TH1D* predJet_etaVar_W_pTge100 = new TH1D("BHCJet_etaVar_W_pTge100","BHCJet_etaVar_W_pTge100",50,0,1);
		//38 - !W jet pair eta variance
		TH1D* predJet_etaVar_notWjj = new TH1D("BHCJet_etaVar_notWjj","BHCJet_etaVar_notWjj",50,0,1);
		//39 - W jet pair phi variance
		TH1D* predJet_phiVar_Wjj = new TH1D("BHCJet_phiVar_Wjj","BHCJet_phiVar_Wjj",50,0,1);
		//40 - W jet phi variance, pT < 100 (resolved)
		TH1D* predJet_phiVar_W_pTl100 = new TH1D("BHCJet_phiVar_W_pTl100","BHCJet_phiVar_W_pTl100",50,0,1);
		//41 - W jet phi variance, pT >= 100 (boosted)
		TH1D* predJet_phiVar_W_pTge100 = new TH1D("BHCJet_phiVar_W_pTge100","BHCJet_phiVar_W_pTge100",50,0,1);
		//42 - !W jet pair phi variance
		TH1D* predJet_phiVar_notWjj = new TH1D("BHCJet_phiVar_notWjj","BHCJet_phiVar_notWjj",50,0,1);
		//43 - W jet pair time variance
		TH1D* predJet_timeVar_Wjj = new TH1D("BHCJet_timeVar_Wjj","BHCJet_timeVar_Wjj",50,0,1);
		//44 - W jet time variance, pT < 100 (resolved)
		TH1D* predJet_timeVar_W_pTl100 = new TH1D("BHCJet_timeVar_W_pTl100","BHCJet_timeVar_W_pTl100",50,0,1);
		//45 - W jet time variance, pT >= 100 (boosted)
		TH1D* predJet_timeVar_W_pTge100 = new TH1D("BHCJet_timeVar_W_pTge100","BHCJet_timeVar_W_pTge100",50,0,1);
		//46 - !W jet pair time variance
		TH1D* predJet_timeVar_notWjj = new TH1D("BHCJet_timeVar_notWjj","BHCJet_timeVar_notWjj",50,0,1);
		//47 - eta sigma for jet
		TH1D* predJet_EtaVar = new TH1D("BHCJet_EtaVar","BHCJet_EtaVar",25,0.,1.);
		//48 - phi sigma for jet
		TH1D* predJet_PhiVar = new TH1D("BHCJet_PhiVar","BHCJet_PhiVar",25,0., 1.);
		//49 - time sigma for jet
		TH1D* predJet_TimeVar = new TH1D("BHCJet_TimeVar","BHCJet_TimeVar",25,0.,1.);
		//50 - eta-phi covariance for jet
		TH1D* predJet_etaPhiCov = new TH1D("BHCJet_etaPhiCov","BHCJet_etaPhiCov",25,-1.,1.);
		//51 - time-eta covariance for jet
		TH1D* predJet_timeEtaCov = new TH1D("BHCJet_timeEtaCov","BHCJet_timeEtaCov",25,-1.,1.);
		//52 - time-phi covariance for jet
		TH1D* predJet_timePhiCov = new TH1D("BHCJet_timePhiCov","BHCJet_timePhiCov",25,-0.75,0.75);
		//53 - dr between predicted jets for W candidate
		TH1D* predJet_WjjDr = new TH1D("BHCJet_WjjDr","BHCJet dR between W candidate jets",25,0,3.5);
		//54 - dr between reco jets for W candidate
		TH1D* recoJet_WjjDr = new TH1D("recoJet_WjjDr","recoJet dR between W candidate jets",25,0,3.5);
		//55 - pt of W candidate from pred jets
		TH1D* predJet_Wpt = new TH1D("BHCJet_Wpt","BHCJet W candidate pt",25,0,250.);
		//56 - pt of W candidate from reco jets
		TH1D* recoJet_Wpt = new TH1D("recoJet_Wpt","recoJet W candidate pt",25,0,250.);
		//57 - reco dr match to gen b's
		TH1D* recoJet_dR_b = new TH1D("recoJet_dR_b","recoJet_dR_b",25,0,1.);
		//58 - bhc dr match to gen b's
		TH1D* BHCJet_dR_b = new TH1D("BHCJet_dR_b","BHCJet_dR_b",25,0,1.);
		//59 - reco dr match to gen q's/g's
		TH1D* recoJet_dR_qg = new TH1D("recoJet_dR_qg","recoJet_dR_qg",25,0,1.);
		//60 - bhc dr match to gen q's/g's
		TH1D* BHCJet_dR_qg = new TH1D("BHCJet_dR_qg","BHCJet_dR_qg",25,0,1.);
		//61 - reco dr match to gen l's  
		TH1D* recoJet_dR_lep = new TH1D("recoJet_dR_lep","recoJet_dR_lep",25,0,1.);
		//62 - bhc dr match to gen l's  
		TH1D* BHCJet_dR_lep = new TH1D("BHCJet_dR_lep","BHCJet_dR_lep",25,0,1.);
		//63 - reco dr match to gen b's
		TH1D* recoJet_genOvRecoE_b = new TH1D("recoJet_genOvRecoE_b","recoJet_genOvRecoE_b",25,0,5.);
		//64 - bhc dr match to gen b's
		TH1D* BHCJet_genOvRecoE_b = new TH1D("BHCJet_genOvRecoE_b","BHCJet_genOvRecoE_b",25,0,5.);
		//65 - reco dr match to gen q's/g's
		TH1D* recoJet_genOvRecoE_qg = new TH1D("recoJet_genOvRecoE_qg","recoJet_genOvRecoE_qg",25,0,5.);
		//66 - bhc dr match to gen q's/g's
		TH1D* BHCJet_genOvRecoE_qg = new TH1D("BHCJet_genOvRecoE_qg","BHCJet_genOvRecoE_qg",25,0,5.);
		//67 - reco dr match to gen l's  
		TH1D* recoJet_genOvRecoE_lep = new TH1D("recoJet_genOvRecoE_lep","recoJet_genOvRecoE_lep",25,0,5.);
		//68 - bhc dr match to gen l's  
		TH1D* BHCJet_genOvRecoE_lep = new TH1D("BHCJet_genOvRecoE_lep","BHCJet_genOvRecoE_lep",25,0,5.);
		//69 - types of quarks in reco jets - d, u, s, c, b (sign agnostic)
		TH1D* qType_recoJet = new TH1D("qType_recoJet","qType_recoJet",5,1,6);
		//70 - types of quarks in BHC jets - d, u, s, c, b (sign agnostic)
		TH1D* qType_BHCJet = new TH1D("qType_BHCJet","qType_BHCJet",5,1,6);
		//71 - n reco jets, tt fully had
		TH1D* reco_nJets_fullHad = new TH1D("reco_nJets_fullHad","reco_nJets_fullHad",10,0,10);
		//72 - n reco jets, tt semi lep
		TH1D* reco_nJets_semiLep = new TH1D("reco_nJets_semiLep","reco_nJets_semiLep",10,0,10);
		//73 - n reco jets, tt fully lep
		TH1D* reco_nJets_fullLep = new TH1D("reco_nJets_fullLep","reco_nJets_fullLep",10,0,10);
		//74 - n BHC jets, tt fully had
		TH1D* BHC_nJets_fullHad = new TH1D("BHC_nJets_fullHad","BHC_nJets_fullHad",10,0,10);
		//75 - n BHC jets, tt semi lep
		TH1D* BHC_nJets_semiLep = new TH1D("BHC_nJets_semiLep","BHC_nJets_semiLep",10,0,10);
		//76 - n BHC jets, tt fully lep
		TH1D* BHC_nJets_fullLep = new TH1D("BHC_nJets_fullLep","BHC_nJets_fullLep",10,0,10);


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
		TH2D* recoJetInvMassW_recoJetPairjetSize = new TH2D("recoJetInvMassW_recoJetPairjetSize","recoJetInvMassW_recoJetPairjetSize;reco m_jj;jetSize_jj",50,0,500,50,0,1); 
		//7 - pred jet mass vs pred jet pt
		TH2D* predJetMass_predJetPt = new TH2D("BHCJetMass_predJetPt","BHCJetMass_predJetPt;predJetMass;predJetPt",50,0,250,50,0,250);
		//8 - pred jet mass vs pred jet jetSize
		TH2D* predJetMass_predJetSize = new TH2D("BHCJetMass_predJetSize","BHCJetMass_predJetSize;predJetMass;predJetSize",50,0,250,50,0,1);
		//9 - pred m_jj ~ W mass vs jet pair jetSize 
		TH2D* predJetInvMassW_predJetPairjetSize = new TH2D("BHCJetInvMassW_predJetPairjetSize","BHCJetInvMassW_predJetPairjetSize;pred m_jj;jetSize_jj",50,0,500,50,0,1.); 
		//10 - pred jet pt vs pred jet jetSize
		TH2D* predJetPt_predJetSize = new TH2D("BHCJetPt_predJetSize","BHCJetPt_predJetSize;predJetPt;predJetSize",50,0,250,50,0,1);
		//11 - pred jet n subclusters vs jet size
		TH2D* prednSubclusters_jetSize = new TH2D("BHCnSubclusters_jetSize","BHCnSubclusters_jetSize;nSubclusters;jetsize",10,0,10,50,0,1);
		//12 - reco dr match to gen b's
		TH2D* recoJet_genOvRecoE_dR_b = new TH2D("recoJet_genOvRecoE_dR_b","recoJet_genOvRecoE_dR_b;ratioE;dR",25,0,5,25,0,4);
		//13 - bhc dr match to gen b's
		TH2D* BHCJet_genOvRecoE_dR_b = new TH2D("BHCJet_genOvRecoE_dR_b","BHCJet_genOvRecoE_dR_b;ratioE;dR",25,0,5,25,0,4);
		//14 - reco dr match to gen q's/g's
		TH2D* recoJet_genOvRecoE_dR_qg = new TH2D("recoJet_genOvRecoE_dR_qg","recoJet_genOvRecoE_dR_qg;ratioE;dR",25,0,5,25,0,4);
		//15 - bhc dr match to gen q's/g's
		TH2D* BHCJet_genOvRecoE_dR_qg = new TH2D("BHCJet_genOvRecoE_dR_qg","BHCJet_genOvRecoE_dR_qg;ratioE;dR",25,0,5,25,0,4);
		//16 - reco dr match to gen l's  
		TH2D* recoJet_genOvRecoE_dR_lep = new TH2D("recoJet_genOvRecoE_dR_lep","recoJet_genOvRecoE_dR_lep;ratioE;dR",25,0,5,25,0,4);
		//17 - bhc dr match to gen l's  
		TH2D* BHCJet_genOvRecoE_dR_lep = new TH2D("BHCJet_genOvRecoE_dR_lep","BHCJet_genOvRecoE_dR_lep;ratioE;dR",25,0,5,25,0,4);
		//18 - reco dr jet-quark match vs W energy (b's excluded)
		TH2D* recoJet_dRquark_Wenergy = new TH2D("recoJet_dRquark_Wenergy","recoJet_dRquark_Wenergy;dRquark;Wenergy",25,0,4,0,1000);
		//19 - BHC dr jet-quark match vs W energy (b's excluded)
		TH2D* BHCJet_dRquark_Wenergy = new TH2D("BHCJet_dRquark_Wenergy","BHCJet_dRquark_Wenergy;dRquark;Wenergy",25,0,4,0,1000);

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

		void FindResonances(vector<Jet>& jets, pair<int,int>& wmass_idxs, int& tmass_idx){
			//find W candidates
			double mW = 80.377;
			double diff = 1e6;
			double mij;
			int j1, j2;
			if(jets.size() < 2){
				wmass_idxs = make_pair(-999,-999);
				tmass_idx = -999;
				return;
			}
			diff = fabs(jets[0].invMass(jets[1]) - diff);
			j1 = 0;
			j2 = 1;
			for(int i = 1; i < jets.size(); i++){
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
			Jet wJet = jets[j1];
			wJet.add(jets[j2]);
			diff = fabs(jets[2].invMass(wJet) - mTop);
			tmass_idx = 2;
			for(int i = 2; i < jets.size(); i++){
				if(i != j1 && i != j2){
					mij = jets[i].invMass(wJet);
					if(fabs(mij - mTop) < diff){
						diff = fabs(mij - mTop);
						tmass_idx = i;
					}
				}
			}

		}



	//find gen particle that most closely is dr matched to jet
	void GenMatchJet(vector<Jet>& jets, vector<int>& bestGenMatchIdxs){
		//loop through gen particles
		//dr match to jet
		double bestDr, dr;
		int nGen = _base->genpart_ngenpart;
		int otherJet, thisgenidx, thisJet;
		bestGenMatchIdxs.clear();
		bestGenMatchIdxs = {};

		if(jets.size() < 1) return;

		//vector<double> drs;
		//drs[i][j] = dr for jet i and gen particle j
		vector<vector<double>> drs;
		vector<int> genIdxs;
		for(int j = 0; j < jets.size(); j++){
			drs.push_back({});
			for(int g = 0; g < nGen; g++){
				//skip W's
				drs[j].push_back(999);
				if(fabs(_base->genpart_id->at(g)) == 24) continue;
				//skip neutrinos
				if(fabs(_base->genpart_id->at(g)) == 14 || fabs(_base->genpart_id->at(g)) == 12 || fabs(_base->genpart_id->at(g)) == 16) continue;
				//rough acceptance cut - if gen particles are out of acceptance any associated jets wouldn't be reconstructed
				if(fabs(_base->genpart_eta->at(g)) > 1.5) continue;
			
				dr = dR(_base->genpart_eta->at(g), _base->genpart_phi->at(g), jets[j].eta(), jets[j].phi());
				drs[j][g] = dr;
			}
		}
		vector<int> best_idxs; //one per jet
		for(int j = 0; j < jets.size(); j++){
			for(int g = 0; g < nGen; g++){
				cout << "jet " << j << " and gen particle " << g << " have dr " << drs[j][g] << endl;
			}
			cout << "jet " << j << " has best dr " << *min_element(drs[j].begin(), drs[j].end()) << " at gen particle " << find(drs[j].begin(), drs[j].end(), *min_element(drs[j].begin(), drs[j].end())) - drs[j].begin() << endl;
			double mindr = *min_element(drs[j].begin(), drs[j].end());
			int genidx = find(drs[j].begin(), drs[j].end(), mindr) - drs[j].begin();
			if(mindr == 999) genidx = -1; //no match found (ie no available gen particle for best match)
			best_idxs.push_back(genidx);
			thisgenidx = genidx;
			thisJet = j;
			while(count(best_idxs.begin(), best_idxs.end(), genidx) > 1){
				otherJet = find(best_idxs.begin(), best_idxs.end(), genidx) - best_idxs.begin();
				//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
				//skip otherJet (ie thisJet) in this case
				if(otherJet == thisJet){
					otherJet = find(best_idxs.begin()+otherJet+1, best_idxs.end(), genidx) - best_idxs.begin();
				}
				for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
				cout << " found another match at jet " << otherJet << " with other dr " << drs[otherJet][genidx] << " against this jet " << thisJet << endl;
				//if other dr is less than this dr
				if(drs[otherJet][genidx] < mindr){
					//set this dr to 999 (is invalid), find new min, reset genidx to this index
					drs[thisJet][genidx] = 999;
					mindr = *min_element(drs[thisJet].begin(), drs[thisJet].end());
					if(mindr == 999) genidx = -1;
					else genidx = find(drs[thisJet].begin(), drs[thisJet].end(), mindr) - drs[thisJet].begin();
					best_idxs[thisJet] = genidx;
					cout << " reset gen match of this jet " << thisJet << " to particle " << genidx << " with dr " << mindr << endl;
		
				}
				else{
					//set other dr to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
					drs[otherJet][genidx] = 999;
					thisgenidx = genidx;
					mindr = *min_element(drs[otherJet].begin(), drs[otherJet].end());
					if(mindr == 999) genidx = -1;
					else genidx = find(drs[otherJet].begin(), drs[otherJet].end(), mindr) - drs[otherJet].begin();
					thisJet = otherJet;
					best_idxs[thisJet] = genidx;
					cout << " reset gen match of other jet " << otherJet << " to particle " << genidx << " with dr " << mindr << endl;
				}	
				cout << "genidx is now " << genidx << " with count " << count(best_idxs.begin(), best_idxs.end(), genidx) << " for jet " << thisJet << endl;

			}
			cout << "jet " << j << " has best exclusive gen match with " << best_idxs[j] << "\n" << endl;
		}
		bestGenMatchIdxs = best_idxs;

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
		//top decay info - 0 = fully had, 1 = semi lep, 2 = fully lep
		int _topDecayType;


		double dR(double eta1, double phi1, double eta2, double phi2){
			//phi wraparound
			double dphi = (phi1-phi2);
			dphi = acos(cos(dphi));
			return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
		}
		
};
#endif
