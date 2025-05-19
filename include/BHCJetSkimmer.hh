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
			_smear = false;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_thresh = 1.;
			
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);
						
			_cell = 0;
			_tresCte = 0;
			_tresNoise = 0;
			_tresStoch = 0;
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
			_smear = false;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_thresh = 1.;
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);
			
			_cell = acos(-1)/180;
			_tresCte = 0.1727;//times given in ns//0.133913 * 1e-9;
			_tresStoch = 0.5109;//1.60666 * 1e-9; 
			_tresNoise = 2.106;//0.00691415 * 1e-9;
				
	
			graphs.push_back(nrhs_comptime);
			graphs.push_back(nrhs_comptime_subcl);

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
			_hists1D.push_back(predJet_subClusterEtaSig);
			_hists1D.push_back(predJet_subClusterPhiSig);
			_hists1D.push_back(predJet_subClusterTimeSig);
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
			_hists1D.push_back(reco_nSubclusters);
			_hists1D.push_back(recoJet_subClusterEnergy);
			_hists1D.push_back(recoJet_subClusterTimeCenter);
			_hists1D.push_back(recoJet_subClusterEtaCenter);
			_hists1D.push_back(recoJet_subClusterPhiCenter);
			_hists1D.push_back(recoJet_subClusterEtaSig);
			_hists1D.push_back(recoJet_subClusterPhiSig);
			_hists1D.push_back(recoJet_subClusterTimeSig);
			_hists1D.push_back(recoJet_subClusteretaPhiCov);
			_hists1D.push_back(recoJet_subClustertimeEtaCov);
			_hists1D.push_back(recoJet_subClustertimePhiCov);
			_hists1D.push_back(recoJet_nRhs);
			_hists1D.push_back(recoJet_rhE);
			_hists1D.push_back(AK4Jet_GenP);
			_hists1D.push_back(AK4Jet_GenPt);
			_hists1D.push_back(AK4JetConstituent_GenP);
			_hists1D.push_back(AK4JetConstituent_GenPt);
			_hists1D.push_back(AK4JetConstJetRatio_GenP);
			_hists1D.push_back(AK4JetConstJetRatio_GenPt);
			_hists1D.push_back(AK4Jet_rhTimes);
			_hists1D.push_back(geoEavg_sigmaDeltaTime_adjRhs);
			_hists1D.push_back(AK4Jet_subClusterEtaCenter_rStat);
			_hists1D.push_back(AK4Jet_subClusterPhiCenter_rStat);
			_hists1D.push_back(AK4Jet_subClusterTimeCenter_rStat);
			_hists1D.push_back(AK4Jet_subClusterEtaSig_rStat);
			_hists1D.push_back(AK4Jet_subClusterPhiSig_rStat);
			_hists1D.push_back(AK4Jet_subClusterTimeSig_rStat);
			_hists1D.push_back(recoJet_subClusteretaPhiCovNorm);
			_hists1D.push_back(recoJet_subClustertimeEtaCovNorm);
			_hists1D.push_back(recoJet_subClustertimePhiCovNorm);
			_hists1D.push_back(nGenParticles);
			_hists1D.push_back(genParticle_eta);
			_hists1D.push_back(genParticle_phi);
			_hists1D.push_back(genParticle_time);
			_hists1D.push_back(genParticle_pt);
			_hists1D.push_back(genParticle_mass);
			_hists1D.push_back(genParticle_energy);
			_hists1D.push_back(genAK4Jet_eta);
			_hists1D.push_back(genAK4Jet_phi);
			_hists1D.push_back(genAK4Jet_time);
			_hists1D.push_back(genAK4Jet_pt);
			_hists1D.push_back(genAK4Jet_mass);
			_hists1D.push_back(genAK4Jet_energy);
			_hists1D.push_back(nJet_genAK4Jet);
			_hists1D.push_back(genAK4Jet_nConstituents);
			_hists1D.push_back(genAK4JetParticle_nDiff);
			_hists1D.push_back(genAK4JetParticle_dR);
			_hists1D.push_back(genAK4JetParticle_Eratio);
			_hists1D.push_back(BHCJetParticle_nDiff);
			_hists1D.push_back(BHCJetParticle_dR);
			_hists1D.push_back(BHCJetParticle_Eratio);
			_hists1D.push_back(BHCJet_EtaCenter);
			_hists1D.push_back(BHCJet_PhiCenter);
			_hists1D.push_back(BHCJet_TimeCenter);
			_hists1D.push_back(rhTime);
			_hists1D.push_back(BHCJet_rhEtaSig);
			_hists1D.push_back(BHCJet_rhPhiSig);
			_hists1D.push_back(BHCJet_rhTimeSig);
			_hists1D.push_back(recoAK4Jet_rhEtaSig);
			_hists1D.push_back(recoAK4Jet_rhPhiSig);
			_hists1D.push_back(recoAK4Jet_rhTimeSig);
			_hists1D.push_back(recoAK4Jet_TimeCenter);
			_hists1D.push_back(recoAK4Jet_EtaCenter);
			_hists1D.push_back(recoAK4Jet_PhiCenter);
			_hists1D.push_back(recoAK4Jet_nSubclustersEvt);
			_hists1D.push_back(BHCJet_nSubclustersEvt);
			_hists1D.push_back(BHCJet_drSubclusters);
			_hists1D.push_back(recoAK4Jet_drSubclusters);
			_hists1D.push_back(BHCJet_rhE);
			_hists1D.push_back(recoAK4Jet_rotundity);
			_hists1D.push_back(BHCJet_rotundity);

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
			_hists2D.push_back(AK4Jet_nRhs_nSubclustersJet);
			_hists2D.push_back(AK4Jet_nGenParts_nSubclusters);
			_hists2D.push_back(AK4Jet_genP_genJetP);
			_hists2D.push_back(AK4Jet_genPt_genJetPt);
			_hists2D.push_back(AK4Jet_genP_genPartJetPRatio);
			_hists2D.push_back(AK4Jet_genJetP_genPartJetPRatio);
			_hists2D.push_back(AK4Jet_genPt_genPartJetPtRatio);
			_hists2D.push_back(AK4Jet_genJetPt_genPartJetPtRatio);
			_hists2D.push_back(AK4Jet_genJetPt_nSubclusters);
			_hists2D.push_back(AK4Jet_nGenPartsptge5_nSubclusters);
			_hists2D.push_back(geoEavg_diffDeltaTime_adjRhs);
			_hists2D.push_back(recoJet_subClusteretaPhiCov_timeEtaCov);
			_hists2D.push_back(recoJet_subClusteretaPhiCovNorm_timeEtaCovNorm);
			_hists2D.push_back(recoAK4Jet_nJets_jetSize);
			_hists2D.push_back(BHCJet_nJets_jetSize);
			_hists2D.push_back(BHCJet_nRhs_nSubclustersJet);
			_hists2D.push_back(recoAK4Jet_nSubclustersJet_nGenConstituents);
			_hists2D.push_back(recoAK4Jet_nSubclustersJet_mass);
			_hists2D.push_back(recoAK4Jet_nSubclustersJet_energy);
			_hists2D.push_back(recoAK4Jet_nSubclustersEvt_nJet);
			_hists2D.push_back(BHCJet_nSubclustersJet_mass);
			_hists2D.push_back(BHCJet_nSubclustersJet_energy);
			_hists2D.push_back(BHCJet_nSubclustersEvt_nJet);
			_hists2D.push_back(recoAK4JetnSubclustersJet_BHCJetnSubclustersJet);	
			_hists2D.push_back(recoAK4Jet_jetEnergy_jetMass);
			_hists2D.push_back(BHCJet_jetEnergy_jetMass);
			_hists2D.push_back(recoAK4Jet_nSubclusters_jetSize);

		}
		void SetMinRhE(double r){ _prod->SetMinRhE(r); }
		void SetRecoMinPt(double r){ _prod->SetRecoMinPt(r); }
		void SetRecoMinE(double r){ _prod->SetRecoMinE(r); }
		void SetGenMinPt(double r){ _prod->SetGenMinPt(r); }
		void SetGenMinE(double r){ _prod->SetGenMinE(r); }
		void SetMinNrhs(int r){ _prod->SetMinNrhs(r); }
		void SetMinNGenConsts(int r){ _prod->SetMinNGenConsts(r); }
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else if(i == 2) _strategy = gmmOnly;
			else if(i == 3) _strategy = NlnNonAK4;
			else return; 
		}
	



		void TreesToJets(){
			_predJets.clear();
			vector<JetPoint> rhs;
			double x, y, z, eta, phi, t, theta, px, py, pz;
			BayesPoint vertex({_pvx, _pvy, _pvz});
			int njets_tot = 0;
			for(int i = 0; i < _trees.size(); i++){
				//get points from tree
				PointCollection* pc = _trees[i]->points;
				//at least 2 points (rhs)
				if(pc->GetNPoints() < 2) continue;
				rhs.clear();
				njets_tot++;
				//loop over points
				//cout << "TREE " << i << endl;
				//cout << "tree points " << endl; pc->Print();
				//create new Jet
				Jet predJet(_trees[i]->model, BayesPoint({_pvx, _pvy, _pvz}), _gev, _radius);
				//put pt cut in for predjets of 5 GeV to match reco AK4 definition
				if(predJet.pt() < 5) continue; 
				//add Jet to jets	
				_predJets.push_back(predJet);	
			}
			cout << njets_tot << " pred jets total" << endl;
			//cout << _predJets.size() << " pred jets pt > 20 GeV" << endl;
			for(auto j : _predJets) cout << "pred jet px " << j.px() << " py " << j.py() << " pz " << j.pz() << " E " << j.E() << " m2 " << j.m2() << " mass " << j.mass_rhs() << " eta " << j.eta() << " phi " << j.phi() << endl;
		}

		void FillPredJetHists(){
			int njets = _predJets.size();
			int tmass_idx, id, widx;
			tmass_idx = -999;

			double wmass, topmass, dr, dr_pair, jetsize;
			pair<int, int> wmass_idxs = make_pair(-999, -999);
			Jet w, top, genpart, genW;
			//do gen matching
			vector<int> genMatchIdxs; //one per jet, follows same indexing as jets
			//GenMatchJet(_predJets,genMatchIdxs);
			vector<int> recoMatchIdxs;
			GenericMatchJet(_predJets,_recojets, recoMatchIdxs); //match BHC jets to reco jets
			int nsubs, bestMatchIdx;
			for(int p = 0; p < _procCats.size(); p++){
				//if(p != 0) cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				_procCats[p].hists1D[0][0]->Fill(njets);
				cout << "# pred jets - # reco jets " << njets - (int)_recojets.size() << endl;
				cout << "# pred jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][11]->Fill(njets - (int)_recojets.size());
				//FindResonances(_predJets,wmass_idxs,tmass_idx);
				//W mass
				//cout << "wmass idxs " << wmass_idxs.first << " " << wmass_idxs.second << endl;
				/*
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
				*/
				//top mass
				/*
				if(tmass_idx != -999){
					//cout << "adding w jet with pt " << w.pt() << endl;
					top = _predJets[tmass_idx];
					top.add(w);
					topmass = top.mass();
					_procCats[p].hists1D[0][30]->Fill(topmass);
				}
				*/
				//loop over pt bins
				njets = _predJets.size();
				double pt_thresh = 50;
				for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
					for(int j = 0; j < _predJets.size(); j++){
						if(p != 0) cout << "pred jet #" << j << " phi " << _predJets[j].phi() << " eta " << _predJets[j].eta() << " energy " << _predJets[j].E() <<  " mass " << _predJets[j].mass_rhs() << " nConstituents " << _predJets[j].GetNConstituents() << " nRhs " << _predJets[j].GetNRecHits() << " pt " << _predJets[j].pt() << endl;
						//define pt bins
						//pt == 1 -> [50,inf)
						if(pt == 1 && _predJets[j].pt() < pt_thresh) continue;
						//pt == 2 -> [0,50)
						if(pt == 2 && _predJets[j].pt() >= pt_thresh) continue;


						Matrix jetcov = CalcJetCovMat(_predJets[j]);
						//get 2D matrix for jet size
						Matrix jetcov2D(2,2);
						Get2DMat(jetcov,jetcov2D);	
						vector<double> eigvals;
						vector<Matrix> eigvecs;
						jetcov2D.eigenCalc(eigvals, eigvecs);
						//define jet size as length of major axis
						//also include rotundity
						dr = sqrt(eigvals[1]);//sqrt(sqrt(jetcov.at(0,0))*sqrt(jetcov.at(1,1)));
						double rot = Rotundity(jetcov2D);
						_procCats[p].hists1D[pt][147]->Fill(rot);
						cout << "calc jet size for BHC jet " << j << ": " << dr << endl;
						nsubs = _predJets[j].GetNConstituents();
						_procCats[p].hists1D[pt][1]->Fill(nsubs);
						_procCats[p].hists2D[pt][11]->Fill(nsubs,dr);
						//cout << "pred jet j " << j << " dr " << dr << " n constituents " << _predJets[j].GetNConstituents() << endl;
						Matrix cov = _predJets[j].GetCovariance();
					
						_procCats[p].hists1D[pt][7]->Fill(_predJets[j].e());
						_procCats[p].hists1D[pt][8]->Fill(_predJets[j].pt());
						_procCats[p].hists1D[pt][9]->Fill(_predJets[j].mass_rhs());
						_procCats[p].hists1D[pt][6]->Fill(dr);

						_procCats[p].hists1D[pt][132]->Fill(sqrt(jetcov.at(0,0)));
						_procCats[p].hists1D[pt][133]->Fill(sqrt(jetcov.at(1,1)));
						_procCats[p].hists1D[pt][134]->Fill(sqrt(jetcov.at(2,2)));
							
						_procCats[p].hists1D[pt][47]->Fill(sqrt(cov.at(0,0)));	
						_procCats[p].hists1D[pt][48]->Fill(sqrt(cov.at(1,1)));	
						_procCats[p].hists1D[pt][49]->Fill(sqrt(cov.at(2,2)));	
						_procCats[p].hists1D[pt][50]->Fill(cov.at(0,1));	
						_procCats[p].hists1D[pt][51]->Fill(cov.at(0,2));	
						_procCats[p].hists1D[pt][52]->Fill(cov.at(2,1));	
						
						_procCats[p].hists2D[pt][7]->Fill(_predJets[j].mass(), _predJets[j].pt());
						_procCats[p].hists2D[pt][8]->Fill(_predJets[j].mass(), dr);
						_procCats[p].hists2D[pt][10]->Fill(_predJets[j].pt(), dr);
						
						_procCats[p].hists2D[pt][34]->Fill((double)_predJets.size(), dr);
						_procCats[p].hists2D[pt][35]->Fill(_predJets[j].GetNRecHits(),_predJets[j].GetNConstituents());
						_procCats[p].hists2D[pt][40]->Fill(_predJets[j].GetNConstituents(), _predJets[j].mass_rhs());
						_procCats[p].hists2D[pt][45]->Fill(_predJets[j].e(), _predJets[j].mass_rhs());
						_procCats[p].hists2D[pt][41]->Fill(_predJets[j].GetNConstituents(), _predJets[j].e());	


						Matrix jet_mu, jet_cov;
						_predJets[j].GetClusterParams(jet_mu,jet_cov);	
						_procCats[p].hists1D[pt][128]->Fill(jet_mu.at(0,0));
						_procCats[p].hists1D[pt][129]->Fill(jet_mu.at(1,0));
						_procCats[p].hists1D[pt][130]->Fill(jet_mu.at(2,0));
						
						vector<JetPoint> rhs = _predJets[j].GetJetPoints();
						for(int r = 0; r < rhs.size(); r++){
							_procCats[p].hists1D[pt][145]->Fill(rhs[r].E());
							
						}

					//cout << "pred jet #" << j << " phi1 " << _predJets[j].phi() << " phi " << jet_mu.at(1,0) << " phi std " << _predJets[j].phi_std() << endl;
						//get subcluster information
						vector<Jet> consts = _predJets[j].GetConstituents();
						//for(auto subcl : _predJets[j].GetConstituents()){
						for(int c = 0; c < (int)consts.size(); c++){
							Jet subcl = consts[c];
							Matrix subcl_cov = subcl.GetCovariance();
							_procCats[p].hists1D[pt][12]->Fill(sqrt(subcl_cov.at(0,0)));
							_procCats[p].hists1D[pt][13]->Fill(sqrt(subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][14]->Fill(sqrt(subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][15]->Fill(subcl_cov.at(0,1));
							_procCats[p].hists1D[pt][16]->Fill(subcl_cov.at(0,2));
							_procCats[p].hists1D[pt][17]->Fill(subcl_cov.at(1,2));
							//do for normalized values too
	
							//E_k = norms[k]/_gev;
							_procCats[p].hists1D[pt][2]->Fill(subcl.E());

							//params = model->GetPriorParameters(k);
							//ceta = params["mean"].at(0,0);
							//cphi = params["mean"].at(1,0);
							//ctime = params["mean"].at(2,0);
							//norm += params["pi"].at(0,0);
						//cout << "pred jet subcluster phi " << subcl.phi() << " phi std " << subcl.phi_std() << " phi 02pi " << subcl.phi_02pi() << endl;	
							_procCats[p].hists1D[pt][3]->Fill(subcl.eta());
							_procCats[p].hists1D[pt][4]->Fill(subcl.phi());
							_procCats[p].hists1D[pt][5]->Fill(subcl.time());
							
							//do dr calcs bw all subclusters
							for(int cc = c+1; cc < consts.size(); cc++){
								double dr = dR(consts[c].eta(), consts[c].phi(), consts[cc].eta(), consts[cc].phi());
								_procCats[p].hists1D[pt][143]->Fill(dr);
							}
		
						}


						//if BHC and reco AK4 jets can be matched 1:1, compare # of subclusters
						if(njets - (int)_recojets.size() == 0){
							if(recoMatchIdxs[j] != -1){
								if(_recojets[recoMatchIdxs[j]].GetNConstituents() == 0) cout << "reco jet #" << recoMatchIdxs[j] << " matched to bhc jet #" << j << " has " << _recojets[recoMatchIdxs[j]].GetNConstituents() << " # subclusters" << endl;
								_procCats[p].hists2D[pt][43]->Fill(_recojets[recoMatchIdxs[j]].GetNConstituents(), _predJets[j].GetNConstituents());
							}
						}

						/*
						if(j != wmass_idxs.first && j != wmass_idxs.second){
							_procCats[p].hists1D[0][38]->Fill(cov.at(0,0));
							_procCats[p].hists1D[0][42]->Fill(cov.at(1,1));
							_procCats[p].hists1D[0][46]->Fill(cov.at(2,2));
						}
						
						if(genMatchIdxs[j] != -1) cout << " jet " << j << " is exclusively matched to gen particle " << genMatchIdxs[j] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs[j]), _base->genpart_phi->at(genMatchIdxs[j]), _recojets[j].eta(), _recojets[j].phi()) << endl;
						else{ cout << " jet " << j << " could not be gen matched" << endl; continue; }

						bestMatchIdx = genMatchIdxs[j];
						dr = dR(_base->genpart_eta->at(bestMatchIdx), _base->genpart_phi->at(bestMatchIdx), _predJets[j].eta(), _predJets[j].phi());
						widx = _base->genpart_momIdx->at(bestMatchIdx); //this branch only filled for products of Ws
						//cout << "widx " << widx << endl;
						if(widx != -1) _procCats[p].hists2D[0][19]->Fill(dr,_base->genpart_energy->at(widx));
						//fill gen match hists
						_procCats[p].hists1D[0][70]->Fill(id);
						//b's
						if(fabs(id) == 5){
							_procCats[p].hists1D[0][58]->Fill(dr);
							_procCats[p].hists1D[0][64]->Fill(_base->genpart_energy->at(bestMatchIdx)/_predJets[j].E());
						}
						//g's, q's (not b's)
						else if(fabs(id) == 1 || fabs(id) == 2 || fabs(id) == 3 || fabs(id) == 4){
							_procCats[p].hists1D[0][60]->Fill(dr);
							_procCats[p].hists1D[0][66]->Fill(_base->genpart_energy->at(bestMatchIdx)/_predJets[j].E());

						}
						//leptons
						else if(fabs(id) == 13 || fabs(id) == 14){
							_procCats[p].hists1D[0][62]->Fill(dr);
							_procCats[p].hists1D[0][68]->Fill(_base->genpart_energy->at(bestMatchIdx)/_predJets[j].E());

						}
						else{ }
						*/
					}
					//fill njets based on top decay type
					//0 -> fully had
					//1 -> fully lep
					//2 -> semi lep 
					//if(_topDecayType == 0)
					//	_procCats[p].hists1D[0][74]->Fill(njets);
					//else if(_topDecayType == 1)
					//	_procCats[p].hists1D[0][76]->Fill(njets);
					//else if(_topDecayType == 2)
					//	_procCats[p].hists1D[0][75]->Fill(njets);
					//else{ }
				}
			}
		}


		//fill hists for gen jets and gen particles
		//ONLY GEN PARTICLES CONSIDERED - tops, b's from tops, quarks from W's, direct daughters of b's
		void FillGenHists(){
			for(int p = 0; p < _procCats.size(); p++){
				//fill gen particle hists - needs GetGenParticles() method in JetSimProducer
				int nGenParts = 0; //only count hadron-izable gen particles (ie quarks)
				vector<int> qids = {1,2,3,4,5,6};
				int id;
				_procCats[p].hists1D[0][107]->Fill((double)nGenParts);
				for(int g = 0; g < _genparts.size(); g++){
					_procCats[p].hists1D[0][108]->Fill(_genparts[g].eta());
					_procCats[p].hists1D[0][109]->Fill(_genparts[g].phi());
					_procCats[p].hists1D[0][110]->Fill(_genparts[g].time());
					_procCats[p].hists1D[0][111]->Fill(_genparts[g].pt());
					_procCats[p].hists1D[0][112]->Fill(_genparts[g].mass());
					_procCats[p].hists1D[0][113]->Fill(_genparts[g].E());
					id = _base->genpart_id->at(_genparts[g].GetUserIdx());
					if(find(qids.begin(), qids.end(), id) != qids.end()) nGenParts++;
					
				}
				_procCats[p].hists1D[0][120]->Fill((double)_genjets.size());
				_procCats[p].hists1D[0][122]->Fill((double)_genjets.size() - (double)nGenParts);
				_procCats[p].hists1D[0][125]->Fill((double)_predJets.size() - (double)nGenParts);
				//gen match jets to particles
				vector<int> genMatchIdxs;
				//cout << "gen matching gen jets to particles - start" << endl;
				GenMatchJet(_genjets,genMatchIdxs);
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genjets.size(); j++){
					if(p == 0) cout << "gen jet #" << j << " phi " << _genjets[j].phi() << " eta " << _genjets[j].eta() << " energy " << _genjets[j].E() <<  " mass " << _genjets[j].mass() << " pt " << _genjets[j].pt() << endl;;
					_procCats[p].hists1D[0][114]->Fill(_genjets[j].eta());
					_procCats[p].hists1D[0][115]->Fill(_genjets[j].phi());
					_procCats[p].hists1D[0][116]->Fill(_genjets[j].time());
					_procCats[p].hists1D[0][117]->Fill(_genjets[j].pt());
					_procCats[p].hists1D[0][118]->Fill(_genjets[j].mass());
					_procCats[p].hists1D[0][119]->Fill(_genjets[j].E());
					_procCats[p].hists1D[0][121]->Fill(_base->Jet_genNConstituents->at(_genjets[j].GetUserIdx()));
					//dr bw gen jet and best exclusive gen particle match
					int genmatch = genMatchIdxs[j];
					if(genmatch != -1){
						double gendR = dR(_genjets[j].eta(), _genjets[j].phi(), _genparts[genmatch].eta(), _genparts[genmatch].phi());
						_procCats[p].hists1D[0][123]->Fill(gendR);
						_procCats[p].hists1D[0][124]->Fill(_genjets[j].E()/_genparts[genmatch].E());
					}	

				}	
				/*
				cout << "gen matching pred jets to particles - start" << endl;
				GenMatchJet(_predJets,genMatchIdxs);
				cout << "gen matching pred jets to particles - end" << endl;
				for(int j = 0; j < _predJets.size(); j++){
					int genmatch = genMatchIdxs[j];
					if(genmatch != -1){
						double gendR = dR(_predJets[j].eta(), _predJets[j].phi(), _genparts[genmatch].eta(), _genparts[genmatch].phi());
				cout << "gen dR for BHC jet " << gendR << endl;
						_procCats[p].hists1D[0][126]->Fill(gendR);
						_procCats[p].hists1D[0][127]->Fill(_predJets[j].E()/_genparts[genmatch].E());
					}	

				}
				*/	

			}

		}	
		void FillRecoJetHists(){
			int njets, tmass_idx, widx, genidx, id;
			double wmass, topmass, dr, dr_pair, jetsize;
			pair<int, int> wmass_idxs;
			Jet w, top, genW;
			//do gen matching
			vector<int> genMatchIdxs_jet; //one per jet, follows same indexing as jets
			vector<int> genMatchIdxs_p; //one per jet, follows same indexing as jets
			GenMatchJetToJet(_recojets,genMatchIdxs_jet);
			//skip gen particle matching for now
			//GenMatchJet(_recojets,genMatchIdxs_p);
			//cout << "final best matches" << endl;
			//for(int b = 0; b < genMatchIdxs_p.size(); b++){
			//	if(genMatchIdxs_p[b] != -1) cout << " jet " << b << " is exclusively matched to gen particle " << genMatchIdxs_p[b] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs_p[b]), _base->genpart_phi->at(genMatchIdxs_p[b]), _recojets[b].eta(), _recojets[b].phi()) << endl;
			//	 else cout << " jet " << b << " could not be gen matched" << endl;

			//}
			double pt_thresh = 50;	
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = _recojets.size();
				_procCats[p].hists1D[0][18]->Fill(njets);
				for(int j = 0; j < _recojets.size(); j++){
					for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
						//define pt bins
						//pt == 1 -> [50,inf)
						if(pt == 1 && _recojets[j].pt() < pt_thresh) continue;
						//pt == 2 -> [0,50)
						if(pt == 2 && _recojets[j].pt() >= pt_thresh) continue;
	
			cout << "pt " << pt << " reco jet pt " << _recojets[j].pt() << endl;
						Matrix jetcov = CalcJetCovMat(_recojets[j]);
						//get 2D matrix for jet size
						Matrix jetcov2D(2,2);
						Get2DMat(jetcov,jetcov2D);	
						vector<double> eigvals;
						vector<Matrix> eigvecs;
						jetcov2D.eigenCalc(eigvals, eigvecs);
						//define jet size as length of major axis
						//also include rotundity
						jetsize = sqrt(eigvals[1]);//sqrt(sqrt(jetcov.at(0,0))*sqrt(jetcov.at(1,1)));
						double rot = Rotundity(jetcov2D);
						_procCats[p].hists1D[pt][146]->Fill(rot);
						if(p == 0) cout << "reco jet #" << j << " phi " << _recojets[j].phi() << " eta " << _recojets[j].eta() << " energy " << _recojets[j].E() <<  " mass " << _recojets[j].mass() << " nConstituents " << _recojets[j].GetNConstituents() << " nRhs " << _recojets[j].GetNRecHits() << " pt " << _recojets[j].pt() << " jetsize " << jetsize << endl;
						_procCats[p].hists1D[pt][19]->Fill(jetsize);
						_procCats[p].hists1D[pt][135]->Fill(sqrt(jetcov.at(0,0)));
						_procCats[p].hists1D[pt][136]->Fill(sqrt(jetcov.at(1,1)));
						_procCats[p].hists1D[pt][137]->Fill(sqrt(jetcov.at(2,2)));
						_procCats[p].hists1D[pt][20]->Fill(_recojets[j].e());
						_procCats[p].hists1D[pt][21]->Fill(_recojets[j].pt());
						_procCats[p].hists1D[pt][22]->Fill(_recojets[j].mass_rhs());
cout << "mass hist for pt " << pt << " has " << _procCats[p].hists1D[pt][22]->GetEntries() << " entries - filled with " << _recojets[j].mass_rhs() << endl;
						_procCats[p].hists1D[pt][138]->Fill(_recojets[j].time());
						_procCats[p].hists1D[pt][139]->Fill(_recojets[j].eta());
						_procCats[p].hists1D[pt][140]->Fill(_recojets[j].phi());
						_procCats[p].hists2D[pt][4]->Fill(_recojets[j].mass(), _recojets[j].pt());
						_procCats[p].hists2D[pt][5]->Fill(_recojets[j].mass(), jetsize);
						_procCats[p].hists2D[pt][33]->Fill((double)_recojets.size(), jetsize);
						
						//fill subcluster hists
						_procCats[p].hists1D[pt][77]->Fill(_recojets[j].GetNConstituents());
						_procCats[p].hists2D[pt][20]->Fill(_recojets[j].GetNRecHits(),_recojets[j].GetNConstituents());
						_procCats[p].hists2D[pt][37]->Fill(_recojets[j].GetNConstituents(), _recojets[j].mass_rhs());
						_procCats[p].hists2D[pt][38]->Fill(_recojets[j].GetNConstituents(), _recojets[j].e());
						_procCats[p].hists2D[pt][44]->Fill(_recojets[j].e(), _recojets[j].mass_rhs());
						_procCats[p].hists2D[pt][46]->Fill(_recojets[j].GetNConstituents(), jetsize);
						if(_recojets[j].GetNConstituents() == 0) cout << _recojets[j].GetNConstituents() << " n subcl " << _recojets[j].GetNRecHits() << " n rhs" << endl;
						vector<Jet> consts = _recojets[j].GetConstituents();
						for(int c = 0; c < (int)consts.size(); c++){
							Jet subcl = consts[c];
							_procCats[p].hists1D[pt][78]->Fill(subcl.E());
							_procCats[p].hists1D[pt][79]->Fill(subcl.time());
							_procCats[p].hists1D[pt][80]->Fill(subcl.eta());
							_procCats[p].hists1D[pt][81]->Fill(subcl.phi());
							
							Matrix subcl_cov = subcl.GetCovariance();
							_procCats[p].hists1D[pt][82]->Fill(sqrt(subcl_cov.at(0,0)));
							_procCats[p].hists1D[pt][83]->Fill(sqrt(subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][84]->Fill(sqrt(subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][85]->Fill(subcl_cov.at(0,1));
							_procCats[p].hists1D[pt][86]->Fill(subcl_cov.at(0,2));
							_procCats[p].hists1D[pt][87]->Fill(subcl_cov.at(1,2));
							_procCats[p].hists1D[pt][104]->Fill(subcl_cov.at(0,1)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][105]->Fill(subcl_cov.at(0,2)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][106]->Fill(subcl_cov.at(1,2)/sqrt(subcl_cov.at(1,1)*subcl_cov.at(2,2)));
							_procCats[p].hists2D[pt][31]->Fill(subcl_cov.at(0,1),subcl_cov.at(0,2));
							_procCats[p].hists2D[pt][32]->Fill(subcl_cov.at(0,1)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(1,1)),subcl_cov.at(0,2)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(2,2)));
						
							//do dr calcs bw all subclusters
							for(int cc = c+1; cc < consts.size(); cc++){
								double dr = dR(consts[c].eta(), consts[c].phi(), consts[cc].eta(), consts[cc].phi());
								_procCats[p].hists1D[pt][144]->Fill(dr);
							}
				

						}
						vector<JetPoint> rhs = _recojets[j].GetJetPoints();
						for(int r = 0; r < rhs.size(); r++){
							_procCats[p].hists1D[pt][89]->Fill(rhs[r].E());
							_procCats[p].hists1D[pt][96]->Fill(rhs[r].t());
							
						}
						vector<pair<double,double>> geoEavg_diffT;
						CalcRhTimeDiff(rhs,geoEavg_diffT);
						for(int r = 0; r < geoEavg_diffT.size(); r++){
						      _procCats[p].hists2D[pt][30]->Fill(geoEavg_diffT[r].first, geoEavg_diffT[r].second);
						}
						_procCats[p].hists1D[pt][88]->Fill(_recojets[j].GetNRecHits());
						//if no gen match, skip
						if(genMatchIdxs_jet[j] == -1) continue;
						genidx = genMatchIdxs_jet[j];	
						
						//cout << "gen jet match idx " << genMatchIdxs_jet[j] << " gen jet user idx  " << _genjets[genMatchIdxs_jet[j]].GetUserIdx() << " gen n const " << _base->Jet_genNConstituents->at(_genjets[genMatchIdxs_jet[j]].GetUserIdx()) << endl;
						int genjetidx = _genjets[genMatchIdxs_jet[j]].GetUserIdx();
						_procCats[p].hists2D[pt][21]->Fill(_base->Jet_genNConstituents->at(genjetidx),_recojets[j].GetNConstituents());
						double genpt, pz, jetpt, jetpz, ratio_p;
						jetpt = _base->Jet_genPt->at(genjetidx);
						jetpz = _base->Jet_genPz->at(genjetidx);

						_procCats[p].hists1D[pt][90]->Fill(sqrt(jetpt*jetpt + jetpz*jetpz));	
						_procCats[p].hists1D[pt][91]->Fill(jetpt);	
						_procCats[p].hists2D[pt][28]->Fill(jetpt, _recojets[j].GetNConstituents());
						_procCats[p].hists2D[pt][36]->Fill(_recojets[j].GetNConstituents(), _base->Jet_genNConstituents->at(genjetidx));
						
						int genpartidx = -1;
						int ngenparts_ptge5 = 0;
						for(int g = 0; g < _base->Jet_genNConstituents->at(genjetidx); g++){
							genpartidx = _base->Jet_genConstituentIdxs->at(genjetidx).at(g);
							genpt = _base->genpart_pt->at(genpartidx);
							pz = _base->genpart_pz->at(genpartidx);
							ratio_p = (sqrt(genpt*genpt + pz*pz))/(sqrt(jetpt*jetpt + jetpz*jetpz));	
							_procCats[p].hists1D[pt][92]->Fill(genpt);	
							_procCats[p].hists1D[pt][93]->Fill(sqrt(genpt*genpt + pz*pz));

							_procCats[p].hists1D[pt][94]->Fill( ratio_p );	
							_procCats[p].hists1D[pt][95]->Fill( genpt/jetpt );

							_procCats[p].hists2D[pt][22]->Fill(sqrt(genpt*genpt + pz*pz), sqrt(jetpt*jetpt + jetpz*jetpz));	
							_procCats[p].hists2D[pt][23]->Fill(genpt, jetpt);	
							_procCats[p].hists2D[pt][24]->Fill(sqrt(genpt*genpt + pz*pz), ratio_p);	
							_procCats[p].hists2D[pt][25]->Fill(sqrt(jetpt*jetpt + jetpz*jetpz), ratio_p);	
							_procCats[p].hists2D[pt][26]->Fill(genpt, genpt/jetpt);	
							_procCats[p].hists2D[pt][27]->Fill(jetpt, genpt/jetpt);
	
							if(pt >= 5) ngenparts_ptge5++;	
						}
						_procCats[p].hists2D[pt][29]->Fill(ngenparts_ptge5,_recojets[j].GetNConstituents());

						//skip gen matching for now
						continue;
						//fill gen match histograms
						if(p == 0) cout << "found gen match to reco jet " << j << " with id " << _base->genpart_id->at(genidx) << " and dr " << dR(_base->genpart_eta->at(genidx), _base->genpart_phi->at(genidx), _recojets[j].eta(), _recojets[j].phi()) << " and gen/reco E ratio " << _base->genpart_energy->at(genidx)/_recojets[j].E() << " user idx " << genidx << endl;
						//get matched gen particle here
						dr = dR(_base->genpart_eta->at(genidx), _base->genpart_phi->at(genidx), _recojets[j].eta(), _recojets[j].phi());
						id = _base->genpart_id->at(genidx);
						widx = _base->genpart_momIdx->at(genidx); //this branch only filled for products of Ws
						if(p == 0)cout << "widx " << widx << endl;
						if(widx != -1) _procCats[p].hists2D[pt][18]->Fill(dr,_base->genpart_energy->at(widx));
						//fill gen match hists
						if(fabs(id) != 13 && fabs(id) != 11 && fabs(id) != 15) _procCats[p].hists1D[pt][69]->Fill(fabs(id));
						else _procCats[p].hists1D[pt][69]->Fill(0);
						//b's
						if(fabs(id) == 5){
							_procCats[p].hists1D[pt][57]->Fill(dr);
							_procCats[p].hists1D[pt][63]->Fill(_base->genpart_energy->at(genidx)/_recojets[j].E());
						}
						//g's, q's (not b's)
						else if(fabs(id) == 1 || fabs(id) == 2 || fabs(id) == 3 || fabs(id) == 4){
							_procCats[p].hists1D[pt][59]->Fill(dr);
							_procCats[p].hists1D[pt][65]->Fill(_base->genpart_energy->at(genidx)/_recojets[j].E());

						}
						//leptons
						else if(fabs(id) == 13 || fabs(id) == 14){
							_procCats[p].hists1D[pt][61]->Fill(dr);
							_procCats[p].hists1D[pt][67]->Fill(_base->genpart_energy->at(genidx)/_recojets[j].E());

						}
						else{ }

					}
				}
				//fill njets based on top decay type
				//0 -> fully had
				//1 -> fully lep
				//2 -> semi lep 
				//skip gen matching for now
				continue;	

				if(_topDecayType == 0)
					_procCats[p].hists1D[0][71]->Fill(njets);
				else if(_topDecayType == 1)
					_procCats[p].hists1D[0][73]->Fill(njets);
				else if(_topDecayType == 2)
					_procCats[p].hists1D[0][72]->Fill(njets);
				else{ }

				//cout << "hist name " << _procCats[p].hists1D[0][18]->GetName() << " nentries " << _procCats[p].hists1D[0][18]->GetEntries() << endl;
				//if(p == 0)cout << "# reco jets - # gen jets " << njets - (int)_genjets.size() << endl;
				_procCats[p].hists1D[0][24]->Fill(njets - (int)_genjets.size());


				FindResonances(_recojets,wmass_idxs,tmass_idx);
				if(p == 0)cout << "wmass idxs " << wmass_idxs.first << " " << wmass_idxs.second << " top idx " << tmass_idx << endl;
				if(wmass_idxs.first != -999){
					w = _recojets[wmass_idxs.first];
					w.add(_recojets[wmass_idxs.second]);
					wmass = w.mass();
					dr_pair = dR(_recojets[wmass_idxs.first].eta(), _recojets[wmass_idxs.first].phi(), _recojets[wmass_idxs.second].eta(), _recojets[wmass_idxs.second].phi());
					//_procCats[p].hists1D[0][27]->Fill(wmass);
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
			vector<int> genMatchIdxs_reco; //one per jet, follows same indexing as jets
			vector<int> genMatchIdxs_bhc; //one per jet, follows same indexing as jets
			GenMatchJetToJet(_recojets,genMatchIdxs_reco);
			GenMatchJetToJet(_predJets,genMatchIdxs_bhc);
			int bestIdx;
			//cout << "n reco jets " << _recojets.size() << " n idxs " << genMatchIdxs_reco.size() << " n pred jets " << _predJets.size() << " n pred idxs " << genMatchIdxs_bhc.size() <<  endl;
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "proc " << p << " # reco jets " << _recojets.size() << " # gen jets " << _genjets.size() << endl;
				//reco jets
				for(int j = 0; j < _recojets.size(); j++){
					if(genMatchIdxs_reco[j] != -1) cout << " reco jet " << j << " is exclusively matched to gen jet " << genMatchIdxs_reco[j] << " with dr " << dR(_base->Jet_genEta->at(genMatchIdxs_reco[j]), _base->Jet_genPhi->at(genMatchIdxs_reco[j]), _recojets[j].eta(), _recojets[j].phi()) << endl;
					 else{ cout << " jet " << j << " could not be gen matched" << endl; continue; }
					 bestIdx = genMatchIdxs_reco[j];
					
					_procCats[p].hists2D[0][1]->Fill(_genjets[bestIdx].E(), _recojets[j].pt() - _genjets[bestIdx].pt());
					_procCats[p].hists1D[0][25]->Fill(_recojets[j].pt()/_genjets[bestIdx].pt());
					_procCats[p].hists1D[0][26]->Fill(_recojets[j].e()/_genjets[bestIdx].e());
					_procCats[p].hists2D[0][2]->Fill(_genjets[bestIdx].pt(), _recojets[j].pt());
					_procCats[p].hists2D[0][3]->Fill(_genjets[bestIdx].mass(), _recojets[j].pt()/_genjets[bestIdx].pt());
				}
				//predicted jets
				for(int j = 0; j < _predJets.size(); j++){
					if(genMatchIdxs_bhc[j] != -1) cout << " bhc jet " << j << " is exclusively matched to gen jet " << genMatchIdxs_bhc[j] << " with dr " << dR(_base->Jet_genEta->at(genMatchIdxs_bhc[j]), _base->Jet_genPhi->at(genMatchIdxs_bhc[j]), _predJets[j].eta(), _predJets[j].phi()) << endl;
					 else{ cout << " jet " << j << " could not be gen matched" << endl; continue; }
					 bestIdx = genMatchIdxs_bhc[j];
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
				//E_k = norms[k]/_gev;
				//Etot += E_k;
				//_procCats[p].hists1D[0][2]->Fill(E_k);

				//params = model->GetPriorParameters(k);
				//ceta = params["mean"].at(0,0);
				//cphi = params["mean"].at(1,0);
				//ctime = params["mean"].at(2,0);
				//norm += params["pi"].at(0,0);
				//
				//_procCats[p].hists1D[0][3]->Fill(ceta);
				//_procCats[p].hists1D[0][4]->Fill(cphi);
				//_procCats[p].hists1D[0][5]->Fill(ctime);
		
				//calculate slopes from eigenvectors
				//params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//total cluster energy
			}
		}
	
		//use rhs - space only 
		Matrix CalcJetCovMat(const Jet& jet){
			Matrix recocov = Matrix(3,3);
			vector<JetPoint> rhs = jet.GetJetPoints();
			int nrhs = rhs.size();
			
			double diffeta, diffphi, difftime;
			double wtot = 0;
			for(int r = 0; r < nrhs; r++){
				diffeta = rhs[r].eta() - jet.eta();
				diffphi = rhs[r].phi() - jet.phi();
				diffphi = acos(cos(diffphi));
				difftime = rhs[r].t() - jet.t();					

				wtot += rhs[r].E();

				recocov.SetEntry( recocov.at(0,0) + rhs[r].E()*diffeta*diffeta, 0, 0 );
				recocov.SetEntry( recocov.at(1,0) + rhs[r].E()*diffphi*diffeta, 1, 0 );
				recocov.SetEntry( recocov.at(0,1) + rhs[r].E()*diffeta*diffphi, 0, 1 );
				recocov.SetEntry( recocov.at(2,0) + rhs[r].E()*difftime*diffeta, 2, 0 );
				recocov.SetEntry( recocov.at(0,2) + rhs[r].E()*diffeta*difftime, 0, 2 );
				recocov.SetEntry( recocov.at(1,1) + rhs[r].E()*diffphi*diffphi, 1, 1 );
				recocov.SetEntry( recocov.at(1,2) + rhs[r].E()*diffphi*difftime, 1, 2 );
				recocov.SetEntry( recocov.at(2,1) + difftime*diffphi, 2, 1 );
				recocov.SetEntry( recocov.at(2,2) + difftime*difftime, 2, 2 );
			}
			recocov.mult(recocov,1./(double)wtot);
			//cout << "recocov for jetsize from " << nrhs << " pts" << endl; recocov.Print();
			return recocov;
		}
	
	

		void CalcMMAvgPhiTime(BasePDFMixture* model, double& phi, double& t){
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double pi, ws;
			map<string, Matrix> params;
			for(int k = 0; k < kmax; k++){
				params = model->GetLHPosteriorParameters(k);
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
			for(int pt = 0; pt < _procCats[0].hists1D.size(); pt++){
				int nhists = _procCats[0].hists1D[pt].size();
				string pt_title = "";
				//if(/
				for(int i = 0; i < nhists; i++){
					if(_procCats[0].hists1D[pt][i] == nullptr) continue;
					name = _procCats[0].hists1D[pt][i]->GetName();
					if(_procCats[0].hists1D[pt][i]->GetEntries() == 0 && ((histname.find("sigma") == string::npos && histname.find("mean") == string::npos) && histname.find("profile") == string::npos)){ continue; }
					
					_procCats[0].hists1D[pt][i]->Write();
					//make process breakdown directory - not making these profiles rn
					TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
					//cout << "  making dir " << dir2->GetName() << " name " << name << endl;
					dir2->cd();
					for(int p = 1; p < _procCats.size(); p++){
					//loop over processes
						//if(name.find("profile") != string::npos) continue;
						if(_procCats[p].hists1D[pt][i] == nullptr) continue;
						histname = _procCats[p].hists1D[pt][i]->GetName();
						//cout << "    proc " << _procCats[p].plotName << " hist " << _procCats[p].hists1D[0][i]->GetName() << " " << _procCats[p].hists1D[0][i]->GetTitle() << " entries " << _procCats[p].hists1D[0][i]->GetEntries() << endl;			
						if(_procCats[p].hists1D[pt][i]->GetEntries() == 0 && ((histname.find("sigma") == string::npos && histname.find("mean") == string::npos) && histname.find("profile") == string::npos)){ continue; }
						if(histname.find("profile") != string::npos){
							histname = histname.substr(0,histname.rfind("_"));
							_procCats[p].hists1D[pt][i]->SetName(histname.c_str());
							_procCats[p].hists1D[pt][i]->SetTitle(_procCats[p].plotName.c_str());
						}
						//cout << "  n hists " << _procCats[0].hists1D[0].size() << endl;
						//cout << "i " << i << " p " << p << " writing " << _procCats[p].hists1D[0][i]->GetName() << " " << _procCats[p].hists1D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;;
						_procCats[p].hists1D[pt][i]->Write();

					}
					ofile->cd(); 
				}
				//cout << "2D hists" << endl;
				//write 2D hists
				nhists = _procCats[0].hists2D[pt].size();
				for(int i = 0; i < nhists; i++){
					name = _procCats[0].hists2D[pt][i]->GetName();
					histname = _procCats[0].hists2D[pt][i]->GetName();
					//cout << "writing for " << name << " i " << i << " entries " << _procCats[0].hists2D[0][i]->GetEntries() << endl;
					//write total method histogram outside process directory
					if(_procCats[0].hists2D[pt][i] == nullptr) continue;
					//cout << "passed null" << endl;
					if(_procCats[0].hists2D[pt][i]->GetEntries() == 0 && histname.find("sigma") == string::npos){ continue; }
					//cout << "writing hist " << _procCats[0].hists2D[0][i]->GetName() << endl;
					_procCats[0].hists2D[pt][i]->Write();
					//write method as directory within directory
					TDirectory *dir2 = ofile->mkdir((name+"_procStack").c_str());
					//cout << "  making dir " << dir2->GetName() << endl;
					dir2->cd();
					for(int p = 1; p < _procCats.size(); p++){
						//loop over processes
						if(_procCats[p].hists2D[pt][i] == nullptr) continue;
						if(_procCats[p].hists2D[pt][i]->GetEntries() == 0 && dirname.find("meanRecoGenDeltaT") == string::npos){ continue;}// cout << "Histogram for proc " << _procCats[p].hists2D[pt][i]->GetName() << " not filled." << endl; continue; }
						//check if data can be run
						//cout << "writing " << _procCats[p].hists2D[0][i]->GetName() <<  " " <<  _procCats[p].hists2D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;
						histname = _procCats[p].hists2D[pt][i]->GetName();
						_procCats[p].hists2D[pt][i]->SetTitle(_procCats[p].plotName.c_str());
						_procCats[p].hists2D[pt][i]->Write();
					} 
					ofile->cd(); 
				}
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
					//cout << "hist " << histname << endl;
					//make sure profiles get written
					nprofs = _procCats[p].hists2D[0][i]->GetNbinsX();
					for(int k = 1; k < nprofs+1; k++){
						histname = _procCats[p].hists2D[0][i]->GetName();
						profname = "profile_"+histname;
						//profname.insert(profname.size(),"_bin"+std::to_string(k));
						if(!_procCats[p].plotName.empty()) profname.insert(profname.find("_"+_procCats[p].plotName),"_bin"+std::to_string(k));
						else profname = profname += "_bin"+std::to_string(k);
						nbins = _procCats[p].hists2D[0][i]->GetNbinsY();
						TH1D* prof = new TH1D(profname.c_str(), profname.c_str(), nbins, _procCats[0].hists2D[0][i]->GetYaxis()->GetBinLowEdge(1), _procCats[0].hists2D[0][i]->GetYaxis()->GetBinUpEdge(nbins));
						prof->SetTitle(_procCats[p].plotName.c_str());
			//cout << "p " << p << " profname " << profname << " plotname " << _procCats[p].plotName << " title " << prof->GetTitle() << endl;
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
		//comp time distribution
		TH1D* comptime_subcl = new TH1D("comptime_subcl","comptime_subcl",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime_subcl = new TGraph();
		
		vector<double> xbins = {0, 1, 2, 3, 4, 5, 10, 15, 20, 30,}; //for time resolution
		//predicted jet plots
		vector<double> xbins_recoGenPt = {0, 20, 30, 50, 100};
		//0 - n bhc jets
		TH1D* nClusters = new TH1D("BHCJet_nJets","BHCJet_nJets",15,0,15);
		//1 - n subclusters per bhc jets
		TH1D* nSubclusters = new TH1D("BHCJet_nSubclustersJet","BHCJet_nSubclustersJet",30,0,30);
		//2 - bhc subcluster energy
		TH1D* predJet_subClusterEnergy = new TH1D("BHCJet_subClusterEnergy","BHCJet_subClusterEnergy",50,0,250);
		//3 - bhc subcluster eta center
		TH1D* predJet_subClusterEtaCenter = new TH1D("BHCJet_subClusterEtaCenter","BHCJet_subClusterEtaCenter",25,-3.2,3.2);
		//4 - bhc subcluster phi center
		TH1D* predJet_subClusterPhiCenter = new TH1D("BHCJet_subClusterPhiCenter","BHCJet_subClusterPhiCenter",25,-0.1,6.3);
		//5 - bhc jet subcluster time center
		TH1D* predJet_subClusterTimeCenter = new TH1D("BHCJet_subClusterTimeCenter","BHCJet_subClusterTimeCenter",25,-1,1);
		//6 - bhc jet size
		TH1D* predJet_jetSize = new TH1D("BHCJet_jetSize","BHCJet_jetSize",50,0,1);
		//7 - bhc jet energy
		TH1D* predJet_energy = new TH1D("BHCJet_energy","BHCJet_energy",25,0,500);
		//8 - bhc jet pt
		TH1D* predJet_pt = new TH1D("BHCJet_pt","BHCJet_pt",25,0,500);
		//9 - bhc jet mass
		TH1D* predJet_mass = new TH1D("BHCJet_mass","BHCJet_mass",50,0,180);
		//10 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_predGen = new TH1D("jetGenE_sigmaDeltaPt_predGen","jetGenE_sigmaDeltaPt_predGen",5,0,100);
		//11 - # pred jets - # reco jets
		TH1D* predGen_nJets = new TH1D("BHCRecoAK4_diffNJets","BHCRecoAK4_diffNJets",20,-10,10);
		//for subclusters
		//12 - eta sigma
		TH1D* predJet_subClusterEtaSig = new TH1D("BHCJet_subClusterEtaSig","BHCJet_subClusterEtaSig",50,0.,0.1);
		//13 - phi sigma
		TH1D* predJet_subClusterPhiSig = new TH1D("BHCJet_subClusterPhiSig","BHCJet_subClusterPhiSig",50,0.,0.1);
		//14 - time sigma
		TH1D* predJet_subClusterTimeSig = new TH1D("BHCJet_subClusterTimeSig","BHCJet_subClusterTimeSig",50,0.,5.);
		//15 - etaphi cov
		TH1D* predJet_subClusteretaPhiCov = new TH1D("BHCJet_subClusteretaPhiCov","BHCJet_subClusteretaPhiCov",50,-0.0005,0.0005);
		//16 - time-eta covariance
		TH1D* predJet_subClustertimeEtaCov = new TH1D("BHCJet_subClustertimeEtaCov","BHCJet_subClustertimeEtaCov",50,-0.05,0.05);
		//17 - time-phi covariance
		TH1D* predJet_subClustertimePhiCov = new TH1D("BHCJet_subClustertimePhiCov","BHCJet_subClustertimePhiCov",50,-0.05,0.05);
		//18 - n reco AK4 jets
		TH1D* nRecoJets = new TH1D("recoAK4_nJets","recoAK4_nJets",10,0,10);
		//19 - reco AK4 jet size
		TH1D* recoJet_jetSize = new TH1D("recoAK4Jet_jetSize","recoAK4Jet_jetSize",50,0,1);
		//20 - reco AK4 jet energy
		TH1D* recoJet_energy = new TH1D("recoAK4Jet_energy","recoAK4Jet_energy",25,0,500);
		//21 - reco AK4 jet pt
		TH1D* recoJet_pt = new TH1D("recoAK4Jet_pt","recoAK4Jet_pt",25,0,500);
		//22 - reco AK4 jet mass
		TH1D* recoJet_mass = new TH1D("recoAK4Jet_mass","recoAK4Jet_mass",50,0,180);
		//23 - resolution of difference of pt between reco and gen jets as a function of gen jet energy
		TH1D* jetGenE_sigmaDeltaPt_recoGen = new TH1D("jetGenE_sigmaDeltaPt_recoAK4Gen","jetGenE_sigmaDeltaPtOvJetGenE_recoAK4Gen",4,&xbins_recoGenPt[0]);
		//24 - # reco jets - # gen jets
		TH1D* recoGen_nJets = new TH1D("recoAK4Gen_diffNJets","recoAK4Gen_diffNJets",20,-10,10);
		//25 - reco jet pt/gen jet pt
		TH1D* recoGen_jetPtRatio = new TH1D("recoAK4Gen_jetPtRatio","recoAK4Gen_jetPtRatio",20,0,1.5);
		//26 - reco jet e - gen jet e		
		TH1D* recoGen_jetERatio = new TH1D("recoAK4Gen_jetERatio","recoAK4Gen_jetERatio",20,-10,10);
		//27 - reco jet W invariant mass
		TH1D* recoJet_Wmass = new TH1D("recoAK4Jet_Wmass","recoAK4Jet_Wmass;Invariant mass for best W candidate",50,0,500);
		//28 - reco jet top invariant mass
		TH1D* recoJet_topmass = new TH1D("recoAK4Jet_topmass","recoAK4Jet_topmass;Invariant mass for best top candidate",50,0,500);
		//29 - pred jet W invariant mass
		TH1D* predJet_Wmass = new TH1D("BHCJet_Wmass","BHCJet_Wmass;Invariant mass for best W candidate",50,0,500);
		//30 - pred jet top invariant mass
		TH1D* predJet_topmass = new TH1D("BHCJet_topmass","BHCJet_topmass;Invariant mass for best top candidate",50,0,500);
		//31 - pred jet W invariant mass, pT_jj < 100
		TH1D* predJet_Wmass_pTjjl100 = new TH1D("BHCJet_Wmass_pTjjl100","BHCJet_Wmass_pTjjl100;Invariant mass_pTjjl100 for best W candidate",50,0,500);
		//32 - pred jet W invariant mass, pT_jj >= 100
		TH1D* predJet_Wmass_pTjjge100 = new TH1D("BHCJet_Wmass_pTjjge100","BHCJet_Wmass_pTjjge100;Invariant mass_pTjjge100 for best W candidate",50,0,500);
		//33 - reco jet W invariant mass, pT_jj < 100
		TH1D* recoJet_Wmass_pTjjl100 = new TH1D("recoAK4Jet_Wmass_pTjjl100","recoAK4Jet_Wmass_pTjjl100;Invariant mass_pTjjl100 for best W candidate",50,0,500);
		//34 - reco jet W invariant mass, pT_jj >= 100
		TH1D* recoJet_Wmass_pTjjge100 = new TH1D("recoAK4Jet_Wmass_pTjjge100","recoAK4Jet_Wmass_pTjjge100;Invariant mass_pTjjge100 for best W candidate",50,0,500);
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
		TH1D* predJet_EtaVar = new TH1D("BHCJet_EtaSig","BHCJet_EtaSig",25,0.,1.);
		//48 - phi sigma for jet
		TH1D* predJet_PhiVar = new TH1D("BHCJet_PhiSig","BHCJet_PhiSig",25,0., 1.);
		//49 - time sigma for jet
		TH1D* predJet_TimeVar = new TH1D("BHCJet_TimeSig","BHCJet_TimeSig",25,0.,1.);
		//50 - eta-phi covariance for jet
		TH1D* predJet_etaPhiCov = new TH1D("BHCJet_etaPhiCov","BHCJet_etaPhiCov",25,-0.0005,0.0005);
		//51 - time-eta covariance for jet
		TH1D* predJet_timeEtaCov = new TH1D("BHCJet_timeEtaCov","BHCJet_timeEtaCov",25,-0.2,0.2);
		//52 - time-phi covariance for jet
		TH1D* predJet_timePhiCov = new TH1D("BHCJet_timePhiCov","BHCJet_timePhiCov",25,-0.2,0.2);
		//53 - dr between predicted jets for W candidate
		TH1D* predJet_WjjDr = new TH1D("BHCJet_WjjDr","BHCJet dR between W candidate jets",25,0,3.5);
		//54 - dr between reco jets for W candidate
		TH1D* recoJet_WjjDr = new TH1D("recoAK4Jet_WjjDr","recoJet dR between W candidate jets",25,0,3.5);
		//55 - pt of W candidate from pred jets
		TH1D* predJet_Wpt = new TH1D("BHCJet_Wpt","BHCJet W candidate pt",25,0,250.);
		//56 - pt of W candidate from reco jets
		TH1D* recoJet_Wpt = new TH1D("recoAK4Jet_Wpt","recoJet W candidate pt",25,0,250.);
		//57 - reco dr match to gen b's
		TH1D* recoJet_dR_b = new TH1D("recoAK4Jet_dR_b","recoAK4Jet_dR_b",25,0,1.);
		//58 - bhc dr match to gen b's
		TH1D* BHCJet_dR_b = new TH1D("BHCJet_dR_b","BHCJet_dR_b",25,0,1.);
		//59 - reco dr match to gen q's/g's
		TH1D* recoJet_dR_qg = new TH1D("recoAK4Jet_dR_qg","recoAK4Jet_dR_qg",25,0,1.);
		//60 - bhc dr match to gen q's/g's
		TH1D* BHCJet_dR_qg = new TH1D("BHCJet_dR_qg","BHCJet_dR_qg",25,0,1.);
		//61 - reco dr match to gen l's  
		TH1D* recoJet_dR_lep = new TH1D("recoAK4Jet_dR_lep","recoAK4Jet_dR_lep",25,0,1.);
		//62 - bhc dr match to gen l's  
		TH1D* BHCJet_dR_lep = new TH1D("BHCJet_dR_lep","BHCJet_dR_lep",25,0,1.);
		//63 - reco dr match to gen b's
		TH1D* recoJet_genOvRecoE_b = new TH1D("recoAK4Jet_genOvRecoE_b","recoAK4Jet_genOvRecoE_b",25,0,5.);
		//64 - bhc dr match to gen b's
		TH1D* BHCJet_genOvRecoE_b = new TH1D("BHCJet_genOvRecoE_b","BHCJet_genOvRecoE_b",25,0,5.);
		//65 - reco dr match to gen q's/g's
		TH1D* recoJet_genOvRecoE_qg = new TH1D("recoAK4Jet_genOvRecoE_qg","recoAK4Jet_genOvRecoE_qg",25,0,5.);
		//66 - bhc dr match to gen q's/g's
		TH1D* BHCJet_genOvRecoE_qg = new TH1D("BHCJet_genOvRecoE_qg","BHCJet_genOvRecoE_qg",25,0,5.);
		//67 - reco dr match to gen l's  
		TH1D* recoJet_genOvRecoE_lep = new TH1D("recoAK4Jet_genOvRecoE_lep","recoAK4Jet_genOvRecoE_lep",25,0,5.);
		//68 - bhc dr match to gen l's  
		TH1D* BHCJet_genOvRecoE_lep = new TH1D("BHCJet_genOvRecoE_lep","BHCJet_genOvRecoE_lep",25,0,5.);
		//69 - types of quarks in reco jets - d, u, s, c, b (sign agnostic)
		TH1D* qType_recoJet = new TH1D("qType_recoAK4Jet","qType_recoAK4Jet",6,0,6);
		//70 - types of quarks in BHC jets - d, u, s, c, b (sign agnostic)
		TH1D* qType_BHCJet = new TH1D("qType_BHCJet","qType_BHCJet",6,0,6);
		//71 - n reco jets, tt fully had
		TH1D* reco_nJets_fullHad = new TH1D("recoAK4_nJets_fullHad","recoAK4_nJets_fullHad",10,0,10);
		//72 - n reco jets, tt semi lep
		TH1D* reco_nJets_semiLep = new TH1D("recoAK4_nJets_semiLep","recoAK4_nJets_semiLep",10,0,10);
		//73 - n reco jets, tt fully lep
		TH1D* reco_nJets_fullLep = new TH1D("recoAK4_nJets_fullLep","recoAK4_nJets_fullLep",10,0,10);
		//74 - n BHC jets, tt fully had
		TH1D* BHC_nJets_fullHad = new TH1D("BHC_nJets_fullHad","BHC_nJets_fullHad",10,0,10);
		//75 - n BHC jets, tt semi lep
		TH1D* BHC_nJets_semiLep = new TH1D("BHC_nJets_semiLep","BHC_nJets_semiLep",10,0,10);
		//76 - n BHC jets, tt fully lep
		TH1D* BHC_nJets_fullLep = new TH1D("BHC_nJets_fullLep","BHC_nJets_fullLep",10,0,10);
		//77 - n GMM clusters in reco jets
		TH1D* reco_nSubclusters = new TH1D("recoAK4Jet_nSubclustersJet","recoAK4Jet_nSubclustersJet",30,0,30);
		//78 - energy per GMM cluster from reco jets
		TH1D* recoJet_subClusterEnergy = new TH1D("recoAK4Jet_subClusterEnergy","recoAK4Jet_subClusterEnergy",50,0,250);
		//79 - time center of GMM cluster from reco jets
		TH1D* recoJet_subClusterTimeCenter = new TH1D("recoAK4Jet_subClusterTimeCenter","recoAK4Jet_subClusterTimeCenter",25,-1,1);
		//80 - eta center of GMM cluster from reco jets
		TH1D* recoJet_subClusterEtaCenter = new TH1D("recoAK4Jet_subClusterEtaCenter","recoAK4Jet_subClusterEtaCenter",25,-3.2,3.2);
		//81 - phi center of GMM cluster from reco jets
		TH1D* recoJet_subClusterPhiCenter = new TH1D("recoAK4Jet_subClusterPhiCenter","recoAK4Jet_subClusterPhiCenter",25,-0.1,6.3);
		//82 - eta sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterEtaSig = new TH1D("recoAK4Jet_subClusterEtaSig","recoAK4Jet_subClusterEtaSig",50,0.,0.1);
		//83 - phi sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterPhiSig = new TH1D("recoAK4Jet_subClusterPhiSig","recoAK4Jet_subClusterPhiSig",50,0.,0.1);
		//84 - time sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterTimeSig = new TH1D("recoAK4Jet_subClusterTimeSig","recoAK4Jet_subClusterTimeSig",50,0.,5.);
		//85 - eta-phi covariance of GMM cluster from reco jets
		TH1D* recoJet_subClusteretaPhiCov = new TH1D("recoAK4Jet_subClusteretaPhiCov","recoAK4Jet_subClusteretaPhiCov",50,-0.0005,0.0005);
		//86 - time-eta covariance of GMM cluster from reco jets
		TH1D* recoJet_subClustertimeEtaCov = new TH1D("recoAK4Jet_subClustertimeEtaCov","recoAK4Jet_subClustertimeEtaCov",50,-0.2,0.2);
		//87 - time-phi covariance of GMM cluster from reco jets
		TH1D* recoJet_subClustertimePhiCov = new TH1D("recoAK4Jet_subClustertimePhiCov","recoAK4Jet_subClustertimePhiCov",50,-0.2,0.2);
		//88 - n rhs in reco jets
		TH1D* recoJet_nRhs = new TH1D("recoAK4Jet_nRhs","recoAK4Jet_nRhs",300,0,300);
		//89 - rh energy in reco jets
		TH1D* recoJet_rhE = new TH1D("recoAK4Jet_rhE","recoAK4Jet_rhE",50,0,300);
		//90 - AK4 Jet P 
		TH1D* AK4Jet_GenP = new TH1D("recoAK4Jet_GenP","recoAK4Jet_GenP",50,0,500);
		//91 - AK4 Jet Pt 
		TH1D* AK4Jet_GenPt = new TH1D("recoAK4Jet_GenPt","recoAK4Jet_GenPt",50,0,500);
		//92 - Gen particle Pt 
		TH1D* AK4JetConstituent_GenP = new TH1D("recoAK4JetConstituent_GenP","recoAK4JetConstituent_GenP",50,0,200);
		//93 - Gen particle Pt 
		TH1D* AK4JetConstituent_GenPt = new TH1D("recoAK4JetConstituent_GenPt","recoAK4JetConstituent_GenPt",50,0,200);
		//94 - Gen particle P/Gen jet P 
		TH1D* AK4JetConstJetRatio_GenP = new TH1D("recoAK4JetConstJetRatio_GenP","recoAK4JetConstJetRatio_GenP",50,0,1);
		//95 - Gen particle Pt/Gen jet Pt 
		TH1D* AK4JetConstJetRatio_GenPt = new TH1D("recoAK4JetConstJetRatio_GenPt","recoAK4JetConstJetRatio_GenPt",50,0,1);
		//96 - ak4 jet rh times
		TH1D* AK4Jet_rhTimes = new TH1D("recoAK4Jet_rhTimes","recoAK4Jet_rhTimes",200,-25,25);
		//97 - time resolution for adjacent xtals
		TH1D* geoEavg_sigmaDeltaTime_adjRhs = new TH1D("geoEavg_sigmaDeltaTime_adjRhs","geoEavg_sigmaDeltaTime_adjRhs;geoEavg;sigmaDeltaTime;a.u.",xbins.size()-1,&xbins[0]);
		//98 - eta center for jet - data statistic
		TH1D* AK4Jet_subClusterEtaCenter_rStat = new TH1D("recoAK4Jet_subClusterEtaCenter_rStat","recoAK4Jet_subClusterEtaCenter_rStat",25,-3.2,3.2);
		//99 - phi center for jet - data statistic
		TH1D* AK4Jet_subClusterPhiCenter_rStat = new TH1D("recoAK4Jet_subClusterPhiCenter_rStat","recoAK4Jet_subClusterPhiCenter_rStat",25,-0.1, 6.3);
		//100 - time center for jet - data statistic
		TH1D* AK4Jet_subClusterTimeCenter_rStat = new TH1D("recoAK4Jet_subClusterTimeCenter_rStat","recoAK4Jet_subClusterTimeCenter_rStat",25,-20.,20.);
		//101 - eta sigma for jet - data statistic
		TH1D* AK4Jet_subClusterEtaSig_rStat = new TH1D("recoAK4Jet_subClusterEtaSig_rStat","recoAK4Jet_subClusterEtaSig_rStat",25,0.,0.1);
		//102 - phi sigma for jet - data statistic
		TH1D* AK4Jet_subClusterPhiSig_rStat = new TH1D("recoAK4Jet_subClusterPhiSig_rStat","recoAK4Jet_subClusterPhiSig_rStat",25,0., 0.1);
		//103 - time sigma for jet - data statistic
		TH1D* AK4Jet_subClusterTimeSig_rStat = new TH1D("recoAK4Jet_subClusterTimeSig_rStat","recoAK4Jet_subClusterTimeSig_rStat",25,0.,5.);
		//104 - eta-phi covariance of GMM cluster from reco jets normalized
		TH1D* recoJet_subClusteretaPhiCovNorm = new TH1D("recoAK4Jet_subClusteretaPhiCovNorm","recoAK4Jet_subClusteretaPhiCovNorm",50,-1.,1.);
		//105 - time-eta covariance of GMM cluster from reco jets normalized
		TH1D* recoJet_subClustertimeEtaCovNorm = new TH1D("recoAK4Jet_subClustertimeEtaCovNorm","recoAK4Jet_subClustertimeEtaCovNorm",50,-1.,1.);
		//106 - time-phi covariance of GMM cluster from reco jets normalized
		TH1D* recoJet_subClustertimePhiCovNorm = new TH1D("recoAK4Jet_subClustertimePhiCovNorm","recoAK4Jet_subClustertimePhiCovNorm",50,-1.,1.);
		//107 - # gen partons (t, b, q from W)
		TH1D* nGenParticles = new TH1D("nGenParticles","nGenParticles",20,0,20);
		//108 - gen particle eta at detector
		TH1D* genParticle_eta = new TH1D("genParticle_eta","genParticle_eta",25,-3.2,3.2);
		//109 - gen particle phi at detector
		TH1D* genParticle_phi = new TH1D("genParticle_phi","genParticle_phi",25,-0.2,6.4);
		//110 - gen particle time at detector
		TH1D* genParticle_time = new TH1D("genParticle_time","genParticle_time",25,-10,10);
		//111 - gen particle pt		
		TH1D* genParticle_pt = new TH1D("genParticle_pt","genParticle_pt",25,0,500);
		//112 - gen particle mass
		TH1D* genParticle_mass = new TH1D("genParticle_mass","genParticle_mass",25,0,200);
		//113 - gen particle energy
		TH1D* genParticle_energy = new TH1D("genParticle_energy","genParticle_energy",25,0,500);
		//114 - gen AK4 jet eta at detector
		TH1D* genAK4Jet_eta = new TH1D("genAK4Jet_EtaCenter","genAK4Jet_EtaCenter",25,-3.2,3.2);
		//115 - gen AK4 jet phi at detector
		TH1D* genAK4Jet_phi = new TH1D("genAK4Jet_PhiCenter","genAK4Jet_PhiCenter",25,-0.2,6.4);
		//116 - gen AK4 jet time at detector
		TH1D* genAK4Jet_time = new TH1D("genAK4Jet_TimeCenter","genAK4Jet_TimeCenter",25,-1,1);
		//117 - gen AK4 jet pt		
		TH1D* genAK4Jet_pt = new TH1D("genAK4Jet_pt","genAK4Jet_pt",25,0,500);
		//118 - gen AK4 jet mass
		TH1D* genAK4Jet_mass = new TH1D("genAK4Jet_mass","genAK4Jet_mass",50,0,50);
		//119 - gen AK4 jet energy
		TH1D* genAK4Jet_energy = new TH1D("genAK4Jet_energy","genAK4Jet_energy",25,0,500);
		//120 - # gen AK4 jets
		TH1D* nJet_genAK4Jet = new TH1D("genAK4_nJets","genAK4_nJets",15,0,15);
		//121 - # constituents per gen AK4 jet
		TH1D* genAK4Jet_nConstituents = new TH1D("genAK4Jet_nConstituents","genAK4Jet_nConstituents",50,0,50);
		//122 - # gen jets - # gen particles
		TH1D* genAK4JetParticle_nDiff = new TH1D("genAK4Jet_genParticle_nDiff","genAK4Jet_genParticle_nDiff",20,-10,10);
		//123 - dR bw gen jet and gen particle its exclusively matched to
		TH1D* genAK4JetParticle_dR = new TH1D("genAK4Jet_genParticle_dR","genAK4Jet_genParticle_dR",25,0,1.5);
		//124 - E ratio bw gen jet and gen particle its exclusively matched to - gen jet energy/gen particle energy
		TH1D* genAK4JetParticle_Eratio = new TH1D("genAK4Jet_genParticle_Eratio","genAK4Jet_genParticle_Eratio",25,0,2);
		//125 - # bhc jets - # gen particles
		TH1D* BHCJetParticle_nDiff = new TH1D("BHCJet_genParticle_nDiff","BHCJet_genParticle_nDiff",20,-10,10);
		//126 - dR bw bhc jet and gen particle its exclusively matched to
		TH1D* BHCJetParticle_dR = new TH1D("BHCJet_genParticle_dR","BHCJet_genParticle_dR",25,0,1.5);
		//127 - E ratio bw bhc jet and gen particle its exclusively matched to - bhc jet energy/gen particle energy
		TH1D* BHCJetParticle_Eratio = new TH1D("BHCJet_genParticle_Eratio","BHCJet_genParticle_Eratio",25,0,2);
		//128 - bhc jet eta center
		TH1D* BHCJet_EtaCenter = new TH1D("BHCJet_EtaCenter","BHCJet_EtaCenter",25,-3.2,3.2);
		//129 - bhc jet phi center
		TH1D* BHCJet_PhiCenter = new TH1D("BHCJet_PhiCenter","BHCJet_PhiCenter",25,-0.1,6.3);
		//130 - bhc jet center
		TH1D* BHCJet_TimeCenter = new TH1D("BHCJet_TimeCenter","BHCJet_TimeCenter",25,-1,1);
		//131 - rh time
		TH1D* rhTime = new TH1D("rhTime","rhTime",25,-10,10);
		//132 - BHC jet rh eta sig
		TH1D* BHCJet_rhEtaSig = new TH1D("BHCJet_rhEtaSig","BHC_rhEtaSig",50,0,1);
		//133 - BHC jet rh phi sig
		TH1D* BHCJet_rhPhiSig = new TH1D("BHCJet_rhPhiSig","BHC_rhPhiSig",50,0,1);
		//134 - BHC jet rh time sig
		TH1D* BHCJet_rhTimeSig = new TH1D("BHCJet_rhTimeSig","BHC_rhTimeSig",50,0,5);
		//135 - recoAK4 jet rh eta sig
		TH1D* recoAK4Jet_rhEtaSig = new TH1D("recoAK4Jet_rhEtaSig","recoAK4_rhEtaSig",50,0,1);
		//136 - recoAK4 jet rh phi sig
		TH1D* recoAK4Jet_rhPhiSig = new TH1D("recoAK4Jet_rhPhiSig","recoAK4_rhPhiSig",50,0,1);
		//137 - BHC jet rh time sig
		TH1D* recoAK4Jet_rhTimeSig = new TH1D("recoAK4Jet_rhTimeSig","recoAK4_rhTimeSig",50,0,5);
		//138 - reco AK4 jet center
		TH1D* recoAK4Jet_TimeCenter = new TH1D("recoAK4Jet_TimeCenter","recoAK4Jet_TimeCenter",25,-1,1);
		//139 - reco AK4 jet eta at detector
		TH1D* recoAK4Jet_EtaCenter = new TH1D("recoAK4Jet_EtaCenter","recoAK4Jet_EtaCenter",25,-3.2,3.2);
		//140 - reco AK4 jet phi at detector
		TH1D* recoAK4Jet_PhiCenter = new TH1D("recoAK4Jet_PhiCenter","recoAK4Jet_PhiCenter",25,-0.2,6.4);
		//141 - reco AK4 jet # subclusters in event
		TH1D* recoAK4Jet_nSubclustersEvt = new TH1D("recoAK4Jet_nSubclustersEvt","recoAK4Jet_nSubclustersEvt",30,0,30);
		//142 - bhc jet # subclusters in event
		TH1D* BHCJet_nSubclustersEvt = new TH1D("BHCJet_nSubclustersEvt","BHCJet_nSubclustersEvt",30,0,30);
		//143 - dr bw subclusters in BHC jet
		TH1D* BHCJet_drSubclusters = new TH1D("BHCJet_drSubclusters","BHCJet_drSubclusters",50,0,0.1);
		//144 - dr bw subclusters in reco AK4 jet
		TH1D* recoAK4Jet_drSubclusters = new TH1D("recoAK4Jet_drSubclusters","recoAK4Jet_drSubclusters",50,0,0.1);
		//145 - rh energy in reco jets
		TH1D* BHCJet_rhE = new TH1D("BHCJet_rhE","BHCJet_rhE",50,0,300);
		//146 - reco AK4 rotundity
		TH1D* recoAK4Jet_rotundity = new TH1D("recoAK4Jet_rotundity","recoAK4Jet_rotundity",50,0.4,1.1);
		//147 - BHC rotundity	
		TH1D* BHCJet_rotundity = new TH1D("BHCJet_rotundity","BHCJet_rotundity",50,0.4,1.1);


		//2D plots
		//0 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_predGen = new TH2D("jetGenE_diffDeltaPt_predGen","jetGenE_diffDeltaPt_predGen;jet_{gen} E (GeV);#Delta p_{T}_{pred, gen} (GeV)",5,0,100,50,-50,50);
		//1 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_recoGen = new TH2D("jetGenE_diffDeltaPt_recoGen","jetGenE_diffDeltaPt_recoGen;jet_{gen} E (GeV);#Delta p_{T}_{reco, gen} (GeV)",4,&xbins_recoGenPt[0],50,-50,50);
		//2 - 2D histogram of gen pT vs reco pT
		TH2D* genPt_recoPt = new TH2D("genPt_recoPt","genPt_recoPt;genpt;recopt",50,5,50,50,5,50);
		//3 - gen jet mass vs reco/gen pt
		TH2D* genJetMass_recoGenPtRatio = new TH2D("genJetMass_recoAK4GenPtRatio","genJetMass_recoAK4GenPtRatio",20,0,10,20,0.8,1.2);
		//4 - reco jet mass vs reco jet pt
		TH2D* recoJetMass_recoJetPt = new TH2D("recoAK4JetMass_recoAK4JetPt","recoAK4JetMass_recoAK4JetPt;recoAK4JetMass;recoAK4JetPt",50,0,250,50,0,250);
		//5 - reco jet mass vs reco jet jetSize
		TH2D* recoJetMass_recoJetSize = new TH2D("recoAK4JetMass_recoAK4JetSize","recoAK4JetMass_recoAK4JetSize;recoAK4JetMass;recoAK4JetSize",50,0,250,50,0,1);
		//6 - reco m_jj ~ W mass vs jet pair jetSize
		TH2D* recoJetInvMassW_recoJetPairjetSize = new TH2D("recoAK4JetInvMassW_recoAK4JetPairjetSize","recoAK4JetInvMassW_recoAK4JetPairjetSize;recoAK4 m_jj;jetSize_jj",50,0,500,50,0,1); 
		//7 - pred jet mass vs pred jet pt
		TH2D* predJetMass_predJetPt = new TH2D("BHCJetMass_BHCJetPt","BHCJetMass_BHCJetPt;BHCJetMass;BHCJetPt",50,0,250,50,0,250);
		//8 - pred jet mass vs pred jet jetSize
		TH2D* predJetMass_predJetSize = new TH2D("BHCJetMass_BHCJetSize","BHCJetMass_BHCJetSize;BHCJetMass;BHCJetSize",50,0,250,50,0,1);
		//9 - pred m_jj ~ W mass vs jet pair jetSize 
		TH2D* predJetInvMassW_predJetPairjetSize = new TH2D("BHCJetInvMassW_BHCJetPairjetSize","BHCJetInvMassW_BHCJetPairjetSize;pred m_jj;jetSize_jj",50,0,500,50,0,1.); 
		//10 - pred jet pt vs pred jet jetSize
		TH2D* predJetPt_predJetSize = new TH2D("BHCJetPt_BHCJetSize","BHCJetPt_BHCJetSize;BHCJetPt;BHCJetSize",50,0,250,50,0,1);
		//11 - pred jet n subclusters vs jet size
		TH2D* prednSubclusters_jetSize = new TH2D("BHCJet_nSubclustersJet_jetSize","BHCJet_nSubclustersJet_jetSize;nSubclusters;jetsize",30,0,30,50,0,1);
		//12 - reco dr match to gen b's
		TH2D* recoJet_genOvRecoE_dR_b = new TH2D("recoAK4Jet_genOvRecoE_dR_b","recoAK4Jet_genOvRecoE_dR_b;ratioE;dR",25,0,5,25,0,4);
		//13 - bhc dr match to gen b's
		TH2D* BHCJet_genOvRecoE_dR_b = new TH2D("BHCJet_genOvRecoE_dR_b","BHCJet_genOvRecoE_dR_b;ratioE;dR",25,0,5,25,0,4);
		//14 - reco dr match to gen q's/g's
		TH2D* recoJet_genOvRecoE_dR_qg = new TH2D("recoAK4Jet_genOvRecoE_dR_qg","recoAK4Jet_genOvRecoE_dR_qg;ratioE;dR",25,0,5,25,0,4);
		//15 - bhc dr match to gen q's/g's
		TH2D* BHCJet_genOvRecoE_dR_qg = new TH2D("BHCJet_genOvRecoE_dR_qg","BHCJet_genOvRecoE_dR_qg;ratioE;dR",25,0,5,25,0,4);
		//16 - reco dr match to gen l's  
		TH2D* recoJet_genOvRecoE_dR_lep = new TH2D("recoAK4Jet_genOvRecoE_dR_lep","recoAK4Jet_genOvRecoE_dR_lep;ratioE;dR",25,0,5,25,0,4);
		//17 - bhc dr match to gen l's  
		TH2D* BHCJet_genOvRecoE_dR_lep = new TH2D("BHCJet_genOvRecoE_dR_lep","BHCJet_genOvRecoE_dR_lep;ratioE;dR",25,0,5,25,0,4);
		//18 - reco dr jet-quark match vs W energy (b's excluded)
		TH2D* recoJet_dRquark_Wenergy = new TH2D("recoAK4Jet_dRquark_Wenergy","recoAK4Jet_dRquark_Wenergy;dRquark;Wenergy",25,0,4,25,0,1000);
		//19 - BHC dr jet-quark match vs W energy (b's excluded)
		TH2D* BHCJet_dRquark_Wenergy = new TH2D("BHCJet_dRquark_Wenergy","BHCJet_dRquark_Wenergy;dRquark;Wenergy",25,0,4,25,0,1000);
		//20 - # rhs vs # subclusters for AK4 jets
		TH2D* AK4Jet_nRhs_nSubclustersJet = new TH2D("recoAK4Jet_nRhs_nSubclustersJet","recoAK4Jet_nRhs_nSubclustersJet;nRhs;nSubclustersJet;a.u.",300,0,300,30,0,30);
		//21 - # gen particles from gen-matched jet vs # subclusters for AK4 jets
		TH2D* AK4Jet_nGenParts_nSubclusters = new TH2D("recoAK4Jet_nGenParts_nSubclusters","recoAK4Jet_nGenParts_nSubclusters;nGenParts;nSubclusters;a.u.",40,0,40,30,0,30);
		//22 - gen particle p vs gen jet p for AK4 jets
		TH2D* AK4Jet_genP_genJetP = new TH2D("recoAK4Jet_genP_genJetP","recoAK4Jet_genP_genJetP;genP;genJetP;a.u.",50,0,200,50,0,500);
		//23 - gen particle pt vs gen jet pt for AK4 jets
		TH2D* AK4Jet_genPt_genJetPt = new TH2D("recoAK4Jet_genPt_genJetPt","recoAK4Jet_genPt_genJetPt;genPt;genJetPt;a.u.",50,0,200,50,0,500);
		//24 - gen particle p vs gen particle p/gen jet p for AK4 jets
		TH2D* AK4Jet_genP_genPartJetPRatio = new TH2D("recoAK4Jet_genP_genPartJetPRatio","recoAK4Jet_genP_genPartJetPRatio;genP;genPartJetPRatio;a.u.",50,0,200,50,0,1);
		//25 - gen jet p vs gen particle p/gen jet p for AK4 jets
		TH2D* AK4Jet_genJetP_genPartJetPRatio = new TH2D("recoAK4Jet_genJetP_genPartJetPRatio","recoAK4Jet_genJetP_genPartJetPRatio;genJetP;genPartJetPRatio;a.u.",50,0,500,50,0,1);
		//26 - gen particle pt vs gen particle pt/gen jet pt for AK4 jets
		TH2D* AK4Jet_genPt_genPartJetPtRatio = new TH2D("recoAK4Jet_genPt_genPartJetPtRatio","recoAK4Jet_genPt_genPartJetPtRatio;genPt;genPartJetPtRatio;a.u.",50,0,200,50,0,1);
		//27 - gen jet pt vs gen particle pt/gen jet pt for AK4 jets
		TH2D* AK4Jet_genJetPt_genPartJetPtRatio = new TH2D("recoAK4Jet_genJetPt_genPartJetPtRatio","recoAK4Jet_genJetPt_genPartJetPtRatio;genJetPt;genPartJetPtRatio;a.u.",50,0,500,50,0,1);
		//28 - gen jet pt vs # subclusters
		TH2D* AK4Jet_genJetPt_nSubclusters = new TH2D("recoAK4Jet_genJetPt_nSubclusters","recoAK4Jet_genJetPt_nSubclusters;genJetPt;nSubclusters",50,0,500,30,0,30);	
		//29 - # gen particles w/ pt > 5 gev from gen-matched jet vs # subclusters for AK4 jets
		TH2D* AK4Jet_nGenPartsptge5_nSubclusters = new TH2D("recoAK4Jet_nGenPartsptge5_nSubclusters","recoAK4Jet_nGenPartsptge5_nSubclusters;nGenPartsptge5;nSubclusters;a.u.",20,0,20,30,0,30);
		//30 - geo energy avg vs difference in time for adjacent crystals in same obj w/in 10% energy 
		TH2D* geoEavg_diffDeltaTime_adjRhs = new TH2D("geoEavg_diffDeltaTime_adjRhs","geoEavg_diffDeltaTime_adjRhs;geoEavg;diffDeltaTime;a.u.",xbins.size()-1,&xbins[0],25,-5,5);
		//31 - eta-phi cov vs time-eta cov 
		TH2D* recoJet_subClusteretaPhiCov_timeEtaCov = new TH2D("recoAK4Jet_subClusteretaPhiCov_timeEtaCov","recoAK4Jet_subClusteretaPhiCov_timeEtaCov;etaPhiCov;timeEtaCov",20,-0.0005,0.0005,20,-0.2,0.2);
		//32 - eta-phi cov norm vs time-eta cov norm 
		TH2D* recoJet_subClusteretaPhiCovNorm_timeEtaCovNorm = new TH2D("recoAK4Jet_subClusteretaPhiCovNorm_timeEtaCovNorm","recoAK4Jet_subClusteretaPhiCovNorm_timeEtaCovNorm;etaPhiCovNorm;timeEtaCovNorm",20,-1.,1.,20, -1, 1.);
		//33 - reco AK4 jet multiplicity vs jet size
		TH2D* recoAK4Jet_nJets_jetSize = new TH2D("recoAK4Jet_nJets_jetSize","recoAK4Jet_nJets_jetSize;nJets;jetSize",10,0,10,50,0,1);
		//34 - BHC jet multiplicity vs jet size
		TH2D* BHCJet_nJets_jetSize = new TH2D("BHCJet_nJets_jetSize","BHCJet_nJets_jetSize;nJets;jetSize",15,0,15,50,0,1);
		//35 - # rhs vs # subclusters for BHC jets
		TH2D* BHCJet_nRhs_nSubclustersJet = new TH2D("BHCJet_nRhs_nSubclustersJet","BHCJet_nRhs_nSubclustersJet;nRhs;nSubclustersJet;a.u.",300,0,300,30,0,30);
		//36 - # subclusters vs # constituents in gen jets for gen-matched AK4 jets
		TH2D* recoAK4Jet_nSubclustersJet_nGenConstituents = new TH2D("recoAK4Jet_nSubclustersJet_nGenConstituents","recoAK4Jet_nSubclustersJet_nGenConstituents;nSubclustersJet;nGenConstituents",30,0,30,100,0,100);
		//37 - # subclusters vs jet mass for reco AK4 jets
		TH2D* recoAK4Jet_nSubclustersJet_mass = new TH2D("recoAK4Jet_nSubclustersJet_mass","recoAK4Jet_nSubclustersJet_mass;nSubclustersJet;mass",30,0,30,50,0,180);
		//38 - # subclusters vs jet energy for reco AK4 jets
		TH2D* recoAK4Jet_nSubclustersJet_energy = new TH2D("recoAK4Jet_nSubclustersJet_energy","recoAK4Jet_nSubclustersJet_energy;nSubclusters;energy",30,0,30,50,0,2000);
		//39 - # subclusters/evt vs # subclusters/jet for reco AK4 jets
		TH2D* recoAK4Jet_nSubclustersEvt_nJet = new TH2D("recoAK4Jet_nSubclustersEvt_nJet","recoAK4Jet_nSubclustersEvt_nJet;nSubclustersEvt;nJet",30,0,30,10,0,10);
		//40 - # subclusters vs jet mass for BHC jets
		TH2D* BHCJet_nSubclustersJet_mass = new TH2D("BHCJet_nSubclustersJet_mass","BHCJet_nSubclustersJet_mass;nSubclustersJet;mass",30,0,30,50,0,180);
		//41 - # subclusters vs jet energy for BHC jets
		TH2D* BHCJet_nSubclustersJet_energy = new TH2D("BHCJet_nSubclustersJet_energy","BHCJet_nSubclustersJet_energy;nSubclusters;energy",30,0,30,50,0,2000);
		//42 - # subclusters/evt vs # subclusters/jet for BHC jets
		TH2D* BHCJet_nSubclustersEvt_nJet = new TH2D("BHCJet_nSubclustersEvt_nJet","BHCJet_nSubclustersEvt_nJet;nSubclustersEvt;nJet",30,0,30,10,0,10);
		//43 - # subclusters in reco AK4 jet and # subclusters in dR matched BHC jet (if matching can be 1:1 ie # BHC jets = # reco AK4 jets)
		TH2D* recoAK4JetnSubclustersJet_BHCJetnSubclustersJet = new TH2D("recoAK4JetnSubclustersJet_BHCJetnSubclustersJet","recoAK4JetnSubclustersJet_BHCJetnSubclustersJet;recoAK4JetnSubclustersJet;BHCJetnSubclustersJet",30,0,30,30,0,30);
		//44 - reco AK4 jet energy vs jet mass
		TH2D* recoAK4Jet_jetEnergy_jetMass = new TH2D("recoAK4Jet_jetEnergy_jetMass","recoAK4Jet_jetEnergy_jetMass;jetEnergy;jetMass",50,0,2000,50,0,180);
		//45 - BHC jet energy vs mass
		TH2D* BHCJet_jetEnergy_jetMass = new TH2D("BHCJet_jetEnergy_jetMass","BHCJet_jetEnergy_jetMass;jetEnergy;jetMass",50,0,2000,50,0,180);
		//46 - reco AK4 jets # subclusters vs jet size
		TH2D* recoAK4Jet_nSubclusters_jetSize = new TH2D("recoAK4Jet_nSubclustersJet_jetSize","recoAK4Jet_nSubclustersJet_jetSize;nSubclusters;jetsize",30,0,30,50,0,1);


		void SetSmear(bool t){ _smear = t; }
		double _cell, _tresCte, _tresNoise, _tresStoch;
		void SetMeasErrParams(double spatial, double tresCte, double tresStoch, double tresNoise){ _cell = spatial; _tresCte = tresCte; _tresStoch = tresStoch; _tresNoise = tresNoise; 
	cout << "Using tres_cte = " << _tresCte << " ns, tres_stoch = " << _tresStoch << " ns and tres_noise = " << _tresNoise << endl;
 }
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

	//find gen jet that most closely is dr matched to jet
	//void GenericMatchJet(vector<Jet>& jets, vector<int>& bestGenMatchIdxs){
	void GenericMatchJet(vector<Jet>& injets, vector<Jet>& matchjets, vector<int>& bestMatchIdxs){
		//loop through gen particles
		//dr match to jet
		double bestDr, dr;
		//int nGen = _genjets.size();
		int nMatch = matchjets.size();
		bestMatchIdxs.clear();
		bestMatchIdxs = {};
	
		//no gen particles to match, set all jets to unmatched	
		if(nMatch < 1){
			for(auto j : injets) bestMatchIdxs.push_back(-1);
			return;
		}
		if(injets.size() < 1) return;

		//vector<double> drs;
		//drs[i][j] = dr for jet i and gen particle j
		vector<vector<double>> drs;
		vector<int> idxs;
		for(int j = 0; j < injets.size(); j++){
			drs.push_back({});
			for(int g = 0; g < nMatch; g++){
				drs[j].push_back(999);
				//rough acceptance cut - if gen jets are out of acceptance any associated jets wouldn't be reconstructed
				//if(fabs(_base->Jet_genEta->at(g)) > 1.5) continue;
				if(fabs(matchjets[g].eta()) > 1.5) continue;
			
				dr = dR(matchjets[g].eta(), matchjets[g].phi(), injets[j].eta(), injets[j].phi());
				drs[j][g] = dr;
			}
		}
		vector<int> best_idxs; //one per jet
		int otherJet, thismatchidx, thisJet;
		//go back through jets can check to see if there are overlapping matches
		for(int j = 0; j < injets.size(); j++){
			//for(int g = 0; g < nMatch; g++){
			//	cout << "jet " << j << " and match jet " << g << " have dr " << drs[j][g] << endl;
			//}
			//cout << "jet " << j << " has best dr " << *min_element(drs[j].begin(), drs[j].end()) << " at match jet " << find(drs[j].begin(), drs[j].end(), *min_element(drs[j].begin(), drs[j].end())) - drs[j].begin() << endl;
			double mindr = *min_element(drs[j].begin(), drs[j].end());
			int matchidx = find(drs[j].begin(), drs[j].end(), mindr) - drs[j].begin();
			if(mindr == 999) matchidx = -1; //no match found (ie no available match jet for best match)
			best_idxs.push_back(matchidx);
			thismatchidx = matchidx;
			thisJet = j;
			//if other jets have the same genidx matched, go through and disambiguate until there is only 1 instance of genidx
			while(count(best_idxs.begin(), best_idxs.end(), matchidx) > 1 && matchidx != -1){
				otherJet = find(best_idxs.begin(), best_idxs.end(), matchidx) - best_idxs.begin();
				//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
				//skip otherJet (ie thisJet) in this case and look at all other jets
				if(otherJet == thisJet){
					otherJet = find(best_idxs.begin()+otherJet+1, best_idxs.end(), matchidx) - best_idxs.begin();
				}
				//for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
				//cout << " found another match at jet " << otherJet << " with other dr " << drs[otherJet][genidx] << " against this jet " << thisJet << endl;
				//if other dr is less than current mindr
				if(drs[otherJet][matchidx] < mindr){
					//set this dr to 999 (is invalid), find new min for this jet, reset genidx to this index
					drs[thisJet][matchidx] = 999;
					mindr = *min_element(drs[thisJet].begin(), drs[thisJet].end());
					if(mindr == 999) matchidx = -1;
					else matchidx = find(drs[thisJet].begin(), drs[thisJet].end(), mindr) - drs[thisJet].begin();
					best_idxs[thisJet] = matchidx;
					//cout << " reset gen match of this jet " << thisJet << " to gen jet " << genidx << " with dr " << mindr << endl;
		
				}
				//if this dr is less than (or equal to) current mindr
				else{
					//set other dr to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
					drs[otherJet][matchidx] = 999;
					thismatchidx = matchidx;
					mindr = *min_element(drs[otherJet].begin(), drs[otherJet].end());
					if(mindr == 999) matchidx = -1;
					else matchidx = find(drs[otherJet].begin(), drs[otherJet].end(), mindr) - drs[otherJet].begin();
					thisJet = otherJet;
					best_idxs[thisJet] = matchidx;
					//cout << " reset match of other jet " << otherJet << " to  jet " << idx << " with dr " << mindr << endl;
				}	
				//cout << "matchidx is now " << matchidx << " with count " << count(best_idxs.begin(), best_idxs.end(), matchidx) << " for jet " << thisJet << endl;

			}
			//cout << "jet " << j << " has best exclusive match with " << best_idxs[j] << "\n" << endl;
		}
		bestMatchIdxs = best_idxs;

	}



	//find gen jet that most closely is dr matched to jet
	void GenMatchJetToJet(vector<Jet>& jets, vector<int>& bestGenMatchIdxs){
		//loop through gen particles
		//dr match to jet
		double bestDr, dr;
		int nGen = _genjets.size();
		//cout << "nGen jets " << nGen << endl;
		int otherJet, thisgenidx, thisJet;
		bestGenMatchIdxs.clear();
		bestGenMatchIdxs = {};
	
		//no gen particles to match, set all jets to unmatched	
		if(nGen < 1){
			for(auto j : jets) bestGenMatchIdxs.push_back(-1);
			return;
		}
		if(jets.size() < 1) return;

		//vector<double> drs;
		//drs[i][j] = dr for jet i and gen particle j
		vector<vector<double>> drs;
		vector<int> genIdxs;
		for(int j = 0; j < jets.size(); j++){
			drs.push_back({});
			for(int g = 0; g < nGen; g++){
				drs[j].push_back(999);
				//rough acceptance cut - if gen jets are out of acceptance any associated jets wouldn't be reconstructed
				//if(fabs(_base->Jet_genEta->at(g)) > 1.5) continue;
				if(fabs(_genjets[g].eta()) > 1.5) continue;
			
				dr = dR(_genjets[g].eta(), _genjets[g].phi(), jets[j].eta(), jets[j].phi());
				drs[j][g] = dr;
			}
		}
		vector<int> best_idxs; //one per jet
		//go back through jets can check to see if there are overlapping matches
		for(int j = 0; j < jets.size(); j++){
			//for(int g = 0; g < nGen; g++){
			//	cout << "jet " << j << " and gen jet " << g << " have dr " << drs[j][g] << endl;
			//}
			//cout << "jet " << j << " has best dr " << *min_element(drs[j].begin(), drs[j].end()) << " at gen jet " << find(drs[j].begin(), drs[j].end(), *min_element(drs[j].begin(), drs[j].end())) - drs[j].begin() << endl;
			double mindr = *min_element(drs[j].begin(), drs[j].end());
			int genidx = find(drs[j].begin(), drs[j].end(), mindr) - drs[j].begin();
			if(mindr == 999) genidx = -1; //no match found (ie no available gen jet for best match)
			best_idxs.push_back(genidx);
			thisgenidx = genidx;
			thisJet = j;
			//if other jets have the same genidx matched, go through and disambiguate until there is only 1 instance of genidx
			while(count(best_idxs.begin(), best_idxs.end(), genidx) > 1 && genidx != -1){
				otherJet = find(best_idxs.begin(), best_idxs.end(), genidx) - best_idxs.begin();
				//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
				//skip otherJet (ie thisJet) in this case and look at all other jets
				if(otherJet == thisJet){
					otherJet = find(best_idxs.begin()+otherJet+1, best_idxs.end(), genidx) - best_idxs.begin();
				}
				//for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
				//cout << " found another match at jet " << otherJet << " with other dr " << drs[otherJet][genidx] << " against this jet " << thisJet << endl;
				//if other dr is less than current mindr
				if(drs[otherJet][genidx] < mindr){
					//set this dr to 999 (is invalid), find new min for this jet, reset genidx to this index
					drs[thisJet][genidx] = 999;
					mindr = *min_element(drs[thisJet].begin(), drs[thisJet].end());
					if(mindr == 999) genidx = -1;
					else genidx = find(drs[thisJet].begin(), drs[thisJet].end(), mindr) - drs[thisJet].begin();
					best_idxs[thisJet] = genidx;
					//cout << " reset gen match of this jet " << thisJet << " to gen jet " << genidx << " with dr " << mindr << endl;
		
				}
				//if this dr is less than (or equal to) current mindr
				else{
					//set other dr to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
					drs[otherJet][genidx] = 999;
					thisgenidx = genidx;
					mindr = *min_element(drs[otherJet].begin(), drs[otherJet].end());
					if(mindr == 999) genidx = -1;
					else genidx = find(drs[otherJet].begin(), drs[otherJet].end(), mindr) - drs[otherJet].begin();
					thisJet = otherJet;
					best_idxs[thisJet] = genidx;
					//cout << " reset gen match of other jet " << otherJet << " to gen jet " << genidx << " with dr " << mindr << endl;
				}	
				//cout << "genidx is now " << genidx << " with count " << count(best_idxs.begin(), best_idxs.end(), genidx) << " for jet " << thisJet << endl;

			}
			//cout << "jet " << j << " has best exclusive gen match with " << best_idxs[j] << "\n" << endl;
		}
		bestGenMatchIdxs = best_idxs;

	}


	//find gen particle that most closely is dr matched to jet
	void GenMatchJet(vector<Jet>& jets, vector<int>& bestGenMatchIdxs){
		//loop through gen particles
		//dr match to jet
		double bestDr, dr;
		int nGen = _genparts.size();//_base->genpart_ngenpart;
		//cout << "nGen " << nGen << endl;
		int otherJet, thisgenidx, thisJet;
		bestGenMatchIdxs.clear();
		bestGenMatchIdxs = {};
	
		//no gen particles to match, set all jets to unmatched	
		if(nGen < 1){
			for(auto j : jets) bestGenMatchIdxs.push_back(-1);
			return;
		}
		if(jets.size() < 1) return;

		//vector<double> drs;
		//drs[i][j] = dr for jet i and gen particle j
		vector<vector<double>> drs;
		vector<int> genIdxs;
		vector<int> qids = {1,2,3,4,5,6};
		for(int j = 0; j < jets.size(); j++){
			drs.push_back({});
			for(int g = 0; g < nGen; g++){
				drs[j].push_back(999);
				////skip W's
				//if(fabs(_base->genpart_id->at(g)) == 24) continue;
				////skip neutrinos
				int truegenidx = _genparts[g].GetUserIdx();
				if(fabs(_base->genpart_id->at(truegenidx)) == 14 || fabs(_base->genpart_id->at(truegenidx)) == 12 || fabs(_base->genpart_id->at(truegenidx)) == 16) continue;
				//only match to quarks
				///if(find(qids.begin(),qids.end(),fabs(_base->genpart_id->at(g))) == qids.end()) continue;
				dr = dR(_base->genpart_eta->at(truegenidx), _base->genpart_phi->at(truegenidx), jets[j].eta(), jets[j].phi());
				drs[j][g] = dr;
			}
		}
		vector<int> best_idxs; //one per jet
		//go back through jets can check to see if there are overlapping matches
		for(int j = 0; j < jets.size(); j++){
			//for(int g = 0; g < nGen; g++){
			//	cout << "jet " << j << " and gen particle " << g << " (id: " << _base->genpart_id->at(_genparts[g].GetUserIdx()) << ") have dr " << drs[j][g] << endl;
			//}
			//cout << "jet " << j << " has best dr " << *min_element(drs[j].begin(), drs[j].end()) << " at gen particle " << find(drs[j].begin(), drs[j].end(), *min_element(drs[j].begin(), drs[j].end())) - drs[j].begin() << endl;
			double mindr = *min_element(drs[j].begin(), drs[j].end());
			int genidx = find(drs[j].begin(), drs[j].end(), mindr) - drs[j].begin();
			if(mindr == 999) genidx = -1; //no match found (ie no available gen particle for best match)
			best_idxs.push_back(genidx);
			thisgenidx = genidx;
			thisJet = j;
			//if other jets have the same genidx matched, go through and disambiguate until there is only 1 instance of genidx
			while(count(best_idxs.begin(), best_idxs.end(), genidx) > 1 && genidx != -1){
				otherJet = find(best_idxs.begin(), best_idxs.end(), genidx) - best_idxs.begin();
				//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
				//skip otherJet (ie thisJet) in this case and look at all other jets
				if(otherJet == thisJet){
					otherJet = find(best_idxs.begin()+otherJet+1, best_idxs.end(), genidx) - best_idxs.begin();
				}
				//for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
				//cout << " found another match at jet " << otherJet << " with other dr " << drs[otherJet][genidx] << " against this jet " << thisJet << endl;
				//if other dr is less than current mindr
				if(drs[otherJet][genidx] < mindr){
					//set this dr to 999 (is invalid), find new min for this jet, reset genidx to this index
					drs[thisJet][genidx] = 999;
					mindr = *min_element(drs[thisJet].begin(), drs[thisJet].end());
					if(mindr == 999) genidx = -1;
					else genidx = find(drs[thisJet].begin(), drs[thisJet].end(), mindr) - drs[thisJet].begin();
					best_idxs[thisJet] = genidx; 
					//cout << " reset gen match of this jet " << thisJet << " to particle " << genidx << " with dr " << mindr << endl;
		
				}
				//if this dr is less than (or equal to) current mindr
				else{
					//set other dr to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
					drs[otherJet][genidx] = 999;
					thisgenidx = genidx;
					mindr = *min_element(drs[otherJet].begin(), drs[otherJet].end());
					if(mindr == 999) genidx = -1;
					else genidx = find(drs[otherJet].begin(), drs[otherJet].end(), mindr) - drs[otherJet].begin();
					thisJet = otherJet;
					best_idxs[thisJet] = genidx; 
					//cout << " reset gen match of other jet " << otherJet << " to particle " << genidx << " with dr " << mindr << endl;
				}	
				//cout << "genidx is now " << genidx << " with count " << count(best_idxs.begin(), best_idxs.end(), genidx) << " for jet " << thisJet << endl;

			}
			//cout << "jet " << j << " has best exclusive gen match with " << best_idxs[j] << " particle with id " << _base->genpart_id->at(_genparts[best_idxs[j]].GetUserIdx()) << "\n" << endl;
		}
		bestGenMatchIdxs = best_idxs;

	}


	void SetAlpha(double a){_alpha = a;}
	void SetSubclusterAlpha(double a){_emAlpha = a; }
	void SetThreshold(double t){ _thresh = t; }
	void SetPriorParameters(map<string, Matrix> params){_prior_params = params;} 		

	void CalcRhTimeDiff(vector<JetPoint>& rhs, vector<pair<double, double>>& geoEavg_diffT){
		geoEavg_diffT.clear();
		double cell = 1;//acos(-1)/180;
		double deta, dphi; //deta, dphi < cell ==> adjacent
		int maxE, lessE;
		double geoEavg, diffT;
		double minE = 1;
		int r_ieta, r_iphi, rr_ieta, rr_iphi;
		for(int r = 0; r < rhs.size(); r++){
			if(rhs[r].E() < minE) continue;	
			for(int rr = r+1; rr < rhs.size(); rr++){
				if(rhs[rr].E() < minE) continue;
				r_ieta = rhs[r].rhId() / 1000;	
				rr_ieta = rhs[rr].rhId() / 1000;	
				r_iphi = rhs[r].rhId() % 1000;	
				rr_iphi = rhs[rr].rhId() % 1000;	
				deta = fabs(r_ieta - rr_ieta);
				dphi = fabs(r_iphi - rr_iphi);
				maxE = rhs[r].E() > rhs[rr].E() ? r : rr;
				lessE = maxE == r ? rr : r;
				geoEavg = sqrt(rhs[r].E() * rhs[rr].E());
				diffT = rhs[maxE].t() - rhs[lessE].t();
				
				if(deta <= cell && dphi <= cell && 0.9*maxE <= lessE){
					//cout << "deta " << deta << " dphi " << dphi << " diffE/maxE " << diffE/maxE << " rhs[r].eta " << rhs[r].eta() << " rhs[rr].eta " << rhs[rr].eta() << " rhs[r].phi " << rhs[r].phi() << " rhs[rr].phi " << rhs[rr].phi() << " rhs[r].E " << rhs[r].E() << " rhs[rr].E " << rhs[rr].E() << " geoEavg " << geoEavg << " diffT " << diffT << endl;
					geoEavg_diffT.push_back(make_pair(geoEavg,diffT));
				}

			}

		}

	}

	private:
		string _oname;
		vector<TH1D*> _hists1D;
		vector<TH2D*> _hists2D;
		vector<TGraph*> graphs;
		vector<Jet> _phos; //photons for event
		vector<procCat> _procCats;
		vector<node*> _trees;
		vector<Jet> _predJets, _genjets, _recojets, _genparts;
		bool _smear;
		enum Strategy{
			//Delauney strategy - NlnN time - for 2pi cylinder
			NlnN = 0,
			//traditional strategy - N^2 time
			N2 = 1,
			//gmm only on reco'd jets
			gmmOnly = 2,
			//NlnN using rhs from reco AK4 jets
			NlnNonAK4 = 3,
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

		

		map<string, Matrix> _prior_params;

		double dR(double eta1, double phi1, double eta2, double phi2){
			//phi wraparound
			double dphi = (phi1-phi2);
			dphi = acos(cos(dphi));
			return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
		}


		void Get2DMat(const Matrix& inmat, Matrix& outmat){
			if(!outmat.square()) return;
			if(outmat.GetDims()[0] != 2) return;
			outmat.reset();
			outmat.SetEntry(inmat.at(0,0),0,0);	
			outmat.SetEntry(inmat.at(0,1),0,1);	
			outmat.SetEntry(inmat.at(1,0),1,0);	
			outmat.SetEntry(inmat.at(1,1),1,1);
		}
		double Rotundity(Matrix& inmat){
			vector<Matrix> eigenvecs;
			vector<double> eigenvals;
			inmat.eigenCalc(eigenvals, eigenvecs);
			int maxd = inmat.GetDims()[0] - 1;
			double rot = 0;
			for(int i = 0; i < (int)eigenvals.size(); i++) rot += eigenvals[i];
			rot = eigenvals[maxd]/rot;
			//if(rot < 0.5 || rot > 1) cout << "rot: " << rot << endl;
			return rot;
		}
		
};
#endif
