#ifndef BHCJETSKIMMER_HH
#define BHCJETSKIMMER_HH
#include "JetSimProducer.hh"
#include "BaseSkimmer.hh"
#include "BaseTree.hh"
#include "TGraph.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TExec.h"
#include "TLine.h"
#include <set>
#include <Math/Vector4D.h>

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
			_evt2disp = 0;
			_evt2disp_z = 0;

			_minTopPt = 0;
			_minTopE = 0;
			_minWPt = 0;		

			_nGhosts = 0;
			_check_merges = false;
	
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

			_strategy = NlnN;
			_sel = def;

			_infile = nullptr;
		}

		virtual ~BHCJetSkimmer(){ }

		//BHCJetSkimmer(TFile* file){
		BHCJetSkimmer(string file){
			if(gSystem->AccessPathName(file.c_str())){ cout << "Error: file " << file << " doesn't exist." << endl; return; }
			InitHists();
			_infile = TFile::Open(file.c_str());
			_prod = new JetSimProducer(_infile);
			_strategy = NlnN;
			_sel = def;

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
			_minTopPt = 0;
			_minWPt = 0;		
			_minTopE = 0;
			_evt2disp = 0;
			_evt2disp_z = 0;
			
			_nGhosts = 0;
			_check_merges = false;
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
			

		}
		void InitHists(){
	
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
			_hists1D.push_back(recoGen_jetPtRatio);
			_hists1D.push_back(recoGen_jetERatio);
			_hists1D.push_back(predJet_EtaVar);
			_hists1D.push_back(predJet_PhiVar);
			_hists1D.push_back(predJet_TimeVar);
			_hists1D.push_back(predJet_etaPhiCov);
			_hists1D.push_back(predJet_timeEtaCov);
			_hists1D.push_back(predJet_timePhiCov);
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
			_hists1D.push_back(genAK4JetTop_dR);
			_hists1D.push_back(genAK4JetTop_Eratio);
			_hists1D.push_back(BHCJetParticle_nDiff);
			_hists1D.push_back(BHCJetTop_dR);
			_hists1D.push_back(BHCJetTop_Eratio);
			_hists1D.push_back(BHCJet_EtaCenter);
			_hists1D.push_back(BHCJet_PhiCenter);
			_hists1D.push_back(BHCJet_TimeCenter);
			_hists1D.push_back(rhTime);
			_hists1D.push_back(recoAK4Jet_rhEtaSig);
			_hists1D.push_back(recoAK4Jet_rhPhiSig);
			_hists1D.push_back(recoAK4Jet_rhTimeSig);
			_hists1D.push_back(recoAK4Jet_TimeCenter);
			_hists1D.push_back(recoAK4Jet_EtaCenter);
			_hists1D.push_back(recoAK4Jet_PhiCenter);
			_hists1D.push_back(recoAK4Jet_nSubclustersEvt);
			_hists1D.push_back(BHCJet_drSubclusters);
			_hists1D.push_back(recoAK4Jet_drSubclusters);
			_hists1D.push_back(BHCJet_rhE);
			_hists1D.push_back(recoAK4Jet_rotundity);
			_hists1D.push_back(BHCJet_rotundity);
			_hists1D.push_back(BHCJet_nRhs);
			_hists1D.push_back(genAK15Jet_eta);
			_hists1D.push_back(genAK15Jet_phi);
			_hists1D.push_back(genAK15Jet_time);
			_hists1D.push_back(genAK15Jet_pt);
			_hists1D.push_back(genAK15Jet_mass);
			_hists1D.push_back(genAK15Jet_energy);
			_hists1D.push_back(nJet_genAK15Jet);
			_hists1D.push_back(genAK15Jet_nConstituents);
			_hists1D.push_back(genAK15JetParticle_nDiff);
			_hists1D.push_back(genAK15JetTop_dR);
			_hists1D.push_back(genAK15JetTop_Eratio);
			_hists1D.push_back(nRecoAK8Jets);
			_hists1D.push_back(recoAK8Jet_eta);
			_hists1D.push_back(recoAK8Jet_phi);
			_hists1D.push_back(recoAK8Jet_time);
			_hists1D.push_back(recoAK8Jet_pt);
			_hists1D.push_back(recoAK8Jet_mass);
			_hists1D.push_back(recoAK8Jet_energy);
			_hists1D.push_back(recoAK8Jet_nConstituents);
			_hists1D.push_back(recoAK8JetTop_dR);
			_hists1D.push_back(recoAK8JetTop_Eratio);
			_hists1D.push_back(recoAK8Jet_jetSize);
			_hists1D.push_back(nRecoAK15Jets);
			_hists1D.push_back(recoAK15Jet_eta);
			_hists1D.push_back(recoAK15Jet_phi);
			_hists1D.push_back(recoAK15Jet_time);
			_hists1D.push_back(recoAK15Jet_pt);
			_hists1D.push_back(recoAK15Jet_mass);
			_hists1D.push_back(recoAK15Jet_energy);
			_hists1D.push_back(recoAK15Jet_nConstituents);
			_hists1D.push_back(recoAK15JetTop_dR);
			_hists1D.push_back(recoAK15JetTop_Eratio);
			_hists1D.push_back(recoAK15Jet_jetSize);
			_hists1D.push_back(nGenAK8Jets);
			_hists1D.push_back(genAK8Jet_eta);
			_hists1D.push_back(genAK8Jet_phi);
			_hists1D.push_back(genAK8Jet_time);
			_hists1D.push_back(genAK8Jet_pt);
			_hists1D.push_back(genAK8Jet_mass);
			_hists1D.push_back(genAK8Jet_energy);
			_hists1D.push_back(genAK8Jet_nConstituents);
			_hists1D.push_back(genAK8JetTop_dR);
			_hists1D.push_back(genAK8JetTop_Eratio);
			_hists1D.push_back(genAK4JetW_dR);
			_hists1D.push_back(genAK4JetW_Eratio);
			_hists1D.push_back(genAK15JetW_dR);
			_hists1D.push_back(genAK15JetW_Eratio);
			_hists1D.push_back(BHCJetW_dR);
			_hists1D.push_back(BHCJetW_Eratio);
			_hists1D.push_back(BHCJetW_nSubclusters);
			_hists1D.push_back(BHCJetW_subClusterEnergy);
                	_hists1D.push_back(BHCJet_subClusterMass);
                	_hists1D.push_back(BHCJet_subClusterEffnRhs);
                	_hists1D.push_back(recoAK4Jet_subClusterMass);
                	_hists1D.push_back(recoAK4Jet_subClusterEffnRhs);
			_hists1D.push_back(BHCJetW_subClusterMass);
			_hists1D.push_back(BHCJetW_subClusterLeadInvMass);
			_hists1D.push_back(BHCJet_nGhosts);
			_hists1D.push_back(BHCJet_ghostSubClusterEnergy);
			_hists1D.push_back(BHCJet_ghostSubClusterEffnRhs);
			_hists1D.push_back(BHCJetW_subclParton_dR);
			_hists1D.push_back(BHCJetW_subclParton_Eratio);
			_hists1D.push_back(BHCJetTop_nSubclusters);
			_hists1D.push_back(BHCJetTop_subClusterMass);
			_hists1D.push_back(BHCJetTop_subClusterLeadInvMass);
			_hists1D.push_back(BHCJetq_dR);
			_hists1D.push_back(BHCJetq_Eratio);
			_hists1D.push_back(BHCJetq_nSubclusters);
			_hists1D.push_back(BHCJetq_subClusterMass);
			_hists1D.push_back(BHCJetW_subclEtaCenter);
			_hists1D.push_back(BHCJetW_subclPhiCenter);
			_hists1D.push_back(BHCJetW_subclTimeCenter);
			_hists1D.push_back(BHCJetW_subClusterEtaSig);
			_hists1D.push_back(BHCJetW_subClusterPhiSig);
			_hists1D.push_back(BHCJetW_subClusterTimeSig);
			_hists1D.push_back(BHCJetW_subClusteretaPhiCov);
			_hists1D.push_back(BHCJetW_subClustertimeEtaCov);
			_hists1D.push_back(BHCJetW_subClustertimePhiCov);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclPt);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclPt);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclSize);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclSize);
			_hists1D.push_back(recoAK8JetW_dR);
			_hists1D.push_back(recoAK8JetW_Eratio);
			_hists1D.push_back(recoAK8JetGluon_dR);
			_hists1D.push_back(recoAK8JetGluon_Eratio);
			_hists1D.push_back(BHCJetGluon_dR);
			_hists1D.push_back(BHCJetGluon_Eratio);
			_hists1D.push_back(BHCJetW_highMass_nSubclustersJet);
			_hists1D.push_back(recoAK8Jetq_dR);
			_hists1D.push_back(recoAK8Jetq_Eratio);
			_hists1D.push_back(recoAK15JetW_dR);
			_hists1D.push_back(recoAK15JetW_Eratio);
			_hists1D.push_back(recoAK15JetGluon_dR);
			_hists1D.push_back(recoAK15JetGluon_Eratio);
			_hists1D.push_back(recoAK15Jetq_dR);
			_hists1D.push_back(recoAK15Jetq_Eratio);
			_hists1D.push_back(recoAK4JetGluon_dR);
			_hists1D.push_back(recoAK4JetGluon_Eratio);
			_hists1D.push_back(recoAK4Jetq_dR);
			_hists1D.push_back(recoAK4Jetq_Eratio);
			_hists1D.push_back(recoAK4JetW_dR);
			_hists1D.push_back(recoAK4JetW_Eratio);
			_hists1D.push_back(BHCJetq_ge2Subcls_partonMatchSubclPt);
			_hists1D.push_back(BHCJetq_ge2Subcls_partonNoMatchSubclPt);
			_hists1D.push_back(BHCJetq_ge2Subcls_partonMatchSubclSize);
			_hists1D.push_back(BHCJetq_ge2Subcls_partonNoMatchSubclSize);
			_hists1D.push_back(BHCJetW_lowMass_nSubclustersJet);
			_hists1D.push_back(BHCJetW_Wmass_nSubclustersJet);
			_hists1D.push_back(BHCJetGluon_nSubclusters);
			_hists1D.push_back(BHCJetq_subclParton_dR);
			_hists1D.push_back(BHCJetq_subclParton_Eratio);
			_hists1D.push_back(BHCJetGluon_subclParton_dR);
			_hists1D.push_back(BHCJetGluon_subclParton_Eratio);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclPtOvJetPt);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclPtOvJetPt);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclSizeOvJetSize);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize);

			_hists2D.push_back(jetGenE_diffDeltaPt_recoGen);
			_hists2D.push_back(genPt_recoPt);
			_hists2D.push_back(genJetMass_recoGenPtRatio);
			_hists2D.push_back(recoJetMass_recoJetPt);
			_hists2D.push_back(recoJetMass_recoJetSize);
			_hists2D.push_back(predJetMass_predJetPt);
			_hists2D.push_back(predJetPt_predJetSize);
			_hists2D.push_back(prednSubclusters_jetSize);
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
			_hists2D.push_back(recoAK4nSubclusters_genAK4nConstituents);
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
			_hists2D.push_back(recoAK4Jet_jetEnergy_jetSize);
			_hists2D.push_back(BHCJet_jetEnergy_jetSize);
                	_hists2D.push_back(recoAK4JetPt_genAK4JetPt);
                	_hists2D.push_back(recoAK4JetE_genAK4JetE);
                	_hists2D.push_back(recoAK4JetMass_genAK4JetMass);
                	_hists2D.push_back(recoAK4JetEtaCenter_genAK4JetEtaCenter);
                	_hists2D.push_back(recoAK4JetPhiCenter_genAK4JetPhiCenter);
                	_hists2D.push_back(recoAK4JetTimeCenter_genAK4JetTimeCenter);
                	_hists2D.push_back(BHCJetPt_genAK4JetPt);
                	_hists2D.push_back(BHCJetE_genAK4JetE);
                	_hists2D.push_back(BHCJetMass_genAK4JetMass);
                	_hists2D.push_back(BHCJetEtaCenter_genAK4JetEtaCenter);
                	_hists2D.push_back(BHCJetPhiCenter_genAK4JetPhiCenter);
                	_hists2D.push_back(BHCJetTimeCenter_genAK4JetTimeCenter);
                	_hists2D.push_back(BHCnSubclusters_genAK4nConstituents);
                	_hists2D.push_back(BHCJetPt_genAK15JetPt);
                	_hists2D.push_back(BHCJetE_genAK15JetE);
                	_hists2D.push_back(BHCJetMass_genAK15JetMass);
                	_hists2D.push_back(BHCJetEtaCenter_genAK15JetEtaCenter);
                	_hists2D.push_back(BHCJetPhiCenter_genAK15JetPhiCenter);
                	_hists2D.push_back(BHCJetTimeCenter_genAK15JetTimeCenter);
                	_hists2D.push_back(BHCnSubclusters_genAK15nConstituents);
			_hists2D.push_back(BHCJetPt_genTopPt);
			_hists2D.push_back(BHCJetE_genTopE);
			_hists2D.push_back(BHCJetMass_genTopMass);
			_hists2D.push_back(BHCJetEtaCenter_genTopEtaCenter);
			_hists2D.push_back(BHCJetPhiCenter_genTopPhiCenter);
			_hists2D.push_back(recoAK4JetPt_genTopPt);
			_hists2D.push_back(recoAK4JetE_genTopE);
			_hists2D.push_back(recoAK4JetMass_genTopMass);
			_hists2D.push_back(recoAK4JetEtaCenter_genTopEtaCenter);
			_hists2D.push_back(recoAK4JetPhiCenter_genTopPhiCenter);
			_hists2D.push_back(recoAK15JetPt_genTopPt);
			_hists2D.push_back(recoAK15JetE_genTopE);
			_hists2D.push_back(recoAK15JetMass_genTopMass);
			_hists2D.push_back(recoAK15JetEtaCenter_genTopEtaCenter);
			_hists2D.push_back(recoAK15JetPhiCenter_genTopPhiCenter);
			_hists2D.push_back(genAK4JetPt_genTopPt);
			_hists2D.push_back(genAK4JetE_genTopE);
			_hists2D.push_back(genAK4JetMass_genTopMass);
			_hists2D.push_back(genAK4JetEtaCenter_genTopEtaCenter);
			_hists2D.push_back(genAK4JetPhiCenter_genTopPhiCenter);
			_hists2D.push_back(genAK15JetPt_genTopPt);
			_hists2D.push_back(genAK15JetE_genTopE);
			_hists2D.push_back(genAK15JetMass_genTopMass);
			_hists2D.push_back(genAK15JetEtaCenter_genTopEtaCenter);
			_hists2D.push_back(genAK15JetPhiCenter_genTopPhiCenter);
			_hists2D.push_back(BHCJetPt_genWPt);
			_hists2D.push_back(BHCJetE_genWE);
			_hists2D.push_back(BHCJetMass_genWMass);
			_hists2D.push_back(BHCJetEtaCenter_genWEtaCenter);
			_hists2D.push_back(BHCJetPhiCenter_genWPhiCenter);
			_hists2D.push_back(recoAK4JetPt_genWPt);
			_hists2D.push_back(recoAK4JetE_genWE);
			_hists2D.push_back(recoAK4JetMass_genWMass);
			_hists2D.push_back(recoAK4JetEtaCenter_genWEtaCenter);
			_hists2D.push_back(recoAK4JetPhiCenter_genWPhiCenter);
			_hists2D.push_back(recoAK15JetPt_genWPt);
			_hists2D.push_back(recoAK15JetE_genWE);
			_hists2D.push_back(recoAK15JetMass_genWMass);
			_hists2D.push_back(recoAK15JetEtaCenter_genWEtaCenter);
			_hists2D.push_back(recoAK15JetPhiCenter_genWPhiCenter);
			_hists2D.push_back(genAK4JetPt_genWPt);
			_hists2D.push_back(genAK4JetE_genWE);
			_hists2D.push_back(genAK4JetMass_genWMass);
			_hists2D.push_back(genAK4JetEtaCenter_genWEtaCenter);
			_hists2D.push_back(genAK4JetPhiCenter_genWPhiCenter);
			_hists2D.push_back(genAK15JetPt_genWPt);
			_hists2D.push_back(genAK15JetE_genWE);
			_hists2D.push_back(genAK15JetMass_genWMass);
			_hists2D.push_back(genAK15JetEtaCenter_genWEtaCenter);
			_hists2D.push_back(genAK15JetPhiCenter_genWPhiCenter);
                	_hists2D.push_back(BHCJet_subclusterEnergy_subclusterMass);
                	_hists2D.push_back(BHCJet_subclusterEnergy_subclusterEffnRhs);
                	_hists2D.push_back(BHCJet_subclusterMass_subclusterEffnRhs);
                	_hists2D.push_back(recoAK4Jet_subclusterEnergy_subclusterMass);
                	_hists2D.push_back(recoAK4Jet_subclusterEnergy_subclusterEffnRhs);
                	_hists2D.push_back(recoAK4Jet_subclusterMass_subclusterEffnRhs);
			_hists2D.push_back(BHCJet_subclusterEnergy_subclusterdRToJet);
			_hists2D.push_back(BHCJet_subclusterEffnRhs_subclusterdRToJet);
			_hists2D.push_back(BHCJetW_nSubclustersJet_mass);
			_hists2D.push_back(BHCJetW_dRGenPartons_nSubclustersJet);
			_hists2D.push_back(BHCJetW_dRGenPartons_jetSize);
			_hists2D.push_back(BHCJetW_EratioSubclGenPart_nSubclustersJet);
			_hists2D.push_back(BHCJetTop_nSubclustersJet_mass);
			_hists2D.push_back(BHCJetW_openAng_nJets);
			_hists2D.push_back(BHCJetW_openAng_nSubclustersJet);
			_hists2D.push_back(BHCJetW_openAng_subclMass);
			_hists2D.push_back(BHCJetW_1subcl_dRGenPartons_jetSize);
			_hists2D.push_back(BHCJetW_ge2subcl_dRGenPartons_jetSize);
			_hists2D.push_back(BHCJetW_1subcl_dRGenPartons_avgPartonEnergy);
			_hists2D.push_back(BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy);
			_hists2D.push_back(BHCJetW_ge2subcl_subclEnergy_subclLeadIdx);
			_hists2D.push_back(BHCJetW_EratioJetGenW_nSubclustersJet);
			_hists2D.push_back(BHCJet_jetMass_jetSize);
			_hists2D.push_back(EvtDisplay_etaCell_phiCell);
			_hists2D.push_back(recoAK8JetMass_recoAK8JetSize);
			_hists2D.push_back(recoAK15JetMass_recoAK15JetSize);
			_hists2D.push_back(recoAK8JetMass_BHCJetMass_matched);
			_hists2D.push_back(recoAK8JetW_dRGenPartons_jetSize);
			_hists2D.push_back(BHCJet_subclTimeSig_subclSize);
			_hists2D.push_back(BHCJet_subclTimeSig_subclPt);
			_hists2D.push_back(recoAK15JetW_dRGenPartons_jetSize);
			_hists2D.push_back(recoAK4JetW_dRGenPartons_jetSize);
			_hists2D.push_back(BHCJetW_highMass_partonMatchSubclPt_subclTime);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchSubclPt_subclTime);
			_hists2D.push_back(BHCJetW_highMass_partonMatchSubclSize_subclTime);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchSubclSize_subclTime);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclMass_RelSubclPt);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclMass_RelSubclPt);
			_hists2D.push_back(BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclPt);
			_hists2D.push_back(BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclPt);
			_hists2D.push_back(BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclSize);
			_hists2D.push_back(BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclSize);
			_hists2D.push_back(BHCJetW_highMassOvPt_partonMatchSubclMassOvPt_JetMassOvPt);
			_hists2D.push_back(BHCJetW_highMassOvPt_partonNoMatchSubclMassOvPt_JetMassOvPt);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclSizeSpatioTemp_RelSubclPt);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclSizeSpatioTemp_RelSubclPt);

			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_W);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_W2);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_gluon);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_q1);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_q2);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_b1);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_b2);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_top1);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_top2);
		}
		void SetMinRhE(double r){ _prod->SetMinRhE(r); }
		void SetRecoMinPt(double r){ _prod->SetRecoMinPt(r); }
		void SetRecoMinE(double r){ _prod->SetRecoMinE(r); }
		void SetGenMinPt(double r){ _prod->SetGenMinPt(r); }
		void SetGenMinE(double r){ _prod->SetGenMinE(r); }
		void SetMinNrhs(int r){ _prod->SetMinNrhs(r); }
		void SetMinNGenConsts(int r){ _prod->SetMinNGenConsts(r); }
		void SetGenTopMinPt(double r){ _minTopPt = r; 
			cout << "Minimum gen top pt: " << _minTopPt << endl;}
		void SetGenWMinPt(double r){ _minWPt = r; 
			cout << "Minimum gen W pt: " << _minWPt << endl;}
		bool _check_merges;
		void CheckMerges(bool t){ _check_merges = t; if(_check_merges) cout << "Checking merges" << endl; }
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else if(i == 2) _strategy = gmmOnly;
			else if(i == 3) _strategy = NlnNonAK4;
			else return; 
		}
		void SetEventSelection(int i){
			//default selection (generic hadronic)
			if(i == 0) _sel = def;
			else if(i == 1) _sel = singW;//single W
			else if(i == 2) _sel = boostTop;//ttbar
			else if(i == 3) _sel = QCDdijets;//QCD
			else return; 
		}
		int _nGhosts;
		void SetNGhosts(int t){ _nGhosts = t; cout << "Adding " << _nGhosts << " ghosts to each BHC merging step." << endl;}
		//set evt display hist name according to event #
		//z is what to put on z axis (0 : energy, 1 : time)
		void SetEvent2Display(int e, int z = 0){
			_evt2disp = e;
			_evt2disp_z = z;	
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
			//cout << "pre pt cut - pred jet px " << predJet.px() << " py " << predJet.py() << " pz " << predJet.pz() << " pt " << predJet.pt() << " E " << predJet.E() << " m2 " << predJet.m2() << " mass " << predJet.mass() << " eta " << predJet.eta() << " phi " << predJet.phi() << endl;
				//put pt cut in for predjets of 5 GeV to match reco AK4 definition
				if(predJet.pt() < 5) continue; 
				//add Jet to jets	
				_predJets.push_back(predJet);	
				//fill hists for mixture model
				for(int p = 0; p < _procCats.size(); p++){
					for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
							//define pt bins
							//pt == 1 -> [_pt_thresh,inf)
							if(pt == 1 && predJet.pt() < _pt_thresh) continue;
							//pt == 2 -> [0,_pt_thresh)
							if(pt == 2 && predJet.pt() >= _pt_thresh) continue;
							_procCats[p].hists1D[pt][148]->Fill(_trees[i]->model->GetNGhosts());
							vector<double> norms;
							_trees[i]->model->GetNorms(norms);
							Matrix r_post = _trees[i]->model->GetPosterior();
							for(int k = 0; k < _trees[i]->model->GetNClusters(); k++){
								auto params = _trees[i]->model->GetLHPosteriorParameters(k);
								bool g = (bool)params["ghost"].at(0,0);
								if(!g) continue;
								_procCats[p].hists1D[pt][149]->Fill(norms[k]/_gev);
								double eff_rhs = 0;
								for(int n = 0; n < _trees[i]->model->GetData()->GetNPoints(); n++){
									eff_rhs += r_post.at(n,k)/_trees[i]->model->GetData()->at(n).w();
								}
								_procCats[p].hists1D[pt][150]->Fill(eff_rhs);
							}
					}
				}
			}
			//sort predicted jets
			sort(_predJets.begin(), _predJets.end(), ptsort);
			cout << _predJets.size()  << " pred jets total" << endl;
			//cout << _predJets.size() << " pred jets pt > 20 GeV" << endl;
			for(auto j : _predJets) cout << "pred jet px " << j.px() << " py " << j.py() << " pz " << j.pz() << " E " << j.E() << " m2 " << j.m2() << " mass " << j.mass() << " eta " << j.eta() << " phi " << j.phi() << " pt " << j.pt() << endl;
		}

		void FillPredJetHists(){
			int njets = _predJets.size();
			int tmass_idx, id, widx;
			tmass_idx = -999;

			double wmass, topmass, dr, dr_pair, jetsize, jetsize_full;
			pair<int, int> wmass_idxs = make_pair(-999, -999);
			Jet w, top, genpart, genW;
			//do gen matching
			vector<int> genAK4MatchIdxs, genAK15MatchIdxs; //one per jet, follows same indexing as jets
			vector<int> recoMatchIdxs(_predJets.size(),-1);
			vector<int> recoAK8MatchIdxs(_predJets.size(), -1);
			GenericMatchJet(_predJets,_recoAK4jets, recoMatchIdxs); //match BHC jets to reco jets
			GenericMatchJet(_predJets,_recoAK4jets, recoAK8MatchIdxs); //match BHC jets to reco AK8 jets
			GenericMatchJet(_predJets,_genAK4jets, genAK4MatchIdxs); //match BHC jets to gen AK4 jets
			GenericMatchJet(_predJets,_genAK15jets, genAK15MatchIdxs); //match BHC jets to gen AK15 jets
			vector<int> genTopMatchIdxs(_predJets.size(),-1);
			vector<int> genWMatchIdxs(_predJets.size(),-1);
			vector<int> genGluonMatchIdxs(_predJets.size(),-1);
			vector<int> genqMatchIdxs(_predJets.size(),-1);
			double openAng = -1;
			if(_sel == boostTop){
				GenericMatchJet(_predJets,_genTop, genTopMatchIdxs); //match BHC jets to good gen tops
			}
			else if(_sel == singW && _genq.size() > 1){
				GenericMatchJet(_predJets,_genW, genWMatchIdxs); //match BHC jets to good gen Ws
				GenericMatchJet(_predJets,_genglu, genGluonMatchIdxs); //match BHC jets to good gen gluons
				GenericMatchJet(_predJets,_genq, genqMatchIdxs); //match BHC jets to good gen qs
				//calculate polar opening angle bw gen partons from W
				vector<double> v1, v2;

				GetXYZVec(_genq[0],v1);
				GetXYZVec(_genq[1],v2);

				double dotprod = dot(v1,v2);
				openAng = acos(dotprod);
				//cout << "dotprod " << dotprod << " openAng " << openAng << endl;
			}
			else if(_sel == QCDdijets){
				GenericMatchJet(_predJets,_genq, genqMatchIdxs); //match BHC jets to good gen Ws
			}
			else{ }
	
			int nsubs, bestMatchIdx;
			for(int p = 0; p < _procCats.size(); p++){
				//if(p != 0) cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				_procCats[p].hists1D[0][0]->Fill(njets);
				cout << "# pred jets - # reco AK4 jets " << njets - (int)_recoAK4jets.size() << endl;
				cout << "# pred jets - # gen AK4 jets " << njets - (int)_genAK4jets.size() << endl;
				_procCats[p].hists1D[0][10]->Fill(njets - (int)_recoAK4jets.size());
				//count # jets over pt thresh for lead/not lead
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _predJets.size(); j++){
					if(_predJets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_predJets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][0]->Fill(njets_lead);
				_procCats[p].hists1D[2][0]->Fill(njets_notlead);
				if(openAng != -1){
					_procCats[p].hists2D[0][119]->Fill(openAng,njets);

					_procCats[p].hists2D[1][119]->Fill(openAng,njets_lead);
					_procCats[p].hists2D[2][119]->Fill(openAng,njets_notlead);
					
				}
				//loop over pt bins
				njets = _predJets.size();
				for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
					for(int j = 0; j < _predJets.size(); j++){
						//define pt bins
						//pt == 1 -> [_pt_thresh,inf)
						if(pt == 1 && _predJets[j].pt() < _pt_thresh) continue;
						//pt == 2 -> [0,_pt_thresh)
						if(pt == 2 && _predJets[j].pt() >= _pt_thresh) continue;

						vector<JetPoint> rhs = _predJets[j].GetJetPoints();
						Matrix jet_mu, jet_cov;
						_predJets[j].GetClusterParams(jet_mu,jet_cov);	
						//also include rotundity
						jetsize = CalcSize(jet_cov);
						jetsize_full = CalcSize(jet_cov, true);
						//define jetsize bins
						//pt == 1 -> [0,0.2) (AK4 level)
						//if(pt == 1 && jetsize >= 0.2) continue;
						////pt == 2 -> [0.2,inf) (AK15 level)
						//if(pt == 2 && jetsize < 0.2) continue;


						if(p != 0 && pt == 0){ 
							cout << "pred jet #" << j << " phi " << _predJets[j].phi() << " eta " << _predJets[j].eta() << " energy " << _predJets[j].E() <<  " mass " << _predJets[j].mass() << " nConstituents " << _predJets[j].GetNConstituents() << " nRhs " << _predJets[j].GetNRecHits() << " pt " << _predJets[j].pt() << " jetsize " << jetsize << " eta var " << jet_cov.at(0,0) << " phi var " << jet_cov.at(1,1) << endl;
						}
			
						double rot = Rotundity(jet_cov);
						_procCats[p].hists1D[pt][89]->Fill(rot);
						
						nsubs = _predJets[j].GetNConstituents();
						_procCats[p].hists1D[pt][1]->Fill(nsubs);
						if(openAng != -1) _procCats[p].hists2D[pt][120]->Fill(openAng,nsubs);
						//cout << "pred jet j " << j << " dr " << dr << " n constituents " << _predJets[j].GetNConstituents() << endl;
						Matrix cov = _predJets[j].GetCovariance();
					
						_procCats[p].hists1D[pt][7]->Fill(_predJets[j].e());
						_procCats[p].hists1D[pt][8]->Fill(_predJets[j].pt());
						_procCats[p].hists1D[pt][9]->Fill(_predJets[j].mass());
						_procCats[p].hists1D[pt][6]->Fill(jetsize);
							
						_procCats[p].hists1D[pt][24]->Fill(sqrt(cov.at(0,0)));	
						_procCats[p].hists1D[pt][25]->Fill(sqrt(cov.at(1,1)));	
						_procCats[p].hists1D[pt][26]->Fill(sqrt(cov.at(2,2)));	
						_procCats[p].hists1D[pt][27]->Fill(cov.at(0,1));	
						_procCats[p].hists1D[pt][28]->Fill(cov.at(0,2));	
						_procCats[p].hists1D[pt][29]->Fill(cov.at(2,1));	

						_procCats[p].hists2D[pt][5]->Fill(_predJets[j].mass(), _predJets[j].pt());
						_procCats[p].hists2D[pt][6]->Fill(_predJets[j].pt(), jetsize);
						
						_procCats[p].hists2D[pt][22]->Fill((double)_predJets.size(), jetsize);
						_procCats[p].hists2D[pt][23]->Fill(_predJets[j].GetNRecHits(),_predJets[j].GetNConstituents());
						_procCats[p].hists2D[pt][28]->Fill(_predJets[j].GetNConstituents(), _predJets[j].mass());
						_procCats[p].hists2D[pt][21]->Fill(_predJets[j].e(), _predJets[j].mass());
						_procCats[p].hists2D[pt][35]->Fill(_predJets[j].e(), jetsize);
						_procCats[p].hists2D[pt][29]->Fill(_predJets[j].GetNConstituents(), _predJets[j].e());	
						_procCats[p].hists2D[pt][7]->Fill(_predJets[j].GetNConstituents(), jetsize);	
						_procCats[p].hists2D[pt][128]->Fill(_predJets[j].mass(), jetsize);	

						_procCats[p].hists1D[pt][74]->Fill(jet_mu.at(0,0));
						_procCats[p].hists1D[pt][75]->Fill(_predJets[j].phi()); //so phi is [0,2pi] but equivalent to jet_mu.at(1,0)
						_procCats[p].hists1D[pt][76]->Fill(jet_mu.at(2,0));
					
						
						_procCats[p].hists1D[pt][90]->Fill((int)rhs.size());
						for(int r = 0; r < rhs.size(); r++){
							_procCats[p].hists1D[pt][87]->Fill(rhs[r].E());
							
						}

					//cout << "pred jet #" << j << " phi1 " << _predJets[j].phi() << " phi " << jet_mu.at(1,0) << " phi std " << _predJets[j].phi_std() << endl;
						//get subcluster information
						vector<Jet> consts = _predJets[j].GetConstituents();
						for(int c = 0; c < (int)consts.size(); c++){
							Jet subcl = consts[c];
							Matrix subcl_cov = subcl.GetCovariance();
							_procCats[p].hists1D[pt][11]->Fill(sqrt(subcl_cov.at(0,0)));
							_procCats[p].hists1D[pt][12]->Fill(sqrt(subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][13]->Fill(sqrt(subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][14]->Fill(subcl_cov.at(0,1));
							_procCats[p].hists1D[pt][15]->Fill(subcl_cov.at(0,2));
							_procCats[p].hists1D[pt][16]->Fill(subcl_cov.at(1,2));
							//do for normalized values too
	
							//E_k = norms[k]/_gev;
							_procCats[p].hists1D[pt][2]->Fill(subcl.E());
							if(openAng != -1) _procCats[p].hists2D[pt][121]->Fill(openAng,subcl.m());
							_procCats[p].hists1D[pt][3]->Fill(subcl.eta());
							_procCats[p].hists1D[pt][4]->Fill(subcl.phi());
							_procCats[p].hists1D[pt][5]->Fill(subcl.time());
							
							_procCats[p].hists1D[pt][142]->Fill(subcl.m());
							_procCats[p].hists2D[pt][106]->Fill(subcl.E(), subcl.m());
							//effective # of rhs is sum of weights bc weights are each the responsibility of rechit n to this subcl
							vector<double> ws;
							subcl.GetWeights(ws);
							double effnRhs = 0;
							for(int n = 0; n < ws.size(); n++){
								effnRhs += ws[n]; 
							}
							_procCats[p].hists1D[pt][143]->Fill(effnRhs);
							_procCats[p].hists2D[pt][107]->Fill(subcl.E(), effnRhs);
							_procCats[p].hists2D[pt][108]->Fill(subcl.m(), effnRhs);
							if(p == 0 && pt == 0) cout << " subcl #" << c << " has " << effnRhs << " # of effective rechits of " << ws.size() << " total rhs and mass " << subcl.m() << " and energy " << subcl.E() << " px " << subcl.px() << " py " << subcl.py() << " pz " << subcl.pz() << endl;
							
							//do dr calcs bw all subclusters
							for(int cc = c+1; cc < consts.size(); cc++){
								double dr = dR(consts[c].eta(), consts[c].phi(), consts[cc].eta(), consts[cc].phi());
								_procCats[p].hists1D[pt][85]->Fill(dr);
							}
							//dr to jet center
							double jet_dr = dR(consts[c].eta(), consts[c].phi(), jet_mu.at(0,0), jet_mu.at(1,0));
							_procCats[p].hists2D[pt][112]->Fill(consts[c].E(), jet_dr);
							_procCats[p].hists2D[pt][113]->Fill(effnRhs, jet_dr);
						}


						if(p == 0 && pt == 0) cout << "pred jet #" << j << " matched to reco AK4 jet #" << recoMatchIdxs[j] << " gen AK4 jet #" << genAK4MatchIdxs[j] << " and gen AK15 jet #" << genAK15MatchIdxs[j] << " and gen top #" << genTopMatchIdxs[j] << " and gen W #" << genWMatchIdxs[j] << endl;
						//exclusively dr-matched to reco jet
						if(recoMatchIdxs[j] != -1){
							if(_recoAK4jets[recoMatchIdxs[j]].GetNConstituents() == 0) cout << "reco jet #" << recoMatchIdxs[j] << " matched to bhc jet #" << j << " has " << _recoAK4jets[recoMatchIdxs[j]].GetNConstituents() << " # subclusters" << endl;
							_procCats[p].hists2D[pt][31]->Fill(_recoAK4jets[recoMatchIdxs[j]].GetNConstituents(), _predJets[j].GetNConstituents());
						}
						//do reco AK8 matching hists
						if(recoAK8MatchIdxs[j] != -1){
							int matchidx = recoAK8MatchIdxs[j];
							_procCats[p].hists2D[pt][132]->Fill(_recoAK8jets[matchidx].m(), _predJets[j].m());
						}
	
						//do gen AK4 matching hists
						if(genAK4MatchIdxs[j] != -1){
							int genAK4jetidx = _genAK4jets[genAK4MatchIdxs[j]].GetUserIdx();
							_procCats[p].hists2D[pt][42]->Fill(_predJets[j].pt(),_genAK4jets[genAK4MatchIdxs[j]].pt());
							_procCats[p].hists2D[pt][43]->Fill(_predJets[j].E(),_genAK4jets[genAK4MatchIdxs[j]].E());
							_procCats[p].hists2D[pt][44]->Fill(_predJets[j].m(),_genAK4jets[genAK4MatchIdxs[j]].m());
							_procCats[p].hists2D[pt][45]->Fill(_predJets[j].eta(),_genAK4jets[genAK4MatchIdxs[j]].eta());
							_procCats[p].hists2D[pt][46]->Fill(_predJets[j].phi(),_genAK4jets[genAK4MatchIdxs[j]].phi());
							_procCats[p].hists2D[pt][47]->Fill(_predJets[j].time(),_genAK4jets[genAK4MatchIdxs[j]].time());
							_procCats[p].hists2D[pt][48]->Fill(_predJets[j].GetNConstituents(), _base->AK4Jet_genNConstituents->at(genAK4jetidx));
							
						}			
						//do gen AK15 matching hists
						if(genAK15MatchIdxs[j] != -1){
							int genAK15jetidx = _genAK15jets[genAK15MatchIdxs[j]].GetUserIdx();
							_procCats[p].hists2D[pt][49]->Fill(_predJets[j].pt(),_genAK15jets[genAK15MatchIdxs[j]].pt());
							_procCats[p].hists2D[pt][50]->Fill(_predJets[j].E(),_genAK15jets[genAK15MatchIdxs[j]].E());
							_procCats[p].hists2D[pt][51]->Fill(_predJets[j].m(),_genAK15jets[genAK15MatchIdxs[j]].m());
							_procCats[p].hists2D[pt][52]->Fill(_predJets[j].eta(),_genAK15jets[genAK15MatchIdxs[j]].eta());
							_procCats[p].hists2D[pt][53]->Fill(_predJets[j].phi(),_genAK15jets[genAK15MatchIdxs[j]].phi());
							_procCats[p].hists2D[pt][54]->Fill(_predJets[j].time(),_genAK15jets[genAK15MatchIdxs[j]].time());
							_procCats[p].hists2D[pt][55]->Fill(_predJets[j].GetNConstituents(), _base->AK15Jet_genNConstituents->at(genAK15jetidx));

							//phi debugging
							double genphi = _genAK15jets[genAK15MatchIdxs[j]].phi_02pi();
							double bhcphi = _predJets[j].phi();
							double true_dphi = acos(cos(genphi - bhcphi));
							double dphi = fabs(genphi - bhcphi);
							if(true_dphi < 0.1 && dphi > 6.){
								cout << "BAD PHI MATCH HERE!!!!! for bhc jet # " << j << " with phi " << bhcphi << " and gen AK15 gen # " << genAK15MatchIdxs[j] << " with phi " << genphi << " evt " << _base->event << endl;
							}

						}			
						//do gen top matching hists
						if(genTopMatchIdxs[j] != -1){
							int gentopidx = genTopMatchIdxs[j];
							double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genTop[gentopidx].eta(), _genTop[gentopidx].phi());
							double eratio = _predJets[j].E()/_genTop[gentopidx].E();
							if(p == 0 && pt == 0) cout << "BHC jet #" << j << " with energy " << _predJets[j].E() << " matched to top " << gentopidx << " with eta " << _genTop[gentopidx].eta() << " phi " << _genTop[gentopidx].phi() << " and energy " << _genTop[gentopidx].E() << " matched with dr " << dr << " and Eratio " << eratio << endl;
							_procCats[p].hists2D[pt][56]->Fill(_predJets[j].pt(),_genTop[gentopidx].pt());
							_procCats[p].hists2D[pt][57]->Fill(_predJets[j].E(),_genTop[gentopidx].E());
							_procCats[p].hists2D[pt][58]->Fill(_predJets[j].m(),_genTop[gentopidx].m());
							_procCats[p].hists2D[pt][59]->Fill(_predJets[j].eta(),_genTop[gentopidx].eta());
							_procCats[p].hists2D[pt][60]->Fill(_predJets[j].phi(),_genTop[gentopidx].phi());
							
							_procCats[p].hists1D[pt][72]->Fill(dr);
							_procCats[p].hists1D[pt][73]->Fill(eratio);
							_procCats[p].hists1D[pt][153]->Fill(_predJets[j].GetNConstituents());
							//sort by energy
							consts = _predJets[j].GetConstituents();
							sort(consts.begin(), consts.end(), Esort_jet);
							for(int c = 0; c < consts.size(); c++){ 
								_procCats[p].hists1D[pt][154]->Fill(consts[c].m());
								_procCats[p].hists2D[pt][118]->Fill((int)consts.size(),consts[c].m());
							}
							//get invariant mass of leading two subclusters
							if(consts.size() > 1){
								double invmass = consts[0].invMass(consts[1]);
								_procCats[p].hists1D[pt][155]->Fill(invmass);
							
								//do gen matching of 2 leading subclusters to partons from W decays
								Jet leadcl = consts[0];
								Jet subleadcl = consts[1];
								vector<Jet> subcls = {leadcl, subleadcl};

							}
						}			
						//do gen W matching hists
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];
							double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							double eratio = _predJets[j].E()/_genW[genWidx].E();
							if(p == 0 && pt == 0) cout << "BHC jet #" << j << " with eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " and energy " << _predJets[j].E() << " matched to W " << genWidx << " with eta " << _genW[genWidx].eta() << " phi " << _genW[genWidx].phi() << " and energy " << _genW[genWidx].E() << " matched with dr " << dr << " and Eratio " << eratio << endl;
							//put loose req on W match
							//if(dr < 0.1 && (Eratio < 1.5 && Eratio > 0.5))
							//cout << "passed gen match reqs" << endl;
							_procCats[p].hists2D[pt][81]->Fill(_predJets[j].pt(),_genW[genWMatchIdxs[j]].pt());
							_procCats[p].hists2D[pt][82]->Fill(_predJets[j].E(),_genW[genWMatchIdxs[j]].E());
							_procCats[p].hists2D[pt][83]->Fill(_predJets[j].m(),_genW[genWMatchIdxs[j]].m());
							_procCats[p].hists2D[pt][84]->Fill(_predJets[j].eta(),_genW[genWMatchIdxs[j]].eta());
							_procCats[p].hists2D[pt][85]->Fill(_predJets[j].phi(),_genW[genWMatchIdxs[j]].phi());
							
							_procCats[p].hists1D[pt][138]->Fill(dr);
							_procCats[p].hists1D[pt][139]->Fill(eratio);
							consts = _predJets[j].GetConstituents();
							if(consts.size() > 1) _procCats[p].hists2D[pt][127]->Fill(_predJets[j].E()/_genW[genWidx].E(),consts.size());
							
							//finding partons/subclusters
							_procCats[p].hists1D[pt][140]->Fill(_predJets[j].GetNConstituents());
							//get gen partons from W decay 
							vector<int> genLeadMatchIdxs(2,-1);
							vector<Jet> Wpartons;
							int ggenWidx = _genW[genWidx].GetUserIdx();
							for(int g = 0; g < _genparts.size(); g++){
								int genidx = _genparts[g].GetUserIdx();
								if(_base->genpart_momIdx->at(genidx) != ggenWidx) continue;
								Wpartons.push_back(_genparts[g]);
							}
							if(Wpartons.size() != 2){
								cout << "Error: " << Wpartons.size() << " daughter particles found for W " << genWidx << " skipping hist filling" << endl;
								continue;
							}
							double gendR = dR(Wpartons[0].eta(), Wpartons[0].phi(), Wpartons[1].eta(), Wpartons[1].phi());
							if(p == 0 && pt == 0) cout << "gen dr " << gendR << endl;
							_procCats[p].hists2D[pt][115]->Fill(gendR, _predJets[j].GetNConstituents());
							_procCats[p].hists2D[pt][116]->Fill(gendR, jetsize);

							double avgPartE = (Wpartons[0].E() + Wpartons[1].E())/2.;
cout << "avgPart E " << avgPartE << endl;
					
							if(consts.size() == 1){
								_procCats[p].hists2D[pt][122]->Fill(gendR, jetsize);
								_procCats[p].hists2D[pt][124]->Fill(gendR, avgPartE);
							}
							if(consts.size() > 1){
								_procCats[p].hists2D[pt][123]->Fill(gendR, jetsize);
								_procCats[p].hists2D[pt][125]->Fill(gendR, avgPartE);

							}
							if(_predJets[j].m() > 100) _procCats[p].hists1D[pt][179]->Fill((int)consts.size());
							if(70 < _predJets[j].m() && _predJets[j].m() < 90) _procCats[p].hists1D[pt][199]->Fill((int)consts.size());
							if(_predJets[j].m() < 50) _procCats[p].hists1D[pt][198]->Fill((int)consts.size());
							//sort by energy
							sort(consts.begin(), consts.end(), Esort_jet);
							//get invariant mass of leading two subclusters
							if(consts.size() > 1){
								double invmass = consts[0].invMass(consts[1]);
								_procCats[p].hists1D[pt][147]->Fill(invmass);
								vector<int> subclGenMatchIdx(consts.size(), -1);	
								GenericMatchJet(consts,Wpartons,subclGenMatchIdx); //match subclusters to W partons
								for(int c = 0; c < consts.size(); c++){
									_procCats[p].hists2D[pt][126]->Fill(consts[c].E(),c);
									
									Matrix subcl_mu, subcl_cov;
									consts[c].GetClusterParams(subcl_mu, subcl_cov);
									double subclsize = CalcSize(subcl_cov);
									double subclsize_full = CalcSize(subcl_cov, true);
									int genmatchidx = subclGenMatchIdx[c];
									if(genmatchidx == -1){
										//not matched
										if(_predJets[j].m() > 100){
											_procCats[p].hists1D[pt][170]->Fill(consts[c].pt());
											_procCats[p].hists1D[pt][206]->Fill(consts[c].pt() / _predJets[j].pt());
											_procCats[p].hists1D[pt][172]->Fill(subclsize);
											_procCats[p].hists1D[pt][208]->Fill(subclsize / jetsize);
											
											_procCats[p].hists2D[pt][139]->Fill(consts[c].pt(), consts[c].t());
											_procCats[p].hists2D[pt][141]->Fill(subclsize, consts[c].t());
											_procCats[p].hists2D[pt][143]->Fill(subclsize / jetsize, consts[c].pt() / _predJets[j].pt());
											_procCats[p].hists2D[pt][145]->Fill(consts[c].m() / _predJets[j].m(), consts[c].pt() / _predJets[j].pt());
											//relative time variance vs relative pt
											_procCats[p].hists2D[pt][147]->Fill(subcl_cov.at(2,2) / jet_cov.at(2,2), consts[c].pt() / _predJets[j].pt());
											//relative time variance vs relative size
											_procCats[p].hists2D[pt][149]->Fill(subcl_cov.at(2,2) / jet_cov.at(2,2), subclsize / jetsize);
											//subcl mass / pt vs jet mass / pt
											_procCats[p].hists2D[pt][151]->Fill(consts[c].m() / consts[c].pt() , _predJets[j].m() / _predJets[j].pt());
											//full (space + time) subcl size / full jet size vs subcl pt / jet pt
											_procCats[p].hists2D[pt][153]->Fill(subclsize_full / jetsize_full, consts[c].pt() / _predJets[j].pt());
										}
				
									}
									else{
										//matched
										if(_predJets[j].m() > 100){
											_procCats[p].hists1D[pt][169]->Fill(consts[c].pt());
											_procCats[p].hists1D[pt][205]->Fill(consts[c].pt()/_predJets[j].pt());
											_procCats[p].hists1D[pt][171]->Fill(subclsize);
											_procCats[p].hists1D[pt][207]->Fill(subclsize / jetsize);

											_procCats[p].hists2D[pt][138]->Fill(consts[c].pt(), consts[c].t());
											_procCats[p].hists2D[pt][140]->Fill(subclsize, consts[c].t());
											_procCats[p].hists2D[pt][142]->Fill(subclsize / jetsize, consts[c].pt() / _predJets[j].pt());
											_procCats[p].hists2D[pt][144]->Fill(consts[c].m() / _predJets[j].m(), consts[c].pt() / _predJets[j].pt());
											//relative time variance vs relative pt
											_procCats[p].hists2D[pt][146]->Fill(subcl_cov.at(2,2) / jet_cov.at(2,2), consts[c].pt() / _predJets[j].pt());
											//relative time variance vs relative size
											_procCats[p].hists2D[pt][148]->Fill(subcl_cov.at(2,2) / jet_cov.at(2,2), subclsize / jetsize);
											//subcl mass / pt vs jet mass / pt
											_procCats[p].hists2D[pt][150]->Fill(consts[c].m() / consts[c].pt() , _predJets[j].m() / _predJets[j].pt());
											//full (space + time) subcl size / full jet size vs subcl pt / jet pt
											_procCats[p].hists2D[pt][152]->Fill(subclsize_full / jetsize_full, consts[c].pt() / _predJets[j].pt());
										}
									}
								}	

								//do gen matching of 2 leading subclusters to partons from W decays
								Jet leadcl = consts[0];
								Jet subleadcl = consts[1];
								vector<Jet> subcls = {leadcl, subleadcl};

								//check that there are two daughters with light quark ids
								//only matching lead subclusters to W daughter partons for now
								GenericMatchJet(subcls,Wpartons,genLeadMatchIdxs); //match subclusters to W partons
								for(int c = 0; c < subcls.size(); c++){
									int genmatchidx = genLeadMatchIdxs[c];
									if(genmatchidx == -1) continue;
									double gen_clDr = dR(subcls[c].eta(), subcls[c].phi(), Wpartons[genmatchidx].eta(), Wpartons[genmatchidx].phi());
									double genEr = subcls[c].E() / Wpartons[genmatchidx].E();
									_procCats[p].hists1D[pt][151]->Fill(gen_clDr);
									_procCats[p].hists1D[pt][152]->Fill(genEr);
									_procCats[p].hists2D[pt][117]->Fill(genEr,_predJets[j].GetNConstituents());

									_procCats[p].hists1D[pt][160]->Fill(subcls[c].eta());
									_procCats[p].hists1D[pt][161]->Fill(subcls[c].phi());
									_procCats[p].hists1D[pt][162]->Fill(subcls[c].time());
									Matrix subcl_cov = subcls[c].GetCovariance();
									_procCats[p].hists1D[pt][163]->Fill(sqrt(subcl_cov.at(0,0)));
									_procCats[p].hists1D[pt][164]->Fill(sqrt(subcl_cov.at(1,1)));
									_procCats[p].hists1D[pt][165]->Fill(sqrt(subcl_cov.at(2,2)));
									_procCats[p].hists1D[pt][166]->Fill(subcl_cov.at(0,1));
									_procCats[p].hists1D[pt][167]->Fill(subcl_cov.at(0,2));
									_procCats[p].hists1D[pt][168]->Fill(subcl_cov.at(2,1));

								}
							}
							for(int c = 0; c < (int)consts.size(); c++){
								_procCats[p].hists1D[pt][141]->Fill(consts[c].E());
								_procCats[p].hists1D[pt][146]->Fill(consts[c].m());
							//cout << "consts #" << c << " mass " << consts[c].m() << endl;
								_procCats[p].hists2D[pt][114]->Fill(_predJets[j].GetNConstituents(),consts[c].m());
							}
							
						}			
						//do gen q matching hists
						if(genqMatchIdxs[j] != -1){
							int genqidx = genqMatchIdxs[j];
							double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genq[genqidx].eta(), _genq[genqidx].phi());
							double eratio = _predJets[j].E()/_genq[genqidx].E();
							if(p == 0 && pt == 0) cout << "BHC jet #" << j << " with eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " and energy " << _predJets[j].E() << " matched to q " << genqidx << " with eta " << _genq[genqidx].eta() << " phi " << _genq[genqidx].phi() << " and energy " << _genq[genqidx].E() << " matched with dr " << dr << " and Eratio " << eratio << endl;
							//put loose req on q match
							//if(dr < 0.1 && (Eratio < 1.5 && Eratio > 0.5))
							//cout << "passed gen match reqs" << endl;
							_procCats[p].hists1D[pt][156]->Fill(dr);
							_procCats[p].hists1D[pt][157]->Fill(eratio);
							//finding partons/subclusters
							_procCats[p].hists1D[pt][158]->Fill(_predJets[j].GetNConstituents());

							consts = _predJets[j].GetConstituents();
							//sort by energy
							sort(consts.begin(), consts.end(), Esort_jet);
							for(int c = 0; c < (int)consts.size(); c++){
								_procCats[p].hists1D[pt][159]->Fill(consts[c].m());
							}
							vector<int> subclGenMatchIdx(consts.size(), -1);	
							GenericMatchJet(consts,_genq,subclGenMatchIdx); //match subclusters to hard light quarks
							for(int c = 0; c < consts.size(); c++){
								int genmatchidx = subclGenMatchIdx[c];
								//1+ subcl jets only
								Matrix subcl_cov = consts[c].GetCovariance();
								double subclsize = CalcSize(subcl_cov);
								if((int)consts.size() > 1){
									//not matched
									if(genmatchidx == -1){
										_procCats[p].hists1D[pt][195]->Fill(consts[c].pt());
										_procCats[p].hists1D[pt][197]->Fill(subclsize);
										continue;
									}
									//matched
									_procCats[p].hists1D[pt][194]->Fill(consts[c].pt());
									_procCats[p].hists1D[pt][196]->Fill(subclsize);
								}
							}	
							//do gen matching of leading subcluster to hard parton
							Jet leadcl = consts[0];
							vector<Jet> subcls = {leadcl};
							vector<int> genLeadMatchIdxs(1,-1);
							//TODO: need to match subcl to THIS (ie genGluonMatchIdxs[j]) gen gluon
							GenericMatchJet(subcls,_genq,genLeadMatchIdxs); //match subclusters to quarks
							for(int c = 0; c < subcls.size(); c++){
								int genmatchidx = genLeadMatchIdxs[c];
								if(genmatchidx == -1) continue;
								double gen_clDr = dR(subcls[c].eta(), subcls[c].phi(), _genq[genmatchidx].eta(), _genq[genmatchidx].phi());
								double genEr = subcls[c].E() / _genq[genmatchidx].E();
								_procCats[p].hists1D[pt][201]->Fill(gen_clDr);
								_procCats[p].hists1D[pt][202]->Fill(genEr);

							}

						}			
						//do gen gluon matching hists
						if(genGluonMatchIdxs[j] != -1){
							int genGluonidx = genGluonMatchIdxs[j];
							double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genglu[genGluonidx].eta(), _genglu[genGluonidx].phi());
							double eratio = _predJets[j].E()/_genglu[genGluonidx].E();
							if(p == 0 && pt == 0) cout << "BHC jet #" << j << " with eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " and energy " << _predJets[j].E() << " matched to gluon " << genGluonidx << " with eta " << _genglu[genGluonidx].eta() << " phi " << _genglu[genGluonidx].phi() << " and energy " << _genglu[genGluonidx].E() << " matched with dr " << dr << " and Eratio " << eratio << endl;
							_procCats[p].hists1D[pt][177]->Fill(dr);
							_procCats[p].hists1D[pt][178]->Fill(eratio);
							_procCats[p].hists1D[pt][200]->Fill(_predJets[j].GetNConstituents());
							//do gen matching of leading subcluster to hard parton
							Jet leadcl = consts[0];
							vector<Jet> subcls = {leadcl};
							vector<int> genLeadMatchIdxs(1,-1);
							//TODO: need to match subcl to THIS (ie genGluonMatchIdxs[j]) gen gluon
							GenericMatchJet(subcls,_genglu,genLeadMatchIdxs); //match subclusters to gluons
							for(int c = 0; c < subcls.size(); c++){
								int genmatchidx = genLeadMatchIdxs[c];
								if(genmatchidx == -1) continue;
								double gen_clDr = dR(subcls[c].eta(), subcls[c].phi(), _genglu[genmatchidx].eta(), _genglu[genmatchidx].phi());
								double genEr = subcls[c].E() / _genglu[genmatchidx].E();
								_procCats[p].hists1D[pt][203]->Fill(gen_clDr);
								_procCats[p].hists1D[pt][204]->Fill(genEr);

							}
						}

					}
				}
			}
		}


		//fill hists for gen jets and gen particles
		//ONLY GEN PARTICLES CONSIDERED - tops, b's from tops, quarks from W's, direct daughters of b's
		void FillGenParticleHists(){
			for(int p = 0; p < _procCats.size(); p++){
				//fill gen particle hists - needs GetGenParticles() method in JetSimProducer
				int nGenParts = 0; //only count hadron-izable gen particles (ie quarks)
				vector<int> qids = {1,2,3,4,5,6};
				int id;
				_procCats[p].hists1D[0][53]->Fill((double)nGenParts);
				_procCats[p].hists1D[0][68]->Fill((double)_genAK4jets.size() - (double)nGenParts);
				_procCats[p].hists1D[0][99]->Fill((double)_genAK15jets.size() - (double)nGenParts);
				_procCats[p].hists1D[0][71]->Fill((double)_predJets.size() - (double)nGenParts);
				for(int g = 0; g < _genparts.size(); g++){
					_procCats[p].hists1D[0][54]->Fill(_genparts[g].eta());
					_procCats[p].hists1D[0][55]->Fill(_genparts[g].phi());
					_procCats[p].hists1D[0][56]->Fill(_genparts[g].time());
					_procCats[p].hists1D[0][57]->Fill(_genparts[g].pt());
					_procCats[p].hists1D[0][58]->Fill(_genparts[g].mass());
					_procCats[p].hists1D[0][59]->Fill(_genparts[g].E());
					id = _base->genpart_id->at(_genparts[g].GetUserIdx());
					if(find(qids.begin(), qids.end(), id) != qids.end()) nGenParts++;
					
				}
			}
		}
	
		void FillGenJetHists(){
			for(int p = 0; p < _procCats.size(); p++){
				_procCats[p].hists1D[0][66]->Fill((double)_genAK4jets.size());
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _genAK4jets.size(); j++){
					if(_genAK4jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_genAK4jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][66]->Fill(njets_lead);
				_procCats[p].hists1D[2][66]->Fill(njets_notlead);
				//gen match jets to particles
				vector<int> genTopMatchIdxs;
				GenericMatchJet(_genAK4jets,_genparts,genTopMatchIdxs, 6); //match gen AK4 jets to gen tops
				vector<int> genWMatchIdxs;
				GenericMatchJet(_genAK4jets,_genparts,genWMatchIdxs, 24); //match gen AK4 jets to gen Ws
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genAK4jets.size(); j++){
					if(p == 0) cout << "gen AK4 jet #" << j << " phi " << _genAK4jets[j].phi() << " eta " << _genAK4jets[j].eta() << " energy " << _genAK4jets[j].E() <<  " mass " << _genAK4jets[j].mass() << " pt " << _genAK4jets[j].pt() << endl;
					_procCats[p].hists1D[0][60]->Fill(_genAK4jets[j].eta());
					_procCats[p].hists1D[0][61]->Fill(_genAK4jets[j].phi());
					_procCats[p].hists1D[0][62]->Fill(_genAK4jets[j].time());
					_procCats[p].hists1D[0][63]->Fill(_genAK4jets[j].pt());
					_procCats[p].hists1D[0][64]->Fill(_genAK4jets[j].mass());
					_procCats[p].hists1D[0][65]->Fill(_genAK4jets[j].E());
					_procCats[p].hists1D[0][67]->Fill(_base->AK4Jet_genNConstituents->at(_genAK4jets[j].GetUserIdx()));
					//dr bw gen jet and best exclusive gen top match
					int genTopMatch = genTopMatchIdxs[j];
					if(genTopMatch != -1){
					if(p == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << " with id " << _base->genpart_id->at(_genparts[genTopMatch].GetUserIdx()) << " and mass " << _genparts[genTopMatch].mass() << " and energy ratio " << _genAK4jets[j].E() / _genparts[genTopMatch].E() << endl;
						double gendR = dR(_genAK4jets[j].eta(), _genAK4jets[j].phi(), _genparts[genTopMatch].eta(), _genparts[genTopMatch].phi());
						_procCats[p].hists1D[0][69]->Fill(gendR);
						_procCats[p].hists1D[0][70]->Fill(_genAK4jets[j].E()/_genparts[genTopMatch].E());

						_procCats[p].hists2D[0][71]->Fill(_genAK4jets[j].pt(),_genparts[genTopMatch].pt());
						_procCats[p].hists2D[0][72]->Fill(_genAK4jets[j].E(),_genparts[genTopMatch].E());
						_procCats[p].hists2D[0][73]->Fill(_genAK4jets[j].m(),_genparts[genTopMatch].m());
						_procCats[p].hists2D[0][74]->Fill(_genAK4jets[j].eta(),_genparts[genTopMatch].eta());
						_procCats[p].hists2D[0][75]->Fill(_genAK4jets[j].phi(),_genparts[genTopMatch].phi());
					}	
					//dr bw gen jet and best exclusive gen W match
					int genWMatch = genWMatchIdxs[j];
					if(genWMatch != -1){
					if(p == 0) cout << " matched to gen W #" << genWMatchIdxs[j] << " with id " << _base->genpart_id->at(_genparts[genWMatch].GetUserIdx()) << " and mass " << _genparts[genWMatch].mass() << " and energy ratio " << _genAK4jets[j].E()/_genparts[genWMatch].E() << endl;
						double gendR = dR(_genAK4jets[j].eta(), _genAK4jets[j].phi(), _genparts[genWMatch].eta(), _genparts[genWMatch].phi());
						_procCats[p].hists1D[0][134]->Fill(gendR);
						_procCats[p].hists1D[0][135]->Fill(_genAK4jets[j].E()/_genparts[genWMatch].E());

						_procCats[p].hists2D[0][96]->Fill(_genAK4jets[j].pt(),_genparts[genWMatch].pt());
						_procCats[p].hists2D[0][97]->Fill(_genAK4jets[j].E(),_genparts[genWMatch].E());
						_procCats[p].hists2D[0][98]->Fill(_genAK4jets[j].m(),_genparts[genWMatch].m());
						_procCats[p].hists2D[0][99]->Fill(_genAK4jets[j].eta(),_genparts[genWMatch].eta());
						_procCats[p].hists2D[0][100]->Fill(_genAK4jets[j].phi(),_genparts[genWMatch].phi());
					}	

				}
				_procCats[p].hists1D[0][124]->Fill((double)_genAK8jets.size());
				njets_lead = 0;
				njets_notlead = 0;
				for(int j = 0; j < _genAK8jets.size(); j++){
					if(_genAK8jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_genAK8jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][124]->Fill(njets_lead);
				_procCats[p].hists1D[2][124]->Fill(njets_notlead);
				//cout << "gen matching gen jets to particles - start" << endl;
				GenericMatchJet(_genAK8jets,_genparts,genTopMatchIdxs,6);
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genAK8jets.size(); j++){
					if(p == 0) cout << "gen AK8 jet #" << j << " phi " << _genAK8jets[j].phi() << " eta " << _genAK8jets[j].eta() << " energy " << _genAK8jets[j].E() <<  " mass " << _genAK8jets[j].mass() << " pt " << _genAK8jets[j].pt() << endl;
					_procCats[p].hists1D[0][125]->Fill(_genAK8jets[j].eta());
					_procCats[p].hists1D[0][126]->Fill(_genAK8jets[j].phi());
					_procCats[p].hists1D[0][127]->Fill(_genAK8jets[j].time());
					_procCats[p].hists1D[0][128]->Fill(_genAK8jets[j].pt());
					_procCats[p].hists1D[0][129]->Fill(_genAK8jets[j].mass());
					_procCats[p].hists1D[0][130]->Fill(_genAK8jets[j].E());
					_procCats[p].hists1D[0][131]->Fill(_base->AK8Jet_genNConstituents->at(_genAK8jets[j].GetUserIdx()));
					//dr bw gen jet and best exclusive gen top match
					//if(p == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << endl;
					int genTopMatch = genTopMatchIdxs[j];
					if(genTopMatch != -1){
						double gendR = dR(_genAK8jets[j].eta(), _genAK8jets[j].phi(), _genparts[genTopMatch].eta(), _genparts[genTopMatch].phi());
						_procCats[p].hists1D[0][132]->Fill(gendR);
						_procCats[p].hists1D[0][133]->Fill(_genAK8jets[j].E()/_genparts[genTopMatch].E());
						
					}	

				}
				
				_procCats[p].hists1D[0][97]->Fill((double)_genAK15jets.size());
				njets_lead = 0;
				njets_notlead = 0;
				for(int j = 0; j < _genAK15jets.size(); j++){
					if(_genAK15jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_genAK15jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][97]->Fill(njets_lead);
				_procCats[p].hists1D[2][97]->Fill(njets_notlead);
				//cout << "gen matching gen jets to particles - start" << endl;
				GenericMatchJet(_genAK15jets, _genparts, genTopMatchIdxs,6);
				GenericMatchJet(_genAK15jets,_genparts,genWMatchIdxs, 24); //match gen AK4 jets to gen Ws
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genAK15jets.size(); j++){
					if(p == 0) cout << "gen AK15 jet #" << j << " phi " << _genAK15jets[j].phi() << " eta " << _genAK15jets[j].eta() << " energy " << _genAK15jets[j].E() <<  " mass " << _genAK15jets[j].mass() << " pt " << _genAK15jets[j].pt() << endl;
					_procCats[p].hists1D[0][91]->Fill(_genAK15jets[j].eta());
					_procCats[p].hists1D[0][92]->Fill(_genAK15jets[j].phi());
					_procCats[p].hists1D[0][93]->Fill(_genAK15jets[j].time());
					_procCats[p].hists1D[0][94]->Fill(_genAK15jets[j].pt());
					_procCats[p].hists1D[0][95]->Fill(_genAK15jets[j].mass());
					_procCats[p].hists1D[0][96]->Fill(_genAK15jets[j].E());
					_procCats[p].hists1D[0][98]->Fill(_base->AK15Jet_genNConstituents->at(_genAK15jets[j].GetUserIdx()));
					//dr bw gen jet and best exclusive gen top match
					int genTopMatch = genTopMatchIdxs[j];
					//if(p == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << endl;
					if(genTopMatch != -1){
						double gendR = dR(_genAK15jets[j].eta(), _genAK15jets[j].phi(), _genparts[genTopMatch].eta(), _genparts[genTopMatch].phi());
						_procCats[p].hists1D[0][100]->Fill(gendR);
						_procCats[p].hists1D[0][101]->Fill(_genAK15jets[j].E()/_genparts[genTopMatch].E());
						
						_procCats[p].hists2D[0][76]->Fill(_genAK15jets[j].pt(),_genparts[genTopMatch].pt());
						_procCats[p].hists2D[0][77]->Fill(_genAK15jets[j].E(),_genparts[genTopMatch].E());
						_procCats[p].hists2D[0][78]->Fill(_genAK15jets[j].m(),_genparts[genTopMatch].m());
						_procCats[p].hists2D[0][79]->Fill(_genAK15jets[j].eta(),_genparts[genTopMatch].eta());
						_procCats[p].hists2D[0][80]->Fill(_genAK15jets[j].phi(),_genparts[genTopMatch].phi());
						
					}	
					
					//dr bw gen jet and best exclusive gen W match
					int genWMatch = genWMatchIdxs[j];
					if(p == 0) cout << " matched to gen W #" << genWMatchIdxs[j] << endl;
					if(genWMatch != -1){
						double gendR = dR(_genAK15jets[j].eta(), _genAK15jets[j].phi(), _genparts[genWMatch].eta(), _genparts[genWMatch].phi());
						_procCats[p].hists1D[0][136]->Fill(gendR);
						_procCats[p].hists1D[0][137]->Fill(_genAK15jets[j].E()/_genparts[genWMatch].E());
						
						_procCats[p].hists2D[0][101]->Fill(_genAK15jets[j].pt(),_genparts[genWMatch].pt());
						_procCats[p].hists2D[0][102]->Fill(_genAK15jets[j].E(),_genparts[genWMatch].E());
						_procCats[p].hists2D[0][103]->Fill(_genAK15jets[j].m(),_genparts[genWMatch].m());
						_procCats[p].hists2D[0][104]->Fill(_genAK15jets[j].eta(),_genparts[genWMatch].eta());
						_procCats[p].hists2D[0][105]->Fill(_genAK15jets[j].phi(),_genparts[genWMatch].phi());
						
					}	


				}
				

			}
		}


		void FillGenHists(){
			FillGenParticleHists();
			FillGenJetHists();
		}


		void FillRecoJetHists(){
			FillRecoAK4JetHists();
			FillRecoAK8JetHists();
			FillRecoAK15JetHists();
		}

		//AK4 reco jets only
		void FillRecoAK4JetHists(){
			int njets, tmass_idx, widx, genidx, id;
			double wmass, topmass, dr, dr_pair, jetsize;
			pair<int, int> wmass_idxs;
			Jet w, top, genW;
			//do gen matching
			vector<int> genMatchIdxs_jet; //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genAK4jets,genMatchIdxs_jet);
			vector<int> genTopMatchIdxs; //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genparts,genTopMatchIdxs, 6);
			vector<int> genWMatchIdxs; //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genparts,genWMatchIdxs, 24);
			vector<int> genqMatchIdxs; //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genq,genqMatchIdxs);
			vector<int> genGluonMatchIdxs; //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genglu,genGluonMatchIdxs);
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				njets = _recoAK4jets.size();
				_procCats[p].hists1D[0][17]->Fill(njets);
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _recoAK4jets.size(); j++){
					if(_recoAK4jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_recoAK4jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][17]->Fill(njets_lead);
				_procCats[p].hists1D[2][17]->Fill(njets_notlead);
				for(int j = 0; j < _recoAK4jets.size(); j++){
					for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
						//define pt bins
						//pt == 1 -> [_pt_thresh,inf)
						if(pt == 1 && _recoAK4jets[j].pt() < _pt_thresh) continue;
						//pt == 2 -> [0,_pt_thresh)
						if(pt == 2 && _recoAK4jets[j].pt() >= _pt_thresh) continue;
	
						Matrix jet_cov = _recoAK4jets[j].GetCovariance();
						jetsize = CalcSize(jet_cov);
						
						//define jetsize bins
						//pt == 1 -> [0,0.2) (AK4 level)
						//if(pt == 1 && jetsize >= 0.2) continue;
						////pt == 2 -> [0.2,inf) (AK15 level)
						//if(pt == 2 && jetsize < 0.2) continue;
						
						double rot = Rotundity(jet_cov);
						_procCats[p].hists1D[pt][88]->Fill(rot);
						if(pt == 0 && p == 0) cout << "reco AK4 jet #" << j << " phi " << _recoAK4jets[j].phi() << " eta " << _recoAK4jets[j].eta() << " energy " << _recoAK4jets[j].E() <<  " mass " << _recoAK4jets[j].mass() << " nConstituents " << _recoAK4jets[j].GetNConstituents() << " nRhs " << _recoAK4jets[j].GetNRecHits() << " pt " << _recoAK4jets[j].pt() << " jetsize " << jetsize << " px " << _recoAK4jets[j].px() << " py " << _recoAK4jets[j].py() << " pz " << _recoAK4jets[j].pz() << " mass from rhs " << _recoAK4jets[j].mass_rhs() << endl;
						_procCats[p].hists1D[pt][18]->Fill(jetsize);
						_procCats[p].hists1D[pt][78]->Fill(sqrt(jet_cov.at(0,0)));
						_procCats[p].hists1D[pt][79]->Fill(sqrt(jet_cov.at(1,1)));
						_procCats[p].hists1D[pt][80]->Fill(sqrt(jet_cov.at(2,2)));
						_procCats[p].hists1D[pt][19]->Fill(_recoAK4jets[j].e());
						_procCats[p].hists1D[pt][20]->Fill(_recoAK4jets[j].pt());
						_procCats[p].hists1D[pt][21]->Fill(_recoAK4jets[j].mass());
						_procCats[p].hists1D[pt][81]->Fill(_recoAK4jets[j].time());
						_procCats[p].hists1D[pt][82]->Fill(_recoAK4jets[j].eta());
						_procCats[p].hists1D[pt][83]->Fill(_recoAK4jets[j].phi());
						_procCats[p].hists2D[pt][3]->Fill(_recoAK4jets[j].mass(), _recoAK4jets[j].pt());
						_procCats[p].hists2D[pt][4]->Fill(_recoAK4jets[j].mass(), jetsize);
						_procCats[p].hists2D[pt][21]->Fill((double)_recoAK4jets.size(), jetsize);
						
						//fill subcluster hists
						_procCats[p].hists1D[pt][30]->Fill(_recoAK4jets[j].GetNConstituents());
						_procCats[p].hists2D[pt][8]->Fill(_recoAK4jets[j].GetNRecHits(),_recoAK4jets[j].GetNConstituents());
						_procCats[p].hists2D[pt][25]->Fill(_recoAK4jets[j].GetNConstituents(), _recoAK4jets[j].mass());
						_procCats[p].hists2D[pt][26]->Fill(_recoAK4jets[j].GetNConstituents(), _recoAK4jets[j].e());
						_procCats[p].hists2D[pt][32]->Fill(_recoAK4jets[j].e(), _recoAK4jets[j].mass());
						_procCats[p].hists2D[pt][22]->Fill(_recoAK4jets[j].GetNConstituents(), jetsize);
						_procCats[p].hists2D[pt][35]->Fill(_recoAK4jets[j].e(), jetsize);
						_procCats[p].hists2D[pt][33]->Fill(_recoAK4jets[j].GetNConstituents(), jetsize);	
						if(_recoAK4jets[j].GetNConstituents() == 0) cout << _recoAK4jets[j].GetNConstituents() << " n subcl " << _recoAK4jets[j].GetNRecHits() << " n rhs" << endl;
						//fill subcluster hists
						vector<Jet> consts = _recoAK4jets[j].GetConstituents();
						double rnk = 0;
						for(int c = 0; c < (int)consts.size(); c++){
							Jet subcl = consts[c];
							_procCats[p].hists1D[pt][31]->Fill(subcl.E());
							_procCats[p].hists1D[pt][32]->Fill(subcl.time());
							_procCats[p].hists1D[pt][33]->Fill(subcl.eta());
							_procCats[p].hists1D[pt][34]->Fill(subcl.phi());
							
							Matrix subcl_cov = subcl.GetCovariance();
							_procCats[p].hists1D[pt][35]->Fill(sqrt(subcl_cov.at(0,0)));
							_procCats[p].hists1D[pt][36]->Fill(sqrt(subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][37]->Fill(sqrt(subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][38]->Fill(subcl_cov.at(0,1));
							_procCats[p].hists1D[pt][39]->Fill(subcl_cov.at(0,2));
							_procCats[p].hists1D[pt][40]->Fill(subcl_cov.at(1,2));
							_procCats[p].hists1D[pt][50]->Fill(subcl_cov.at(0,1)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][51]->Fill(subcl_cov.at(0,2)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][52]->Fill(subcl_cov.at(1,2)/sqrt(subcl_cov.at(1,1)*subcl_cov.at(2,2)));
							_procCats[p].hists2D[pt][19]->Fill(subcl_cov.at(0,1),subcl_cov.at(0,2));
							_procCats[p].hists2D[pt][20]->Fill(subcl_cov.at(0,1)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(1,1)),subcl_cov.at(0,2)/sqrt(subcl_cov.at(0,0)*subcl_cov.at(2,2)));
					
							_procCats[p].hists1D[pt][144]->Fill(subcl.m());
							_procCats[p].hists2D[pt][109]->Fill(subcl.E(), subcl.m());
							//effective # of rhs is sum of weights bc weights are each the responsibility of rechit n to this subcl
							vector<double> ws;
							subcl.GetWeights(ws);
							double effnRhs = 0;
							for(int n = 0; n < ws.size(); n++){
								effnRhs += ws[n]; 
							}
							if(p == 0 && pt == 0) cout << " subcl #" << c << " has " << effnRhs << " # of effective rechits of " << ws.size() << " total rhs and mass " << subcl.m() << " and energy " << subcl.E() << " px " << subcl.px() << " py " << subcl.py() << " pz " << subcl.pz() << endl;
							_procCats[p].hists1D[pt][145]->Fill(effnRhs);
							_procCats[p].hists2D[pt][110]->Fill(subcl.E(), effnRhs);
							_procCats[p].hists2D[pt][111]->Fill(subcl.m(), effnRhs);

	
							//do dr calcs bw all subclusters
							for(int cc = c+1; cc < consts.size(); cc++){
								double dr = dR(consts[c].eta(), consts[c].phi(), consts[cc].eta(), consts[cc].phi());
								_procCats[p].hists1D[pt][86]->Fill(dr);
							}
				

						}
						vector<Jet> rhs;
						_recoAK4jets[j].GetJets(rhs);
						for(int r = 0; r < rhs.size(); r++){
							_procCats[p].hists1D[pt][42]->Fill(rhs[r].E());
							_procCats[p].hists1D[pt][49]->Fill(rhs[r].t());
						}
						vector<pair<double,double>> geoEavg_diffT;
						vector<JetPoint> rhs_pt = _recoAK4jets[j].GetJetPoints();
						CalcRhTimeDiff(rhs_pt,geoEavg_diffT);
						for(int r = 0; r < geoEavg_diffT.size(); r++){
						      _procCats[p].hists2D[pt][18]->Fill(geoEavg_diffT[r].first, geoEavg_diffT[r].second);
						}
						_procCats[p].hists1D[pt][41]->Fill(_recoAK4jets[j].GetNRecHits());
						//if no gen match, skip
						if(p == 0 && pt == 0) cout << "reco AK4 jet #" << j << " gen top idx " << genTopMatchIdxs[j] << " gen gluon idx " << genGluonMatchIdxs[j] << " gen q idx " << genqMatchIdxs[j] << " gen w idx " << genWMatchIdxs[j] << endl; 
						if(genMatchIdxs_jet[j] != -1){
							genidx = genMatchIdxs_jet[j];	
							//if(p == 0 && pt == 0) cout << " matched to gen AK4 jet #" << genMatchIdxs_jet[j] << endl;
							int genjetidx = _genAK4jets[genMatchIdxs_jet[j]].GetUserIdx();
							_procCats[p].hists2D[pt][9]->Fill(_base->AK4Jet_genNConstituents->at(genjetidx),_recoAK4jets[j].GetNConstituents());
							double genpt, pz, jetpt, jetpz, ratio_p;
							jetpt = _genAK4jets[genidx].pt();//_base->Jet_genPt->at(genjetidx);
							jetpz = _genAK4jets[genidx].pz();//_base->Jet_genPz->at(genjetidx);

							_procCats[p].hists1D[pt][43]->Fill(sqrt(jetpt*jetpt + jetpz*jetpz));	
							_procCats[p].hists1D[pt][44]->Fill(jetpt);	
							_procCats[p].hists2D[pt][16]->Fill(jetpt, _recoAK4jets[j].GetNConstituents());
							_procCats[p].hists2D[pt][24]->Fill(_recoAK4jets[j].GetNConstituents(), _base->AK4Jet_genNConstituents->at(genjetidx));
				

							_procCats[p].hists2D[pt][36]->Fill(_recoAK4jets[j].pt(),_genAK4jets[genidx].pt());
							_procCats[p].hists2D[pt][37]->Fill(_recoAK4jets[j].E(),_genAK4jets[genidx].E());
							_procCats[p].hists2D[pt][38]->Fill(_recoAK4jets[j].m(),_genAK4jets[genidx].m());
							_procCats[p].hists2D[pt][39]->Fill(_recoAK4jets[j].eta(),_genAK4jets[genidx].eta());
							_procCats[p].hists2D[pt][40]->Fill(_recoAK4jets[j].phi(),_genAK4jets[genidx].phi());
							_procCats[p].hists2D[pt][41]->Fill(_recoAK4jets[j].time(),_genAK4jets[genidx].time());
							if(_sel != QCDdijets){
								int genpartidx = -1;
								int ngenparts_ptge5 = 0;
								for(int g = 0; g < _base->AK4Jet_genNConstituents->at(genjetidx); g++){
									genpartidx = _base->AK4Jet_genConstituentIdxs->at(genjetidx).at(g);
									if(genpartidx >= _base->genpart_pt->size()) break;
									genpt = _base->genpart_pt->at(genpartidx);
									pz = _base->genpart_pz->at(genpartidx);
									ratio_p = (sqrt(genpt*genpt + pz*pz))/(sqrt(jetpt*jetpt + jetpz*jetpz));	
									_procCats[p].hists1D[pt][45]->Fill(genpt);	
									_procCats[p].hists1D[pt][46]->Fill(sqrt(genpt*genpt + pz*pz));
									_procCats[p].hists1D[pt][47]->Fill( ratio_p );	
									_procCats[p].hists1D[pt][48]->Fill( genpt/jetpt );

									_procCats[p].hists2D[pt][10]->Fill(sqrt(genpt*genpt + pz*pz), sqrt(jetpt*jetpt + jetpz*jetpz));	
									_procCats[p].hists2D[pt][11]->Fill(genpt, jetpt);	
									_procCats[p].hists2D[pt][12]->Fill(sqrt(genpt*genpt + pz*pz), ratio_p);	
									_procCats[p].hists2D[pt][13]->Fill(sqrt(jetpt*jetpt + jetpz*jetpz), ratio_p);	
									_procCats[p].hists2D[pt][14]->Fill(genpt, genpt/jetpt);	
									_procCats[p].hists2D[pt][15]->Fill(jetpt, genpt/jetpt);
	
									if(pt >= 5) ngenparts_ptge5++;	
								}
								_procCats[p].hists2D[pt][17]->Fill(ngenparts_ptge5,_recoAK4jets[j].GetNConstituents());
							}
						}
						//if no gen top match, skip
						//if(p == 0 && pt == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << endl;
						if(genTopMatchIdxs[j] != -1){
							int gentopidx = genTopMatchIdxs[j];
							_procCats[p].hists2D[pt][61]->Fill(_recoAK4jets[j].pt(), _genTop[gentopidx].pt());
							_procCats[p].hists2D[pt][62]->Fill(_recoAK4jets[j].E(), _genTop[gentopidx].E());
							_procCats[p].hists2D[pt][63]->Fill(_recoAK4jets[j].m(), _genTop[gentopidx].m());
							_procCats[p].hists2D[pt][64]->Fill(_recoAK4jets[j].eta(), _genTop[gentopidx].eta());
							_procCats[p].hists2D[pt][65]->Fill(_recoAK4jets[j].phi(), _genTop[gentopidx].phi());
		
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genTop[gentopidx].eta(), _genTop[gentopidx].phi());
							_procCats[p].hists1D[pt][123]->Fill(dr);
							_procCats[p].hists1D[pt][124]->Fill(_recoAK4jets[j].E()/_genTop[gentopidx].E());
						}
						if(genGluonMatchIdxs[j] != -1){
							int genGluonidx = genGluonMatchIdxs[j];
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genglu[genGluonidx].eta(), _genglu[genGluonidx].phi());
							double eratio = _recoAK4jets[j].E() / _genglu[genGluonidx].E();
							_procCats[p].hists1D[pt][188]->Fill(dr);
							_procCats[p].hists1D[pt][189]->Fill(eratio);
						}
						if(genqMatchIdxs[j] != -1){
							int genqidx = genqMatchIdxs[j];
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genq[genqidx].eta(), _genq[genqidx].phi());
							double eratio = _recoAK4jets[j].E() / _genq[genqidx].E();
							_procCats[p].hists1D[pt][190]->Fill(dr);
							_procCats[p].hists1D[pt][191]->Fill(eratio);
						}
						//if no gen W match, skip
						//if(p == 0 && pt == 0) cout << " matched to gen W #" << genWMatchIdxs[j] << endl;
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];
							_procCats[p].hists2D[pt][86]->Fill(_recoAK4jets[j].pt(), _genW[genWidx].pt());
							_procCats[p].hists2D[pt][87]->Fill(_recoAK4jets[j].E(), _genW[genWidx].E());
							_procCats[p].hists2D[pt][88]->Fill(_recoAK4jets[j].m(), _genW[genWidx].m());
							_procCats[p].hists2D[pt][89]->Fill(_recoAK4jets[j].eta(), _genW[genWidx].eta());
							_procCats[p].hists2D[pt][90]->Fill(_recoAK4jets[j].phi(), _genW[genWidx].phi());
		
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							_procCats[p].hists1D[pt][192]->Fill(dr);
							_procCats[p].hists1D[pt][193]->Fill(_recoAK4jets[j].E()/_genW[genWidx].E());
							//get gen partons from W decay 
							vector<int> genLeadMatchIdxs(2,-1);
							vector<Jet> Wpartons;
							int ggenWidx = _genW[genWidx].GetUserIdx();
							for(int g = 0; g < _genparts.size(); g++){
								int genidx = _genparts[g].GetUserIdx();
								if(_base->genpart_momIdx->at(genidx) != ggenWidx) continue;
								Wpartons.push_back(_genparts[g]);
							}
							if(Wpartons.size() != 2){
								cout << "Error: " << Wpartons.size() << " daughter particles found for W " << genWidx << " skipping hist filling" << endl;
								continue;
							}
							double gendR = dR(Wpartons[0].eta(), Wpartons[0].phi(), Wpartons[1].eta(), Wpartons[1].phi());
							_procCats[p].hists2D[pt][137]->Fill(gendR, jetsize);	
						}
					}
				}
			}
		}

		//_recoAK8 hists 
		void FillRecoAK8JetHists(){
			//do gen matching
			vector<int> genTopMatchIdxs; 
			GenericMatchJet(_recoAK8jets,_genTop,genTopMatchIdxs);
			vector<int> genWMatchIdxs; 
			GenericMatchJet(_recoAK8jets,_genW,genWMatchIdxs);
			vector<int> genGluonMatchIdxs; 
			GenericMatchJet(_recoAK8jets,_genglu,genGluonMatchIdxs);
			vector<int> genqMatchIdxs; 
			GenericMatchJet(_recoAK8jets,_genq,genqMatchIdxs);
			//cout << "final best matches" << endl;
			//for(int b = 0; b < genMatchIdxs_p.size(); b++){
			//	if(genMatchIdxs_p[b] != -1) cout << " jet " << b << " is exclusively matched to gen particle " << genMatchIdxs_p[b] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs_p[b]), _base->genpart_phi->at(genMatchIdxs_p[b]), _recoAK4jets[b].eta(), _recoAK4jets[b].phi()) << endl;
			//	 else cout << " jet " << b << " could not be gen matched" << endl;

			//}
			//else return; //not defined for other jet sizes and nominal FillRecoJetHists() fills for AK4
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				int njets = _recoAK8jets.size();
				_procCats[p].hists1D[0][102]->Fill(njets);
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _recoAK8jets.size(); j++){
					if(_recoAK8jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_recoAK8jets[j].pt() >= _pt_thresh) njets_lead++;
				}
				_procCats[p].hists1D[1][102]->Fill(njets_lead);
				_procCats[p].hists1D[2][102]->Fill(njets_notlead);
				for(int j = 0; j < _recoAK8jets.size(); j++){
					for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
						//define pt bins
						//pt == 1 -> [_pt_thresh,inf)
						if(pt == 1 && _recoAK8jets[j].pt() < _pt_thresh) continue;
						//pt == 2 -> [0,_pt_thresh)
						if(pt == 2 && _recoAK8jets[j].pt() >= _pt_thresh) continue;
	
						Matrix jet_cov = _recoAK8jets[j].GetCovariance();
						double jetsize = CalcSize(jet_cov);
						if(pt == 0 && p == 0) cout << "reco AK8 jet #" << j << " phi " << _recoAK8jets[j].phi() << " eta " << _recoAK8jets[j].eta() << " energy " << _recoAK8jets[j].E() <<  " mass " << _recoAK8jets[j].mass() << " nConstituents " << _recoAK8jets[j].GetNConstituents() << " nRhs " << _recoAK8jets[j].GetNRecHits() << " pt " << _recoAK8jets[j].pt() << " jetsize " << jetsize << endl;
						_procCats[p].hists1D[pt][103]->Fill(_recoAK8jets[j].eta());
						_procCats[p].hists1D[pt][104]->Fill(_recoAK8jets[j].phi_02pi());
						_procCats[p].hists1D[pt][105]->Fill(_recoAK8jets[j].time());
						_procCats[p].hists1D[pt][106]->Fill(_recoAK8jets[j].pt());
						_procCats[p].hists1D[pt][107]->Fill(_recoAK8jets[j].mass());
						_procCats[p].hists1D[pt][108]->Fill(_recoAK8jets[j].e());
						_procCats[p].hists1D[pt][109]->Fill(_recoAK8jets[j].GetNConstituents());
						_procCats[p].hists1D[pt][112]->Fill(jetsize);
						
						_procCats[p].hists2D[pt][130]->Fill(_recoAK8jets[j].m(), jetsize);

						//if no gen top match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen top #" << genTopMatchIdxs[j] << endl;
						if(p == 0 && pt == 0) cout << "reco AK8 jet #" << j << " gen top idx " << genTopMatchIdxs[j] << " gen gluon idx " << genGluonMatchIdxs[j] << " gen q idx " << genqMatchIdxs[j] << " gen w idx " << genWMatchIdxs[j] << endl; 
						if(genTopMatchIdxs[j] != -1){
							int genidx = genTopMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genTop[genidx].eta(), _genTop[genidx].phi());
							double eratio = _recoAK8jets[j].E()/_genTop[genidx].E();
							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to top " << genidx << " with dr " << dr << " and eratio " << eratio << endl;

							_procCats[p].hists1D[pt][110]->Fill(dr);
							_procCats[p].hists1D[pt][111]->Fill(eratio);

						}
						if(genGluonMatchIdxs[j] != -1){
							int genidx = genGluonMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genglu[genidx].eta(), _genglu[genidx].phi());
							double eratio = _recoAK8jets[j].E()/_genglu[genidx].E();
							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to gluon " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][175]->Fill(dr);
							_procCats[p].hists1D[pt][176]->Fill(eratio);


						}
						if(genqMatchIdxs[j] != -1){
							int genidx = genqMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genq[genidx].eta(), _genq[genidx].phi());
							double eratio = _recoAK8jets[j].E()/_genq[genidx].E();

							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to q " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][180]->Fill(dr);
							_procCats[p].hists1D[pt][181]->Fill(eratio);

	
						}

						//if no gen W match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen W #" << genWMatchIdxs[j] << endl;
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							double eratio = _recoAK8jets[j].E()/_genW[genWidx].E();
							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to W " << genWidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][173]->Fill(dr);
							_procCats[p].hists1D[pt][174]->Fill(eratio);


							//get gen partons from W decay 
							vector<int> genLeadMatchIdxs(2,-1);
							vector<Jet> Wpartons;
							int ggenWidx = _genW[genWidx].GetUserIdx();
							for(int g = 0; g < _genparts.size(); g++){
								int genidx = _genparts[g].GetUserIdx();
								if(_base->genpart_momIdx->at(genidx) != ggenWidx) continue;
								Wpartons.push_back(_genparts[g]);
							}
							if(Wpartons.size() != 2){
								cout << "Error: " << Wpartons.size() << " daughter particles found for W " << genWidx << " skipping hist filling" << endl;
								continue;
							}
							double gendR = dR(Wpartons[0].eta(), Wpartons[0].phi(), Wpartons[1].eta(), Wpartons[1].phi());
							_procCats[p].hists2D[pt][133]->Fill(gendR, jetsize);	

						}
					}	
				}
			}
		}

		//_recoAK15 hists 
		void FillRecoAK15JetHists(){
			//do gen matching
			vector<int> genTopMatchIdxs; 
			GenericMatchJet(_recoAK15jets,_genTop,genTopMatchIdxs);
			vector<int> genWMatchIdxs; 
			GenericMatchJet(_recoAK15jets,_genW,genWMatchIdxs);
			vector<int> genGluonMatchIdxs; 
			GenericMatchJet(_recoAK15jets,_genglu,genGluonMatchIdxs);
			vector<int> genqMatchIdxs; 
			GenericMatchJet(_recoAK15jets,_genq,genqMatchIdxs);
			//cout << "final best matches" << endl;
			//for(int b = 0; b < genMatchIdxs_p.size(); b++){
			//	if(genMatchIdxs_p[b] != -1) cout << " jet " << b << " is exclusively matched to gen particle " << genMatchIdxs_p[b] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs_p[b]), _base->genpart_phi->at(genMatchIdxs_p[b]), _recoAK4jets[b].eta(), _recoAK4jets[b].phi()) << endl;
			//	 else cout << " jet " << b << " could not be gen matched" << endl;

			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				int njets = _recoAK15jets.size();
				_procCats[p].hists1D[0][113]->Fill(njets);
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _recoAK15jets.size(); j++){
					if(_recoAK15jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_recoAK15jets[j].pt() >= _pt_thresh) njets_lead++;
				}
				_procCats[p].hists1D[1][113]->Fill(njets_lead);
				_procCats[p].hists1D[2][113]->Fill(njets_notlead);
				for(int j = 0; j < _recoAK15jets.size(); j++){
					for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
						//define pt bins
						//pt == 1 -> [_pt_thresh,inf)
						if(pt == 1 && _recoAK15jets[j].pt() < _pt_thresh) continue;
						//pt == 2 -> [0,_pt_thresh)
						if(pt == 2 && _recoAK15jets[j].pt() >= _pt_thresh) continue;
	
						Matrix jet_cov = _recoAK15jets[j].GetCovariance();
						double jetsize = CalcSize(jet_cov);
						if(pt == 0 && p == 0) cout << "reco AK15 jet #" << j << " phi " << _recoAK15jets[j].phi() << " eta " << _recoAK15jets[j].eta() << " energy " << _recoAK15jets[j].E() <<  " mass " << _recoAK15jets[j].mass() << " nConstituents " << _recoAK15jets[j].GetNConstituents() << " nRhs " << _recoAK15jets[j].GetNRecHits() << " pt " << _recoAK15jets[j].pt() << " jetsize " << jetsize << endl;
						_procCats[p].hists1D[pt][114]->Fill(_recoAK15jets[j].eta());
						_procCats[p].hists1D[pt][115]->Fill(_recoAK15jets[j].phi_02pi());
						_procCats[p].hists1D[pt][116]->Fill(_recoAK15jets[j].time());
						_procCats[p].hists1D[pt][117]->Fill(_recoAK15jets[j].pt());
						_procCats[p].hists1D[pt][118]->Fill(_recoAK15jets[j].mass());
						_procCats[p].hists1D[pt][119]->Fill(_recoAK15jets[j].e());
						_procCats[p].hists1D[pt][120]->Fill(_recoAK15jets[j].GetNConstituents());
						_procCats[p].hists1D[pt][123]->Fill(jetsize);
						
						_procCats[p].hists2D[pt][131]->Fill(_recoAK15jets[j].m(), jetsize);

						//if no gen top match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen top #" << genTopMatchIdxs[j] << endl;
						if(p == 0 && pt == 0) cout << "reco AK15 jet #" << j << " gen top idx " << genTopMatchIdxs[j] << " gen gluon idx " << genGluonMatchIdxs[j] << " gen q idx " << genqMatchIdxs[j] << " gen w idx " << genWMatchIdxs[j] << endl; 
						if(genTopMatchIdxs[j] != -1){
							int genidx = genTopMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genTop[genidx].eta(), _genTop[genidx].phi());
							double eratio = _recoAK15jets[j].E()/_genTop[genidx].E();
							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to top " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	

							_procCats[p].hists1D[pt][121]->Fill(dr);
							_procCats[p].hists1D[pt][122]->Fill(eratio);

							_procCats[p].hists2D[pt][66]->Fill(_recoAK15jets[j].pt(), _genTop[genidx].pt());	
							_procCats[p].hists2D[pt][67]->Fill(_recoAK15jets[j].E(), _genTop[genidx].E());		
							_procCats[p].hists2D[pt][68]->Fill(_recoAK15jets[j].m(), _genTop[genidx].m());	
							_procCats[p].hists2D[pt][69]->Fill(_recoAK15jets[j].eta(), _genTop[genidx].eta());	
							_procCats[p].hists2D[pt][70]->Fill(_recoAK15jets[j].phi(), _genTop[genidx].phi());
							
						
						}
						if(genGluonMatchIdxs[j] != -1){
							int genidx = genGluonMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genglu[genidx].eta(), _genglu[genidx].phi());
							double eratio = _recoAK15jets[j].E()/_genglu[genidx].E();
							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to gluon " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][184]->Fill(dr);
							_procCats[p].hists1D[pt][185]->Fill(eratio);


						}
						if(genqMatchIdxs[j] != -1){
							int genidx = genqMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genq[genidx].eta(), _genq[genidx].phi());
							double eratio = _recoAK15jets[j].E()/_genq[genidx].E();

							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to q " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][186]->Fill(dr);
							_procCats[p].hists1D[pt][187]->Fill(eratio);

	
						}

						//if no gen W match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen W #" << genWMatchIdxs[j] << endl;
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							double eratio = _recoAK15jets[j].E()/_genW[genWidx].E();
							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to W " << genWidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][182]->Fill(dr);
							_procCats[p].hists1D[pt][183]->Fill(eratio);
							_procCats[p].hists2D[pt][91]->Fill(_recoAK15jets[j].pt(), _genW[genWidx].pt());	
							_procCats[p].hists2D[pt][92]->Fill(_recoAK15jets[j].E(), _genW[genWidx].E());		
							_procCats[p].hists2D[pt][93]->Fill(_recoAK15jets[j].m(), _genW[genWidx].m());	
							_procCats[p].hists2D[pt][94]->Fill(_recoAK15jets[j].eta(), _genW[genWidx].eta());	
							_procCats[p].hists2D[pt][95]->Fill(_recoAK15jets[j].phi(), _genW[genWidx].phi());	


							//get gen partons from W decay 
							vector<int> genLeadMatchIdxs(2,-1);
							vector<Jet> Wpartons;
							int ggenWidx = _genW[genWidx].GetUserIdx();
							for(int g = 0; g < _genparts.size(); g++){
								int genidx = _genparts[g].GetUserIdx();
								if(_base->genpart_momIdx->at(genidx) != ggenWidx) continue;
								Wpartons.push_back(_genparts[g]);
							}
							if(Wpartons.size() != 2){
								cout << "Error: " << Wpartons.size() << " daughter particles found for W " << genWidx << " skipping hist filling" << endl;
								continue;
							}
							double gendR = dR(Wpartons[0].eta(), Wpartons[0].phi(), Wpartons[1].eta(), Wpartons[1].phi());
							_procCats[p].hists2D[pt][136]->Fill(gendR, jetsize);	

						}
					}	
				}
			}
		}

		void FillResolutionHists(){
			//need to gen match jets to find difference in pt bw reco - gen
			//dr match (maybe do dE match later)
			vector<int> genMatchIdxs_reco; //one per jet, follows same indexing as jets
			vector<int> genMatchIdxs_bhc; //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genAK4jets, genMatchIdxs_reco); //match BHC jets to reco jets
			int bestIdx;
			//cout << "n reco jets " << _recoAK4jets.size() << " n idxs " << genMatchIdxs_reco.size() << " n pred jets " << _predJets.size() << " n pred idxs " << genMatchIdxs_bhc.size() <<  endl;
			for(int p = 0; p < _procCats.size(); p++){
				//reco jets
				for(int j = 0; j < _recoAK4jets.size(); j++){
					if(genMatchIdxs_reco[j] != -1) cout << " reco jet " << j << " is exclusively matched to gen jet " << genMatchIdxs_reco[j] << " with dr " << dR(_genAK4jets[genMatchIdxs_reco[j]].eta(), _genAK4jets[genMatchIdxs_reco[j]].phi(), _recoAK4jets[j].eta(), _recoAK4jets[j].phi()) << endl;
					 else{ cout << " jet " << j << " could not be gen matched" << endl; continue; }
					 bestIdx = genMatchIdxs_reco[j];
					
					_procCats[p].hists2D[0][0]->Fill(_genAK4jets[bestIdx].E(), _recoAK4jets[j].pt() - _genAK4jets[bestIdx].pt());
					_procCats[p].hists1D[0][22]->Fill(_recoAK4jets[j].pt()/_genAK4jets[bestIdx].pt());
					_procCats[p].hists1D[0][23]->Fill(_recoAK4jets[j].e()/_genAK4jets[bestIdx].e());
					_procCats[p].hists2D[0][1]->Fill(_genAK4jets[bestIdx].pt(), _recoAK4jets[j].pt());
					_procCats[p].hists2D[0][2]->Fill(_genAK4jets[bestIdx].mass(), _recoAK4jets[j].pt()/_genAK4jets[bestIdx].pt());
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
			else if(sample.find("Wgluon") != string::npos){
				procCat Wgluon(_hists1D, _hists2D, "Wgluon", "single W^{#pm} + g",leadsep);
				Wgluon.ids = {-999};
				_procCats.push_back(Wgluon);
			}
			else if(sample.find("singleW") != string::npos){
				procCat singleW(_hists1D, _hists2D, "singleW", "single W^{#pm}",leadsep);
				singleW.ids = {-999};
				_procCats.push_back(singleW);
			}
			else if(sample.find("singleGamma") != string::npos){
				procCat singleGamma(_hists1D, _hists2D, "singleGamma", "single #gamma",leadsep);
				singleGamma.ids = {-999};
				_procCats.push_back(singleGamma);
			}
			else if(sample.find("singleGluon") != string::npos){
				procCat singleGluon(_hists1D, _hists2D, "singleGluon", "single gluon",leadsep);
				singleGluon.ids = {-999};
				_procCats.push_back(singleGluon);
			}
			else return;

		}
		//draw ellispes + tmarkers to canvas and write canvas to file
		//center_coords = center of event display in [eta_c, phi_c] per gen object s.t. center_coords[i] = BayesPoint(eta_c, phi_c) for object i
		//window_width = width of event display in [deta, dphi] per gen object s.t. window_width[i] = BayesPoint(deta, dphi) for object i
		//void WriteEventDisplays(TFile* ofile, vector<BayesPoint> center_coords = {}, vector<BayesPoint> window_width = {}){
		void WriteEventDisplays(TFile* ofile, map<string,BayesPoint> center_coords = {}, map<string,BayesPoint> window_width = {}){
			ofile->cd();
			//write overall event display
			TCanvas* cv = new TCanvas("evtdisp","evtdisp");
			_procCats[0].hists2D[0][129]->Draw("colz");
			//plot jets
			for(int j = 0; j < _predJets.size(); j++){
				_jellipses[j].Draw();
				_jcenters[j].Draw();
				
				for(int k = 0; k < _predJets[j].GetNConstituents(); k++){
					_subclellipses[j][k].Draw();
					_subclcenters[j][k].Draw();
				}
			}
			//plot gen particles
			for(int m = 0; m < _plot_particles.size(); m++){
				_plot_particles[m].Draw();
			}
			cv->Write();
		
			vector<string> names;
			for(auto it = center_coords.begin(); it != center_coords.end(); it++){
				names.push_back(it->first);
			}
	
			//write object specific plots
			for(int h = 0; h < _evtdisps_obj.size(); h++){
				//get match string for center + width
				string name = _evtdisps_obj[h]->GetName();
				string objmatch = "EvtDisplay_etaCell_phiCell_";
				name = name.substr(objmatch.size());
				//skip hists that aren't filled
				if(find(names.begin(), names.end(), name) == names.end()) continue;

				BayesPoint center = center_coords[name];
				double eta_max = 0;
				double phi_max = 0;
				double eta_min = 999;
				double phi_min = 999;
				BayesPoint width = BayesPoint({eta_max, phi_max});// - set by rhs drawn to get everything in frame = window_width[name];
cout << "drawing hist #" << h << " of " << _evtdisps_obj.size() << " with name " << name << endl;
				for(int j = 0; j < _predJets.size(); j++){
					BayesPoint ell_center({_predJets[j].eta(), _predJets[j].phi()});
					ell_center.Translate(center.at(0),0);
					ell_center.CircularTranslate(center.at(1),1);
					//
					//double ell_maj_r = _jellipses[j].GetR1();
					//double ell_min_r = _jellipses[j].GetR2();
					//double theta = _jellipses[j].GetTheta();
					////put in rad
					//theta *= acos(-1)/180;
					//double r_eta = ell_maj_r*cos(theta);
					//double r_phi = ell_maj_r*sin(theta);
					////if full ellipse cannot be drawn in window, skip
					double dr = dR(ell_center.at(0), ell_center.at(1), 0., 0.);
					cout << "jet #" << j << " window width eta " << window_width[name].at(0) << " phi " << window_width[name].at(1) << " this jet center eta " << ell_center.at(0) << " phi " << ell_center.at(1) << endl;
					if(fabs(ell_center.at(0)) > window_width[name].at(0)) continue; //out of frame in eta
					if(fabs(ell_center.at(1)) > window_width[name].at(1)) continue; //out of frame in phi
					//if(fabs(r_eta) > fabs(width.at(0))) continue;
					//if(fabs(r_phi) > fabs(width.at(1))) continue;

					vector<JetPoint> rhs = _predJets[j].GetJetPoints();	
					for(auto rh : rhs){
						double w;
						if(_evt2disp_z == 0)
							w = rh.E();
						else if(_evt2disp_z == 1){
							w = rh.t();
						}
						else
							w = rh.E(); //default energy weighted
						//center according to main gen particle
						BayesPoint rh_pt({rh.eta(), rh.phi()}); //save as BayesPoint to do correct circular translation to (0,0)
						rh_pt.SetWeight(w);
						//translate into local eta, phi coords
						rh_pt.Translate(center.at(0),0);
						rh_pt.CircularTranslate(center.at(1),1);
						if(rh_pt.at(0) > eta_max)
							eta_max = rh_pt.at(0);
						if(rh_pt.at(1) > phi_max)
							phi_max = rh_pt.at(1);
						if(rh_pt.at(0) < eta_min)
							eta_min = rh_pt.at(0);
						if(rh_pt.at(1) < phi_min)
							phi_min = rh_pt.at(1);
					
						_evtdisps_obj[h]->Fill(rh_pt.at(0), rh_pt.at(1), rh_pt.w());
					}
				}
//TODO - check eta/phi_max/min values (see drawn canvases)
cout << "eta_max " << eta_max << " eta_min " << eta_min << " phi_max " << phi_max << " phi_min " << phi_min << endl;
				width.SetValue(max(fabs(eta_max), fabs(eta_min)), 0);
				width.SetValue(max(fabs(phi_max), fabs(phi_min)), 1);
				if(_evtdisps_obj[h]->GetEntries() == 0) continue; //don't draw if not filled for this particle gen obj
				TCanvas* cv_obj = new TCanvas(_evtdisps_obj[h]->GetName(),_evtdisps_obj[h]->GetTitle());
				cv_obj->cd();
				_evtdisps_obj[h]->Draw("colz");
				_evtdisps_obj[h]->GetXaxis()->SetRangeUser(-width.at(0), width.at(0));
				_evtdisps_obj[h]->GetYaxis()->SetRangeUser(-width.at(1), width.at(1));
cout << "hist for " << name << " integral " << _evtdisps_obj[h]->Integral() << " entries " << _evtdisps_obj[h]->GetEntries() << endl;
				//do for gen particles too 
				for(int m = 0; m < _plot_particles.size(); m++){
					BayesPoint m_center({_plot_particles[m].GetX(), _plot_particles[m].GetY()});
					m_center.Translate(center.at(0),0);
					m_center.CircularTranslate(center.at(1),1);
				
					double dr = dR(m_center.at(0), m_center.at(1), 0., 0.);
					if(dr > sqrt(width.at(0)*width.at(0) + width.at(1)*width.at(1))) continue;
					_plot_particles[m].DrawMarker(m_center.at(0), m_center.at(1));

				}
				for(int j = 0; j < _predJets.size(); j++){
					BayesPoint ell_center({_jellipses[j].GetX1(), _jellipses[j].GetY1()});
					ell_center.Translate(center.at(0),0);
					ell_center.CircularTranslate(center.at(1),1);
					if(fabs(ell_center.at(0)) > width.at(0)) continue; //out of frame in eta
					if(fabs(ell_center.at(1)) > width.at(1)) continue; //out of frame in phi
					
					double ell_center_eta = _jellipses[j].GetX1(); 
					double ell_center_phi = _jellipses[j].GetY1();
					double ell_maj_r = _jellipses[j].GetR1();
					double ell_min_r = _jellipses[j].GetR2();
					double theta = _jellipses[j].GetTheta();
					//put in rad
					theta *= acos(-1)/180;
					double r_eta = ell_maj_r*cos(theta);
					double r_phi = ell_maj_r*sin(theta);
					//if full ellipse cannot be drawn in window, skip
					double dr = dR(ell_center.at(0), ell_center.at(1), 0., 0.);
					if(dr > sqrt(width.at(0)*width.at(0) + width.at(1)*width.at(1))) continue;
					if(fabs(r_eta) > fabs(width.at(0))) continue;
					if(fabs(r_phi) > fabs(width.at(1))) continue;

		
					_jellipses[j].DrawEllipse(ell_center.at(0), ell_center.at(1), ell_maj_r, ell_min_r, 0, 360, _jellipses[j].GetTheta());
				

					
					_jcenters[j].DrawMarker(ell_center.at(0), ell_center.at(1));
				
					for(int k = 0; k < _predJets[j].GetNConstituents(); k++){
						BayesPoint subcl_center({_subclellipses[j][k].GetX1(), _subclellipses[j][k].GetY1()});
						subcl_center.Translate(center.at(0),0);
						subcl_center.CircularTranslate(center.at(1),1);
						ell_center_eta = _subclellipses[j][k].GetX1(); 
						ell_center_phi = _subclellipses[j][k].GetY1();
						ell_maj_r = _subclellipses[j][k].GetR1();
						ell_min_r = _subclellipses[j][k].GetR2();
						theta = _subclellipses[j][k].GetTheta();
						
						_subclellipses[j][k].DrawEllipse(subcl_center.at(0), subcl_center.at(1), ell_maj_r, ell_min_r, 0, 360, theta);
					
						_subclcenters[j][k].DrawMarker(subcl_center.at(0), subcl_center.at(1));
					}

				}

				/*

				PointCollection ell_centers, ell_centers_og, m_centers, m_centers_og;
				//center all points
				ell_centers.Translate(center.at(0),0);
				ell_centers.CircularTranslate(center.at(1),1);
				m_centers.Translate(center.at(0),0);
				m_centers.CircularTranslate(center.at(1),1);
		
				for(int el = 0; el < _ellipses.size(); el++){
					_ellipses[el].SetX1(ell_centers.at(el).at(0));
					_ellipses[el].SetY1(ell_centers.at(el).at(1));
					double ell_center_eta = _ellipses[el].GetX1(); 
					double ell_center_phi = _ellipses[el].GetY1();
					double ell_maj_r = _ellipses[el].GetR1();
					double theta = _ellipses[el].GetTheta();
					//put in rad
					theta *= acos(-1)/180;
					double r_eta = ell_maj_r*cos(theta);
					double r_phi = ell_maj_r*sin(theta);

cout << "ell center " << ell_center_eta << ", " << ell_center_phi << " eta plot window " << width.at(0) << " r_eta " << fabs(r_eta) << " phi plot window " << width.at(1) << " r_phi " << fabs(r_phi) << " ell_maj_r " << ell_maj_r << " theta " << theta << endl;
					//if full ellipse cannot be drawn in window, skip
					//if(fabs(r_eta) > fabs(width.at(0))) continue;
					//if(fabs(r_phi) > fabs(width.at(1))) continue;

					double newwidth_eta = width.at(0);
					double newwidth_phi = width.at(1);
					//maybe just set according to r and not r_eta?
					if(fabs(r_eta) > width.at(0)){
						newwidth_eta = ell_center_eta + r_eta;
					}
					if(fabs(r_phi) > width.at(1)){
						newwidth_phi = ell_center_phi + r_phi;
					}
					double newwidth = max(newwidth_eta, newwidth_phi);
cout << "newwidth " << newwidth << endl;
					width.SetValue(newwidth,0);
					width.SetValue(newwidth,1);

				}



				//center plot
				//double xhi = center.at(0) + width.at(0);
				double xhi = width.at(0);
				//xhi = min(xhi, _procCats[0].hists2D[0][129]->GetXaxis()->GetXmax()); 
				//double xlo = center.at(0) - width.at(0); 
				double xlo = -width.at(0); 
				//xlo = max(xlo, _procCats[0].hists2D[0][129]->GetXaxis()->GetXmin()); 
				_evtdisps_obj[h]->GetXaxis()->SetRangeUser(xlo, xhi);
			
				//double yhi = center.at(1) + width.at(1); 
				double yhi = width.at(1); 
				//yhi = min(yhi, _procCats[0].hists2D[0][129]->GetYaxis()->GetXmax()); 
				//double ylo = center.at(1) - width.at(1);
				double ylo = -width.at(1);
				//ylo = min(ylo, _procCats[0].hists2D[0][129]->GetYaxis()->GetXmax()); 
				_evtdisps_obj[h]->GetYaxis()->SetRangeUser(ylo, yhi);
		cout << "plot center for " << name << endl; center.Print();	
				*/	
				/*
cout << "name " << name << " 1 - ylo " << ylo << " yhi " << yhi << " hist y axis upper bound " << _evtdisps_obj[h]->GetYaxis()->GetXmax() << " this hist y upper bound " << _evtdisps_obj[h]->GetYaxis()->GetXmax() << endl;
				//TODO: fix grid
				//DrawGrid(cv_obj,xlo,xhi,ylo,yhi);
   				//gPad->Update();
				//draw all jet + subcl ellipses
				for(int el = 0; el < _ellipses.size(); el++){
					//_ellipses[el].SetX1(ell_centers.at(el).at(0));
					//_ellipses[el].SetY1(ell_centers.at(el).at(1));
					////skip ellipses that are more than width away from center
					//double ell_center_eta = _ellipses[el].GetX1(); 
					//double ell_center_phi = _ellipses[el].GetY1();
					//double ell_maj_r = _ellipses[el].GetR1();
					//double theta = _ellipses[el].GetTheta();
					////put in rad
					//theta *= acos(-1)/180;
					//double r_eta = ell_maj_r*cos(theta);
					//double r_phi = ell_maj_r*sin(theta);

//cout << "ell center " << ell_center_eta << ", " << ell_center_phi << " eta plot window " << width.at(0) << " r_eta " << fabs(r_eta) << " phi plot window " << width.at(1) << " r_phi " << fabs(r_phi) << " ell_maj_r " << ell_maj_r << " theta " << theta << endl;
					//if full ellipse cannot be drawn in window, skip
					//if(fabs(r_eta) > fabs(width.at(0))) continue;
					//if(fabs(r_phi) > fabs(width.at(1))) continue;

					//_ellipses[el].Draw();
				}
				//draw all particles
				for(int m = 0; m < _plot_particles.size(); m++){
				cout << "m " << m << " original - obj center eta " << _plot_particles[m].GetX() << " obj center phi " << _plot_particles[m].GetY() << endl;	
					_plot_particles[m].SetX(m_centers.at(m).at(0));
					_plot_particles[m].SetY(m_centers.at(m).at(1));
					//skip ellipses that are more than width away from center
					double obj_center_eta = _plot_particles[m].GetX(); 
					double obj_center_phi = _plot_particles[m].GetY();
				cout << "centered - obj center eta " << _plot_particles[m].GetX() << " obj center phi " << _plot_particles[m].GetY() << endl;	
					double dr = dR(obj_center_eta, obj_center_phi, 0., 0.);
					if(dr > sqrt(width.at(0)*width.at(0) + width.at(1)*width.at(1))) continue;
					_plot_particles[m].Draw();
				}
				*/
				cv_obj->Write();

			}


		}

		void WriteOutput(TFile* ofile, map<string,BayesPoint> center_coords = {}, map<string,BayesPoint> window_width = {}){
			WriteEmptyProfiles(ofile);
			WriteStackHists(ofile);
			WriteHists(ofile);
			WriteEventDisplays(ofile, center_coords, window_width);
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
				if(_hists1D[i] == nullptr) continue;
				if(_hists1D[i]->GetEntries() == 0) continue;
				_hists1D[i]->Write();
			}
			for(int i = 0; i < (int)_hists2D.size(); i++){
				//name = hists2D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDR2DHist(hists2D[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				if(_hists2D[i] == nullptr) continue;
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
						if(_procCats[p].hists2D[pt][i]->GetEntries() == 0 && dirname.find("meanRecoGenDeltaT") == string::npos) continue;
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
		TH1D* nClusters = new TH1D("BHCJet_nJets","BHCJet_nJets",30,0,30);
		//1 - n subclusters per bhc jets
		TH1D* nSubclusters = new TH1D("BHCJet_nSubclustersJet","BHCJet_nSubclustersJet",30,0,30);
		//2 - bhc subcluster energy
		TH1D* predJet_subClusterEnergy = new TH1D("BHCJet_subclusterEnergy","BHCJet_subclusterEnergy",50,0,500);
		//3 - bhc subcluster eta center
		TH1D* predJet_subClusterEtaCenter = new TH1D("BHCJet_subclusterEtaCenter","BHCJet_subclusterEtaCenter",25,-3.2,3.2);
		//4 - bhc subcluster phi center
		TH1D* predJet_subClusterPhiCenter = new TH1D("BHCJet_subclusterPhiCenter","BHCJet_subclusterPhiCenter",25,0.,8*atan(1));
		//5 - bhc jet subcluster time center
		TH1D* predJet_subClusterTimeCenter = new TH1D("BHCJet_subclusterTimeCenter","BHCJet_subclusterTimeCenter",25,-1,1);
		//6 - bhc jet size
		TH1D* predJet_jetSize = new TH1D("BHCJet_jetSize","BHCJet_jetSize",50,0,2.);
		//7 - bhc jet energy
		TH1D* predJet_energy = new TH1D("BHCJet_energy","BHCJet_energy",25,0,500);
		//8 - bhc jet pt
		TH1D* predJet_pt = new TH1D("BHCJet_pt","BHCJet_pt",25,0,1000);
		//9 - bhc jet mass
		TH1D* predJet_mass = new TH1D("BHCJet_mass","BHCJet_mass",50,0,250);
		//break!!!!!! -1 of all below
		//10 - # pred jets - # reco jets
		TH1D* predGen_nJets = new TH1D("BHCRecoAK4_diffNJets","BHCRecoAK4_diffNJets",20,-10,10);
		//11 - eta sigma
		TH1D* predJet_subClusterEtaSig = new TH1D("BHCJet_subclusterEtaSig","BHCJet_subclusterEtaSig",50,0.,1.);
		//12 - phi sigma
		TH1D* predJet_subClusterPhiSig = new TH1D("BHCJet_subclusterPhiSig","BHCJet_subclusterPhiSig",50,0.,1.);
		//13 - time sigma
		//TODO: make sure there are no negative entries - redo hist bounds after this is checked
		TH1D* predJet_subClusterTimeSig = new TH1D("BHCJet_subclusterTimeSig","BHCJet_subclusterTimeSig",50,-1.,15.);
		//14 - etaphi cov
		TH1D* predJet_subClusteretaPhiCov = new TH1D("BHCJet_subclusterEtaPhiCov","BHCJet_subclusterEtaPhiCov",50,-0.0005,0.0005);
		//15 - time-eta covariance
		TH1D* predJet_subClustertimeEtaCov = new TH1D("BHCJet_subclusterTimeEtaCov","BHCJet_subclusterTimeEtaCov",50,-0.05,0.05);
		//16 - time-phi covariance
		TH1D* predJet_subClustertimePhiCov = new TH1D("BHCJet_subclusterTimePhiCov","BHCJet_subclusterTimePhiCov",50,-0.05,0.05);
		//17 - n reco AK4 jets
		TH1D* nRecoJets = new TH1D("recoAK4_nJets","recoAK4_nJets",30,0,30);
		//18 - reco AK4 jet size
		TH1D* recoJet_jetSize = new TH1D("recoAK4Jet_jetSize","recoAK4Jet_jetSize",50,0,2.);
		//19 - reco AK4 jet energy
		TH1D* recoJet_energy = new TH1D("recoAK4Jet_energy","recoAK4Jet_energy",25,0,500);
		//20 - reco AK4 jet pt
		TH1D* recoJet_pt = new TH1D("recoAK4Jet_pt","recoAK4Jet_pt",25,0,1000);
		//21 - reco AK4 jet mass
		TH1D* recoJet_mass = new TH1D("recoAK4Jet_mass","recoAK4Jet_mass",50,0,250);
		//22 - reco jet pt/gen jet pt
		TH1D* recoGen_jetPtRatio = new TH1D("recoAK4Gen_jetPtRatio","recoAK4Gen_jetPtRatio",20,0,1.5);
		//23 - reco jet e - gen jet e		
		TH1D* recoGen_jetERatio = new TH1D("recoAK4Gen_jetERatio","recoAK4Gen_jetERatio",20,-10,10);
		//24 - eta sigma for jet
		TH1D* predJet_EtaVar = new TH1D("BHCJet_EtaSig","BHCJet_EtaSig",50,0.,1.);
		//25 - phi sigma for jet
		TH1D* predJet_PhiVar = new TH1D("BHCJet_PhiSig","BHCJet_PhiSig",50,0., 1.);
		//26 - time sigma for jet
		TH1D* predJet_TimeVar = new TH1D("BHCJet_TimeSig","BHCJet_TimeSig",50,0.,10.);
		//27 - eta-phi covariance for jet
		TH1D* predJet_etaPhiCov = new TH1D("BHCJet_etaPhiCov","BHCJet_etaPhiCov",25,-0.0005,0.0005);
		//28 - time-eta covariance for jet
		TH1D* predJet_timeEtaCov = new TH1D("BHCJet_timeEtaCov","BHCJet_timeEtaCov",25,-0.2,0.2);
		//29 - time-phi covariance for jet
		TH1D* predJet_timePhiCov = new TH1D("BHCJet_timePhiCov","BHCJet_timePhiCov",25,-0.2,0.2);
		//30 - n GMM clusters in reco jets
		TH1D* reco_nSubclusters = new TH1D("recoAK4Jet_nSubclustersJet","recoAK4Jet_nSubclustersJet",30,0,30);
		//31 - energy per GMM cluster from reco jets
		TH1D* recoJet_subClusterEnergy = new TH1D("recoAK4Jet_subclusterEnergy","recoAK4Jet_subclusterEnergy",50,0,500);
		//32 - time center of GMM cluster from reco jets
		TH1D* recoJet_subClusterTimeCenter = new TH1D("recoAK4Jet_subclusterTimeCenter","recoAK4Jet_subclusterTimeCenter",25,-1,1);
		//33 - eta center of GMM cluster from reco jets
		TH1D* recoJet_subClusterEtaCenter = new TH1D("recoAK4Jet_subclusterEtaCenter","recoAK4Jet_subclusterEtaCenter",25,-3.2,3.2);
		//34 - phi center of GMM cluster from reco jets
		TH1D* recoJet_subClusterPhiCenter = new TH1D("recoAK4Jet_subclusterPhiCenter","recoAK4Jet_subclusterPhiCenter",25,0.,8*atan(1));
		//35 - eta sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterEtaSig = new TH1D("recoAK4Jet_subclusterEtaSig","recoAK4Jet_subclusterEtaSig",50,0.,0.5);
		//36 - phi sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterPhiSig = new TH1D("recoAK4Jet_subclusterPhiSig","recoAK4Jet_subclusterPhiSig",50,0.,0.5);
		//37 - time sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterTimeSig = new TH1D("recoAK4Jet_subclusterTimeSig","recoAK4Jet_subclusterTimeSig",50,0.,10.);
		//38 - eta-phi covariance of GMM cluster from reco jets
		TH1D* recoJet_subClusteretaPhiCov = new TH1D("recoAK4Jet_subclusterEtaPhiCov","recoAK4Jet_subclusterEtaPhiCov",50,-0.0005,0.0005);
		//39 - time-eta covariance of GMM cluster from reco jets
		TH1D* recoJet_subClustertimeEtaCov = new TH1D("recoAK4Jet_subclusterTimeEtaCov","recoAK4Jet_subclusterTimeEtaCov",50,-0.05,0.05);
		//40 - time-phi covariance of GMM cluster from reco jets
		TH1D* recoJet_subClustertimePhiCov = new TH1D("recoAK4Jet_subclusterTimePhiCov","recoAK4Jet_subclusterTimePhiCov",50,-0.05,0.05);
		//41 - n rhs in reco jets
		TH1D* recoJet_nRhs = new TH1D("recoAK4Jet_nRhs","recoAK4Jet_nRhs",300,0,300);
		//42 - rh energy in reco jets
		TH1D* recoJet_rhE = new TH1D("recoAK4Jet_rhE","recoAK4Jet_rhE",50,0,300);
		//43 - AK4 Jet P 
		TH1D* AK4Jet_GenP = new TH1D("recoAK4Jet_GenP","recoAK4Jet_GenP",50,0,500);
		//44 - AK4 Jet Pt 
		TH1D* AK4Jet_GenPt = new TH1D("recoAK4Jet_GenPt","recoAK4Jet_GenPt",50,0,500);
		//45 - Gen particle Pt 
		TH1D* AK4JetConstituent_GenP = new TH1D("recoAK4JetConstituent_GenP","recoAK4JetConstituent_GenP",50,0,200);
		//46 - Gen particle Pt 
		TH1D* AK4JetConstituent_GenPt = new TH1D("recoAK4JetConstituent_GenPt","recoAK4JetConstituent_GenPt",50,0,200);
		//47 - Gen particle P/Gen jet P 
		TH1D* AK4JetConstJetRatio_GenP = new TH1D("recoAK4JetConstJetRatio_GenP","recoAK4JetConstJetRatio_GenP",50,0,1);
		//48 - Gen particle Pt/Gen jet Pt 
		TH1D* AK4JetConstJetRatio_GenPt = new TH1D("recoAK4JetConstJetRatio_GenPt","recoAK4JetConstJetRatio_GenPt",50,0,1);
		//49 - ak4 jet rh times
		TH1D* AK4Jet_rhTimes = new TH1D("recoAK4Jet_rhTimes","recoAK4Jet_rhTimes",200,-25,25);
		//50 - eta-phi covariance of GMM cluster from reco jets normalized
		TH1D* recoJet_subClusteretaPhiCovNorm = new TH1D("recoAK4Jet_subClusteretaPhiCovNorm","recoAK4Jet_subClusteretaPhiCovNorm",50,-1.,1.);
		//51 - time-eta covariance of GMM cluster from reco jets normalized
		TH1D* recoJet_subClustertimeEtaCovNorm = new TH1D("recoAK4Jet_subClustertimeEtaCovNorm","recoAK4Jet_subClustertimeEtaCovNorm",50,-1.,1.);
		//52 - time-phi covariance of GMM cluster from reco jets normalized
		TH1D* recoJet_subClustertimePhiCovNorm = new TH1D("recoAK4Jet_subClustertimePhiCovNorm","recoAK4Jet_subClustertimePhiCovNorm",50,-1.,1.);
		//53 - # gen partons (t, b, q from W)
		TH1D* nGenParticles = new TH1D("nGenParticles","nGenParticles",20,0,20);
		//54 - gen particle eta at detector
		TH1D* genParticle_eta = new TH1D("genParticle_eta","genParticle_eta",25,-3.2,3.2);
		//55 - 62 - gen particle phi at detector
		TH1D* genParticle_phi = new TH1D("genParticle_phi","genParticle_phi",25,0.,8*atan(1));
		//56 - 63 - gen particle time at detector
		TH1D* genParticle_time = new TH1D("genParticle_time","genParticle_time",25,-10,10);
		//57 - gen particle pt		
		TH1D* genParticle_pt = new TH1D("genParticle_pt","genParticle_pt",25,0,500);
		//58 - gen particle mass
		TH1D* genParticle_mass = new TH1D("genParticle_mass","genParticle_mass",25,0,250);
		//59 - gen particle energy
		TH1D* genParticle_energy = new TH1D("genParticle_energy","genParticle_energy",50,0,500);
		//60 - gen AK4 jet eta at detector
		TH1D* genAK4Jet_eta = new TH1D("genAK4Jet_EtaCenter","genAK4Jet_EtaCenter",25,-3.2,3.2);
		//61 - gen AK4 jet phi at detector
		TH1D* genAK4Jet_phi = new TH1D("genAK4Jet_PhiCenter","genAK4Jet_PhiCenter",25,0.,8*atan(1));
		//62 - gen AK4 jet time at detector
		TH1D* genAK4Jet_time = new TH1D("genAK4Jet_TimeCenter","genAK4Jet_TimeCenter",25,-1,1);
		//63 - gen AK4 jet pt		
		TH1D* genAK4Jet_pt = new TH1D("genAK4Jet_pt","genAK4Jet_pt",25,0,1000);
		//64 - gen AK4 jet mass
		TH1D* genAK4Jet_mass = new TH1D("genAK4Jet_mass","genAK4Jet_mass",50,0,50);
		//65 - gen AK4 jet energy
		TH1D* genAK4Jet_energy = new TH1D("genAK4Jet_energy","genAK4Jet_energy",25,0,500);
		//66 - # gen AK4 jets
		TH1D* nJet_genAK4Jet = new TH1D("genAK4_nJets","genAK4_nJets",30,0,30);
		//67 - # constituents per gen AK4 jet
		TH1D* genAK4Jet_nConstituents = new TH1D("genAK4Jet_nConstituents","genAK4Jet_nConstituents",50,0,50);
		//68 - # gen jets - # gen particles
		TH1D* genAK4JetParticle_nDiff = new TH1D("genAK4Jet_genParticle_nDiff","genAK4Jet_genParticle_nDiff",20,-10,10);
		//69 - dR bw gen jet and gen top its exclusively matched to
		TH1D* genAK4JetTop_dR = new TH1D("genAK4Jet_genTop_dR","genAK4Jet_genTop_dR",25,0,1.5);
		//70 - E ratio bw gen jet and gen top its exclusively matched to - gen jet energy/gen top energy
		TH1D* genAK4JetTop_Eratio = new TH1D("genAK4Jet_genTop_Eratio","genAK4Jet_genTop_Eratio",25,0,2);
		//71 - # bhc jets - # gen particles
		TH1D* BHCJetParticle_nDiff = new TH1D("BHCJet_genParticle_nDiff","BHCJet_genParticle_nDiff",20,-10,10);
		//72 - dR bw bhc jet and gen top its exclusively matched to
		TH1D* BHCJetTop_dR = new TH1D("BHCJet_genTop_dR","BHCJet_genTop_dR",25,0,1.5);
		//73 - E ratio bw bhc jet and gen top its exclusively matched to - bhc jet energy/gen top energy
		TH1D* BHCJetTop_Eratio = new TH1D("BHCJet_genTop_Eratio","BHCJet_genTop_Eratio",25,0,2);
		//74 - bhc jet eta center
		TH1D* BHCJet_EtaCenter = new TH1D("BHCJet_EtaCenter","BHCJet_EtaCenter",25,-3.2,3.2);
		//75 - bhc jet phi center
		TH1D* BHCJet_PhiCenter = new TH1D("BHCJet_PhiCenter","BHCJet_PhiCenter",25,0.,8*atan(1));
		//76 - bhc jet center
		TH1D* BHCJet_TimeCenter = new TH1D("BHCJet_TimeCenter","BHCJet_TimeCenter",25,-1,1);
		//77 - rh time
		TH1D* rhTime = new TH1D("rhTime","rhTime",25,-10,10);
		//78 - recoAK4 jet rh eta sig
		TH1D* recoAK4Jet_rhEtaSig = new TH1D("recoAK4Jet_EtaSig","recoAK4_EtaSig",50,0,1);
		//79 - recoAK4 jet rh phi sig
		TH1D* recoAK4Jet_rhPhiSig = new TH1D("recoAK4Jet_PhiSig","recoAK4_PhiSig",50,0,1);
		//80 - BHC jet rh time sig
		TH1D* recoAK4Jet_rhTimeSig = new TH1D("recoAK4Jet_TimeSig","recoAK4_TimeSig",50,0,15.);
		//81 - reco AK4 jet center
		TH1D* recoAK4Jet_TimeCenter = new TH1D("recoAK4Jet_TimeCenter","recoAK4Jet_TimeCenter",25,-1,1);
		//82 - reco AK4 jet eta at detector
		TH1D* recoAK4Jet_EtaCenter = new TH1D("recoAK4Jet_EtaCenter","recoAK4Jet_EtaCenter",25,-3.2,3.2);
		//83 - reco AK4 jet phi at detector
		TH1D* recoAK4Jet_PhiCenter = new TH1D("recoAK4Jet_PhiCenter","recoAK4Jet_PhiCenter",25,0.,8*atan(1));
		//84 - reco AK4 jet # subclusters in event
		TH1D* recoAK4Jet_nSubclustersEvt = new TH1D("recoAK4Jet_nSubclustersEvt","recoAK4Jet_nSubclustersEvt",30,0,30);
		//85 - dr bw subclusters in BHC jet
		TH1D* BHCJet_drSubclusters = new TH1D("BHCJet_drSubclusters","BHCJet_drSubclusters",50,0,0.1);
		//86 - dr bw subclusters in reco AK4 jet
		TH1D* recoAK4Jet_drSubclusters = new TH1D("recoAK4Jet_drSubclusters","recoAK4Jet_drSubclusters",50,0,0.1);
		//87 - rh energy in reco jets
		TH1D* BHCJet_rhE = new TH1D("BHCJet_rhE","BHCJet_rhE",50,0,300);
		//88 - reco AK4 rotundity
		TH1D* recoAK4Jet_rotundity = new TH1D("recoAK4Jet_rotundity","recoAK4Jet_rotundity",50,0.4,1.1);
		//89 - BHC rotundity	
		TH1D* BHCJet_rotundity = new TH1D("BHCJet_rotundity","BHCJet_rotundity",50,0.4,1.1);
		//90 - BHC jet n rhs
		TH1D* BHCJet_nRhs = new TH1D("BHCJet_nRhs","BHCJet_nRhs",300,0,300);
		//91 - gen AK15 jet eta at detector
		TH1D* genAK15Jet_eta = new TH1D("genAK15Jet_EtaCenter","genAK15Jet_EtaCenter",25,-3.2,3.2);
		//92 - gen AK15 jet phi at detector
		TH1D* genAK15Jet_phi = new TH1D("genAK15Jet_PhiCenter","genAK15Jet_PhiCenter",25,0.,8*atan(1));
		//93 - gen AK15 jet time at detector
		TH1D* genAK15Jet_time = new TH1D("genAK15Jet_TimeCenter","genAK15Jet_TimeCenter",25,-1,1);
		//94 - gen AK15 jet pt		
		TH1D* genAK15Jet_pt = new TH1D("genAK15Jet_pt","genAK15Jet_pt",25,0,1000);
		//95 - gen AK15 jet mass
		TH1D* genAK15Jet_mass = new TH1D("genAK15Jet_mass","genAK15Jet_mass",50,0,50);
		//96 - gen AK15 jet energy
		TH1D* genAK15Jet_energy = new TH1D("genAK15Jet_energy","genAK15Jet_energy",25,0,500);
		//97 - # gen AK15 jets
		TH1D* nJet_genAK15Jet = new TH1D("genAK15_nJets","genAK15_nJets",30,0,30);
		//98 - # constituents per gen AK15 jet
		TH1D* genAK15Jet_nConstituents = new TH1D("genAK15Jet_nConstituents","genAK15Jet_nConstituents",50,0,50);
		//99 - # gen jets - # gen particles
		TH1D* genAK15JetParticle_nDiff = new TH1D("genAK15Jet_genParticle_nDiff","genAK15Jet_genParticle_nDiff",20,-10,10);
		//100 - dR bw gen jet and gen top its exclusively matched to
		TH1D* genAK15JetTop_dR = new TH1D("genAK15Jet_genTop_dR","genAK15Jet_genTop_dR",25,0,1.5);
		//101 - E ratio bw gen jet and gen top its exclusively matched to - gen jet energy/gen top energy
		TH1D* genAK15JetTop_Eratio = new TH1D("genAK15Jet_genTop_Eratio","genAK15Jet_genTop_Eratio",25,0,2);
		//102 - n reco AK8 jets
		TH1D* nRecoAK8Jets = new TH1D("recoAK8_nJets","recoAK8_nJets",30,0,30);
		//103 - reco AK8 jet eta at detector
		TH1D* recoAK8Jet_eta = new TH1D("recoAK8Jet_EtaCenter","recoAK8Jet_EtaCenter",25,-3.2,3.2);
		//104 - reco AK8 jet phi at detector
		TH1D* recoAK8Jet_phi = new TH1D("recoAK8Jet_PhiCenter","recoAK8Jet_PhiCenter",25,0.,8*atan(1));
		//105 - reco AK8 jet time at detector
		TH1D* recoAK8Jet_time = new TH1D("recoAK8Jet_TimeCenter","recoAK8Jet_TimeCenter",25,-1,1);
		//106 - reco AK8 jet pt		
		TH1D* recoAK8Jet_pt = new TH1D("recoAK8Jet_pt","recoAK8Jet_pt",25,0,1000);
		//107 - reco AK8 jet mass
		TH1D* recoAK8Jet_mass = new TH1D("recoAK8Jet_mass","recoAK8Jet_mass",50,0,250);
		//108 - reco AK8 jet energy
		TH1D* recoAK8Jet_energy = new TH1D("recoAK8Jet_energy","recoAK8Jet_energy",25,0,500);
		//109 - # constituents per reco AK8 jet
		TH1D* recoAK8Jet_nConstituents = new TH1D("recoAK8Jet_nSubclustersJet","recoAK8Jet_nSubclustersJet",30,0,30);
		//110 - dR bw reco AK8 jet and gen top its exclusively matched to
		TH1D* recoAK8JetTop_dR = new TH1D("recoAK8Jet_genTop_dR","recoAK8Jet_genTop_dR",25,0,1.5);
		//111 - E ratio bw reco AK8 jet and gen top its exclusively matched to - reco AK8 jet energy/gen top energy
		TH1D* recoAK8JetTop_Eratio = new TH1D("recoAK8Jet_genTop_Eratio","recoAK8Jet_genTop_Eratio",25,0,2);
		//112 - reco AK8 jet size
		TH1D* recoAK8Jet_jetSize = new TH1D("recoAK8Jet_jetSize","recoAK8Jet_jetSize",50,0,2.);
		//113 - n reco AK15 jets
		TH1D* nRecoAK15Jets = new TH1D("recoAK15_nJets","recoAK15_nJets",30,0,30);
		//114 - gen AK15 jet eta at detector
		TH1D* recoAK15Jet_eta = new TH1D("recoAK15Jet_EtaCenter","recoAK15Jet_EtaCenter",25,-3.2,3.2);
		//115 - reco AK15 jet phi at detector
		TH1D* recoAK15Jet_phi = new TH1D("recoAK15Jet_PhiCenter","recoAK15Jet_PhiCenter",25,0.,8*atan(1));
		//116 - reco AK15 jet time at detector
		TH1D* recoAK15Jet_time = new TH1D("recoAK15Jet_TimeCenter","recoAK15Jet_TimeCenter",25,-1,1);
		//117 - reco AK15 jet pt		
		TH1D* recoAK15Jet_pt = new TH1D("recoAK15Jet_pt","recoAK15Jet_pt",25,0,1000);
		//118 - reco AK15 jet mass
		TH1D* recoAK15Jet_mass = new TH1D("recoAK15Jet_mass","recoAK15Jet_mass",50,0,250);
		//119 - reco AK15 jet energy
		TH1D* recoAK15Jet_energy = new TH1D("recoAK15Jet_energy","recoAK15Jet_energy",25,0,500);
		//120 - # constituents per reco AK15 jet
		TH1D* recoAK15Jet_nConstituents = new TH1D("recoAK15Jet_nSubclustersJet","recoAK15Jet_nSubclustersJet",30,0,30);
		//121 - dR bw reco AK15 jet and gen top its exclusively matched to
		TH1D* recoAK15JetTop_dR = new TH1D("recoAK15Jet_genTop_dR","recoAK15Jet_genTop_dR",25,0,1.5);
		//122 - E ratio bw reco AK15 jet and gen top its exclusively matched to - reco AK15 jet energy/gen top energy
		TH1D* recoAK15JetTop_Eratio = new TH1D("recoAK15Jet_genTop_Eratio","recoAK15Jet_genTop_Eratio",25,0,2);
		//123 - reco AK15 jet size
		TH1D* recoAK15Jet_jetSize = new TH1D("recoAK15Jet_jetSize","recoAK15Jet_jetSize",50,0,2.);
		//124 - n gen AK8 jets
		TH1D* nGenAK8Jets = new TH1D("genAK8_nJets","genAK8_nJets",30,0,30);
		//125 - gen AK15 jet eta at detector
		TH1D* genAK8Jet_eta = new TH1D("genAK8Jet_EtaCenter","genAK8Jet_EtaCenter",25,-3.2,3.2);
		//126 - gen AK8 jet phi at detector
		TH1D* genAK8Jet_phi = new TH1D("genAK8Jet_PhiCenter","genAK8Jet_PhiCenter",25,0.,8*atan(1));
		//127 - gen AK8 jet time at detector
		TH1D* genAK8Jet_time = new TH1D("genAK8Jet_TimeCenter","genAK8Jet_TimeCenter",25,-1,1);
		//128 - gen AK8 jet pt		
		TH1D* genAK8Jet_pt = new TH1D("genAK8Jet_pt","genAK8Jet_pt",25,0,1000);
		//129 - gen AK8 jet mass
		TH1D* genAK8Jet_mass = new TH1D("genAK8Jet_mass","genAK8Jet_mass",50,0,50);
		//130 - gen AK8 jet energy
		TH1D* genAK8Jet_energy = new TH1D("genAK8Jet_energy","genAK8Jet_energy",25,0,500);
		//131 - # constituents per gen AK8 jet
		TH1D* genAK8Jet_nConstituents = new TH1D("genAK8Jet_nConstituents","genAK8Jet_nConstituents",50,0,50);
		//132 - dR bw gen AK8 jet and gen top its exclusively matched to
		TH1D* genAK8JetTop_dR = new TH1D("genAK8Jet_genTop_dR","genAK8Jet_genTop_dR",25,0,1.5);
		//133 - E ratio bw gen AK8 jet and gen top its exclusively matched to - gen AK8 jet energy/gen top energy
		TH1D* genAK8JetTop_Eratio = new TH1D("genAK8Jet_genTop_Eratio","genAK8Jet_genTop_Eratio",25,0,2);
		//134 - dR bw gen AK4 jet and gen top its exclusively matched to
		TH1D* genAK4JetW_dR = new TH1D("genAK4Jet_genW_dR","genAK4Jet_genW_dR",25,0,1.5);
		//135 - E ratio bw gen AK4 jet and gen top its exclusively matched to - gen AK4 jet energy/gen top energy
		TH1D* genAK4JetW_Eratio = new TH1D("genAK4Jet_genW_Eratio","genAK4Jet_genW_Eratio",25,0,2.5);
		//136 - dR bw gen AK15 jet and gen top its exclusively matched to
		TH1D* genAK15JetW_dR = new TH1D("genAK15Jet_genW_dR","genAK15Jet_genW_dR",25,0,1.5);
		//137 - E ratio bw gen AK15 jet and gen top its exclusively matched to - gen AK15 jet energy/gen top energy
		TH1D* genAK15JetW_Eratio = new TH1D("genAK15Jet_genW_Eratio","genAK15Jet_genW_Eratio",25,0,2.5);
		//138 - dR bw BHC jet and gen W its exclusively matched to
		TH1D* BHCJetW_dR = new TH1D("BHCJet_genW_dR","BHCJet_genW_dR",25,0,1.5);
		//139 - E ratio bw BHC jet and gen W its exclusively matched to - BHC jet energy/gen W energy
		TH1D* BHCJetW_Eratio = new TH1D("BHCJet_genW_Eratio","BHCJet_genW_Eratio",25,0,2.5);
		//140 - # subclusters in BHC jets matched to Ws
		TH1D* BHCJetW_nSubclusters = new TH1D("BHCJetW_nSubclusters","BHCJetW_nSubclusters",10,0,10);
		//141 - subcluster energy in BHC jets matched to Ws
		TH1D* BHCJetW_subClusterEnergy = new TH1D("BHCJetW_subclusterEnergy","BHCJetW_subclusterEnergy",25,0,500);
		//142 - BHC jet subcluster mass
		TH1D* BHCJet_subClusterMass = new TH1D("BHCJet_subclusterMass","BHCJet_subclusterMass",25,0,150);
		//143 - BHC jet subcluster # effective rechits
		TH1D* BHCJet_subClusterEffnRhs = new TH1D("BHCJet_subclusterEffnRhs","BHCJet_subclusterEffnRhs",25,0,150);
		//144 - reco AK4 jet subcluster mass
		TH1D* recoAK4Jet_subClusterMass = new TH1D("recoAK4Jet_subclusterMass","recoAK4Jet_subclusterMass",25,0,100);
		//145 - reco AK4 jet subcluster # effective rechits
		TH1D* recoAK4Jet_subClusterEffnRhs = new TH1D("recoAK4Jet_subclusterEffnRhs","recoAK4Jet_subclusterEffnRhs",25,0,150);
		//146 - subcluster mass in BHC jets matched to Ws
		TH1D* BHCJetW_subClusterMass = new TH1D("BHCJetW_subclusterMass","BHCJetW_subclusterMass",25,0,200);
		//147 - BHC jets - invariant mass of lead two subclusters (for jets with at least 2 subclusters)
		TH1D* BHCJetW_subClusterLeadInvMass = new TH1D("BHCJetW_subclusterLeadInvMass","BHCJetW_subclusterLeadInvMass",25,0,200);
		//148 - BHC jets - # ghost leftover after clustering
		TH1D* BHCJet_nGhosts = new TH1D("BHCJet_nGhosts","BHCJet_nGhosts",10,0,10);		
		//149 - BHC jets - ghost subcl energy
		TH1D* BHCJet_ghostSubClusterEnergy = new TH1D("BHCJet_ghostSubclusterEnergy","BHCJet_ghostSubclusterEnergy",25,0,500);
		//150 - BHC jets - ghost subcl eff # rhs
		TH1D* BHCJet_ghostSubClusterEffnRhs = new TH1D("BHCJet_ghostSubclusterEffnRhs","BHCJet_ghostSubclusterEffnRhs",25,0,200);
		//151 - BHC jets - gen-matched W - Eratio (reco/gen) of gen partons in W decay and 2 lead subclusters in BHC jet
		TH1D* BHCJetW_subclParton_dR = new TH1D("BHCJetW_subclParton_dR","BHCJetW_subclParton_dR",25,0,2);
		//152 - BHC jets - gen-matched W - Eratio (reco/gen) of gen partons in W decay and 2 lead subclusters in BHC jet
		TH1D* BHCJetW_subclParton_Eratio = new TH1D("BHCJetW_subclParton_Eratio","BHCJetW_subclParton_Eratio",25,0,2);
		//153 - # subclusters in BHC jets matched to Tops
		TH1D* BHCJetTop_nSubclusters = new TH1D("BHCJetTop_nSubclusters","BHCJetTop_nSubclusters",10,0,10);
		//154 - subcluster mass in BHC jets matched to Tops
		TH1D* BHCJetTop_subClusterMass = new TH1D("BHCJetTop_subclusterMass","BHCJetTop_subclusterMass",25,0,200);
		//155 - BHC jets - invariant mass of lead two subclusters (for jets with at least 2 subclusters)
		TH1D* BHCJetTop_subClusterLeadInvMass = new TH1D("BHCJetTop_subclusterLeadInvMass","BHCJetTop_subclusterLeadInvMass",25,0,200);
		//156 - dR bw BHC jet and gen q its exclusively matched to
		TH1D* BHCJetq_dR = new TH1D("BHCJet_genq_dR","BHCJet_genq_dR",25,0,1.5);
		//157 - E ratio bw BHC jet and gen q its exclusively matched to - BHC jet energy/gen q energy
		TH1D* BHCJetq_Eratio = new TH1D("BHCJet_genq_Eratio","BHCJet_genq_Eratio",25,0,2.5);
		//158 - # subclusters in BHC jets matched to qs
		TH1D* BHCJetq_nSubclusters = new TH1D("BHCJetq_nSubclusters","BHCJetq_nSubclusters",10,0,10);
		//159 - subcluster mass in BHC jets matched to qs
		TH1D* BHCJetq_subClusterMass = new TH1D("BHCJetq_subclusterMass","BHCJetq_subclusterMass",25,0,500);
		//160 - BHC jets - gen-matched to W - subcluster eta center		
		TH1D* BHCJetW_subclEtaCenter = new TH1D("BHCJetW_subclEtaCenter","BHCJetW_subclEtaCenter",25,-3.2,3.2);
		//161 - BHC jets - gen-matched to W - subcluster phi center		
		TH1D* BHCJetW_subclPhiCenter = new TH1D("BHCJetW_subclPhiCenter","BHCJetW_subclPhiCenter",25,-3.2,3.2);
		//162 - BHC jets - gen-matched to W - subcluster time center		
		TH1D* BHCJetW_subclTimeCenter = new TH1D("BHCJetW_subclTimeCenter","BHCJetW_subclTimeCenter",25,-3.2,3.2);
		//163 - BHC jets - gen-matched to W - eta sigma of GMM cluster 
		TH1D* BHCJetW_subClusterEtaSig = new TH1D("BHCJetW_subclusterEtaSig","BHCJetW_subclusterEtaSig",50,0.,0.5);
		//164 - BHC jets - gen-matched to W - phi sigma of GMM cluster 
		TH1D* BHCJetW_subClusterPhiSig = new TH1D("BHCJetW_subclusterPhiSig","BHCJetW_subclusterPhiSig",50,0.,0.5);
		//165 - BHC jets - gen-matched to W - time sigma of GMM cluster 
		TH1D* BHCJetW_subClusterTimeSig = new TH1D("BHCJetW_subclusterTimeSig","BHCJetW_subclusterTimeSig",50,0.,5.);
		//166 - BHC jets - gen-matched to W - eta-phi covariance of GMM cluster 
		TH1D* BHCJetW_subClusteretaPhiCov = new TH1D("BHCJetW_subclusterEtaPhiCov","BHCJetW_subclusterEtaPhiCov",50,-0.0005,0.0005);
		//167 - BHC jets - gen-matched to W - time-eta covariance of GMM cluster 
		TH1D* BHCJetW_subClustertimeEtaCov = new TH1D("BHCJetW_subclusterTimeEtaCov","BHCJetW_subclusterTimeEtaCov",50,-0.05,0.05);
		//168 - 182 - BHC jets - gen-matched to W - time-phi covariance of GMM cluster 
		TH1D* BHCJetW_subClustertimePhiCov = new TH1D("BHCJetW_subclusterTimePhiCov","BHCJetW_subclusterTimePhiCov",50,-0.05,0.05);
		//169 - high mass + W-matched BHC jets - pt of subclusters gen-matched to W partons
		TH1D* BHCJetW_highMass_partonMatchSubclPt = new TH1D("BHCJetW_highMass_partonMatchSubclPt","BHCJetW_highMass_partonMatchSubclPt;SubclPt",25,0,500);
		//170 - high mass + W-matched BHC jets - pt of subclusters NOT gen-matched W partons
		TH1D* BHCJetW_highMass_partonNoMatchSubclPt = new TH1D("BHCJetW_highMass_partonNoMatchSubclPt","BHCJetW_highMass_partonNoMatchSubclPt;SubclPt",25,0,500);
		//171 - high mass + W-matched BHC jets - subclSize of subclusters gen-matched to W partons
		TH1D* BHCJetW_highMass_partonMatchSubclSize = new TH1D("BHCJetW_highMass_partonMatchSubclSize","BHCJetW_highMass_partonMatchSubclSize;SubclSize",25,0,1.5);
		//172 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched W partons
		TH1D* BHCJetW_highMass_partonNoMatchSubclSize = new TH1D("BHCJetW_highMass_partonNoMatchSubclSize","BHCJetW_highMass_partonNoMatchSubclSize;SubclSize",25,0,1.5);
		//173 - dR bw reco AK8 jet and gen W its exclusively matched to
		TH1D* recoAK8JetW_dR = new TH1D("recoAK8Jet_genW_dR","recoAK8Jet_genW_dR",25,0,1.5);
		//174 - E ratio bw reco AK8 jet and gen W its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK8JetW_Eratio = new TH1D("recoAK8Jet_genW_Eratio","recoAK8Jet_genW_Eratio",25,0,2.5);
		//175 - dR bw reco AK8 jet and gen gluon its exclusively matched to
		TH1D* recoAK8JetGluon_dR = new TH1D("recoAK8Jet_genGluon_dR","recoAK8Jet_genGluon_dR",25,0,1.5);
		//176 - E ratio bw reco AK8 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK8JetGluon_Eratio = new TH1D("recoAK8Jet_genGluon_Eratio","recoAK8Jet_genGluon_Eratio",25,0,2);
		//177 - dR bw bhc jet and gen gluon its exclusively matched to
		TH1D* BHCJetGluon_dR = new TH1D("BHCJet_genGluon_dR","BHCJet_genGluon_dR",25,0,1.5);
		//178 - E ratio bw bhc jet and gen gluon its exclusively matched to - bhc jet energy/gen top energy
		TH1D* BHCJetGluon_Eratio = new TH1D("BHCJet_genGluon_Eratio","BHCJet_genGluon_Eratio",25,0,2);
		//179 - high mass + W-matched BHC jets - # subclusters
		TH1D* BHCJetW_highMass_nSubclustersJet = new TH1D("BHCJetW_highMass_nSubclustersJet","BHCJetW_highMass_nSubclustersJet",10,0,10);
		//180 - dR bw reco AK8 jet and gen q its exclusively matched to
		TH1D* recoAK8Jetq_dR = new TH1D("recoAK8Jet_genq_dR","recoAK8Jet_genq_dR",25,0,1.5);
		//181 - E ratio bw reco AK8 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK8Jetq_Eratio = new TH1D("recoAK8Jet_genq_Eratio","recoAK8Jet_genq_Eratio",25,0,2.5);
		//182 - dR bw reco AK15 jet and gen W its exclusively matched to
		TH1D* recoAK15JetW_dR = new TH1D("recoAK15Jet_genW_dR","recoAK15Jet_genW_dR",25,0,1.5);
		//183 - E ratio bw reco AK15 jet and gen W its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK15JetW_Eratio = new TH1D("recoAK15Jet_genW_Eratio","recoAK15Jet_genW_Eratio",25,0,2.5);
		//184 - dR bw reco AK15 jet and gen gluon its exclusively matched to
		TH1D* recoAK15JetGluon_dR = new TH1D("recoAK15Jet_genGluon_dR","recoAK15Jet_genGluon_dR",25,0,1.5);
		//185 - E ratio bw reco AK15 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK15JetGluon_Eratio = new TH1D("recoAK15Jet_genGluon_Eratio","recoAK15Jet_genGluon_Eratio",25,0,2);
		//186 - dR bw reco AK15 jet and gen q its exclusively matched to
		TH1D* recoAK15Jetq_dR = new TH1D("recoAK15Jet_genq_dR","recoAK15Jet_genq_dR",25,0,1.5);
		//187 - E ratio bw reco AK15 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK15Jetq_Eratio = new TH1D("recoAK15Jet_genq_Eratio","recoAK15Jet_genq_Eratio",25,0,2.5);
		//188 - dR bw reco AK4 jet and gen gluon its exclusively matched to
		TH1D* recoAK4JetGluon_dR = new TH1D("recoAK4Jet_genGluon_dR","recoAK4Jet_genGluon_dR",25,0,1.5);
		//189 - E ratio bw reco AK4 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK4JetGluon_Eratio = new TH1D("recoAK4Jet_genGluon_Eratio","recoAK4Jet_genGluon_Eratio",25,0,2);
		//190 - dR bw reco AK4 jet and gen q its exclusively matched to
		TH1D* recoAK4Jetq_dR = new TH1D("recoAK4Jet_genq_dR","recoAK4Jet_genq_dR",25,0,1.5);
		//191 - E ratio bw reco AK4 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK4Jetq_Eratio = new TH1D("recoAK4Jet_genq_Eratio","recoAK4Jet_genq_Eratio",25,0,2.5);
		//192 - dR bw reco AK4 jet and gen W its exclusively matched to
		TH1D* recoAK4JetW_dR = new TH1D("recoAK4Jet_genW_dR","recoAK4Jet_genW_dR",25,0,1.5);
		//193 - E ratio bw reco AK4 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK4JetW_Eratio = new TH1D("recoAK4Jet_genW_Eratio","recoAK4Jet_genW_Eratio",25,0,2.5);
		//194 - high mass + q-matched BHC jets - pt of subclusters gen-matched to W partons
		TH1D* BHCJetq_ge2Subcls_partonMatchSubclPt = new TH1D("BHCJetq_ge2Subcls_partonMatchSubclPt","BHCJetq_ge2Subcls_partonMatchSubclPt;SubclPt",25,0,500);
		//195 - high mass + q-matched BHC jets - pt of subclusters NOT gen-matched W partons
		TH1D* BHCJetq_ge2Subcls_partonNoMatchSubclPt = new TH1D("BHCJetq_ge2Subcls_partonNoMatchSubclPt","BHCJetq_ge2Subcls_partonNoMatchSubclPt;SubclPt",25,0,500);
		//196 - high mass + q-matched BHC jets - subclSize of subclusters gen-matched to W partons
		TH1D* BHCJetq_ge2Subcls_partonMatchSubclSize = new TH1D("BHCJetq_ge2Subcls_partonMatchSubclSize","BHCJetq_ge2Subcls_partonMatchSubclSize;SubclSize",25,0,1.5);
		//197 - high mass + q-matched BHC jets - subclSize of subclusters NOT gen-matched W partons
		TH1D* BHCJetq_ge2Subcls_partonNoMatchSubclSize = new TH1D("BHCJetq_ge2Subcls_partonNoMatchSubclSize","BHCJetq_ge2Subcls_partonNoMatchSubclSize;SubclSize",25,0,1.5);
		//198 - low mass + W-matched BHC jets - # subclusters
		TH1D* BHCJetW_lowMass_nSubclustersJet = new TH1D("BHCJetW_lowMass_nSubclustersJet","BHCJetW_lowMass_nSubclustersJet",10,0,10);
		//199 - high mass + W-matched BHC jets - # subclusters
		TH1D* BHCJetW_Wmass_nSubclustersJet = new TH1D("BHCJetW_Wmass_nSubclustersJet","BHCJetW_Wmass_nSubclustersJet",10,0,10);
		//200 - # subclusters in BHC jets matched to Gluons
		TH1D* BHCJetGluon_nSubclusters = new TH1D("BHCJetGluon_nSubclusters","BHCJetGluon_nSubclusters",10,0,10);
		//201 - BHC jets - gen-matched q - Eratio (reco/gen) of gen q to lead subcluster in BHC jet
		TH1D* BHCJetq_subclParton_dR = new TH1D("BHCJetq_subclParton_dR","BHCJetq_subclParton_dR",25,0,2);
		//202 - BHC jets - gen-matched q - Eratio (reco/gen) of gen q to lead subcluster in BHC jet
		TH1D* BHCJetq_subclParton_Eratio = new TH1D("BHCJetq_subclParton_Eratio","BHCJetq_subclParton_Eratio",25,0,2);
		//203 - BHC jets - gen-matched gluon - Eratio (reco/gen) of gen gluon to lead subcluster in BHC jet
		TH1D* BHCJetGluon_subclParton_dR = new TH1D("BHCJetGluon_subclParton_dR","BHCJetGluon_subclParton_dR",25,0,2);
		//204 - BHC jets - gen-matched gluon - Eratio (reco/gen) of gen gluon to lead subcluster in BHC jet
		TH1D* BHCJetGluon_subclParton_Eratio = new TH1D("BHCJetGluon_subclParton_Eratio","BHCJetGluon_subclParton_Eratio",25,0,2);
		//205 - high mass + W-matched BHC jets - pt of subclusters gen-matched to W partons / pt jet
		TH1D* BHCJetW_highMass_partonMatchSubclPtOvJetPt = new TH1D("BHCJetW_highMass_partonMatchSubclPtOvJetPt","BHCJetW_highMass_partonMatchSubclPtOvJetPt;SubclPtOvJetPt",25,0,1.2);
		//206 - high mass + W-matched BHC jets - pt of subclusters NOT gen-matched W partons / pt jet
		TH1D* BHCJetW_highMass_partonNoMatchSubclPtOvJetPt = new TH1D("BHCJetW_highMass_partonNoMatchSubclPtOvJetPt","BHCJetW_highMass_partonNoMatchSubclPtOvJetPt;SubclPtOvJetPt",25,0,1.2);
		//207 - high mass + W-matched BHC jets - subclSize of subclusters gen-matched to W partons / size jet
		TH1D* BHCJetW_highMass_partonMatchSubclSizeOvJetSize = new TH1D("BHCJetW_highMass_partonMatchSubclSizeOvJetSize","BHCJetW_highMass_partonMatchSubclSizeOvJetSize;SubclSizeOvJetSize",25,0,5.);
		//208 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched W partons / size jet
		TH1D* BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize = new TH1D("BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize","BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize;SubclSizeOvJetSize",25,0,5.);


		////////////////////////////////////////////////////////////////////////
		////////////////////////////////2D plots////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		
		//0 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_recoGen = new TH2D("jetGenE_diffDeltaPt_recoGen","jetGenE_diffDeltaPt_recoGen;jet_{gen} E (GeV);#Delta p_{T}_{reco, gen} (GeV)",4,&xbins_recoGenPt[0],50,-50,50);
		//1 - 2D histogram of gen pT vs reco pT
		TH2D* genPt_recoPt = new TH2D("genPt_recoPt","genPt_recoPt;genpt;recopt",50,5,50,50,5,50);
		//2 - gen jet mass vs reco/gen pt
		TH2D* genJetMass_recoGenPtRatio = new TH2D("genJetMass_recoAK4GenPtRatio","genJetMass_recoAK4GenPtRatio",20,0,10,20,0.8,1.2);
		//3 - reco jet mass vs reco jet pt
		TH2D* recoJetMass_recoJetPt = new TH2D("recoAK4Jet_jetMass_jetPt","recoAK4Jet_jetMass_jetPt;recoAK4JetMass;recoAK4JetPt",50,0,250,50,0,250);
		//4 - reco jet mass vs reco jet jetSize
		TH2D* recoJetMass_recoJetSize = new TH2D("recoAK4Jet_jetMass_jetSize","recoAK4Jet_jetMass_jetSize;recoAK4JetMass;recoAK4JetSize",50,0,250,50,0,2.);
		//5 - pred jet mass vs pred jet pt
		TH2D* predJetMass_predJetPt = new TH2D("BHCJet_jetMass_jetPt","BHCJet_jetMass_jetPt;BHCJetMass;BHCJetPt",50,0,250,50,0,250);
		//6 - pred jet pt vs pred jet jetSize
		TH2D* predJetPt_predJetSize = new TH2D("BHCJetPt_BHCJetSize","BHCJetPt_BHCJetSize;BHCJetPt;BHCJetSize",50,0,250,50,0,2);
		//7 - pred jet n subclusters vs jet size
		TH2D* prednSubclusters_jetSize = new TH2D("BHCJet_nSubclustersJet_jetSize","BHCJet_nSubclustersJet_jetSize;nSubclusters;jetsize",30,0,30,50,0,2);
		//8 - # rhs vs # subclusters for AK4 jets
		TH2D* AK4Jet_nRhs_nSubclustersJet = new TH2D("recoAK4Jet_nRhs_nSubclustersJet","recoAK4Jet_nRhs_nSubclustersJet;nRhs;nSubclustersJet;a.u.",300,0,300,30,0,30);
		//9 - # gen particles from gen-matched jet vs # subclusters for AK4 jets
		TH2D* AK4Jet_nGenParts_nSubclusters = new TH2D("recoAK4Jet_nGenParts_nSubclusters","recoAK4Jet_nGenParts_nSubclusters;nGenParts;nSubclusters;a.u.",40,0,40,30,0,30);
		//10 - gen particle p vs gen jet p for AK4 jets
		TH2D* AK4Jet_genP_genJetP = new TH2D("recoAK4Jet_genP_genJetP","recoAK4Jet_genP_genJetP;genP;genJetP;a.u.",50,0,200,50,0,500);
		//11 - gen particle pt vs gen jet pt for AK4 jets
		TH2D* AK4Jet_genPt_genJetPt = new TH2D("recoAK4Jet_genPt_genJetPt","recoAK4Jet_genPt_genJetPt;genPt;genJetPt;a.u.",50,0,200,50,0,500);
		//12 - gen particle p vs gen particle p/gen jet p for AK4 jets
		TH2D* AK4Jet_genP_genPartJetPRatio = new TH2D("recoAK4Jet_genP_genPartJetPRatio","recoAK4Jet_genP_genPartJetPRatio;genP;genPartJetPRatio;a.u.",50,0,200,50,0,1);
		//13 - gen jet p vs gen particle p/gen jet p for AK4 jets
		TH2D* AK4Jet_genJetP_genPartJetPRatio = new TH2D("recoAK4Jet_genJetP_genPartJetPRatio","recoAK4Jet_genJetP_genPartJetPRatio;genJetP;genPartJetPRatio;a.u.",50,0,500,50,0,1);
		//14 - gen particle pt vs gen particle pt/gen jet pt for AK4 jets
		TH2D* AK4Jet_genPt_genPartJetPtRatio = new TH2D("recoAK4Jet_genPt_genPartJetPtRatio","recoAK4Jet_genPt_genPartJetPtRatio;genPt;genPartJetPtRatio;a.u.",50,0,200,50,0,1);
		//15 - gen jet pt vs gen particle pt/gen jet pt for AK4 jets
		TH2D* AK4Jet_genJetPt_genPartJetPtRatio = new TH2D("recoAK4Jet_genJetPt_genPartJetPtRatio","recoAK4Jet_genJetPt_genPartJetPtRatio;genJetPt;genPartJetPtRatio;a.u.",50,0,500,50,0,1);
		//16 - gen jet pt vs # subclusters
		TH2D* AK4Jet_genJetPt_nSubclusters = new TH2D("recoAK4Jet_genJetPt_nSubclusters","recoAK4Jet_genJetPt_nSubclusters;genJetPt;nSubclusters",50,0,500,30,0,30);	
		//17 - # gen particles w/ pt > 5 gev from gen-matched jet vs # subclusters for AK4 jets
		TH2D* AK4Jet_nGenPartsptge5_nSubclusters = new TH2D("recoAK4Jet_nGenPartsptge5_nSubclusters","recoAK4Jet_nGenPartsptge5_nSubclusters;nGenPartsptge5;nSubclusters;a.u.",20,0,20,30,0,30);
		//18 - geo energy avg vs difference in time for adjacent crystals in same obj w/in 10% energy 
		TH2D* geoEavg_diffDeltaTime_adjRhs = new TH2D("geoEavg_diffDeltaTime_adjRhs","geoEavg_diffDeltaTime_adjRhs;geoEavg;diffDeltaTime;a.u.",xbins.size()-1,&xbins[0],25,-5,5);
		//19 - eta-phi cov vs time-eta cov 
		TH2D* recoJet_subClusteretaPhiCov_timeEtaCov = new TH2D("recoAK4Jet_subClusteretaPhiCov_timeEtaCov","recoAK4Jet_subClusteretaPhiCov_timeEtaCov;etaPhiCov;timeEtaCov",20,-0.0005,0.0005,20,-0.2,0.2);
		//20 - 32 - eta-phi cov norm vs time-eta cov norm 
		TH2D* recoJet_subClusteretaPhiCovNorm_timeEtaCovNorm = new TH2D("recoAK4Jet_subClusteretaPhiCovNorm_timeEtaCovNorm","recoAK4Jet_subClusteretaPhiCovNorm_timeEtaCovNorm;etaPhiCovNorm;timeEtaCovNorm",20,-1.,1.,20, -1, 1.);
		//21 - 33 - reco AK4 jet multiplicity vs jet size
		TH2D* recoAK4Jet_nJets_jetSize = new TH2D("recoAK4Jet_nJets_jetSize","recoAK4Jet_nJets_jetSize;nJets;jetSize",30,0,30,50,0,2);
		//22 - 34 - BHC jet multiplicity vs jet size
		TH2D* BHCJet_nJets_jetSize = new TH2D("BHCJet_nJets_jetSize","BHCJet_nJets_jetSize;nJets;jetSize",30,0,30,50,0,2);
		//23 - 36 - # subclusters vs # constituents in gen jets for gen-matched AK4 jets
		TH2D* recoAK4nSubclusters_genAK4nConstituents = new TH2D("recoAK4nSubclusters_genAK4nConstituents","recoAK4nSubclusters_genAK4nConstituents;recoAK4nSubclusters;genAK4nConstituents",30,0,30,30,0,30);
		//24 - 37 - # subclusters vs jet mass for reco AK4 jets
		TH2D* recoAK4Jet_nSubclustersJet_mass = new TH2D("recoAK4Jet_nSubclustersJet_mass","recoAK4Jet_nSubclustersJet_mass;nSubclustersJet;mass",30,0,30,50,0,250);
		//25 - 38 - # subclusters vs jet energy for reco AK4 jets
		TH2D* recoAK4Jet_nSubclustersJet_energy = new TH2D("recoAK4Jet_nSubclustersJet_energy","recoAK4Jet_nSubclustersJet_energy;nSubclustersJet;energy",30,0,30,50,0,2000);
		//26 - 39 - # subclusters/evt vs # subclusters/jet for reco AK4 jets
		TH2D* recoAK4Jet_nSubclustersEvt_nJet = new TH2D("recoAK4Jet_nSubclustersEvt_nJet","recoAK4Jet_nSubclustersEvt_nJet;nSubclustersEvt;nJet",30,0,30,10,0,10);
		//27 - 40 - # subclusters vs jet mass for BHC jets
		TH2D* BHCJet_nSubclustersJet_mass = new TH2D("BHCJet_nSubclustersJet_mass","BHCJet_nSubclustersJet_mass;nSubclustersJet;mass",30,0,30,50,0,250);
		//28 - 41 - # subclusters vs jet energy for BHC jets
		TH2D* BHCJet_nSubclustersJet_energy = new TH2D("BHCJet_nSubclustersJet_energy","BHCJet_nSubclustersJet_energy;nSubclusters;energy",30,0,30,50,0,2000);
		//29 - 42 - # subclusters/evt vs # subclusters/jet for BHC jets
		TH2D* BHCJet_nSubclustersEvt_nJet = new TH2D("BHCJet_nSubclustersEvt_nJet","BHCJet_nSubclustersEvt_nJet;nSubclustersEvt;nJet",30,0,30,10,0,10);
		//30 - 43 - # subclusters in reco AK4 jet and # subclusters in dR matched BHC jet (if matching can be 1:1 ie # BHC jets = # reco AK4 jets)
		TH2D* recoAK4JetnSubclustersJet_BHCJetnSubclustersJet = new TH2D("recoAK4JetnSubclustersJet_BHCJetnSubclustersJet","recoAK4JetnSubclustersJet_BHCJetnSubclustersJet;recoAK4JetnSubclustersJet;BHCJetnSubclustersJet",30,0,30,30,0,30);
		//31 - 44 - reco AK4 jet energy vs jet mass
		TH2D* recoAK4Jet_jetEnergy_jetMass = new TH2D("recoAK4Jet_jetEnergy_jetMass","recoAK4Jet_jetEnergy_jetMass;jetEnergy;jetMass",50,0,2000,50,0,250);
		//32 - 45 - BHC jet energy vs mass
		TH2D* BHCJet_jetEnergy_jetMass = new TH2D("BHCJet_jetEnergy_jetMass","BHCJet_jetEnergy_jetMass;jetEnergy;jetMass",50,0,2000,50,0,250);
		//33 - 46 - reco AK4 jets # subclusters vs jet size
		TH2D* recoAK4Jet_nSubclusters_jetSize = new TH2D("recoAK4Jet_nSubclustersJet_jetSize","recoAK4Jet_nSubclustersJet_jetSize;nSubclustersJet;jetsize",30,0,30,50,0,2);
		//34 - 47 - reco AK4 jet energy vs jet size
		TH2D* recoAK4Jet_jetEnergy_jetSize = new TH2D("recoAK4Jet_jetEnergy_jetSize","recoAK4Jet_jetEnergy_jetSize;jetEnergy;jetSize",50,0,500,50,0,2);
		//35 - 48 - BHC jet energy vs jet size
		TH2D* BHCJet_jetEnergy_jetSize = new TH2D("BHCJet_jetEnergy_jetSize","BHCJet_jetEnergy_jetSize;jetEnergy;jetSize",50,0,500,50,0,2);
		//36 - 49 - gen jet pt vs reco jet pt (AK4)
		TH2D* recoAK4JetPt_genAK4JetPt = new TH2D("recoAK4JetPt_genAK4JetPt","recoAK4JetPt_genAK4JetPt;recoAK4JetPt;genAK4JetPt",25,0,500,25,0,500);
		//37 - 50 - gen jet E vs reco jet E (AK4)
		TH2D* recoAK4JetE_genAK4JetE = new TH2D("recoAK4JetE_genAK4JetE","recoAK4JetE_genAK4JetE;recoAK4JetE;genAK4JetE",25,0,2000,25,0,2000);
		//38 - 51 - gen jet mass vs reco jet mass (AK4)
		TH2D* recoAK4JetMass_genAK4JetMass = new TH2D("recoAK4JetMass_genAK4JetMass","recoAK4JetMass_genAK4JetMass;recoAK4JetMass;genAK4JetMass",25,0,250,25,0,250);
		//39 - 52 - gen jet eta vs reco jet eta (AK4)
		TH2D* recoAK4JetEtaCenter_genAK4JetEtaCenter = new TH2D("recoAK4JetEtaCenter_genAK4JetEtaCenter","recoAK4JetEtaCenter_genAK4JetEtaCenter;recoAK4JetEtaCenter;genAK4JetEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//40 - 53 - gen jet phi vs reco jet phi (AK4)
		TH2D* recoAK4JetPhiCenter_genAK4JetPhiCenter = new TH2D("recoAK4JetPhiCenter_genAK4JetPhiCenter","recoAK4JetPhiCenter_genAK4JetPhiCenter;recoAK4JetPhiCenter;genAK4JetPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//41 - 54 - gen jet time vs reco jet time (AK4)
		TH2D* recoAK4JetTimeCenter_genAK4JetTimeCenter = new TH2D("recoAK4JetTimeCenter_genAK4JetTimeCenter","recoAK4JetTimeCenter_genAK4JetTimeCenter;recoAK4JetTimeCenter;genAK4JetTimeCenter",25,-1,1,25,-1,1.);
		//42 - 55 - gen AK4 jet pt vs BHC jet pt 
		TH2D* BHCJetPt_genAK4JetPt = new TH2D("BHCJetPt_genAK4JetPt","BHCJetPt_genAK4JetPt;BHCJetPt;genAK4JetPt",25,0,500,25,0,500);
		//43 - 56 - gen AK4 jet E vs BHC jet E 
		TH2D* BHCJetE_genAK4JetE = new TH2D("BHCJetE_genAK4JetE","BHCJetE_genAK4JetE;BHCJetE;genAK4JetE",25,0,2000,25,0,2000);
		//44 - 57 - gen AK4 jet mass vs BHC jet mass 
		TH2D* BHCJetMass_genAK4JetMass = new TH2D("BHCJetMass_genAK4JetMass","BHCJetMass_genAK4JetMass;BHCJetMass;genAK4JetMass",25,0,250,25,0,250);
		//45 - 58 - gen AK4 jet eta vs BHC jet eta 
		TH2D* BHCJetEtaCenter_genAK4JetEtaCenter = new TH2D("BHCJetEtaCenter_genAK4JetEtaCenter","BHCJetEtaCenter_genAK4JetEtaCenter;BHCJetEtaCenter;genAK4JetEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//46 - 59 - gen AK4 jet phi vs BHC jet phi 
		TH2D* BHCJetPhiCenter_genAK4JetPhiCenter = new TH2D("BHCJetPhiCenter_genAK4JetPhiCenter","BHCJetPhiCenter_genAK4JetPhiCenter;BHCJetPhiCenter;genAK4JetPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//47 - 60 - gen AK4 jet time vs BHC jet time 
		TH2D* BHCJetTimeCenter_genAK4JetTimeCenter = new TH2D("BHCJetTimeCenter_genAK4JetTimeCenter","BHCJetTimeCenter_genAK4JetTimeCenter;BHCJetTimeCenter;genAK4JetTimeCenter",25,-1,1,25,-1,1.);
		//48 - 61 - gen AK4 jet # constituents vs BHC jet # subclusters
		TH2D* BHCnSubclusters_genAK4nConstituents = new TH2D("BHCnSubclusters_genAK4nConstituents","BHCnSubclusters_genAK4nConstituents;BHCnSubclusters;genAK4nConstituents",30,0,30,30,0,30);
		//49 - 62 - gen AK15 jet pt vs BHC jet pt 
		TH2D* BHCJetPt_genAK15JetPt = new TH2D("BHCJetPt_genAK15JetPt","BHCJetPt_genAK15JetPt;BHCJetPt;genAK15JetPt",25,0,500,25,0,500);
		//50 - 63 - gen AK15 jet E vs BHC jet E 
		TH2D* BHCJetE_genAK15JetE = new TH2D("BHCJetE_genAK15JetE","BHCJetE_genAK15JetE;BHCJetE;genAK15JetE",25,0,2000,25,0,2000);
		//51 - 64 - gen AK15 jet mass vs BHC jet mass 
		TH2D* BHCJetMass_genAK15JetMass = new TH2D("BHCJetMass_genAK15JetMass","BHCJetMass_genAK15JetMass;BHCJetMass;genAK15JetMass",25,0,250,25,0,250);
		//52 - 65 - gen AK15 jet eta vs BHC jet eta 
		TH2D* BHCJetEtaCenter_genAK15JetEtaCenter = new TH2D("BHCJetEtaCenter_genAK15JetEtaCenter","BHCJetEtaCenter_genAK15JetEtaCenter;BHCJetEtaCenter;genAK15JetEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//53 - 66 - gen AK15 jet phi vs BHC jet phi 
		TH2D* BHCJetPhiCenter_genAK15JetPhiCenter = new TH2D("BHCJetPhiCenter_genAK15JetPhiCenter","BHCJetPhiCenter_genAK15JetPhiCenter;BHCJetPhiCenter;genAK15JetPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//54 - 67 - gen AK15 jet time vs BHC jet time 
		TH2D* BHCJetTimeCenter_genAK15JetTimeCenter = new TH2D("BHCJetTimeCenter_genAK15JetTimeCenter","BHCJetTimeCenter_genAK15JetTimeCenter;BHCJetTimeCenter;genAK15JetTimeCenter",25,-1,1,25,-1,1.);
		//55 - 68 - gen AK15 jet # constituents vs BHC jet # subclusters 
		TH2D* BHCnSubclusters_genAK15nConstituents = new TH2D("BHCnSubclusters_genAK15nConstituents","BHCnSubclusters_genAK15nConstituents;BHCnSubclusters;genAK15nConstituents",30,0,30,30,0,30);
		//56 - 69 - BHC jet pt vs gen top pt 
		TH2D* BHCJetPt_genTopPt = new TH2D("BHCJetPt_genTopPt","BHCJetPt_genTopPt;BHCJetPt;genTopPt",25,0,2000,25,0,2000);
		//57 - 70 - BHC jet E vs gen top pt
		TH2D* BHCJetE_genTopE = new TH2D("BHCJetE_genTopE","BHCJetE_genTopE;BHCJetE;genTopE",25,0,2000,25,0,2000);
		//58 - 71 - BHC jet mass vs gen top jet mass 
		TH2D* BHCJetMass_genTopMass = new TH2D("BHCJetMass_genTopMass","BHCJetMass_genTopMass;BHCJetMass;genTopMass",25,0,250,25,0,250);
		//59 - 72 - BHC jet eta vs gen top eta 
		TH2D* BHCJetEtaCenter_genTopEtaCenter = new TH2D("BHCJetEtaCenter_genTopEtaCenter","BHCJetEtaCenter_genTopEtaCenter;BHCJetEtaCenter;genTopEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//60 - 73 - BHC jet phi vs gen top phi 
		TH2D* BHCJetPhiCenter_genTopPhiCenter = new TH2D("BHCJetPhiCenter_genTopPhiCenter","BHCJetPhiCenter_genTopPhiCenter;BHCJetPhiCenter;genTopPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//61 - 74 - recoAK4 jet pt vs gen top pt 
		TH2D* recoAK4JetPt_genTopPt = new TH2D("recoAK4JetPt_genTopPt","recoAK4JetPt_genTopPt;recoAK4JetPt;genTopPt",25,0,2000,25,0,2000);
		//62 - 75 - recoAK4 jet E vs gen top pt
		TH2D* recoAK4JetE_genTopE = new TH2D("recoAK4JetE_genTopE","recoAK4JetE_genTopE;recoAK4JetE;genTopE",25,0,2000,25,0,2000);
		//63 - 76 - recoAK4 jet mass vs gen top jet mass 
		TH2D* recoAK4JetMass_genTopMass = new TH2D("recoAK4JetMass_genTopMass","recoAK4JetMass_genTopMass;recoAK4JetMass;genTopMass",25,0,250,25,0,250);
		//64 - 77 - recoAK4 jet eta vs gen top eta 
		TH2D* recoAK4JetEtaCenter_genTopEtaCenter = new TH2D("recoAK4JetEtaCenter_genTopEtaCenter","recoAK4JetEtaCenter_genTopEtaCenter;recoAK4JetEtaCenter;genTopEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//65 - 78 - recoAK4 jet phi vs gen top phi 
		TH2D* recoAK4JetPhiCenter_genTopPhiCenter = new TH2D("recoAK4JetPhiCenter_genTopPhiCenter","recoAK4JetPhiCenter_genTopPhiCenter;recoAK4JetPhiCenter;genTopPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//66 - 79 - recoAK15 jet pt vs gen top pt 
		TH2D* recoAK15JetPt_genTopPt = new TH2D("recoAK15JetPt_genTopPt","recoAK15JetPt_genTopPt;recoAK15JetPt;genTopPt",25,0,2000,25,0,2000);
		//67 - 80 - recoAK15 jet E vs gen top pt
		TH2D* recoAK15JetE_genTopE = new TH2D("recoAK15JetE_genTopE","recoAK15JetE_genTopE;recoAK15JetE;genTopE",25,0,2000,25,0,2000);
		//68 - 81 - recoAK15 jet mass vs gen top jet mass 
		TH2D* recoAK15JetMass_genTopMass = new TH2D("recoAK15JetMass_genTopMass","recoAK15JetMass_genTopMass;recoAK15JetMass;genTopMass",25,0,250,25,0,250);
		//69 - 82 - recoAK15 jet eta vs gen top eta 
		TH2D* recoAK15JetEtaCenter_genTopEtaCenter = new TH2D("recoAK15JetEtaCenter_genTopEtaCenter","recoAK15JetEtaCenter_genTopEtaCenter;recoAK15JetEtaCenter;genTopEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//70 - 83 - recoAK15 jet phi vs gen top phi 
		TH2D* recoAK15JetPhiCenter_genTopPhiCenter = new TH2D("recoAK15JetPhiCenter_genTopPhiCenter","recoAK15JetPhiCenter_genTopPhiCenter;recoAK15JetPhiCenter;genTopPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//71 - 84 - genAK4 jet pt vs gen top pt 
		TH2D* genAK4JetPt_genTopPt = new TH2D("genAK4JetPt_genTopPt","genAK4JetPt_genTopPt;genAK4JetPt;genTopPt",25,0,2000,25,0,2000);
		//72 - 85 - genAK4 jet E vs gen top pt
		TH2D* genAK4JetE_genTopE = new TH2D("genAK4JetE_genTopE","genAK4JetE_genTopE;genAK4JetE;genTopE",25,0,2000,25,0,2000);
		//73 - genAK4 jet mass vs gen top jet mass 
		TH2D* genAK4JetMass_genTopMass = new TH2D("genAK4JetMass_genTopMass","genAK4JetMass_genTopMass;genAK4JetMass;genTopMass",25,0,250,25,0,250);
		//74 - 87 - genAK4 jet eta vs gen top eta 
		TH2D* genAK4JetEtaCenter_genTopEtaCenter = new TH2D("genAK4JetEtaCenter_genTopEtaCenter","genAK4JetEtaCenter_genTopEtaCenter;genAK4JetEtaCenter;genTopEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//75 - 88 - genAK4 jet phi vs gen top phi 
		TH2D* genAK4JetPhiCenter_genTopPhiCenter = new TH2D("genAK4JetPhiCenter_genTopPhiCenter","genAK4JetPhiCenter_genTopPhiCenter;genAK4JetPhiCenter;genTopPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//76 - 89 - genAK15 jet pt vs gen top pt 
		TH2D* genAK15JetPt_genTopPt = new TH2D("genAK15JetPt_genTopPt","genAK15JetPt_genTopPt;genAK15JetPt;genTopPt",25,0,2000,25,0,2000);
		//77 - 90 - genAK15 jet E vs gen top pt
		TH2D* genAK15JetE_genTopE = new TH2D("genAK15JetE_genTopE","genAK15JetE_genTopE;genAK15JetE;genTopE",25,0,2000,25,0,2000);
		//78 - 91 - genAK15 jet mass vs gen top jet mass 
		TH2D* genAK15JetMass_genTopMass = new TH2D("genAK15JetMass_genTopMass","genAK15JetMass_genTopMass;genAK15JetMass;genTopMass",25,0,250,25,0,250);
		//79 - 92 - genAK15 jet eta vs gen top eta 
		TH2D* genAK15JetEtaCenter_genTopEtaCenter = new TH2D("genAK15JetEtaCenter_genTopEtaCenter","genAK15JetEtaCenter_genTopEtaCenter;genAK15JetEtaCenter;genTopEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//80 - 93 - genAK15 jet phi vs gen top phi 
		TH2D* genAK15JetPhiCenter_genTopPhiCenter = new TH2D("genAK15JetPhiCenter_genTopPhiCenter","genAK15JetPhiCenter_genTopPhiCenter;genAK15JetPhiCenter;genTopPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));

		//81 - 94 - BHC jet pt vs gen W pt 
		TH2D* BHCJetPt_genWPt = new TH2D("BHCJetPt_genWPt","BHCJetPt_genWPt;BHCJetPt;genWPt",25,0,500,25,0,500);
		//82 - 95 - BHC jet E vs gen W pt
		TH2D* BHCJetE_genWE = new TH2D("BHCJetE_genWE","BHCJetE_genWE;BHCJetE;genWE",25,0,2000,25,0,2000);
		//83 - 96 - BHC jet mass vs gen W jet mass 
		TH2D* BHCJetMass_genWMass = new TH2D("BHCJetMass_genWMass","BHCJetMass_genWMass;BHCJetMass;genWMass",25,0,250,25,0,250);
		//84 - 97 - BHC jet eta vs gen W eta 
		TH2D* BHCJetEtaCenter_genWEtaCenter = new TH2D("BHCJetEtaCenter_genWEtaCenter","BHCJetEtaCenter_genWEtaCenter;BHCJetEtaCenter;genWEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//85 - 98 - BHC jet phi vs gen W phi 
		TH2D* BHCJetPhiCenter_genWPhiCenter = new TH2D("BHCJetPhiCenter_genWPhiCenter","BHCJetPhiCenter_genWPhiCenter;BHCJetPhiCenter;genWPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//86 - recoAK4 jet pt vs gen W pt 
		TH2D* recoAK4JetPt_genWPt = new TH2D("recoAK4JetPt_genWPt","recoAK4JetPt_genWPt;recoAK4JetPt;genWPt",25,0,500,25,0,500);
		//87 - 100 - recoAK4 jet E vs gen W pt
		TH2D* recoAK4JetE_genWE = new TH2D("recoAK4JetE_genWE","recoAK4JetE_genWE;recoAK4JetE;genWE",25,0,2000,25,0,2000);
		//88 - 101 - recoAK4 jet mass vs gen W jet mass 
		TH2D* recoAK4JetMass_genWMass = new TH2D("recoAK4JetMass_genWMass","recoAK4JetMass_genWMass;recoAK4JetMass;genWMass",25,0,250,25,0,250);
		//89 - 102 - recoAK4 jet eta vs gen W eta 
		TH2D* recoAK4JetEtaCenter_genWEtaCenter = new TH2D("recoAK4JetEtaCenter_genWEtaCenter","recoAK4JetEtaCenter_genWEtaCenter;recoAK4JetEtaCenter;genWEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//90 - 103 - recoAK4 jet phi vs gen W phi 
		TH2D* recoAK4JetPhiCenter_genWPhiCenter = new TH2D("recoAK4JetPhiCenter_genWPhiCenter","recoAK4JetPhiCenter_genWPhiCenter;recoAK4JetPhiCenter;genWPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//91 - 104 - recoAK15 jet pt vs gen W pt 
		TH2D* recoAK15JetPt_genWPt = new TH2D("recoAK15JetPt_genWPt","recoAK15JetPt_genWPt;recoAK15JetPt;genWPt",25,0,500,25,0,500);
		//92 - 105 - recoAK15 jet E vs gen W pt
		TH2D* recoAK15JetE_genWE = new TH2D("recoAK15JetE_genWE","recoAK15JetE_genWE;recoAK15JetE;genWE",25,0,2000,25,0,2000);
		//93 - 106 - recoAK15 jet mass vs gen W jet mass 
		TH2D* recoAK15JetMass_genWMass = new TH2D("recoAK15JetMass_genWMass","recoAK15JetMass_genWMass;recoAK15JetMass;genWMass",25,0,250,25,0,250);
		//94 - 107 - recoAK15 jet eta vs gen W eta 
		TH2D* recoAK15JetEtaCenter_genWEtaCenter = new TH2D("recoAK15JetEtaCenter_genWEtaCenter","recoAK15JetEtaCenter_genWEtaCenter;recoAK15JetEtaCenter;genWEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//95 - 108 - recoAK15 jet phi vs gen W phi 
		TH2D* recoAK15JetPhiCenter_genWPhiCenter = new TH2D("recoAK15JetPhiCenter_genWPhiCenter","recoAK15JetPhiCenter_genWPhiCenter;recoAK15JetPhiCenter;genWPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//96 - 109 - genAK4 jet pt vs gen W pt 
		TH2D* genAK4JetPt_genWPt = new TH2D("genAK4JetPt_genWPt","genAK4JetPt_genWPt;genAK4JetPt;genWPt",25,0,500,25,0,500);
		//97 - 110 - genAK4 jet E vs gen W pt
		TH2D* genAK4JetE_genWE = new TH2D("genAK4JetE_genWE","genAK4JetE_genWE;genAK4JetE;genWE",25,0,2000,25,0,2000);
		//98 - 111 - genAK4 jet mass vs gen W jet mass 
		TH2D* genAK4JetMass_genWMass = new TH2D("genAK4JetMass_genWMass","genAK4JetMass_genWMass;genAK4JetMass;genWMass",25,0,250,25,0,250);
		//99 - 112 - genAK4 jet eta vs gen W eta 
		TH2D* genAK4JetEtaCenter_genWEtaCenter = new TH2D("genAK4JetEtaCenter_genWEtaCenter","genAK4JetEtaCenter_genWEtaCenter;genAK4JetEtaCenter;genWEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//100 - 113 - genAK4 jet phi vs gen W phi 
		TH2D* genAK4JetPhiCenter_genWPhiCenter = new TH2D("genAK4JetPhiCenter_genWPhiCenter","genAK4JetPhiCenter_genWPhiCenter;genAK4JetPhiCenter;genWPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//101 - 114 - genAK15 jet pt vs gen W pt 
		TH2D* genAK15JetPt_genWPt = new TH2D("genAK15JetPt_genWPt","genAK15JetPt_genWPt;genAK15JetPt;genWPt",25,0,500,25,0,500);
		//102 - 115 - genAK15 jet E vs gen W pt
		TH2D* genAK15JetE_genWE = new TH2D("genAK15JetE_genWE","genAK15JetE_genWE;genAK15JetE;genWE",25,0,2000,25,0,2000);
		//103 - 116 - genAK15 jet mass vs gen W jet mass 
		TH2D* genAK15JetMass_genWMass = new TH2D("genAK15JetMass_genWMass","genAK15JetMass_genWMass;genAK15JetMass;genWMass",25,0,250,25,0,250);
		//104 - 117 - genAK15 jet eta vs gen W eta 
		TH2D* genAK15JetEtaCenter_genWEtaCenter = new TH2D("genAK15JetEtaCenter_genWEtaCenter","genAK15JetEtaCenter_genWEtaCenter;genAK15JetEtaCenter;genWEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//105 - 118 - genAK15 jet phi vs gen W phi 
		TH2D* genAK15JetPhiCenter_genWPhiCenter = new TH2D("genAK15JetPhiCenter_genWPhiCenter","genAK15JetPhiCenter_genWPhiCenter;genAK15JetPhiCenter;genWPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//106 - 119 - BHC jet subcluster energy vs subcluster mass
		TH2D* BHCJet_subclusterEnergy_subclusterMass = new TH2D("BHCJet_subclusterEnergy_subclusterMass","BHCJet_subclusterEnergy_subclusterMass;subclusterEnergy;subclusterMass",25,0,500,25,0,100);
		//107 - 120 - BHC jet subcluster energy vs # effective rec hits
		TH2D* BHCJet_subclusterEnergy_subclusterEffnRhs = new TH2D("BHCJet_subclusterEnergy_subclusterEffnRhs","BHCJet_subclusterEnergy_subclusterEffnRhs;subclusterEnergy;subclusterEffnRhs",25,0,500,25,0,500);
		//108 - BHC jet subcluster mass vs # effective rechits 	
		TH2D* BHCJet_subclusterMass_subclusterEffnRhs = new TH2D("BHCJet_subclusterMass_subclusterEffnRhs","BHCJet_subclusterMass_subclusterEffnRhs;subclusterMass;subclusterEffnRhs",25,0,100,25,0,200);
		//109 - reco AK4 jet subcluster energy vs subcluster mass
		TH2D* recoAK4Jet_subclusterEnergy_subclusterMass = new TH2D("recoAK4Jet_subclusterEnergy_subclusterMass","recoAK4Jet_subclusterEnergy_subclusterMass;subclusterEnergy;subclusterMass",25,0,500,25,0,100);
		//110 - recoAK4 jet subcluster energy vs # effective rec hits
		TH2D* recoAK4Jet_subclusterEnergy_subclusterEffnRhs = new TH2D("recoAK4Jet_subclusterEnergy_subclusterEffnRhs","recoAK4Jet_subclusterEnergy_subclusterEffnRhs;subclusterEnergy;subclusterEffnRhs",25,0,500,25,0,200);
		//111 - recoAK4 jet subcluster mass vs # effective rechits 	
		TH2D* recoAK4Jet_subclusterMass_subclusterEffnRhs = new TH2D("recoAK4Jet_subclusterMass_subclusterEffnRhs","recoAK4Jet_subclusterMass_subclusterEffnRhs;subclusterMass;subclusterEffnRhs",25,0,100,25,0,200);
		//112 - subcluster energy vs dR of subcluster to jet center
		TH2D* BHCJet_subclusterEnergy_subclusterdRToJet = new TH2D("BHCJet_subclusterEnergy_subclusterdRToJet","BHCJet_subclusterEnergy_subclusterdRToJet;subclusterEnergy;subclusterdRToJet",25,0,500,25,0,0.5);
		//113 - subcluster eff # rhs vs dR of subcluster to jet center
		TH2D* BHCJet_subclusterEffnRhs_subclusterdRToJet = new TH2D("BHCJet_subclusterEffnRhs_subclusterdRToJet","BHCJet_subclusterEffnRhs_subclusterdRToJet;subclusterEffnRhs;subclusterdRToJet",25,0,100,25,0,0.5);
		//114 - BHC jets gen-matched to Ws - subcluster mass vs # subclusters/jet
		TH2D* BHCJetW_nSubclustersJet_mass = new TH2D("BHCJetW_nSubclustersJet_mass","BHCJetW_nSubclustersJet_mass;nSubclustersJet;mass",10,0,10,50,0,250);
		//115 - BHC jets gen-matched to Ws - dR bw gen partons of W vs # subclusters/jet
		TH2D* BHCJetW_dRGenPartons_nSubclustersJet = new TH2D("BHCJetW_dRGenPartons_nSubclustersJet","BHCJetW_dRGenPartons_nSubclustersJet;dRGenPartons;nSubclustersJet",50,0,2.,10,0,10);
		//116 - BHC jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* BHCJetW_dRGenPartons_jetSize = new TH2D("BHCJetW_dRGenPartons_jetSize","BHCJetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//117 - BHC jets gen-matched to Ws - Eratio of subcl E/gen parton E vs # subclusters/jet
		TH2D* BHCJetW_EratioSubclGenPart_nSubclustersJet = new TH2D("BHCJetW_EratioSubclGenPart_nSubclustersJet","BHCJetW_EratioSubclGenPart_nSubclustersJet;EratioSubclGenPart;nSubclustersJet",50,0,2.,10,0,10);
		//118 - BHC jets gen-matched to Tops - subcluster mass vs # subclusters/jet
		TH2D* BHCJetTop_nSubclustersJet_mass = new TH2D("BHCJetTop_nSubclustersJet_mass","BHCJetTop_nSubclustersJet_mass;nSubclustersJet;subclMass",10,0,10,50,0,250);
		//119 - # bhc jets as a function of opening angle of partons from W
		TH2D* BHCJetW_openAng_nJets = new TH2D("BHCJetW_openAng_nJets","BHCJetW_openAng_nJets;openAng;nJets",25,0,3.3,30,0,30);
		//120 - # bhc subclusters/jet as a function of opening angle of partons from W
		TH2D* BHCJetW_openAng_nSubclustersJet = new TH2D("BHCJetW_openAng_nSubclustersJet","BHCJetW_openAng_nSubclustersJet;openAng;nSubclustersJet",25,0,3.3,10,0,10);
		//121 - bhc subcluster mass as a function of opening angle of partons from W
		TH2D* BHCJetW_openAng_subclMass = new TH2D("BHCJetW_openAng_subclMass","BHCJetW_openAng_subclMass;openAng;subclMass",25,0,3.3,50,0,200);
		//122 - BHC jets gen-matched to Ws with exactly 1 subcluster - dR bw gen partons of W vs jet size
		TH2D* BHCJetW_1subcl_dRGenPartons_jetSize = new TH2D("BHCJetW_1subcl_dRGenPartons_jetSize","BHCJetW_1subcl_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//123 - BHC jets gen-matched to Ws with 2+ subclusters - dR bw gen partons of W vs jet size
		TH2D* BHCJetW_ge2subcl_dRGenPartons_jetSize = new TH2D("BHCJetW_ge2subcl_dRGenPartons_jetSize","BHCJetW_ge2subcl_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//124 - BHC jets gen-matched to Ws with exactly 1 subcluster - dR bw gen partons of W vs avg. parton energy
		TH2D* BHCJetW_1subcl_dRGenPartons_avgPartonEnergy = new TH2D("BHCJetW_1subcl_dRGenPartons_avgPartonEnergy","BHCJetW_1subcl_dRGenPartons_avgPartonEnergy;dRGenPartons;avgPartonEnergy",50,0,2.,50,0,1000.);
		//125 - BHC jets gen-matched to Ws with 2+ subclusters - dR bw gen partons of W vs avg. parton energy
		TH2D* BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy = new TH2D("BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy","BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy;dRGenPartons;avgPartonEnergy",50,0,2.,50,0,1000.);
		//126 - BHC jets gen-matched to Ws with 2+ subclusters - subcluster energy vs lead index of subcluster (ie 0 = lead, 1 = sublead, etc)
		TH2D* BHCJetW_ge2subcl_subclEnergy_subclLeadIdx = new TH2D("BHCJetW_ge2subcl_subclEnergy_subclLeadIdx","BHCJetW_ge2subcl_subclEnergy_subclLeadIdx;subclEnergy;subclLeadIdx",5,0,5,25,0,500);
		//127 - BHC jets gen-matched to Ws - Eratio of jet E/gen W E vs # subclusters/jet
		TH2D* BHCJetW_EratioJetGenW_nSubclustersJet = new TH2D("BHCJetW_ge2subcl_EratioJetGenW_nSubclustersJet","BHCJetW_EratioJetGenW_nSubclustersJet;EratioJetGenW;nSubclustersJet",50,0,2.,10,0,10);
		//128 - BHC jets - jet mass vs jet size
		TH2D* BHCJet_jetMass_jetSize = new TH2D("BHCJet_jetMass_jetSize","BHCJet_jetMass_jetSize;jetMass;jetSize",50,0,250.,50,0,2.);
		//129 - eta-phi event display of rechits for specified _evt2disp with cell energy on the z axis (overall event)
		TH2D* EvtDisplay_etaCell_phiCell = new TH2D("EvtDisplay_etaCell_phiCell","EvtDisplay_etaCell_phiCell;eta;phi;energy",344,-3,3,360,0,8*atan(1));
		//130 - reco AK8 jet mass vs reco jet pt
		TH2D* recoAK8JetMass_recoAK8JetSize = new TH2D("recoAK8Jet_jetMass_jetSize","recoAK8Jet_jetMass_jetSize;recoAK8JetMass;recoAK8JetSize",50,0,250,50,0,2.);
		//131 - reco AK15 jet mass vs reco jet pt
		TH2D* recoAK15JetMass_recoAK15JetSize = new TH2D("recoAK15Jet_jetMass_jetSize","recoAK15Jet_jetMass_jetSize;recoAK15JetMass;recoAK15JetSize",50,0,250,50,0,2.);
		//132 - matched reco AK8 jet mass vs BHC jet mass
		TH2D* recoAK8JetMass_BHCJetMass_matched = new TH2D("recoAK8JetMass_BHCJetMass_BHCtoAK8matched","recoAK8JetMass_BHCJetMass_BHCtoAK8matched;recoAK8 Jet Mass;BHC Jet Mass",25,0,250,25,0,250);
		//133 - recoAK8 jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* recoAK8JetW_dRGenPartons_jetSize = new TH2D("recoAK8JetW_dRGenPartons_jetSize","recoAK8JetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//134 - BHC jet subcluster time sig vs subcluster size
		TH2D* BHCJet_subclTimeSig_subclSize = new TH2D("BHCJet_subclTimeSig_subclSize","BHCJet_subclTimeSig_subclSize;subclTimeSig;subclSize",50,0,15,50,0,1.5);
		//135 - BHC jet subcluster time sig vs subcluster pt 
		TH2D* BHCJet_subclTimeSig_subclPt = new TH2D("BHCJet_subclTimeSig_subclPt","BHCJet_subclTimeSig_subclPt;subclTimeSig;subclPt",50,0,15,50,0,500);
		//136 - recoAK15 jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* recoAK15JetW_dRGenPartons_jetSize = new TH2D("recoAK15JetW_dRGenPartons_jetSize","recoAK15JetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//137 - recoAK4 jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* recoAK4JetW_dRGenPartons_jetSize = new TH2D("recoAK4JetW_dRGenPartons_jetSize","recoAK4JetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//138 - high mass + W-matched BHC jets - pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchSubclPt_subclTime = new TH2D("BHCJetW_highMass_partonMatchSubclPt_subclTime","BHCJetW_highMass_partonMatchSubclPt_subclTime;SubclPt;SubclTime",25,0,500,25,-3,3);
		//139 - high mass + W-matched BHC jets - pt of subclusters NOT gen-matched W partons / pt jet vs subcl time
		TH2D* BHCJetW_highMass_partonNoMatchSubclPt_subclTime = new TH2D("BHCJetW_highMass_partonNoMatchSubclPt_subclTime","BHCJetW_highMass_partonNoMatchSubclPt_subclTime;SubclPt;SubclTime",25,0,500,25,-3,3);
		//140 - high mass + W-matched BHC jets - subclSize of subclusters gen-matched to W partons / size jet vs subcl time
		TH2D* BHCJetW_highMass_partonMatchSubclSize_subclTime = new TH2D("BHCJetW_highMass_partonMatchSubclSize_subclTime","BHCJetW_highMass_partonMatchSubclSize_subclTime;SubclSize;SubclTime",25,0,1.5,25,-3,3);
		//141 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched W partons / size jet vs subcl time
		TH2D* BHCJetW_highMass_partonNoMatchSubclSize_subclTime = new TH2D("BHCJetW_highMass_partonNoMatchSubclSize_subclTime","BHCJetW_highMass_partonNoMatchSubclSize_subclTime;SubclSize;SubclTime",25,0,1.5,25,-3,3);
		//142 - high mass + W-matched BHC jets - subclSize of subclsuters gen-matched to W partons / size jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt = new TH2D("BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt","BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt;RelSubclSize;RelSubclPt",25,0,5,25,0,1.2);
		//143 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched to W partons / size jet vs pt of subclusters NOT gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt","BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt;RelSubclSize;RelSubclPt",25,0,5,25,0,1.2);
		//144 - high mass + W-matched BHC jets - subclMass of subclusters gen-matched to W partons / mass jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclMass_RelSubclPt = new TH2D("BHCJetW_highMass_partonMatchRelSubclMass_RelSubclPt","BHCJetW_highMass_partonMatchRelSubclMass_RelSubclPt;RelSubclMass;RelSubclPt",25,0,1.2,25,0,1.2);
		//145 - high mass + W-matched BHC jets - subclMass of subclusters NOT gen-matched to W partons / mass jet vs pt of subclusters NOT gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclMass_RelSubclPt = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclMass_RelSubclPt","BHCJetW_highMass_partonNoMatchRelSubclMass_RelSubclPt;RelSubclMass;RelSubclPt",25,0,1.2,25,0,1.2);
		//146 - high mass + W-matched BHC jets - subcl time var of subclusters gen-matched to W partons / time var jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclPt = new TH2D("BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclPt","BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclPt;RelSubclTimeVar;RelSubclPt",25,0,5,25,0,1.2);
		//147 - high mass + W-matched BHC jets - subcl time var of subclusters NOT gen-matched to W partons / time var jet vs pt of subclusters NOT gen-matched to W partons / pt jet
		TH2D* BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclPt = new TH2D("BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclPt","BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclPt;RelSubclTimeVar;RelSubclPt",25,0,5,25,0,1.2);
		//148 - high mass + W-matched BHC jets - subcl time var of subclusters gen-matched to W partons / time var jet vs size of subclusters gen-matched to W partons / size jet
		TH2D* BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclSize = new TH2D("BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclSize","BHCJetW_highTimeVar_partonMatchRelSubclTimeVar_RelSubclSize;RelSubclTimeVar;RelSubclSize",25,0,5,25,0,5);
		//149 - high mass + W-matched BHC jets - subcl time var of subclusters NOT gen-matched to W partons / time var jet vs size of subclusters NOT gen-matched to W partons / size jet
		TH2D* BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclSize = new TH2D("BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclSize","BHCJetW_highTimeVar_partonNoMatchRelSubclTimeVar_RelSubclSize;RelSubclTimeVar;RelSubclSize",25,0,5,25,0,5);
		//150 - high mass + W-matched BHC jets - matched subcl mass / pt vs jet mass / pt 
		TH2D* BHCJetW_highMassOvPt_partonMatchSubclMassOvPt_JetMassOvPt = new TH2D("BHCJetW_highMassOvPt_partonMatchSubclMassOvPt_JetMassOvPt","BHCJetW_highMassOvPt_partonMatchSubclMassOvPt_JetMassOvPt;SubclMassOvPt;JetMassOvPt",25,0,50,25,0,50);
		//151 - high mass + W-matched BHC jets - NOT matched subcl mass / pt vs jet mass / pt 
		TH2D* BHCJetW_highMassOvPt_partonNoMatchSubclMassOvPt_JetMassOvPt = new TH2D("BHCJetW_highMassOvPt_partonNoMatchSubclMassOvPt_JetMassOvPt","BHCJetW_highMassOvPt_partonNoMatchSubclMassOvPt_JetMassOvPt;SubclMassOvPt;JetMassOvPt",25,0,50,25,0,50);
		//152 - high mass + W-matched BHC jets - subclSize (space + time) of subclsuters gen-matched to W partons / size (space + time) jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclSizeSpatioTemp_RelSubclPt = new TH2D("BHCJetW_highMass_partonMatchRelSubclSizeSpatioTemp_RelSubclPt","BHCJetW_highMass_partonMatchRelSubclSizeSpatioTemp_RelSubclPt;RelSubclSizeSpatioTemp;RelSubclPt",25,0,5,25,0,1.2);
		//153 - high mass + W-matched BHC jets - subclSize (space + time) of subclusters NOT gen-matched to W partons / size (space + time) jet vs pt of subclusters NOT gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclSizeSpatioTemp_RelSubclPt = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclSizeSpatioTemp_RelSubclPt","BHCJetW_highMass_partonNoMatchRelSubclSizeSpatioTemp_RelSubclPt;RelSubclSizeSpatioTemp;RelSubclPt",25,0,5,25,0,1.2);





		//event display histograms (per gen object)
		vector<TH2D*> _evtdisps_obj;
		//these should be centered on their respective gen particles
		//0 - single W, W+gluon, first W in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_W = new TH2D("EvtDisplay_etaCell_phiCell_W","EvtDisplay_etaCell_phiCell_W;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//1 - second W in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_W2 = new TH2D("EvtDisplay_etaCell_phiCell_W2","EvtDisplay_etaCell_phiCell_W2;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//2 - gluon in W+gluon
		TH2D* EvtDisplay_etaCell_phiCell_gluon = new TH2D("EvtDisplay_etaCell_phiCell_gluon","EvtDisplay_etaCell_phiCell_gluon;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//3 - q1 in QCD dijets
		TH2D* EvtDisplay_etaCell_phiCell_q1 = new TH2D("EvtDisplay_etaCell_phiCell_q1","EvtDisplay_etaCell_phiCell_q1;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//4 - q2 in QCD dijets
		TH2D* EvtDisplay_etaCell_phiCell_q2 = new TH2D("EvtDisplay_etaCell_phiCell_q2","EvtDisplay_etaCell_phiCell_q2;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//5 - b1 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_b1 = new TH2D("EvtDisplay_etaCell_phiCell_b1","EvtDisplay_etaCell_phiCell_b1;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//6 - b2 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_b2 = new TH2D("EvtDisplay_etaCell_phiCell_b2","EvtDisplay_etaCell_phiCell_b2;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//7 - top1 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_top1 = new TH2D("EvtDisplay_etaCell_phiCell_top1","EvtDisplay_etaCell_phiCell_top1;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//8 - top2 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_top2 = new TH2D("EvtDisplay_etaCell_phiCell_top2","EvtDisplay_etaCell_phiCell_top2;#eta;#phi;Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));


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

	//find gen object that most closely is dr matched to jet
	//when gen parts are used, need to specify id to match to (ie W - 24, top - 6, any parton/quark/gluon - 0)
	void GenericMatchJet(vector<Jet>& injets, vector<Jet>& matchjets, vector<int>& bestMatchIdxs, int id = -1){
		//loop through gen particles
		//dr match to jet
		double bestDr, dr;
		vector<Jet> matchobjs;
		//map 'local' match obj idxs to 'global' matchjets idxs
		map<int,int> matchIdxToObjIdx;
		//if id is specified, make matchobjs the genparts that satisfy this criteria
		if(id != -1){
			vector<int> qids = {1,2,3,4,5,21};
			for(int g = 0; g < matchjets.size(); g++){
				int idx = matchjets[g].GetUserIdx();
				int this_id = fabs(_base->genpart_id->at(idx));
				if(id != 0 && this_id != id) continue;
				if(id == 0){
					if(find(qids.begin(), qids.end(), this_id) == qids.end()) continue;
				}
				matchobjs.push_back(matchjets[g]);
				matchIdxToObjIdx[matchobjs.size()-1] = g;
			}
			//add no match option
			matchIdxToObjIdx[-1] = -1;
		}
		else{ matchobjs = matchjets; }

		int nMatch = matchobjs.size();
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
			
				dr = dR(matchobjs[g].eta(), matchobjs[g].phi(), injets[j].eta(), injets[j].phi());
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
		//map 'local' match obj index to 'global' all obj index
		if(!matchIdxToObjIdx.empty()){
			for(int i = 0; i < bestMatchIdxs.size(); i++){
				//cout << "# " << i << " best match idx " << bestMatchIdxs[i] << " best obj idx " << matchIdxToObjIdx[bestMatchIdxs[i]] << endl;
				bestMatchIdxs[i] = matchIdxToObjIdx[bestMatchIdxs[i]];
			}
		}
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
		TFile* _infile;
		string _oname;
		vector<TH1D*> _hists1D;
		vector<TH2D*> _hists2D;
		vector<TGraph*> graphs;
		vector<TEllipse> _ellipses;
		vector<TMarker> _plot_particles;
		map<int, TEllipse> _jellipses;
		map<int, TMarker> _jcenters;
		map<int,map<int,TEllipse>> _subclellipses;
		map<int,map<int,TMarker>> _subclcenters;

		//do for subclusters
		//map<

		vector<Jet> _phos; //photons for event
		vector<procCat> _procCats;
		vector<node*> _trees;
		vector<Jet> _predJets, _genAK4jets, _recoAK4jets, _genparts, _genAK15jets, _genAK8jets, _recoAK8jets, _recoAK15jets;
		double _minTopPt, _minTopE, _minWPt;
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
		//event selection type
		enum EvtSel{
			//default selection (generic hadronic)
			def = 0,	
			//boosted W selection
			singW = 1,
			//boosted top selection
			boostTop = 2,
			//QCD dijets selection
			QCDdijets = 3
		};
		EvtSel _sel;
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

		double _pt_thresh = 175.;


		//good gen objects available for matching 
		vector<Jet> _genW, _genb, _genTop, _genq, _genglu;
		int _evt2disp, _evt2disp_z;

		map<string, Matrix> _prior_params;

		double dR(double eta1, double phi1, double eta2, double phi2){
			//phi wraparound
			double dphi = (phi1-phi2);
			dphi = acos(cos(dphi));
			return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
		}

		double dot(vector<double> v1, vector<double> v2){
			if(v1.size() != v2.size()){
				cout << "BHCJetSkimmer::dot - vectors given are not the same size" << endl;
				return -1;
			}
			//normalize vecs
			double mag1 = 0;
			for(auto v : v1) mag1 += v*v;
			mag1 = sqrt(mag1);
			for(int i = 0; i < v1.size(); i++)
				v1[i] = v1[i]/mag1;
			
			double mag2 = 0;
			for(auto v : v2) mag2 += v*v;
			mag2 = sqrt(mag2);
			for(int i = 0; i < v2.size(); i++)
				v2[i] = v2[i]/mag2;
			//cout << "mag1 " << mag1 << " mag2 " << mag2 << endl;
			//do dot prod
			double dotprod = 0;
			for(int i = 0; i < v1.size(); i++){
				//cout << "i " << i << " v1 " << v1[i] << " v2 " << v2[i] << endl;
				dotprod += v1[i]*v2[i];
			}
			return dotprod;
		}


		void GetXYZVec(Jet j, vector<double>& v){
			double x = _radius*cos(j.phi());
			double y = _radius*sin(j.phi());
			double theta = 2*atan2(1,exp(j.eta()));
			double z = _radius/tan(theta);
		
			v.clear();
			v.push_back(x);
			v.push_back(y);
			v.push_back(z);
		}

		double CalcSize(Matrix& cov, bool time = false){
			if(cov.GetDims()[0] != 3 || cov.GetDims()[1] != 3){
				cout << "Error: can't calculate size for matrix of size " << cov.GetDims()[0] << " x " << cov.GetDims()[1] << endl;
				return -1;
			}
			vector<double> eigvals;
			vector<Matrix> eigvecs;
			if(time){
				cov.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[2]);
			}
			else{
				Matrix cov2D(2,2);
				Get2DMat(cov,cov2D);	
				cov2D.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[1]);
			}
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
		double Rotundity(Matrix& cov){
			if(cov.GetDims()[0] != 3 || cov.GetDims()[1] != 3){
				cout << "Error: can't calculate rotundity for matrix of size " << cov.GetDims()[0] << " x " << cov.GetDims()[1] << endl;
				return -1;
			}
			Matrix cov2D(2,2);
			Get2DMat(cov,cov2D);	
			vector<Matrix> eigenvecs;
			vector<double> eigenvals;
			cov.eigenCalc(eigenvals, eigenvecs);
			int maxd = cov.GetDims()[0] - 1;
			double rot = 0;
			for(int i = 0; i < (int)eigenvals.size(); i++) rot += eigenvals[i];
			rot = eigenvals[maxd]/rot;
			//if(rot < 0.5 || rot > 1) cout << "rot: " << rot << endl;
			return rot;
		}




		TEllipse PlotEll(const Jet& jet){
			Matrix mu, cov;
			jet.GetClusterParams(mu, cov);
			//ellipse center
			double eta = mu.at(0,0); 
			double phi = mu.at(1,0);

			//get 2D matrix for jet size
			Matrix cov2D(2,2);
			Get2DMat(cov,cov2D);	
		
			double var_eta = cov2D.at(0,0);
			double var_phi = cov2D.at(1,1);
			double covar = cov2D.at(0,1);	
			//use analytic form of eigenvalues for 2x2 matrix - so correct r1, r2 can be passed to TEllipse
			double temp_a  = 0.5*(var_eta + var_phi);
        		double temp_b  = 0.5*sqrt((var_eta-var_phi)*(var_eta-var_phi) + 4*covar*covar);
        		double lambda_eta = temp_a + temp_b;
        		double lambda_phi = temp_a - temp_b;
			double x_r = sqrt(lambda_eta);
			double y_r = sqrt(lambda_phi); 
			
			double theta = 0;
			if(covar > 0) theta = atan((lambda_eta - var_eta)/covar);
			else theta = -atan((lambda_eta - var_eta)/-covar);
			theta = 180 * theta/(4*atan(1)); //put to degrees
			TEllipse el = TEllipse(eta, phi, x_r, y_r, 0, 360, theta);	
			return el;
		}

		void DrawGrid(TCanvas* cv, double true_xmin=-3, double true_xmax=3, double true_ymin=0, double true_ymax=8*atan(1)){
   			double xmin = -3;
			double ymin = 0;
			double xmax = 3;
			double ymax = 8*atan(1);
			double nbins_x = 344;
			double nbins_y = 360;
			double dx = (xmax - xmin)/(nbins_x);
			double dy = (ymax - ymin)/(nbins_y);
	cout << "ymin " << true_ymin << " ymax " << true_ymax << endl;	
			int nline = 0;
			for(int i = 0; i < nbins_x; i++){
				double xcoord = xmin + (double)i*dx;
				if(true_xmax != xmax && xcoord > true_xmax) continue;
				if(true_xmin != xmin && xcoord < true_xmin) continue;
				if(nline == 0) cout << "draw line at x " << xcoord << " from y min " << true_ymin << " to " << true_ymax << " with dx " << dx << " i*dx " << (double)i*dx << " i " << i << endl;
   				auto aline = new TLine(xcoord,true_ymin,xcoord,true_ymax);
				aline->SetLineColor(17);
				//aline->SetLineStyle(2);
				//aline->Draw();
				nline++;

			}
			nline = 0;
			for(int i = 0; i < nbins_y; i++){
				double ycoord = ymin + (double)i*dy;
				if(true_ymax != ymax && ycoord > true_ymax) continue;
				if(true_ymin != ymin && ycoord < true_ymin) continue;
				if(nline == 0) cout << "draw line at y " << ycoord << " from x min " << true_xmin << " to " << true_xmax << " with dy " << dy << " i*dy " <<  (double)i*dy << " i " << i << endl;
				//cout << "draw line at x " << xcoord << " xmax " << xmax << " xmin " << xmin << " dx " << dx << endl;
   				auto aline = new TLine(true_xmin,ycoord,true_xmax,ycoord);
				aline->SetLineColor(17);
				//aline->SetLineStyle(2);
				//aline->Draw();
				nline++;

			}
	
		
		}
		
};
#endif
