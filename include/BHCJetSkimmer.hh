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

			_minRelPt = 0.2;
			_maxRelTimeVar = 1.;
			

			_zoom_window = false;
			
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
			
			_zoom_window = false;
			_minRelPt = 0.2;
			_maxRelTimeVar = 1.;
			
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
			_hists1D.push_back(AK4Jet_rhTimes);
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
			_hists1D.push_back(genAK4JetTop_dR);
			_hists1D.push_back(genAK4JetTop_Eratio);
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
			_hists1D.push_back(BHCJet_nRhs);
			_hists1D.push_back(genAK15Jet_eta);
			_hists1D.push_back(genAK15Jet_phi);
			_hists1D.push_back(genAK15Jet_time);
			_hists1D.push_back(genAK15Jet_pt);
			_hists1D.push_back(genAK15Jet_mass);
			_hists1D.push_back(genAK15Jet_energy);
			_hists1D.push_back(nJet_genAK15Jet);
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
			_hists1D.push_back(BHCJetq_highMass_partonMatchSubclPt);
			_hists1D.push_back(BHCJetq_highMass_partonNoMatchSubclPt);
			_hists1D.push_back(BHCJetq_highMass_partonMatchSubclSize);
			_hists1D.push_back(BHCJetq_highMass_partonNoMatchSubclSize);
			_hists1D.push_back(BHCJetW_lowMass_nSubclustersJet);
			_hists1D.push_back(BHCJetW_Wmass_nSubclustersJet);
			_hists1D.push_back(BHCJetGluon_nSubclusters);
			_hists1D.push_back(BHCJetq_subclParton_dR);
			_hists1D.push_back(BHCJetq_subclParton_Eratio);
			_hists1D.push_back(BHCJetGluon_subclParton_dR);
			_hists1D.push_back(BHCJetGluon_subclParton_Eratio);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclEOvJetE);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclEOvJetE);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclSizeOvJetSize);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize);
			_hists1D.push_back(BHCJet_PUdownweighted_genFrameW_mass);
			_hists1D.push_back(BHCJet_PUremoved_mass);
			_hists1D.push_back(BHCJetW_subclParton_PUcleaned_dR);
			_hists1D.push_back(nBHCJetsW_eq2cleanedSubcls);
			_hists1D.push_back(BHCJet_subclTimeCenter_PUlike);
			_hists1D.push_back(BHCJet_subclTimeCenter_PUcleaned);
			_hists1D.push_back(BHCJet_subclTimeSig_PUlike);
			_hists1D.push_back(BHCJet_subclTimeSig_PUcleaned);
			_hists1D.push_back(BHCJet_PUdownweighted_nSubcleq2_mass);
			_hists1D.push_back(BHCJet_PUdownweighted_nSubcllt2_mass);
			_hists1D.push_back(BHCJet_PUdownweighted_nSubclgt2_mass);
			_hists1D.push_back(BHCJetW_highMass_partonMatchSubclTimeSigOvJetTimeSig);
			_hists1D.push_back(BHCJetW_highMass_partonNoMatchSubclTimeSigOvJetTimeSig);
			_hists1D.push_back(BHCJetq_highMass_partonMatchSubclTimeSigOvJetTimeSig);
			_hists1D.push_back(BHCJetq_highMass_partonNoMatchSubclTimeSigOvJetTimeSig);
			_hists1D.push_back(BHCJetq_highMass_partonMatchSubclEOvJetE);
			_hists1D.push_back(BHCJetq_highMass_partonNoMatchSubclEOvJetE);
			_hists1D.push_back(BHCJetTop_highMass_partonMatchSubclTimeSigOvJetTimeSig);
			_hists1D.push_back(BHCJetTop_highMass_partonNoMatchSubclTimeSigOvJetTimeSig);
			_hists1D.push_back(BHCJetTop_highMass_partonMatchSubclEOvJetE);
			_hists1D.push_back(BHCJetTop_highMass_partonNoMatchSubclEOvJetE);
			_hists1D.push_back(BHCJet_PUdownweighted_mass);
			_hists1D.push_back(BHCJet_massle30_nJets);
			_hists1D.push_back(BHCJet_massle30_PUremoved_SubclInvMass);
			_hists1D.push_back(BHCJet_massle30_PUremoved_nSubclusters);
		

			_hists2D.push_back(jetGenE_diffDeltaPt_recoGen);
			_hists2D.push_back(geoEavg_diffDeltaTime_adjRhs);
			_hists2D.push_back(BHCJet_nJets_jetSize);
			_hists2D.push_back(BHCJet_nSubclustersJet_mass);
			_hists2D.push_back(BHCJet_nSubclustersJet_energy);
			_hists2D.push_back(BHCJet_nSubclustersEvt_nJet);
			_hists2D.push_back(recoAK4Jet_jetEnergy_jetMass);
			_hists2D.push_back(BHCJet_jetEnergy_jetMass);
			_hists2D.push_back(recoAK4Jet_nSubclusters_jetSize);
			_hists2D.push_back(recoAK4Jet_jetEnergy_jetSize);
			_hists2D.push_back(BHCJet_jetEnergy_jetSize);
			_hists2D.push_back(BHCJetPt_genTopPt);
			_hists2D.push_back(BHCJetE_genTopE);
			_hists2D.push_back(BHCJetMass_genTopMass);
			_hists2D.push_back(BHCJetEtaCenter_genTopEtaCenter);
			_hists2D.push_back(BHCJetPhiCenter_genTopPhiCenter);
			_hists2D.push_back(BHCJetPt_genWPt);
			_hists2D.push_back(BHCJetE_genWE);
			_hists2D.push_back(BHCJetMass_genWMass);
			_hists2D.push_back(BHCJetEtaCenter_genWEtaCenter);
			_hists2D.push_back(BHCJetPhiCenter_genWPhiCenter);
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
			_hists2D.push_back(recoAK15JetW_dRGenPartons_jetSize);
			_hists2D.push_back(recoAK4JetW_dRGenPartons_jetSize);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclSize);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclSize);
			_hists2D.push_back(BHCJetq_highMass_partonMatchRelSubclTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetq_highMass_partonNoMatchRelSubclTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetTop_highMass_partonMatchRelSubclTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetTop_highMass_partonNoMatchRelSubclTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclEtaSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclEtaSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclPhiSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclPhiSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE);
			_hists2D.push_back(BHCJetTop_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE);
			_hists2D.push_back(BHCJetTop_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE);
			_hists2D.push_back(BHCJetq_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE);
			_hists2D.push_back(BHCJetq_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE_zoom);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE_zoom);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclPhiVar);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclPhiVar);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclTimeVar);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclTimeVar);
			_hists2D.push_back(BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclPhiVar);
			_hists2D.push_back(BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclPhiVar);

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
							_procCats[p].hists1D[pt][130]->Fill(_trees[i]->model->GetNGhosts());
							vector<double> norms;
							_trees[i]->model->GetNorms(norms);
							Matrix r_post = _trees[i]->model->GetPosterior();
							for(int k = 0; k < _trees[i]->model->GetNClusters(); k++){
								auto params = _trees[i]->model->GetLHPosteriorParameters(k);
								bool g = (bool)params["ghost"].at(0,0);
								if(!g) continue;
								_procCats[p].hists1D[pt][131]->Fill(norms[k]/_gev);
								double eff_rhs = 0;
								for(int n = 0; n < _trees[i]->model->GetData()->GetNPoints(); n++){
									eff_rhs += r_post.at(n,k)/_trees[i]->model->GetData()->at(n).w();
								}
								_procCats[p].hists1D[pt][132]->Fill(eff_rhs);
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

			double wmass, topmass, dr, dr_pair, jetsize, jetsize_full, jet_subleadsize;
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
			if(_genTop.size() > 0)
				GenericMatchJet(_predJets,_genTop, genTopMatchIdxs); //match BHC jets to good gen tops
			if(_genq.size() > 0)
				GenericMatchJet(_predJets,_genq, genqMatchIdxs); //match BHC jets to good gen qs
			if(_genW.size() > 0)
				GenericMatchJet(_predJets,_genW, genWMatchIdxs); //match BHC jets to good gen Ws
			if(_genglu.size() > 0)
				GenericMatchJet(_predJets,_genglu, genGluonMatchIdxs); //match BHC jets to good gen gluons
			if(_sel == singW && _genq.size() > 0){
				//calculate polar opening angle bw gen partons from W
				vector<double> v1, v2;

				GetXYZVec(_genq[0],v1);
				GetXYZVec(_genq[1],v2);

				double dotprod = dot(v1,v2);
				openAng = acos(dotprod);
				//cout << "dotprod " << dotprod << " openAng " << openAng << endl;
			}
	
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
					_procCats[p].hists2D[0][34]->Fill(openAng,njets);

					_procCats[p].hists2D[1][34]->Fill(openAng,njets_lead);
					_procCats[p].hists2D[2][34]->Fill(openAng,njets_notlead);
					
				}
				//loop over pt bins
				njets = _predJets.size();
				int njets_eq2CleanedSubcls = 0;
				int njets_eq2CleanedSubcls_lead = 0;
				int njets_eq2CleanedSubcls_notlead = 0;
				for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
					vector<Jet> lowMassJets;
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
						jet_subleadsize = CalcSubleadSize(jet_cov);
						jetsize_full = CalcSize(jet_cov, true);
						//define jetsize bins
						//pt == 1 -> [0,0.2) (AK4 level)
						//if(pt == 1 && jetsize >= 0.2) continue;
						////pt == 2 -> [0.2,inf) (AK15 level)
						//if(pt == 2 && jetsize < 0.2) continue;


						if(p != 0 && pt == 0){ 
							cout << "pred jet #" << j << " phi " << _predJets[j].phi() << " eta " << _predJets[j].eta() << " energy " << _predJets[j].E() <<  " mass " << _predJets[j].mass() << " nConstituents " << _predJets[j].GetNConstituents() << " nRhs " << _predJets[j].GetNRecHits() << " pt " << _predJets[j].pt() << " jetsize " << jetsize << " eta var " << jet_cov.at(0,0) << " phi var " << jet_cov.at(1,1) << endl;
						}
			
						
						nsubs = _predJets[j].GetNConstituents();
						_procCats[p].hists1D[pt][1]->Fill(nsubs);
						if(openAng != -1) _procCats[p].hists2D[pt][35]->Fill(openAng,nsubs);
						//cout << "pred jet j " << j << " dr " << dr << " n constituents " << _predJets[j].GetNConstituents() << endl;
						Matrix cov = _predJets[j].GetCovariance();
					
						_procCats[p].hists1D[pt][7]->Fill(_predJets[j].e());
						_procCats[p].hists1D[pt][8]->Fill(_predJets[j].pt());
						_procCats[p].hists1D[pt][9]->Fill(_predJets[j].mass());
						_procCats[p].hists1D[pt][6]->Fill(jetsize);
							
						_procCats[p].hists1D[pt][22]->Fill(sqrt(cov.at(0,0)));	
						_procCats[p].hists1D[pt][23]->Fill(sqrt(cov.at(1,1)));	
						_procCats[p].hists1D[pt][24]->Fill(sqrt(cov.at(2,2)));	
						_procCats[p].hists1D[pt][25]->Fill(cov.at(0,1));	
						_procCats[p].hists1D[pt][26]->Fill(cov.at(0,2));	
						_procCats[p].hists1D[pt][27]->Fill(cov.at(2,1));	

						
						_procCats[p].hists2D[pt][2]->Fill((double)_predJets.size(), jetsize);
						_procCats[p].hists2D[pt][3]->Fill(_predJets[j].GetNConstituents(), _predJets[j].mass());
						_procCats[p].hists2D[pt][4]->Fill(_predJets[j].GetNConstituents(), _predJets[j].e());
						_procCats[p].hists2D[pt][10]->Fill(_predJets[j].e(), jetsize);
						_procCats[p].hists2D[pt][7]->Fill(_predJets[j].e(), _predJets[j].mass());
						_procCats[p].hists2D[pt][43]->Fill(_predJets[j].mass(), jetsize);	

						_procCats[p].hists1D[pt][60]->Fill(jet_mu.at(0,0));
						_procCats[p].hists1D[pt][61]->Fill(_predJets[j].phi()); //so phi is [0,2pi] but equivalent to jet_mu.at(1,0)
						_procCats[p].hists1D[pt][62]->Fill(jet_mu.at(2,0));
					
						
						_procCats[p].hists1D[pt][74]->Fill((int)rhs.size());
						for(int r = 0; r < rhs.size(); r++){
							_procCats[p].hists1D[pt][73]->Fill(rhs[r].E());
							
						}

					//cout << "pred jet #" << j << " phi1 " << _predJets[j].phi() << " phi " << jet_mu.at(1,0) << " phi std " << _predJets[j].phi_std() << endl;
						//get subcluster information
						vector<Jet> consts = _predJets[j].GetConstituents();
						vector<Jet> consts_puCleaned;
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
							if(openAng != -1) _procCats[p].hists2D[pt][36]->Fill(openAng,subcl.m());
							_procCats[p].hists1D[pt][3]->Fill(subcl.eta());
							_procCats[p].hists1D[pt][4]->Fill(subcl.phi());
							_procCats[p].hists1D[pt][5]->Fill(subcl.time());
							
							_procCats[p].hists1D[pt][124]->Fill(subcl.m());
							_procCats[p].hists2D[pt][21]->Fill(subcl.E(), subcl.m());
							//effective # of rhs is sum of weights bc weights are each the responsibility of rechit n to this subcl
							vector<double> ws;
							subcl.GetWeights(ws);
							double effnRhs = 0;
							for(int n = 0; n < ws.size(); n++){
								effnRhs += ws[n]; 
							}
							_procCats[p].hists1D[pt][125]->Fill(effnRhs);
							_procCats[p].hists2D[pt][22]->Fill(subcl.E(), effnRhs);
							_procCats[p].hists2D[pt][23]->Fill(subcl.m(), effnRhs);
							if(p == 0 && pt == 0) cout << " subcl #" << c << " has " << effnRhs << " # of effective rechits of " << ws.size() << " total rhs and mass " << subcl.m() << " and energy " << subcl.E() << " px " << subcl.px() << " py " << subcl.py() << " pz " << subcl.pz() << endl;
							
							//do dr calcs bw all subclusters
							for(int cc = c+1; cc < consts.size(); cc++){
								double dr = dR(consts[c].eta(), consts[c].phi(), consts[cc].eta(), consts[cc].phi());
								_procCats[p].hists1D[pt][71]->Fill(dr);
							}
							//dr to jet center
							double jet_dr = dR(consts[c].eta(), consts[c].phi(), jet_mu.at(0,0), jet_mu.at(1,0));
							_procCats[p].hists2D[pt][27]->Fill(consts[c].E(), jet_dr);
							_procCats[p].hists2D[pt][28]->Fill(effnRhs, jet_dr);
						
							//PU cleaning subclusters
							double relE = consts[c].E() / _predJets[j].E();
							double relEtaVar = subcl_cov.at(0,0) / cov.at(0,0);
							double relPhiVar = subcl_cov.at(1,1) / cov.at(1,1);
							double relTimeVar = subcl_cov.at(2,2) / cov.at(2,2);
							double relGeoAvg = pow( relEtaVar * relPhiVar * relTimeVar, 1./3.);
							//hard subcl
							if(relGeoAvg - relE < 0){
								_procCats[p].hists1D[pt][196]->Fill(consts[c].t());
								_procCats[p].hists1D[pt][198]->Fill(sqrt(subcl_cov.at(2,2)));
								consts_puCleaned.push_back(consts[c]);
							}
							else{ //PU-like subcl
								_procCats[p].hists1D[pt][195]->Fill(consts[c].t());
								_procCats[p].hists1D[pt][197]->Fill(sqrt(subcl_cov.at(2,2)));
					
							}


						}
						Jet cleanedJet_downweight = _predJets[j].CleanOutPU(false);
						Jet cleanedJet_remove = _predJets[j].CleanOutPU(true);

						_procCats[p].hists1D[pt][212]->Fill(cleanedJet_downweight.mass());
//cout << "original jet mass " << _predJets[j].mass() << " PU removed jet mass " << cleanedJet_remove.mass() << " PU dwnwt mass " << cleanedJet_downweight.mass() << endl;
//cout << "original # subclusters " << _predJets[j].GetNConstituents() << " PU removed # subcls " << cleanedJet_remove.GetNConstituents() << endl;
						//look at low mass jets (resolved)
						if(cleanedJet_remove.mass() < 30){
							if(pt == 1 && cleanedJet_remove.pt() < _pt_thresh) continue;
							if(pt == 2 && cleanedJet_remove.pt() >= _pt_thresh) continue;
							lowMassJets.push_back(cleanedJet_remove);
						}
//cout << "lowMassJets.size() " << lowMassJets.size() << endl;
						//only plot for jets within a "gen frame" of gen particle
						if(_genparts.size() > 0){ 
							double maxDr = 1.;
							double bestdR = 999;
							double dR_dwnwt, dR_remove;
							for(auto part : _genparts){
								dR_dwnwt = dR(part.eta(), part.phi(), cleanedJet_downweight.eta(), cleanedJet_downweight.phi());
								if(dR_dwnwt < bestdR)
									bestdR = dR_dwnwt;
							}
							if(bestdR < maxDr && cleanedJet_downweight.mass() > 0){
								_procCats[p].hists1D[pt][191]->Fill(cleanedJet_downweight.mass());
								if(cleanedJet_remove.GetNConstituents() == 2){
									_procCats[p].hists1D[pt][199]->Fill(cleanedJet_downweight.mass());
								}
								else if(cleanedJet_remove.GetNConstituents() < 2){
									_procCats[p].hists1D[pt][200]->Fill(cleanedJet_downweight.mass());
								}
								else{
									_procCats[p].hists1D[pt][201]->Fill(cleanedJet_downweight.mass());
								}
							}
							bestdR = 999;
							for(auto part : _genparts){
								dR_remove = dR(part.eta(), part.phi(), cleanedJet_downweight.eta(), cleanedJet_downweight.phi());
								if(dR_remove < bestdR)
									bestdR = dR_remove;
							}
							if(bestdR < maxDr && cleanedJet_remove.mass() > 0)	
								_procCats[p].hists1D[pt][192]->Fill(cleanedJet_remove.mass());

						}
						if(p != 0 && pt == 0){ 
							cout << "cleaned (dwnwt) pred jet #" << j << " phi " << cleanedJet_downweight.phi() << " eta " << cleanedJet_downweight.eta() << " energy " << cleanedJet_downweight.E() <<  " mass " << cleanedJet_downweight.mass() << " nConstituents " << cleanedJet_downweight.GetNConstituents() << " nRhs " << cleanedJet_downweight.GetNRecHits() << " pt " << cleanedJet_downweight.pt() << " jetsize " << jetsize << " eta var " << jet_cov.at(0,0) << " phi var " << jet_cov.at(1,1) << endl;
							cout << "cleaned (remove) pred jet #" << j << " phi " << cleanedJet_remove.phi() << " eta " << cleanedJet_remove.eta() << " energy " << cleanedJet_remove.E() <<  " mass " << cleanedJet_remove.mass() << " nConstituents " << cleanedJet_remove.GetNConstituents() << " nRhs " << cleanedJet_remove.GetNRecHits() << " pt " << cleanedJet_remove.pt() << " jetsize " << jetsize << " eta var " << jet_cov.at(0,0) << " phi var " << jet_cov.at(1,1) << endl;
						}

						if(p == 0 && pt == 0) cout << "pred jet #" << j << " matched to reco AK4 jet #" << recoMatchIdxs[j] << " gen AK4 jet #" << genAK4MatchIdxs[j] << " and gen AK15 jet #" << genAK15MatchIdxs[j] << " and gen top #" << genTopMatchIdxs[j] << " and gen W #" << genWMatchIdxs[j] << endl;
						//exclusively dr-matched to reco jet
						if(recoMatchIdxs[j] != -1){
							if(_recoAK4jets[recoMatchIdxs[j]].GetNConstituents() == 0) cout << "reco jet #" << recoMatchIdxs[j] << " matched to bhc jet #" << j << " has " << _recoAK4jets[recoMatchIdxs[j]].GetNConstituents() << " # subclusters" << endl;
						}
						//do reco AK8 matching hists
						if(recoAK8MatchIdxs[j] != -1){
							int matchidx = recoAK8MatchIdxs[j];
							_procCats[p].hists2D[pt][47]->Fill(_recoAK8jets[matchidx].m(), _predJets[j].m());
						}
	
						//do gen AK4 matching hists
						//do gen top matching hists
						if(genTopMatchIdxs[j] != -1){
							int gentopidx = genTopMatchIdxs[j];
							double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genTop[gentopidx].eta(), _genTop[gentopidx].phi());
							double eratio = _predJets[j].E()/_genTop[gentopidx].E();
							if(p == 0 && pt == 0) cout << "BHC jet #" << j << " with energy " << _predJets[j].E() << " matched to top " << gentopidx << " with eta " << _genTop[gentopidx].eta() << " phi " << _genTop[gentopidx].phi() << " and energy " << _genTop[gentopidx].E() << " matched with dr " << dr << " and Eratio " << eratio << endl;
							_procCats[p].hists2D[pt][11]->Fill(_predJets[j].pt(),_genTop[gentopidx].pt());
							_procCats[p].hists2D[pt][12]->Fill(_predJets[j].E(),_genTop[gentopidx].E());
							_procCats[p].hists2D[pt][13]->Fill(_predJets[j].m(),_genTop[gentopidx].m());
							_procCats[p].hists2D[pt][14]->Fill(_predJets[j].eta(),_genTop[gentopidx].eta());
							_procCats[p].hists2D[pt][15]->Fill(_predJets[j].phi(),_genTop[gentopidx].phi());
							
							_procCats[p].hists1D[pt][58]->Fill(dr);
							_procCats[p].hists1D[pt][59]->Fill(eratio);
							_procCats[p].hists1D[pt][135]->Fill(_predJets[j].GetNConstituents());
							//sort by energy
							consts = _predJets[j].GetConstituents();
							sort(consts.begin(), consts.end(), Esort_jet);
							//get invariant mass of leading two subclusters
							vector<int> subclGenMatchIdx(consts.size(), -1);	
							if(consts.size() > 1){
								double invmass = consts[0].invMass(consts[1]);
								_procCats[p].hists1D[pt][137]->Fill(invmass);
							
								//do gen matching of 2 leading subclusters to partons from W decays
								Jet leadcl = consts[0];
								Jet subleadcl = consts[1];
								vector<Jet> subcls = {leadcl, subleadcl};
								//get gen partons from top decay 
								vector<int> genLeadMatchIdxs(3,-1);
								vector<Jet> Tpartons;
								//partons from top are saved in src/ main event loop with boostTop selection
								for(auto g : _genq) Tpartons.push_back(g);
								for(auto g : _genb) Tpartons.push_back(g);
								GenericMatchJet(consts,Tpartons,subclGenMatchIdx); //match subclusters to W partons
								for(int c = 0; c < consts.size(); c++){ 
									Matrix subcl_cov = consts[c].GetCovariance();
									double subclsize = CalcSize(subcl_cov);
									double subcl_subleadsize = CalcSubleadSize(subcl_cov);
									double subclsize_full = CalcSize(subcl_cov, true);
									int genmatchidx = subclGenMatchIdx[c];
									//do 3D cuts in relE-relSize-relTimeSig space
									//relTimeSig - relE <= 0 is nonPU else is PU
									//relTimeSig + relSize <= 1 is nonPU else is PU
									double relTimeSig = sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2));
									double relEtaVar = subcl_cov.at(0,0) / jet_cov.at(0,0);
									double relPhiVar = subcl_cov.at(1,1) / jet_cov.at(1,1);
									double relTimeVar = subcl_cov.at(2,2) / jet_cov.at(2,2);
									double relSize = subclsize / jetsize;
									double relSubleadSize = subcl_subleadsize / jet_subleadsize;
									double relE = consts[c].E() / _predJets[j].E();
									if(_predJets[j].mass() > 200){
										//not-matched
										if(subclGenMatchIdx[c] == -1){
											_procCats[p].hists1D[pt][209]->Fill( sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)));
											_procCats[p].hists1D[pt][211]->Fill(consts[c].E() / _predJets[j].E() );
											_procCats[p].hists2D[pt][60]->Fill( relTimeSig, relE);
											_procCats[p].hists2D[pt][72]->Fill( pow(relTimeVar * relEtaVar * relPhiVar, 1./3.), relE);
										}
										//matched
										else{
											_procCats[p].hists1D[pt][208]->Fill( sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)));
											_procCats[p].hists1D[pt][210]->Fill(consts[c].E() / _predJets[j].E() );
											_procCats[p].hists2D[pt][59]->Fill( relTimeSig, relE);
											_procCats[p].hists2D[pt][71]->Fill( pow(relTimeVar * relEtaVar * relPhiVar, 1./3.), relE);

										}


									}


								}

							}
							for(int c = 0; c < consts.size(); c++){ 
								_procCats[p].hists1D[pt][136]->Fill(consts[c].m());
								_procCats[p].hists2D[pt][33]->Fill((int)consts.size(),consts[c].m());
							

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
							_procCats[p].hists2D[pt][16]->Fill(_predJets[j].pt(),_genW[genWMatchIdxs[j]].pt());
							_procCats[p].hists2D[pt][17]->Fill(_predJets[j].E(),_genW[genWMatchIdxs[j]].E());
							_procCats[p].hists2D[pt][18]->Fill(_predJets[j].m(),_genW[genWMatchIdxs[j]].m());
							_procCats[p].hists2D[pt][19]->Fill(_predJets[j].eta(),_genW[genWMatchIdxs[j]].eta());
							_procCats[p].hists2D[pt][20]->Fill(_predJets[j].phi(),_genW[genWMatchIdxs[j]].phi());
							
							_procCats[p].hists1D[pt][120]->Fill(dr);
							_procCats[p].hists1D[pt][121]->Fill(eratio);
							consts = _predJets[j].GetConstituents();
							if(consts.size() > 1) _procCats[p].hists2D[pt][42]->Fill(_predJets[j].E()/_genW[genWidx].E(),consts.size());
							
							//finding partons/subclusters
							_procCats[p].hists1D[pt][122]->Fill(_predJets[j].GetNConstituents());
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
							_procCats[p].hists2D[pt][30]->Fill(gendR, _predJets[j].GetNConstituents());
							_procCats[p].hists2D[pt][31]->Fill(gendR, jetsize);

							double avgPartE = (Wpartons[0].E() + Wpartons[1].E())/2.;
cout << "avgPart E " << avgPartE << endl;
					
							if(consts.size() == 1){
								_procCats[p].hists2D[pt][37]->Fill(gendR, jetsize);
								_procCats[p].hists2D[pt][39]->Fill(gendR, avgPartE);
							}
							if(consts.size() > 1){
								_procCats[p].hists2D[pt][38]->Fill(gendR, jetsize);
								_procCats[p].hists2D[pt][40]->Fill(gendR, avgPartE);

							}
							if(_predJets[j].m() > 100) _procCats[p].hists1D[pt][161]->Fill((int)consts.size());
							if(70 < _predJets[j].m() && _predJets[j].m() < 90) _procCats[p].hists1D[pt][181]->Fill((int)consts.size());
							if(_predJets[j].m() < 50) _procCats[p].hists1D[pt][180]->Fill((int)consts.size());
							//sort by energy
							sort(consts.begin(), consts.end(), Esort_jet);
							//get invariant mass of leading two subclusters
							if(consts.size() > 1){
								double invmass = consts[0].invMass(consts[1]);
								_procCats[p].hists1D[pt][129]->Fill(invmass);
								vector<int> subclGenMatchIdx(consts.size(), -1);	
								GenericMatchJet(consts,Wpartons,subclGenMatchIdx); //match subclusters to W partons
								for(int c = 0; c < consts.size(); c++){
									_procCats[p].hists2D[pt][41]->Fill(consts[c].E(),c);
									
									Matrix subcl_mu, subcl_cov;
									consts[c].GetClusterParams(subcl_mu, subcl_cov);
									double subclsize = CalcSize(subcl_cov);
									double subcl_subleadsize = CalcSubleadSize(subcl_cov);
									double subclsize_full = CalcSize(subcl_cov, true);
									int genmatchidx = subclGenMatchIdx[c];
									//do 3D cuts in relE-relSize-relTimeSig space
									//relTimeSig - relE <= 0 is nonPU else is PU
									//relTimeSig + relSize <= 1 is nonPU else is PU
									double relTimeSig = sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2));
									double relEtaVar = subcl_cov.at(0,0) / jet_cov.at(0,0);
									double relPhiVar = subcl_cov.at(1,1) / jet_cov.at(1,1);
									double relTimeVar = subcl_cov.at(2,2) / jet_cov.at(2,2);
									double relSize = subclsize / jetsize;
									double relSubleadSize = subcl_subleadsize / jet_subleadsize;
									double relE = consts[c].E() / _predJets[j].E();
									if(genmatchidx == -1){
										//not matched
										if(_predJets[j].m() > 100){
											_procCats[p].hists1D[pt][152]->Fill(consts[c].pt());
											_procCats[p].hists1D[pt][188]->Fill(consts[c].E() / _predJets[j].E());
											_procCats[p].hists1D[pt][154]->Fill(subclsize);
											_procCats[p].hists1D[pt][190]->Fill(subclsize / jetsize);
											_procCats[p].hists1D[pt][203]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)));
											
											_procCats[p].hists2D[pt][52]->Fill(subclsize / jetsize, consts[c].pt() / _predJets[j].pt());
											//relative time variance vs relative pt
											_procCats[p].hists2D[pt][54]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)), consts[c].E() / _predJets[j].E());
											//relative time variance vs relative size
											_procCats[p].hists2D[pt][56]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)), subclsize / jetsize);
											//subcl mass / pt vs jet mass / pt
											//full (space + time) subcl size / full jet size vs subcl pt / jet pt
											_procCats[p].hists2D[pt][62]->Fill(sqrt(subcl_cov.at(0,0)) / sqrt(jet_cov.at(0,0)), relE);
											_procCats[p].hists2D[pt][64]->Fill(sqrt(subcl_cov.at(1,1)) / sqrt(jet_cov.at(1,1)), relE);
											_procCats[p].hists2D[pt][66]->Fill(subcl_cov.at(2,2) / jet_cov.at(2,2), relE);
											_procCats[p].hists2D[pt][68]->Fill(sqrt(relSize * relTimeSig), relE);
											_procCats[p].hists2D[pt][70]->Fill(pow(relEtaVar * relPhiVar * relTimeVar, 1./3.), relE);
																		
										}
				
									}
									else{
										//matched
										if(_predJets[j].m() > 100){
											_procCats[p].hists1D[pt][151]->Fill(consts[c].pt());
											_procCats[p].hists1D[pt][187]->Fill(consts[c].E()/_predJets[j].E());
											_procCats[p].hists1D[pt][153]->Fill(subclsize);
											_procCats[p].hists1D[pt][189]->Fill(subclsize / jetsize);
											_procCats[p].hists1D[pt][202]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)));

											_procCats[p].hists2D[pt][51]->Fill(subclsize / jetsize, consts[c].pt() / _predJets[j].pt());
											//relative time variance vs relative pt
											_procCats[p].hists2D[pt][53]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)), consts[c].E() / _predJets[j].E());
											//relative time variance vs relative size
											_procCats[p].hists2D[pt][55]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)), subclsize / jetsize);
											//subcl mass / pt vs jet mass / pt
											//full (space + time) subcl size / full jet size vs subcl pt / jet pt
											_procCats[p].hists2D[pt][61]->Fill(sqrt(subcl_cov.at(0,0)) / sqrt(jet_cov.at(0,0)), relE);
											_procCats[p].hists2D[pt][63]->Fill(sqrt(subcl_cov.at(1,1)) / sqrt(jet_cov.at(1,1)), relE);
											_procCats[p].hists2D[pt][65]->Fill(subcl_cov.at(2,2) / jet_cov.at(2,2), relE);
											_procCats[p].hists2D[pt][67]->Fill(sqrt(relSize * relTimeSig), relE);
											_procCats[p].hists2D[pt][69]->Fill(pow(relEtaVar * relPhiVar * relTimeVar, 1./3.), relE);
										}
									}

								}	
								vector<int> subclGenMatchIdx_puCleaned(consts_puCleaned.size(), -1);	
								GenericMatchJet(consts_puCleaned,Wpartons,subclGenMatchIdx_puCleaned); //match subclusters to W partons
								for(int c = 0; c < consts_puCleaned.size(); c++){
									int genmatchidx = subclGenMatchIdx_puCleaned[c];
									if(genmatchidx == -1) continue;
									double gen_clDr = dR(consts_puCleaned[c].eta(), consts_puCleaned[c].phi(), Wpartons[genmatchidx].eta(), Wpartons[genmatchidx].phi());
									_procCats[p].hists1D[pt][193]->Fill(gen_clDr);


								}
								if(consts_puCleaned.size() == 2){
									if(pt == 0) njets_eq2CleanedSubcls++;
									if(pt == 1 && cleanedJet_remove.pt() > _pt_thresh) njets_eq2CleanedSubcls_lead++;
									if(pt == 2 && cleanedJet_remove.pt() <= _pt_thresh) njets_eq2CleanedSubcls_notlead++;
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
									_procCats[p].hists1D[pt][133]->Fill(gen_clDr);
									_procCats[p].hists1D[pt][134]->Fill(genEr);
									_procCats[p].hists2D[pt][32]->Fill(genEr,_predJets[j].GetNConstituents());

									_procCats[p].hists1D[pt][142]->Fill(subcls[c].eta());
									_procCats[p].hists1D[pt][143]->Fill(subcls[c].phi());
									_procCats[p].hists1D[pt][144]->Fill(subcls[c].time());
									Matrix subcl_cov = subcls[c].GetCovariance();
									_procCats[p].hists1D[pt][145]->Fill(sqrt(subcl_cov.at(0,0)));
									_procCats[p].hists1D[pt][146]->Fill(sqrt(subcl_cov.at(1,1)));
									_procCats[p].hists1D[pt][147]->Fill(sqrt(subcl_cov.at(2,2)));
									_procCats[p].hists1D[pt][148]->Fill(subcl_cov.at(0,1));
									_procCats[p].hists1D[pt][149]->Fill(subcl_cov.at(0,2));
									_procCats[p].hists1D[pt][150]->Fill(subcl_cov.at(2,1));

								}
							}
							for(int c = 0; c < (int)consts.size(); c++){
								_procCats[p].hists1D[pt][123]->Fill(consts[c].E());
								_procCats[p].hists1D[pt][128]->Fill(consts[c].m());
							//cout << "consts #" << c << " mass " << consts[c].m() << endl;
								_procCats[p].hists2D[pt][29]->Fill(_predJets[j].GetNConstituents(),consts[c].m());
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
							_procCats[p].hists1D[pt][138]->Fill(dr);
							_procCats[p].hists1D[pt][139]->Fill(eratio);
							//finding partons/subclusters
							_procCats[p].hists1D[pt][140]->Fill(_predJets[j].GetNConstituents());

							consts = _predJets[j].GetConstituents();
							//sort by energy
							sort(consts.begin(), consts.end(), Esort_jet);
							for(int c = 0; c < (int)consts.size(); c++){
								_procCats[p].hists1D[pt][141]->Fill(consts[c].m());
							}
							vector<int> subclGenMatchIdx(consts.size(), -1);	
							GenericMatchJet(consts,_genq,subclGenMatchIdx); //match subclusters to hard light quarks
							for(int c = 0; c < consts.size(); c++){
								int genmatchidx = subclGenMatchIdx[c];
								//1+ subcl jets only
								Matrix subcl_cov = consts[c].GetCovariance();
								double subclsize = CalcSize(subcl_cov);
								double subcl_subleadsize = CalcSubleadSize(subcl_cov);
								//do 3D cuts in relE-relSize-relTimeSig space
								//relTimeSig - relE <= 0 is nonPU else is PU
								//relTimeSig + relSize <= 1 is nonPU else is PU
								double relTimeSig = sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2));
								double relEtaVar = subcl_cov.at(0,0) / jet_cov.at(0,0);
								double relPhiVar = subcl_cov.at(1,1) / jet_cov.at(1,1);
								double relTimeVar = subcl_cov.at(2,2) / jet_cov.at(2,2);
								double relSize = subclsize / jetsize;
								double relSubleadSize = subcl_subleadsize / jet_subleadsize;
								double relE = consts[c].E() / _predJets[j].E();
								if((int)consts.size() > 1 && _predJets[j].m() > 50){
									//not matched
									if(genmatchidx == -1){
										_procCats[p].hists1D[pt][177]->Fill(consts[c].pt());
										_procCats[p].hists1D[pt][179]->Fill(subclsize);
										_procCats[p].hists1D[pt][205]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)) );
										_procCats[p].hists1D[pt][207]->Fill(consts[c].E() / _predJets[j].E());
										_procCats[p].hists2D[pt][58]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)) ,consts[c].E() / _predJets[j].E());
										_procCats[p].hists2D[pt][74]->Fill( pow(relTimeVar * relEtaVar * relPhiVar, 1./3.), relE);
										_procCats[p].hists2D[pt][76]->Fill( pow(relTimeVar * relEtaVar * relPhiVar, 1./3.), relE);
										_procCats[p].hists2D[pt][78]->Fill( relEtaVar, relPhiVar);
										_procCats[p].hists2D[pt][80]->Fill( relEtaVar, relTimeVar);
										_procCats[p].hists2D[pt][82]->Fill( relTimeVar, relPhiVar);
									}
									//matched
									else{
										_procCats[p].hists1D[pt][176]->Fill(consts[c].pt());
										_procCats[p].hists1D[pt][178]->Fill(subclsize);
										_procCats[p].hists1D[pt][204]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)) );
										_procCats[p].hists1D[pt][206]->Fill(consts[c].E() / _predJets[j].E());
										_procCats[p].hists2D[pt][57]->Fill(sqrt(subcl_cov.at(2,2)) / sqrt(jet_cov.at(2,2)) ,consts[c].E() / _predJets[j].E());
										_procCats[p].hists2D[pt][73]->Fill( pow(relTimeVar * relEtaVar * relPhiVar, 1./3.), relE);
										_procCats[p].hists2D[pt][75]->Fill( pow(relTimeVar * relEtaVar * relPhiVar, 1./3.), relE);
										_procCats[p].hists2D[pt][77]->Fill( relEtaVar, relPhiVar);
										_procCats[p].hists2D[pt][79]->Fill( relEtaVar, relTimeVar);
										_procCats[p].hists2D[pt][81]->Fill( relTimeVar, relPhiVar);
									}
								}
							}	
							//do gen matching of leading subcluster to hard parton
							Jet leadcl = consts[0];
							vector<Jet> subcls = {leadcl};
							vector<int> genLeadMatchIdxs(1,-1);
							GenericMatchJet(subcls,_genq,genLeadMatchIdxs); //match subclusters to quarks
							for(int c = 0; c < subcls.size(); c++){
								int genmatchidx = genLeadMatchIdxs[c];
								if(genmatchidx == -1) continue;
								double gen_clDr = dR(subcls[c].eta(), subcls[c].phi(), _genq[genmatchidx].eta(), _genq[genmatchidx].phi());
								double genEr = subcls[c].E() / _genq[genmatchidx].E();
								_procCats[p].hists1D[pt][183]->Fill(gen_clDr);
								_procCats[p].hists1D[pt][184]->Fill(genEr);

							}

						}			
						//do gen gluon matching hists
						if(genGluonMatchIdxs[j] != -1){
							int genGluonidx = genGluonMatchIdxs[j];
							double dr = dR(_predJets[j].eta(), _predJets[j].phi(), _genglu[genGluonidx].eta(), _genglu[genGluonidx].phi());
							double eratio = _predJets[j].E()/_genglu[genGluonidx].E();
							if(p == 0 && pt == 0) cout << "BHC jet #" << j << " with eta " << _predJets[j].eta() << " phi " << _predJets[j].phi() << " and energy " << _predJets[j].E() << " matched to gluon " << genGluonidx << " with eta " << _genglu[genGluonidx].eta() << " phi " << _genglu[genGluonidx].phi() << " and energy " << _genglu[genGluonidx].E() << " matched with dr " << dr << " and Eratio " << eratio << endl;
							_procCats[p].hists1D[pt][159]->Fill(dr);
							_procCats[p].hists1D[pt][160]->Fill(eratio);
							_procCats[p].hists1D[pt][182]->Fill(_predJets[j].GetNConstituents());
							//do gen matching of leading subcluster to hard parton
							Jet leadcl = consts[0];
							vector<Jet> subcls = {leadcl};
							vector<int> genLeadMatchIdxs(1,-1);
							GenericMatchJet(subcls,_genglu,genLeadMatchIdxs); //match subclusters to gluons
							for(int c = 0; c < subcls.size(); c++){
								int genmatchidx = genLeadMatchIdxs[c];
								if(genmatchidx == -1) continue;
								double gen_clDr = dR(subcls[c].eta(), subcls[c].phi(), _genglu[genmatchidx].eta(), _genglu[genmatchidx].phi());
								double genEr = subcls[c].E() / _genglu[genmatchidx].E();
								_procCats[p].hists1D[pt][185]->Fill(gen_clDr);
								_procCats[p].hists1D[pt][186]->Fill(genEr);

							}
						}

					}
					_procCats[p].hists1D[pt][213]->Fill((double)lowMassJets.size());	
					if(p == 0) 
						_procCats[p].hists1D[pt][194]->Fill((double)njets/(double)njets_eq2CleanedSubcls);
					if(p == 1) 
						_procCats[p].hists1D[pt][194]->Fill((double)njets_lead/(double)njets_eq2CleanedSubcls_lead);
					if(p == 2) 
						_procCats[p].hists1D[pt][194]->Fill((double)njets_notlead/(double)njets_eq2CleanedSubcls_notlead);
					vector<int> genLowMassIdxs(lowMassJets.size(),-1);
					GenericMatchJet(lowMassJets,_genq,genLowMassIdxs); //match low mass jets to quarks
					//make sure there are at least two jets
					if(lowMassJets.size() > 1){
						//indices for inv mass
						double minDiff = 999;
						double Wmass = 80.4;
						int nsubcls1, nsubcls2;
						double bestInvMass = -1;
						for(int j = 0; j < lowMassJets.size(); j++){
							int genmatchidx = genLowMassIdxs[j];
							if(genmatchidx == -1) continue;
							//cout << "lowMassJet #" << j << " matched to q # " << genmatchidx << " with dR " << dR(lowMassJets[j].eta(), lowMassJets[j].phi(), _genq[genmatchidx].eta(), _genq[genmatchidx].phi()) << endl;
							for(int jj = j+1; jj < lowMassJets.size(); jj++){
								int ggenmatchidx = genLowMassIdxs[jj];
								if(ggenmatchidx == -1) continue;
								double invMass = lowMassJets[j].invMass(lowMassJets[jj]);
								if(fabs(invMass - Wmass) < minDiff){
									minDiff = fabs(invMass - Wmass);
									bestInvMass = invMass;
									nsubcls1 = lowMassJets[j].GetNConstituents();
									nsubcls2 = lowMassJets[jj].GetNConstituents();
								}	
								_procCats[p].hists1D[pt][214]->Fill(bestInvMass);
								_procCats[p].hists1D[pt][215]->Fill(nsubcls1);
								_procCats[p].hists1D[pt][215]->Fill(nsubcls2);
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
				_procCats[p].hists1D[0][42]->Fill((double)nGenParts);
				for(int g = 0; g < _genparts.size(); g++){
					_procCats[p].hists1D[0][43]->Fill(_genparts[g].eta());
					_procCats[p].hists1D[0][44]->Fill(_genparts[g].phi());
					_procCats[p].hists1D[0][45]->Fill(_genparts[g].time());
					_procCats[p].hists1D[0][46]->Fill(_genparts[g].pt());
					_procCats[p].hists1D[0][47]->Fill(_genparts[g].mass());
					_procCats[p].hists1D[0][48]->Fill(_genparts[g].E());
					id = _base->genpart_id->at(_genparts[g].GetUserIdx());
					if(find(qids.begin(), qids.end(), id) != qids.end()) nGenParts++;
					
				}
			}
		}
	
		void FillGenJetHists(){
			for(int p = 0; p < _procCats.size(); p++){
				_procCats[p].hists1D[0][55]->Fill((double)_genAK4jets.size());
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _genAK4jets.size(); j++){
					if(_genAK4jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_genAK4jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][55]->Fill(njets_lead);
				_procCats[p].hists1D[2][55]->Fill(njets_notlead);
				//gen match jets to particles
				vector<int> genTopMatchIdxs;
				GenericMatchJet(_genAK4jets,_genparts,genTopMatchIdxs, 6); //match gen AK4 jets to gen tops
				vector<int> genWMatchIdxs;
				GenericMatchJet(_genAK4jets,_genparts,genWMatchIdxs, 24); //match gen AK4 jets to gen Ws
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genAK4jets.size(); j++){
					if(p == 0) cout << "gen AK4 jet #" << j << " phi " << _genAK4jets[j].phi() << " eta " << _genAK4jets[j].eta() << " energy " << _genAK4jets[j].E() <<  " mass " << _genAK4jets[j].mass() << " pt " << _genAK4jets[j].pt() << endl;
					_procCats[p].hists1D[0][49]->Fill(_genAK4jets[j].eta());
					_procCats[p].hists1D[0][50]->Fill(_genAK4jets[j].phi());
					_procCats[p].hists1D[0][51]->Fill(_genAK4jets[j].time());
					_procCats[p].hists1D[0][52]->Fill(_genAK4jets[j].pt());
					_procCats[p].hists1D[0][53]->Fill(_genAK4jets[j].mass());
					_procCats[p].hists1D[0][54]->Fill(_genAK4jets[j].E());
					//dr bw gen jet and best exclusive gen top match
					int genTopMatch = genTopMatchIdxs[j];
					if(genTopMatch != -1){
					if(p == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << " with id " << _base->genpart_id->at(_genparts[genTopMatch].GetUserIdx()) << " and mass " << _genparts[genTopMatch].mass() << " and energy ratio " << _genAK4jets[j].E() / _genparts[genTopMatch].E() << endl;
						double gendR = dR(_genAK4jets[j].eta(), _genAK4jets[j].phi(), _genparts[genTopMatch].eta(), _genparts[genTopMatch].phi());
						_procCats[p].hists1D[0][56]->Fill(gendR);
						_procCats[p].hists1D[0][57]->Fill(_genAK4jets[j].E()/_genparts[genTopMatch].E());

					}	
					//dr bw gen jet and best exclusive gen W match
					int genWMatch = genWMatchIdxs[j];
					if(genWMatch != -1){
					if(p == 0) cout << " matched to gen W #" << genWMatchIdxs[j] << " with id " << _base->genpart_id->at(_genparts[genWMatch].GetUserIdx()) << " and mass " << _genparts[genWMatch].mass() << " and energy ratio " << _genAK4jets[j].E()/_genparts[genWMatch].E() << endl;
						double gendR = dR(_genAK4jets[j].eta(), _genAK4jets[j].phi(), _genparts[genWMatch].eta(), _genparts[genWMatch].phi());
						_procCats[p].hists1D[0][116]->Fill(gendR);
						_procCats[p].hists1D[0][117]->Fill(_genAK4jets[j].E()/_genparts[genWMatch].E());

					}	

				}
				_procCats[p].hists1D[0][106]->Fill((double)_genAK8jets.size());
				njets_lead = 0;
				njets_notlead = 0;
				for(int j = 0; j < _genAK8jets.size(); j++){
					if(_genAK8jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_genAK8jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][106]->Fill(njets_lead);
				_procCats[p].hists1D[2][106]->Fill(njets_notlead);
				//cout << "gen matching gen jets to particles - start" << endl;
				GenericMatchJet(_genAK8jets,_genparts,genTopMatchIdxs,6);
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genAK8jets.size(); j++){
					if(p == 0) cout << "gen AK8 jet #" << j << " phi " << _genAK8jets[j].phi() << " eta " << _genAK8jets[j].eta() << " energy " << _genAK8jets[j].E() <<  " mass " << _genAK8jets[j].mass() << " pt " << _genAK8jets[j].pt() << endl;
					_procCats[p].hists1D[0][107]->Fill(_genAK8jets[j].eta());
					_procCats[p].hists1D[0][108]->Fill(_genAK8jets[j].phi());
					_procCats[p].hists1D[0][109]->Fill(_genAK8jets[j].time());
					_procCats[p].hists1D[0][110]->Fill(_genAK8jets[j].pt());
					_procCats[p].hists1D[0][111]->Fill(_genAK8jets[j].mass());
					_procCats[p].hists1D[0][112]->Fill(_genAK8jets[j].E());
					_procCats[p].hists1D[0][113]->Fill(_base->AK8Jet_genNConstituents->at(_genAK8jets[j].GetUserIdx()));
					//dr bw gen jet and best exclusive gen top match
					//if(p == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << endl;
					int genTopMatch = genTopMatchIdxs[j];
					if(genTopMatch != -1){
						double gendR = dR(_genAK8jets[j].eta(), _genAK8jets[j].phi(), _genparts[genTopMatch].eta(), _genparts[genTopMatch].phi());
						_procCats[p].hists1D[0][114]->Fill(gendR);
						_procCats[p].hists1D[0][115]->Fill(_genAK8jets[j].E()/_genparts[genTopMatch].E());
						
					}	

				}
				
				_procCats[p].hists1D[0][81]->Fill((double)_genAK15jets.size());
				njets_lead = 0;
				njets_notlead = 0;
				for(int j = 0; j < _genAK15jets.size(); j++){
					if(_genAK15jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_genAK15jets[j].pt() >= _pt_thresh) njets_lead++;

				}
				_procCats[p].hists1D[1][81]->Fill(njets_lead);
				_procCats[p].hists1D[2][81]->Fill(njets_notlead);
				//cout << "gen matching gen jets to particles - start" << endl;
				GenericMatchJet(_genAK15jets, _genparts, genTopMatchIdxs,6);
				GenericMatchJet(_genAK15jets,_genparts,genWMatchIdxs, 24); //match gen AK4 jets to gen Ws
				//cout << "gen matching gen jets to particles - end" << endl;
				for(int j = 0; j < _genAK15jets.size(); j++){
					if(p == 0) cout << "gen AK15 jet #" << j << " phi " << _genAK15jets[j].phi() << " eta " << _genAK15jets[j].eta() << " energy " << _genAK15jets[j].E() <<  " mass " << _genAK15jets[j].mass() << " pt " << _genAK15jets[j].pt() << endl;
					_procCats[p].hists1D[0][75]->Fill(_genAK15jets[j].eta());
					_procCats[p].hists1D[0][76]->Fill(_genAK15jets[j].phi());
					_procCats[p].hists1D[0][77]->Fill(_genAK15jets[j].time());
					_procCats[p].hists1D[0][78]->Fill(_genAK15jets[j].pt());
					_procCats[p].hists1D[0][79]->Fill(_genAK15jets[j].mass());
					_procCats[p].hists1D[0][80]->Fill(_genAK15jets[j].E());
					//dr bw gen jet and best exclusive gen top match
					int genTopMatch = genTopMatchIdxs[j];
					//if(p == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << endl;
					if(genTopMatch != -1){
						double gendR = dR(_genAK15jets[j].eta(), _genAK15jets[j].phi(), _genparts[genTopMatch].eta(), _genparts[genTopMatch].phi());
						_procCats[p].hists1D[0][82]->Fill(gendR);
						_procCats[p].hists1D[0][83]->Fill(_genAK15jets[j].E()/_genparts[genTopMatch].E());
						
						
					}	
					
					//dr bw gen jet and best exclusive gen W match
					int genWMatch = genWMatchIdxs[j];
					if(p == 0) cout << " matched to gen W #" << genWMatchIdxs[j] << endl;
					if(genWMatch != -1){
						double gendR = dR(_genAK15jets[j].eta(), _genAK15jets[j].phi(), _genparts[genWMatch].eta(), _genparts[genWMatch].phi());
						_procCats[p].hists1D[0][118]->Fill(gendR);
						_procCats[p].hists1D[0][119]->Fill(_genAK15jets[j].E()/_genparts[genWMatch].E());
						
						
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
			vector<int> genMatchIdxs_jet(_recoAK4jets.size(),-1); //one per jet, follows same indexing as jets
			GenericMatchJet(_recoAK4jets,_genAK4jets,genMatchIdxs_jet);
			vector<int> genTopMatchIdxs(_recoAK4jets.size(), -1); //one per jet, follows same indexing as jets
			if(_genTop.size() > 0)
				GenericMatchJet(_recoAK4jets,_genTop,genTopMatchIdxs);
			vector<int> genWMatchIdxs(_recoAK4jets.size(), -1); //one per jet, follows same indexing as jets
			if(_genW.size() > 0)
				GenericMatchJet(_recoAK4jets,_genW,genWMatchIdxs);
			vector<int> genqMatchIdxs(_recoAK4jets.size(), -1); //one per jet, follows same indexing as jets
			if(_genq.size() > 0)
				GenericMatchJet(_recoAK4jets,_genq,genqMatchIdxs);
			vector<int> genGluonMatchIdxs(_recoAK4jets.size(), -1); //one per jet, follows same indexing as jets
			if(_genglu.size() > 0)
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
						
						if(pt == 0 && p == 0) cout << "reco AK4 jet #" << j << " phi " << _recoAK4jets[j].phi() << " eta " << _recoAK4jets[j].eta() << " energy " << _recoAK4jets[j].E() <<  " mass " << _recoAK4jets[j].mass() << " nConstituents " << _recoAK4jets[j].GetNConstituents() << " nRhs " << _recoAK4jets[j].GetNRecHits() << " pt " << _recoAK4jets[j].pt() << " jetsize " << jetsize << " px " << _recoAK4jets[j].px() << " py " << _recoAK4jets[j].py() << " pz " << _recoAK4jets[j].pz() << " mass from rhs " << _recoAK4jets[j].mass_rhs() << endl;
						_procCats[p].hists1D[pt][18]->Fill(jetsize);
						_procCats[p].hists1D[pt][64]->Fill(sqrt(jet_cov.at(0,0)));
						_procCats[p].hists1D[pt][65]->Fill(sqrt(jet_cov.at(1,1)));
						_procCats[p].hists1D[pt][66]->Fill(sqrt(jet_cov.at(2,2)));
						_procCats[p].hists1D[pt][19]->Fill(_recoAK4jets[j].e());
						_procCats[p].hists1D[pt][20]->Fill(_recoAK4jets[j].pt());
						_procCats[p].hists1D[pt][21]->Fill(_recoAK4jets[j].mass());
						_procCats[p].hists1D[pt][67]->Fill(_recoAK4jets[j].time());
						_procCats[p].hists1D[pt][68]->Fill(_recoAK4jets[j].eta());
						_procCats[p].hists1D[pt][69]->Fill(_recoAK4jets[j].phi());
						_procCats[p].hists2D[pt][21]->Fill((double)_recoAK4jets.size(), jetsize);
						
						//fill subcluster hists
						_procCats[p].hists1D[pt][28]->Fill(_recoAK4jets[j].GetNConstituents());
						_procCats[p].hists2D[pt][6]->Fill(_recoAK4jets[j].e(), _recoAK4jets[j].mass());
						_procCats[p].hists2D[pt][35]->Fill(_recoAK4jets[j].e(), jetsize);
						_procCats[p].hists2D[pt][8]->Fill(_recoAK4jets[j].GetNConstituents(), jetsize);	
						_procCats[p].hists2D[pt][9]->Fill(_recoAK4jets[j].e(), jetsize);	
						if(_recoAK4jets[j].GetNConstituents() == 0) cout << _recoAK4jets[j].GetNConstituents() << " n subcl " << _recoAK4jets[j].GetNRecHits() << " n rhs" << endl;
						//fill subcluster hists
						vector<Jet> consts = _recoAK4jets[j].GetConstituents();
						double rnk = 0;
						for(int c = 0; c < (int)consts.size(); c++){
							Jet subcl = consts[c];
							_procCats[p].hists1D[pt][29]->Fill(subcl.E());
							_procCats[p].hists1D[pt][30]->Fill(subcl.time());
							_procCats[p].hists1D[pt][31]->Fill(subcl.eta());
							_procCats[p].hists1D[pt][32]->Fill(subcl.phi());
							
							Matrix subcl_cov = subcl.GetCovariance();
							_procCats[p].hists1D[pt][33]->Fill(sqrt(subcl_cov.at(0,0)));
							_procCats[p].hists1D[pt][34]->Fill(sqrt(subcl_cov.at(1,1)));
							_procCats[p].hists1D[pt][35]->Fill(sqrt(subcl_cov.at(2,2)));
							_procCats[p].hists1D[pt][36]->Fill(subcl_cov.at(0,1));
							_procCats[p].hists1D[pt][37]->Fill(subcl_cov.at(0,2));
							_procCats[p].hists1D[pt][38]->Fill(subcl_cov.at(1,2));
					
							_procCats[p].hists1D[pt][126]->Fill(subcl.m());
							_procCats[p].hists2D[pt][24]->Fill(subcl.E(), subcl.m());
							//effective # of rhs is sum of weights bc weights are each the responsibility of rechit n to this subcl
							vector<double> ws;
							subcl.GetWeights(ws);
							double effnRhs = 0;
							for(int n = 0; n < ws.size(); n++){
								effnRhs += ws[n]; 
							}
							if(p == 0 && pt == 0) cout << " subcl #" << c << " has " << effnRhs << " # of effective rechits of " << ws.size() << " total rhs and mass " << subcl.m() << " and energy " << subcl.E() << " px " << subcl.px() << " py " << subcl.py() << " pz " << subcl.pz() << endl;
							_procCats[p].hists1D[pt][127]->Fill(effnRhs);
							_procCats[p].hists2D[pt][25]->Fill(subcl.E(), effnRhs);
							_procCats[p].hists2D[pt][26]->Fill(subcl.m(), effnRhs);

	
							//do dr calcs bw all subclusters
							for(int cc = c+1; cc < consts.size(); cc++){
								double dr = dR(consts[c].eta(), consts[c].phi(), consts[cc].eta(), consts[cc].phi());
								_procCats[p].hists1D[pt][72]->Fill(dr);
							}
				

						}
						vector<Jet> rhs;
						_recoAK4jets[j].GetJets(rhs);
						for(int r = 0; r < rhs.size(); r++){
							_procCats[p].hists1D[pt][40]->Fill(rhs[r].E());
							_procCats[p].hists1D[pt][41]->Fill(rhs[r].t());
						}
						vector<pair<double,double>> geoEavg_diffT;
						vector<JetPoint> rhs_pt = _recoAK4jets[j].GetJetPoints();
						CalcRhTimeDiff(rhs_pt,geoEavg_diffT);
						for(int r = 0; r < geoEavg_diffT.size(); r++){
						      _procCats[p].hists2D[pt][1]->Fill(geoEavg_diffT[r].first, geoEavg_diffT[r].second);
						}
						_procCats[p].hists1D[pt][41]->Fill(_recoAK4jets[j].GetNRecHits());
						//if no gen match, skip
						if(p == 0 && pt == 0) cout << "reco AK4 jet #" << j << " gen top idx " << genTopMatchIdxs[j] << " gen gluon idx " << genGluonMatchIdxs[j] << " gen q idx " << genqMatchIdxs[j] << " gen w idx " << genWMatchIdxs[j] << endl; 
						if(genMatchIdxs_jet[j] != -1){
							genidx = genMatchIdxs_jet[j];	
							//if(p == 0 && pt == 0) cout << " matched to gen AK4 jet #" << genMatchIdxs_jet[j] << endl;
							int genjetidx = _genAK4jets[genMatchIdxs_jet[j]].GetUserIdx();
							double genpt, pz, jetpt, jetpz, ratio_p;
							jetpt = _genAK4jets[genidx].pt();//_base->Jet_genPt->at(genjetidx);
							jetpz = _genAK4jets[genidx].pz();//_base->Jet_genPz->at(genjetidx);

				

							if(_sel != QCDdijets){
								int genpartidx = -1;
								int ngenparts_ptge5 = 0;
								for(int g = 0; g < _base->AK4Jet_genNConstituents->at(genjetidx); g++){
									genpartidx = _base->AK4Jet_genConstituentIdxs->at(genjetidx).at(g);
									if(genpartidx >= _base->genpart_pt->size()) break;
									genpt = _base->genpart_pt->at(genpartidx);
									pz = _base->genpart_pz->at(genpartidx);
									ratio_p = (sqrt(genpt*genpt + pz*pz))/(sqrt(jetpt*jetpt + jetpz*jetpz));	

	
									if(pt >= 5) ngenparts_ptge5++;	
								}
							}
						}
						//if no gen top match, skip
						//if(p == 0 && pt == 0) cout << " matched to gen top #" << genTopMatchIdxs[j] << endl;
						if(genTopMatchIdxs[j] != -1){
							int gentopidx = genTopMatchIdxs[j];
		
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genTop[gentopidx].eta(), _genTop[gentopidx].phi());
						}
						if(genGluonMatchIdxs[j] != -1){
							int genGluonidx = genGluonMatchIdxs[j];
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genglu[genGluonidx].eta(), _genglu[genGluonidx].phi());
							double eratio = _recoAK4jets[j].E() / _genglu[genGluonidx].E();
							_procCats[p].hists1D[pt][170]->Fill(dr);
							_procCats[p].hists1D[pt][171]->Fill(eratio);
						}
						if(genqMatchIdxs[j] != -1){
							int genqidx = genqMatchIdxs[j];
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genq[genqidx].eta(), _genq[genqidx].phi());
							double eratio = _recoAK4jets[j].E() / _genq[genqidx].E();
							_procCats[p].hists1D[pt][172]->Fill(dr);
							_procCats[p].hists1D[pt][173]->Fill(eratio);
						}
						//if no gen W match, skip
						//if(p == 0 && pt == 0) cout << " matched to gen W #" << genWMatchIdxs[j] << endl;
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];
		
							double dr = dR(_recoAK4jets[j].eta(), _recoAK4jets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							_procCats[p].hists1D[pt][174]->Fill(dr);
							_procCats[p].hists1D[pt][175]->Fill(_recoAK4jets[j].E()/_genW[genWidx].E());
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
								cout << "Warning: " << Wpartons.size() << " daughter particles found for W " << genWidx << " skipping hist filling" << endl;
								continue;
							}
							double gendR = dR(Wpartons[0].eta(), Wpartons[0].phi(), Wpartons[1].eta(), Wpartons[1].phi());
							_procCats[p].hists2D[pt][50]->Fill(gendR, jetsize);	
						}
					}
				}
			}
		}

		//_recoAK8 hists 
		void FillRecoAK8JetHists(){
			//do gen matching
			vector<int> genTopMatchIdxs(_recoAK8jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genTop.size() > 0)
				GenericMatchJet(_recoAK8jets,_genTop,genTopMatchIdxs);
			vector<int> genWMatchIdxs(_recoAK8jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genW.size() > 0)
				GenericMatchJet(_recoAK8jets,_genW,genWMatchIdxs);
			vector<int> genqMatchIdxs(_recoAK8jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genq.size() > 0)
				GenericMatchJet(_recoAK8jets,_genq,genqMatchIdxs);
			vector<int> genGluonMatchIdxs(_recoAK8jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genglu.size() > 0)
				GenericMatchJet(_recoAK8jets,_genglu,genGluonMatchIdxs);
			//cout << "final best matches" << endl;
			//for(int b = 0; b < genMatchIdxs_p.size(); b++){
			//	if(genMatchIdxs_p[b] != -1) cout << " jet " << b << " is exclusively matched to gen particle " << genMatchIdxs_p[b] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs_p[b]), _base->genpart_phi->at(genMatchIdxs_p[b]), _recoAK4jets[b].eta(), _recoAK4jets[b].phi()) << endl;
			//	 else cout << " jet " << b << " could not be gen matched" << endl;

			//}
			//else return; //not defined for other jet sizes and nominal FillRecoJetHists() fills for AK4
			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				int njets = _recoAK8jets.size();
				_procCats[p].hists1D[0][84]->Fill(njets);
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _recoAK8jets.size(); j++){
					if(_recoAK8jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_recoAK8jets[j].pt() >= _pt_thresh) njets_lead++;
				}
				_procCats[p].hists1D[1][84]->Fill(njets_lead);
				_procCats[p].hists1D[2][84]->Fill(njets_notlead);
				for(int j = 0; j < _recoAK8jets.size(); j++){
					for(int pt = 0; pt < _procCats[p].hists1D.size(); pt++){
						//define pt bins
						//pt == 1 -> [_pt_thresh,inf)
						if(pt == 1 && _recoAK8jets[j].pt() < _pt_thresh) continue;
						//pt == 2 -> [0,_pt_thresh)
						if(pt == 2 && _recoAK8jets[j].pt() >= _pt_thresh) continue;
	
						Matrix jet_cov = _recoAK8jets[j].GetCovariance();
						double jetsize = CalcSize(jet_cov);
						double jet_subleadsize = CalcSubleadSize(jet_cov);
						if(pt == 0 && p == 0) cout << "reco AK8 jet #" << j << " phi " << _recoAK8jets[j].phi() << " eta " << _recoAK8jets[j].eta() << " energy " << _recoAK8jets[j].E() <<  " mass " << _recoAK8jets[j].mass() << " nConstituents " << _recoAK8jets[j].GetNConstituents() << " nRhs " << _recoAK8jets[j].GetNRecHits() << " pt " << _recoAK8jets[j].pt() << " jetsize " << jetsize << endl;
						_procCats[p].hists1D[pt][85]->Fill(_recoAK8jets[j].eta());
						_procCats[p].hists1D[pt][86]->Fill(_recoAK8jets[j].phi_02pi());
						_procCats[p].hists1D[pt][87]->Fill(_recoAK8jets[j].time());
						_procCats[p].hists1D[pt][88]->Fill(_recoAK8jets[j].pt());
						_procCats[p].hists1D[pt][89]->Fill(_recoAK8jets[j].mass());
						_procCats[p].hists1D[pt][90]->Fill(_recoAK8jets[j].e());
						_procCats[p].hists1D[pt][91]->Fill(_recoAK8jets[j].GetNConstituents());
						_procCats[p].hists1D[pt][94]->Fill(jetsize);
						
						_procCats[p].hists2D[pt][45]->Fill(_recoAK8jets[j].m(), jetsize);

						//if no gen top match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen top #" << genTopMatchIdxs[j] << endl;
						if(p == 0 && pt == 0) cout << "reco AK8 jet #" << j << " gen top idx " << genTopMatchIdxs[j] << " gen gluon idx " << genGluonMatchIdxs[j] << " gen q idx " << genqMatchIdxs[j] << " gen w idx " << genWMatchIdxs[j] << endl; 
						if(genTopMatchIdxs[j] != -1){
							int genidx = genTopMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genTop[genidx].eta(), _genTop[genidx].phi());
							double eratio = _recoAK8jets[j].E()/_genTop[genidx].E();
							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to top " << genidx << " with dr " << dr << " and eratio " << eratio << endl;

							_procCats[p].hists1D[pt][92]->Fill(dr);
							_procCats[p].hists1D[pt][93]->Fill(eratio);

						}
						if(genGluonMatchIdxs[j] != -1){
							int genidx = genGluonMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genglu[genidx].eta(), _genglu[genidx].phi());
							double eratio = _recoAK8jets[j].E()/_genglu[genidx].E();
							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to gluon " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][157]->Fill(dr);
							_procCats[p].hists1D[pt][158]->Fill(eratio);


						}
						if(genqMatchIdxs[j] != -1){
							int genidx = genqMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genq[genidx].eta(), _genq[genidx].phi());
							double eratio = _recoAK8jets[j].E()/_genq[genidx].E();

							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to q " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][162]->Fill(dr);
							_procCats[p].hists1D[pt][163]->Fill(eratio);

	
						}

						//if no gen W match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen W #" << genWMatchIdxs[j] << endl;
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];	
							double dr = dR(_recoAK8jets[j].eta(), _recoAK8jets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							double eratio = _recoAK8jets[j].E()/_genW[genWidx].E();
							if(p == 0 && pt == 0) cout << "AK8 jet " << j << " matched to W " << genWidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][155]->Fill(dr);
							_procCats[p].hists1D[pt][156]->Fill(eratio);


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
							_procCats[p].hists2D[pt][48]->Fill(gendR, jetsize);	

						}
					}	
				}
			}
		}

		//_recoAK15 hists 
		void FillRecoAK15JetHists(){
			//do gen matching
			vector<int> genTopMatchIdxs(_recoAK15jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genTop.size() > 0)
				GenericMatchJet(_recoAK15jets,_genTop,genTopMatchIdxs);
			vector<int> genWMatchIdxs(_recoAK15jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genW.size() > 0)
				GenericMatchJet(_recoAK15jets,_genW,genWMatchIdxs);
			vector<int> genqMatchIdxs(_recoAK15jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genq.size() > 0)
				GenericMatchJet(_recoAK15jets,_genq,genqMatchIdxs);
			vector<int> genGluonMatchIdxs(_recoAK15jets.size(),-1); //one per jet, follows same indexing as jets
			if(_genglu.size() > 0)
				GenericMatchJet(_recoAK15jets,_genglu,genGluonMatchIdxs);
			//cout << "final best matches" << endl;
			//for(int b = 0; b < genMatchIdxs_p.size(); b++){
			//	if(genMatchIdxs_p[b] != -1) cout << " jet " << b << " is exclusively matched to gen particle " << genMatchIdxs_p[b] << " with dr " << dR(_base->genpart_eta->at(genMatchIdxs_p[b]), _base->genpart_phi->at(genMatchIdxs_p[b]), _recoAK4jets[b].eta(), _recoAK4jets[b].phi()) << endl;
			//	 else cout << " jet " << b << " could not be gen matched" << endl;

			for(int p = 0; p < _procCats.size(); p++){
				//cout << "process #" << p << ": " << _procCats[p].plotName << endl;
				int njets = _recoAK15jets.size();
				_procCats[p].hists1D[0][95]->Fill(njets);
				int njets_lead = 0;
				int njets_notlead = 0;
				for(int j = 0; j < _recoAK15jets.size(); j++){
					if(_recoAK15jets[j].pt() < _pt_thresh) njets_notlead++;
					//pt == 2 -> [0,_pt_thresh)
					if(_recoAK15jets[j].pt() >= _pt_thresh) njets_lead++;
				}
				_procCats[p].hists1D[1][95]->Fill(njets_lead);
				_procCats[p].hists1D[2][95]->Fill(njets_notlead);
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
						_procCats[p].hists1D[pt][96]->Fill(_recoAK15jets[j].eta());
						_procCats[p].hists1D[pt][97]->Fill(_recoAK15jets[j].phi_02pi());
						_procCats[p].hists1D[pt][98]->Fill(_recoAK15jets[j].time());
						_procCats[p].hists1D[pt][99]->Fill(_recoAK15jets[j].pt());
						_procCats[p].hists1D[pt][100]->Fill(_recoAK15jets[j].mass());
						_procCats[p].hists1D[pt][101]->Fill(_recoAK15jets[j].e());
						_procCats[p].hists1D[pt][102]->Fill(_recoAK15jets[j].GetNConstituents());
						_procCats[p].hists1D[pt][105]->Fill(jetsize);
						
						_procCats[p].hists2D[pt][46]->Fill(_recoAK15jets[j].m(), jetsize);

						//if no gen top match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen top #" << genTopMatchIdxs[j] << endl;
						if(p == 0 && pt == 0) cout << "reco AK15 jet #" << j << " gen top idx " << genTopMatchIdxs[j] << " gen gluon idx " << genGluonMatchIdxs[j] << " gen q idx " << genqMatchIdxs[j] << " gen w idx " << genWMatchIdxs[j] << endl; 
						if(genTopMatchIdxs[j] != -1){
							int genidx = genTopMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genTop[genidx].eta(), _genTop[genidx].phi());
							double eratio = _recoAK15jets[j].E()/_genTop[genidx].E();
							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to top " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	

							_procCats[p].hists1D[pt][103]->Fill(dr);
							_procCats[p].hists1D[pt][104]->Fill(eratio);

							
						
						}
						if(genGluonMatchIdxs[j] != -1){
							int genidx = genGluonMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genglu[genidx].eta(), _genglu[genidx].phi());
							double eratio = _recoAK15jets[j].E()/_genglu[genidx].E();
							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to gluon " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][166]->Fill(dr);
							_procCats[p].hists1D[pt][167]->Fill(eratio);


						}
						if(genqMatchIdxs[j] != -1){
							int genidx = genqMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genq[genidx].eta(), _genq[genidx].phi());
							double eratio = _recoAK15jets[j].E()/_genq[genidx].E();

							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to q " << genidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][168]->Fill(dr);
							_procCats[p].hists1D[pt][169]->Fill(eratio);

	
						}

						//if no gen W match, skip
						//if(p == 0 && pt == 0) cout << "reco AK" << AK << " jet #" << j << " matched to gen W #" << genWMatchIdxs[j] << endl;
						if(genWMatchIdxs[j] != -1){
							int genWidx = genWMatchIdxs[j];	
							double dr = dR(_recoAK15jets[j].eta(), _recoAK15jets[j].phi(), _genW[genWidx].eta(), _genW[genWidx].phi());
							double eratio = _recoAK15jets[j].E()/_genW[genWidx].E();
							if(p == 0 && pt == 0) cout << "AK15 jet " << j << " matched to W " << genWidx << " with dr " << dr << " and eratio " << eratio << endl;	
							_procCats[p].hists1D[pt][164]->Fill(dr);
							_procCats[p].hists1D[pt][165]->Fill(eratio);


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
							_procCats[p].hists2D[pt][49]->Fill(gendR, jetsize);	

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
			if(_procCats[0].hists2D[0][44]->Integral() == 0){ cout << "skipping all evt disps" << endl; return;} //don't write canvases if this event display isn't filled
			string plot_title;
			if(_sel == singW){
				plot_title = "single W^{#pm}";
				if(_oname.find("Wgluon") != string::npos){
					plot_title = "single W^{#pm}+gluon";
				}
			}
			if(_sel == QCDdijets){
				plot_title = "QCD dijets";
			}
			if(_sel == boostTop){
				plot_title = "t#bar{t}";
			}
			cout << "writing event display hist" << endl;
			ofile->cd();
			//write overall event display
			TCanvas* cv = new TCanvas("evtdisp","evtdisp");
			cv->SetTitle("");
			cv->SetRightMargin(0.15);
			if(_evt2disp_z == 1){ //time
				_procCats[0].hists2D[0][44]->GetZaxis()->SetTitle("time [ns]");
			}
			else if(_evt2disp_z % 2 == 0 && _evt2disp_z > 0){ //responsibility is event numbers where _evt2disp_z / 2 == k of subcl responsibility to display
				_procCats[0].hists2D[0][44]->GetZaxis()->SetTitle("responsibility");
			}
			else{
			}
			_procCats[0].hists2D[0][44]->SetTitle("");
			//do formatting
			_procCats[0].hists2D[0][44]->GetXaxis()->CenterTitle(true);
			_procCats[0].hists2D[0][44]->GetYaxis()->CenterTitle(true);
			_procCats[0].hists2D[0][44]->GetZaxis()->CenterTitle(true);
			_procCats[0].hists2D[0][44]->GetYaxis()->SetLabelFont(132);
			_procCats[0].hists2D[0][44]->GetXaxis()->SetLabelFont(132);
			_procCats[0].hists2D[0][44]->GetZaxis()->SetLabelFont(132);
			_procCats[0].hists2D[0][44]->GetYaxis()->SetTitleFont(132);
			_procCats[0].hists2D[0][44]->GetXaxis()->SetTitleFont(132);
			_procCats[0].hists2D[0][44]->GetZaxis()->SetTitleFont(132);
			_procCats[0].hists2D[0][44]->GetYaxis()->SetTitleSize(0.04);
			_procCats[0].hists2D[0][44]->GetXaxis()->SetTitleOffset(1.05);
			_procCats[0].hists2D[0][44]->GetXaxis()->SetTitleSize(0.04);
			_procCats[0].hists2D[0][44]->GetZaxis()->SetTitleSize(0.04);
			_procCats[0].hists2D[0][44]->Draw("colz1");
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
			//write labels
			string lat_cms = "#font[22]{Pythia 8} #font[132]{event generation, #sqrt{s} = 13 TeV}";
        		TLatex lat;
        		lat.SetNDC();
        		lat.SetTextSize(0.04);
        		lat.SetTextFont(42);
        		lat.DrawLatex(0.03,0.92,lat_cms.c_str());
			TLatex lat2;
        		lat2.SetNDC();
        		lat2.SetTextSize(0.04);
        		lat2.SetTextFont(42);
        		plot_title = "#font[132]{"+plot_title+"}";
			lat2.DrawLatex(0.8,0.92,plot_title.c_str());
			cv->Write();
		
			vector<string> names;
			for(auto it = center_coords.begin(); it != center_coords.end(); it++){
				names.push_back(it->first);
			}
cout << _evtdisps_obj.size() << " # of obj evtdisps" << endl;	
			//write object specific plots
			for(int h = 0; h < _evtdisps_obj.size(); h++){
				//get match string for center + width
				string name = _evtdisps_obj[h]->GetName();
				string objmatch = "EvtDisplay_etaCell_phiCell_";
				name = name.substr(objmatch.size());
				//skip hists that aren't filled
				cout << "hist #" << h << " name of obj evtdisp to be filled " << name << endl;
				if(find(names.begin(), names.end(), name) == names.end()) continue;
		
				/*
				vector<int> jet_idxs; //which jets to draw	
				if(name.find("W") != string::npos && _genW.size() > 0){
					vector<int> genWMatchIdxs(_predJets.size(),-1);
					GenericMatchJet(_predJets,_genW, genWMatchIdxs); //match BHC jets to good gen Ws
					for(auto idx : genWMatchIdxs){
						if(idx != -1) jet_idxs.push_back(idx);
					}
				}
				if(name.find("q") != string::npos && _genq.size() > 0){
					vector<int> genqMatchIdxs(_predJets.size(),-1);
					GenericMatchJet(_predJets,_genq, genqMatchIdxs); //match BHC jets to good gen qs
					for(auto idx : genqMatchIdxs){
						if(idx != -1) jet_idxs.push_back(idx);
					}
				}
				*/
				BayesPoint center = center_coords[name];
				double eta_max = 0;
				double phi_max = 0;
				double eta_min = 999;
				double phi_min = 999;
				BayesPoint width = BayesPoint({eta_max, phi_max});// - set by rhs drawn to get everything in frame = window_width[name];
				BayesPoint max_width = BayesPoint({eta_max, phi_max});// - set by rhs drawn to get everything in frame = window_width[name];
				BayesPoint min_width = BayesPoint({eta_min, phi_min});// - set by rhs drawn to get everything in frame = window_width[name];
cout << "drawing hist #" << h << " of " << _evtdisps_obj.size() << " with name " << name << endl;
				for(int j = 0; j < _predJets.size(); j++){
					//if(find(jet_idxs.begin(), jet_idxs.end(), j) == jet_idxs.end()) continue;
					BayesPoint ell_center({_predJets[j].eta(), _predJets[j].phi()});
					ell_center.Translate(center.at(0),0);
					ell_center.CircularTranslate(center.at(1),1);
					cout << "jet #" << j << " window width eta " << window_width[name].at(0) << " phi " << window_width[name].at(1) << " this jet center eta " << ell_center.at(0) << " phi " << ell_center.at(1) << endl;
					if(fabs(ell_center.at(0)) > window_width[name].at(0)) continue; //out of frame in eta
					if(fabs(ell_center.at(1)) > window_width[name].at(1)) continue; //out of frame in phi

cout << "drawing rhs from jet #" << j << endl;

					vector<JetPoint> rhs;
					if(_evt2disp_z % 2 == 0 && _evt2disp_z > 0){
						int k = _evt2disp_z / 2;
						if(k >= _predJets[j].GetNConstituents()){
							rhs = _predJets[j].GetJetPoints();
						}
						else{
							Jet subcl = _predJets[j].GetConstituent(k);
							rhs = subcl.GetJetPoints();
						}
					}
					else{
						rhs = _predJets[j].GetJetPoints();
					}
					for(auto rh : rhs){
						double w;
						if(_evt2disp_z == 0){ //energy
							w = rh.E();
						}
						else if(_evt2disp_z == 1){ //time
							w = rh.t();
						}
						else if(_evt2disp_z % 2 == 0){ //responsibility is event numbers where _evt2disp_z / 2 == k of subcl responsibility to display
							w = rh.GetWeight();
						}
						else{
							w = rh.E(); //default energy weighted
						}
						if(w == 0) continue;
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
//cout << "at jet #" << j << " rh eta " << rh_pt.at(0) << " rh_phi " << rh_pt.at(1) << " eta_max " << eta_max << " eta_min " << eta_min << " phi_max " << phi_max << " phi_min " << phi_min << endl;
						_evtdisps_obj[h]->Fill(rh_pt.at(0), rh_pt.at(1), rh_pt.w());
					}
				}
//TODO - check eta/phi_max/min values (see drawn canvases)
cout << "eta_max " << eta_max << " eta_min " << eta_min << " phi_max " << phi_max << " phi_min " << phi_min << endl;
				if(_evt2disp_z == 1){ //time
					_evtdisps_obj[h]->GetZaxis()->SetTitle("time [ns]");
				}
				else if(_evt2disp_z % 2 == 0 && _evt2disp_z > 0){ //responsibility is event numbers where _evt2disp_z / 2 == k of subcl responsibility to display
					_evtdisps_obj[h]->GetZaxis()->SetTitle("responsibility");
				}
				else{
				}
				min_width.SetValue(eta_min, 0);
				min_width.SetValue(phi_min, 1);
				max_width.SetValue(eta_max, 0);
				max_width.SetValue(phi_max, 1);
				width.SetValue(max(fabs(eta_max), fabs(eta_min)), 0);
				width.SetValue(max(fabs(phi_max), fabs(phi_min)), 1);
				if(_evtdisps_obj[h]->GetEntries() == 0) continue; //don't draw if not filled for this particle gen obj
				TCanvas* cv_obj = new TCanvas(_evtdisps_obj[h]->GetName(),_evtdisps_obj[h]->GetTitle());
				cv_obj->SetRightMargin(0.15);
				cv_obj->cd();
				_evtdisps_obj[h]->GetXaxis()->SetRangeUser(min_width.at(0), max_width.at(0));
				_evtdisps_obj[h]->GetYaxis()->SetRangeUser(min_width.at(1), max_width.at(1));
				_evtdisps_obj[h]->SetTitle("");
				_evtdisps_obj[h]->GetXaxis()->CenterTitle(true);
				_evtdisps_obj[h]->GetYaxis()->CenterTitle(true);
				_evtdisps_obj[h]->GetZaxis()->CenterTitle(true);
				_evtdisps_obj[h]->GetYaxis()->SetLabelFont(132);
				_evtdisps_obj[h]->GetXaxis()->SetLabelFont(132);
				_evtdisps_obj[h]->GetZaxis()->SetLabelFont(132);
				_evtdisps_obj[h]->GetYaxis()->SetTitleFont(132);
				_evtdisps_obj[h]->GetXaxis()->SetTitleFont(132);
				_evtdisps_obj[h]->GetZaxis()->SetTitleFont(132);
				_evtdisps_obj[h]->GetYaxis()->SetTitleSize(0.04);
				_evtdisps_obj[h]->GetXaxis()->SetTitleOffset(1.05);
				_evtdisps_obj[h]->GetXaxis()->SetTitleSize(0.04);
				_evtdisps_obj[h]->GetZaxis()->SetTitleSize(0.04);
				_evtdisps_obj[h]->Draw("colz1");
cout << "hist for " << name << " integral " << _evtdisps_obj[h]->Integral() << " entries " << _evtdisps_obj[h]->GetEntries() << endl;
				//do for gen particles too 
				for(int m = 0; m < _plot_particles.size(); m++){
					BayesPoint m_center({_plot_particles[m].GetX(), _plot_particles[m].GetY()});
					m_center.Translate(center.at(0),0);
					m_center.CircularTranslate(center.at(1),1);
				
					//double dr = dR(m_center.at(0), m_center.at(1), 0., 0.);
					//if(dr > sqrt(width.at(0)*width.at(0) + width.at(1)*width.at(1))) continue;
					if(m_center.at(0) > max_width.at(0) || m_center.at(0) < min_width.at(0)) continue;
					if(m_center.at(1) > max_width.at(1) || m_center.at(1) < min_width.at(1)) continue;
					_plot_particles[m].DrawMarker(m_center.at(0), m_center.at(1));

				}
				for(int j = 0; j < _predJets.size(); j++){
					//if(find(jet_idxs.begin(), jet_idxs.end(), j) == jet_idxs.end()) continue;
					BayesPoint ell_center({_jellipses[j].GetX1(), _jellipses[j].GetY1()});
					ell_center.Translate(center.at(0),0);
					ell_center.CircularTranslate(center.at(1),1);
					cout << "ellipse (jet) #" << j << " window width eta " << width.at(0) << " phi " << width.at(1) << " this ellipse center eta " << ell_center.at(0) << " phi " << ell_center.at(1) << endl;
					//make sure ellipse center is in window
					if(ell_center.at(0) > max_width.at(0)) continue; //out of frame in eta
					if(ell_center.at(1) > max_width.at(1)) continue; //out of frame in eta
					if(ell_center.at(0) < min_width.at(0)) continue; //out of frame in eta
					if(ell_center.at(1) < min_width.at(1)) continue; //out of frame in eta
					
					double ell_center_eta = ell_center.at(0);//_jellipses[j].GetX1(); 
					double ell_center_phi = ell_center.at(1);//_jellipses[j].GetY1();
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
					
					//if(ell_center.at(0) + r_eta > max_width.at(0)) continue; 
					//if(ell_center.at(0) - r_eta < min_width.at(0)) continue; 

					//if(ell_center.at(1) + r_phi > max_width.at(1)) continue; 
					//if(ell_center.at(1) - r_phi < min_width.at(1)) continue; 
					//if(fabs(r_eta) > fabs(width.at(0))) continue;
					//if(fabs(r_phi) > fabs(width.at(1))) continue;

				cout << "drawing ellipse (jet) #" << j << endl; 

		
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

        			TLatex lat1;
        			lat1.SetNDC();
        			lat1.SetTextSize(0.04);
        			lat1.SetTextFont(42);
        			lat1.DrawLatex(0.03,0.92,lat_cms.c_str());
				TLatex lat3;
        			lat3.SetNDC();
        			lat3.SetTextSize(0.04);
        			lat3.SetTextFont(42);
        			lat3.DrawLatex(0.8,0.92,plot_title.c_str());
				cv_obj->SetTitle("");
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
		TH1D* nClusters = new TH1D("BHCJet_nJets","BHCJet_nJets",10,0,10);
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
		//10 - # pred jets - # reco jets
		TH1D* predGen_nJets = new TH1D("BHCRecoAK4_diffNJets","BHCRecoAK4_diffNJets",20,-10,10);
		//11 - eta sigma
		TH1D* predJet_subClusterEtaSig = new TH1D("BHCJet_subclusterEtaSig","BHCJet_subclusterEtaSig",50,0.,1.);
		//12 - phi sigma
		TH1D* predJet_subClusterPhiSig = new TH1D("BHCJet_subclusterPhiSig","BHCJet_subclusterPhiSig",50,0.,1.);
		//13 - time sigma
		TH1D* predJet_subClusterTimeSig = new TH1D("BHCJet_subclusterTimeSig","BHCJet_subclusterTimeSig",50,-1.,15.);
		//14 - etaphi cov
		TH1D* predJet_subClusteretaPhiCov = new TH1D("BHCJet_subclusterEtaPhiCov","BHCJet_subclusterEtaPhiCov",50,-0.0005,0.0005);
		//15 - time-eta covariance
		TH1D* predJet_subClustertimeEtaCov = new TH1D("BHCJet_subclusterTimeEtaCov","BHCJet_subclusterTimeEtaCov",50,-0.05,0.05);
		//16 - time-phi covariance
		TH1D* predJet_subClustertimePhiCov = new TH1D("BHCJet_subclusterTimePhiCov","BHCJet_subclusterTimePhiCov",50,-0.05,0.05);
		//17 - n reco AK4 jets
		TH1D* nRecoJets = new TH1D("recoAK4Jet_nJets","recoAK4Jet_nJets",10,0,10);
		//18 - reco AK4 jet size
		TH1D* recoJet_jetSize = new TH1D("recoAK4Jet_jetSize","recoAK4Jet_jetSize",50,0,2.);
		//19 - reco AK4 jet energy
		TH1D* recoJet_energy = new TH1D("recoAK4Jet_energy","recoAK4Jet_energy",25,0,500);
		//20 - reco AK4 jet pt
		TH1D* recoJet_pt = new TH1D("recoAK4Jet_pt","recoAK4Jet_pt",25,0,1000);
		//21 - reco AK4 jet mass
		TH1D* recoJet_mass = new TH1D("recoAK4Jet_mass","recoAK4Jet_mass",50,0,250);
		//BREAK!!!!! -2
		//22 - 24 - eta sigma for jet
		TH1D* predJet_EtaVar = new TH1D("BHCJet_EtaSig","BHCJet_EtaSig",50,0.,1.);
		//23 - 25 - phi sigma for jet
		TH1D* predJet_PhiVar = new TH1D("BHCJet_PhiSig","BHCJet_PhiSig",50,0., 1.);
		//24 - 26 - time sigma for jet
		TH1D* predJet_TimeVar = new TH1D("BHCJet_TimeSig","BHCJet_TimeSig",50,0.,10.);
		//25 - 27 - eta-phi covariance for jet
		TH1D* predJet_etaPhiCov = new TH1D("BHCJet_etaPhiCov","BHCJet_etaPhiCov",25,-0.0005,0.0005);
		//26 - 28 - time-eta covariance for jet
		TH1D* predJet_timeEtaCov = new TH1D("BHCJet_timeEtaCov","BHCJet_timeEtaCov",25,-0.2,0.2);
		//27 - 29 - time-phi covariance for jet
		TH1D* predJet_timePhiCov = new TH1D("BHCJet_timePhiCov","BHCJet_timePhiCov",25,-0.2,0.2);
		//28 - 30 - n GMM clusters in reco jets
		TH1D* reco_nSubclusters = new TH1D("recoAK4Jet_nSubclustersJet","recoAK4Jet_nSubclustersJet",30,0,30);
		//29 - 31 - energy per GMM cluster from reco jets
		TH1D* recoJet_subClusterEnergy = new TH1D("recoAK4Jet_subclusterEnergy","recoAK4Jet_subclusterEnergy",50,0,500);
		//30 - 32 - time center of GMM cluster from reco jets
		TH1D* recoJet_subClusterTimeCenter = new TH1D("recoAK4Jet_subclusterTimeCenter","recoAK4Jet_subclusterTimeCenter",25,-1,1);
		//31 - 33 - eta center of GMM cluster from reco jets
		TH1D* recoJet_subClusterEtaCenter = new TH1D("recoAK4Jet_subclusterEtaCenter","recoAK4Jet_subclusterEtaCenter",25,-3.2,3.2);
		//32 - 34 - phi center of GMM cluster from reco jets
		TH1D* recoJet_subClusterPhiCenter = new TH1D("recoAK4Jet_subclusterPhiCenter","recoAK4Jet_subclusterPhiCenter",25,0.,8*atan(1));
		//33 - 35 - eta sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterEtaSig = new TH1D("recoAK4Jet_subclusterEtaSig","recoAK4Jet_subclusterEtaSig",50,0.,0.5);
		//34 - 36 - phi sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterPhiSig = new TH1D("recoAK4Jet_subclusterPhiSig","recoAK4Jet_subclusterPhiSig",50,0.,0.5);
		//35 - 37 - time sigma of GMM cluster from reco jets
		TH1D* recoJet_subClusterTimeSig = new TH1D("recoAK4Jet_subclusterTimeSig","recoAK4Jet_subclusterTimeSig",50,0.,10.);
		//36 - 38 - eta-phi covariance of GMM cluster from reco jets
		TH1D* recoJet_subClusteretaPhiCov = new TH1D("recoAK4Jet_subclusterEtaPhiCov","recoAK4Jet_subclusterEtaPhiCov",50,-0.0005,0.0005);
		//37 - 39 - time-eta covariance of GMM cluster from reco jets
		TH1D* recoJet_subClustertimeEtaCov = new TH1D("recoAK4Jet_subclusterTimeEtaCov","recoAK4Jet_subclusterTimeEtaCov",50,-0.05,0.05);
		//38 - 40 - time-phi covariance of GMM cluster from reco jets
		TH1D* recoJet_subClustertimePhiCov = new TH1D("recoAK4Jet_subclusterTimePhiCov","recoAK4Jet_subclusterTimePhiCov",50,-0.05,0.05);
		//39 - 41 - n rhs in reco jets
		TH1D* recoJet_nRhs = new TH1D("recoAK4Jet_nRhs","recoAK4Jet_nRhs",300,0,300);
		//40 - 42 - rh energy in reco jets
		TH1D* recoJet_rhE = new TH1D("recoAK4Jet_rhE","recoAK4Jet_rhE",50,0,300);
		//BREAK!!!!! -8
		//41 - 49 - ak4 jet rh times
		TH1D* AK4Jet_rhTimes = new TH1D("recoAK4Jet_rhTimes","recoAK4Jet_rhTimes",200,-25,25);
		//BREAK!!!!! -11
		//42 - 53 - # gen partons (t, b, q from W)
		TH1D* nGenParticles = new TH1D("nGenParticles","nGenParticles",20,0,20);
		//43 - 54 - gen particle eta at detector
		TH1D* genParticle_eta = new TH1D("genParticle_eta","genParticle_eta",25,-3.2,3.2);
		//44 - 55 - 62 - gen particle phi at detector
		TH1D* genParticle_phi = new TH1D("genParticle_phi","genParticle_phi",25,0.,8*atan(1));
		//45 - 56 - 63 - gen particle time at detector
		TH1D* genParticle_time = new TH1D("genParticle_time","genParticle_time",25,-10,10);
		//46 - 57 - gen particle pt		
		TH1D* genParticle_pt = new TH1D("genParticle_pt","genParticle_pt",25,0,500);
		//47 - 58 - gen particle mass
		TH1D* genParticle_mass = new TH1D("genParticle_mass","genParticle_mass",25,0,250);
		//48 - 59 - gen particle energy
		TH1D* genParticle_energy = new TH1D("genParticle_energy","genParticle_energy",50,0,500);
		//49 - 60 - gen AK4 jet eta at detector
		TH1D* genAK4Jet_eta = new TH1D("genAK4Jet_EtaCenter","genAK4Jet_EtaCenter",25,-3.2,3.2);
		//50 - 61 - gen AK4 jet phi at detector
		TH1D* genAK4Jet_phi = new TH1D("genAK4Jet_PhiCenter","genAK4Jet_PhiCenter",25,0.,8*atan(1));
		//51 - 62 - gen AK4 jet time at detector
		TH1D* genAK4Jet_time = new TH1D("genAK4Jet_TimeCenter","genAK4Jet_TimeCenter",25,-1,1);
		//52 - 63 - gen AK4 jet pt		
		TH1D* genAK4Jet_pt = new TH1D("genAK4Jet_pt","genAK4Jet_pt",25,0,1000);
		//53 - 64 - gen AK4 jet mass
		TH1D* genAK4Jet_mass = new TH1D("genAK4Jet_mass","genAK4Jet_mass",50,0,250);
		//54 - 65 - gen AK4 jet energy
		TH1D* genAK4Jet_energy = new TH1D("genAK4Jet_energy","genAK4Jet_energy",25,0,500);
		//55 - 66 - # gen AK4 jets
		TH1D* nJet_genAK4Jet = new TH1D("genAK4Jet_nJets","genAK4Jet_nJets",10,0,10);
		//BREAK!!!!! -13
		//56 - 69 - dR bw gen jet and gen top its exclusively matched to
		TH1D* genAK4JetTop_dR = new TH1D("genAK4Jet_genTop_dR","genAK4Jet_genTop_dR",25,0,1.5);
		//57 - 70 - E ratio bw gen jet and gen top its exclusively matched to - gen jet energy/gen top energy
		TH1D* genAK4JetTop_Eratio = new TH1D("genAK4Jet_genTop_Eratio","genAK4Jet_genTop_Eratio",25,0,2);
		//BREAK!!!!! -14
		//58 - 72 - dR bw bhc jet and gen top its exclusively matched to
		TH1D* BHCJetTop_dR = new TH1D("BHCJet_genTop_dR","BHCJet_genTop_dR",25,0,1.5);
		//59 - 73 - E ratio bw bhc jet and gen top its exclusively matched to - bhc jet energy/gen top energy
		TH1D* BHCJetTop_Eratio = new TH1D("BHCJet_genTop_Eratio","BHCJet_genTop_Eratio",25,0,2);
		//60 - 74 - bhc jet eta center
		TH1D* BHCJet_EtaCenter = new TH1D("BHCJet_EtaCenter","BHCJet_EtaCenter",25,-3.2,3.2);
		//61 - 75 - bhc jet phi center
		TH1D* BHCJet_PhiCenter = new TH1D("BHCJet_PhiCenter","BHCJet_PhiCenter",25,0.,8*atan(1));
		//62 - 76 - bhc jet center
		TH1D* BHCJet_TimeCenter = new TH1D("BHCJet_TimeCenter","BHCJet_TimeCenter",25,-1,1);
		//63 - 77 - rh time
		TH1D* rhTime = new TH1D("rhTime","rhTime",25,-10,10);
		//64 - 78 - recoAK4 jet rh eta sig
		TH1D* recoAK4Jet_rhEtaSig = new TH1D("recoAK4Jet_EtaSig","recoAK4_EtaSig",50,0,1);
		//65 - 79 - recoAK4 jet rh phi sig
		TH1D* recoAK4Jet_rhPhiSig = new TH1D("recoAK4Jet_PhiSig","recoAK4_PhiSig",50,0,1);
		//66 - 80 - recoAK4 jet rh time sig
		TH1D* recoAK4Jet_rhTimeSig = new TH1D("recoAK4Jet_TimeSig","recoAK4_TimeSig",50,0,15.);
		//67 - 81 - reco AK4 jet center
		TH1D* recoAK4Jet_TimeCenter = new TH1D("recoAK4Jet_TimeCenter","recoAK4Jet_TimeCenter",25,-1,1);
		//68 - 82 - reco AK4 jet eta at detector
		TH1D* recoAK4Jet_EtaCenter = new TH1D("recoAK4Jet_EtaCenter","recoAK4Jet_EtaCenter",25,-3.2,3.2);
		//69 - 83 - reco AK4 jet phi at detector
		TH1D* recoAK4Jet_PhiCenter = new TH1D("recoAK4Jet_PhiCenter","recoAK4Jet_PhiCenter",25,0.,8*atan(1));
		//70 -84 - reco AK4 jet # subclusters in event
		TH1D* recoAK4Jet_nSubclustersEvt = new TH1D("recoAK4Jet_nSubclustersEvt","recoAK4Jet_nSubclustersEvt",30,0,30);
		//71 - 85 - dr bw subclusters in BHC jet
		TH1D* BHCJet_drSubclusters = new TH1D("BHCJet_drSubclusters","BHCJet_drSubclusters",50,0,0.1);
		//72 - 86 - dr bw subclusters in reco AK4 jet
		TH1D* recoAK4Jet_drSubclusters = new TH1D("recoAK4Jet_drSubclusters","recoAK4Jet_drSubclusters",50,0,0.1);
		//73 - 87 - rh energy in reco jets
		TH1D* BHCJet_rhE = new TH1D("BHCJet_rhE","BHCJet_rhE",50,0,300);
		//BREAK!!!! -16
		//74 - 90 - BHC jet n rhs
		TH1D* BHCJet_nRhs = new TH1D("BHCJet_nRhs","BHCJet_nRhs",300,0,300);
		//75 - 91 - gen AK15 jet eta at detector
		TH1D* genAK15Jet_eta = new TH1D("genAK15Jet_EtaCenter","genAK15Jet_EtaCenter",25,-3.2,3.2);
		//76 - 92 - gen AK15 jet phi at detector
		TH1D* genAK15Jet_phi = new TH1D("genAK15Jet_PhiCenter","genAK15Jet_PhiCenter",25,0.,8*atan(1));
		//77 - 93 - gen AK15 jet time at detector
		TH1D* genAK15Jet_time = new TH1D("genAK15Jet_TimeCenter","genAK15Jet_TimeCenter",25,-1,1);
		//78 - 94 - gen AK15 jet pt		
		TH1D* genAK15Jet_pt = new TH1D("genAK15Jet_pt","genAK15Jet_pt",25,0,1000);
		//79 - 95 - gen AK15 jet mass
		TH1D* genAK15Jet_mass = new TH1D("genAK15Jet_mass","genAK15Jet_mass",50,0,250);
		//80 - 96 - gen AK15 jet energy
		TH1D* genAK15Jet_energy = new TH1D("genAK15Jet_energy","genAK15Jet_energy",25,0,500);
		//81 - 97 - # gen AK15 jets
		TH1D* nJet_genAK15Jet = new TH1D("genAK15Jet_nJets","genAK15Jet_nJets",10,0,10);
		//BREAK!!!! -18
		//82 - 100 - dR bw gen jet and gen top its exclusively matched to
		TH1D* genAK15JetTop_dR = new TH1D("genAK15Jet_genTop_dR","genAK15Jet_genTop_dR",25,0,1.5);
		//83 - 101 - E ratio bw gen jet and gen top its exclusively matched to - gen jet energy/gen top energy
		TH1D* genAK15JetTop_Eratio = new TH1D("genAK15Jet_genTop_Eratio","genAK15Jet_genTop_Eratio",25,0,2);
		//84 - 102 - n reco AK8 jets
		TH1D* nRecoAK8Jets = new TH1D("recoAK8Jet_nJets","recoAK8Jet_nJets",10,0,10);
		//85 - 103 - reco AK8 jet eta at detector
		TH1D* recoAK8Jet_eta = new TH1D("recoAK8Jet_EtaCenter","recoAK8Jet_EtaCenter",25,-3.2,3.2);
		//86 - 104 - reco AK8 jet phi at detector
		TH1D* recoAK8Jet_phi = new TH1D("recoAK8Jet_PhiCenter","recoAK8Jet_PhiCenter",25,0.,8*atan(1));
		//87 - 105 - reco AK8 jet time at detector
		TH1D* recoAK8Jet_time = new TH1D("recoAK8Jet_TimeCenter","recoAK8Jet_TimeCenter",25,-1,1);
		//88 - 106 - reco AK8 jet pt		
		TH1D* recoAK8Jet_pt = new TH1D("recoAK8Jet_pt","recoAK8Jet_pt",25,0,1000);
		//89 - 107 - reco AK8 jet mass
		TH1D* recoAK8Jet_mass = new TH1D("recoAK8Jet_mass","recoAK8Jet_mass",50,0,250);
		//90 - 108 - reco AK8 jet energy
		TH1D* recoAK8Jet_energy = new TH1D("recoAK8Jet_energy","recoAK8Jet_energy",25,0,500);
		//91 - 109 - # constituents per reco AK8 jet
		TH1D* recoAK8Jet_nConstituents = new TH1D("recoAK8Jet_nSubclustersJet","recoAK8Jet_nSubclustersJet",30,0,30);
		//92 - 110 - dR bw reco AK8 jet and gen top its exclusively matched to
		TH1D* recoAK8JetTop_dR = new TH1D("recoAK8Jet_genTop_dR","recoAK8Jet_genTop_dR",25,0,1.5);
		//93 - 111 - E ratio bw reco AK8 jet and gen top its exclusively matched to - reco AK8 jet energy/gen top energy
		TH1D* recoAK8JetTop_Eratio = new TH1D("recoAK8Jet_genTop_Eratio","recoAK8Jet_genTop_Eratio",25,0,2);
		//94 - 112 - reco AK8 jet size
		TH1D* recoAK8Jet_jetSize = new TH1D("recoAK8Jet_jetSize","recoAK8Jet_jetSize",50,0,2.);
		//95 - 113 - n reco AK15 jets
		TH1D* nRecoAK15Jets = new TH1D("recoAK15Jet_nJets","recoAK15Jet_nJets",10,0,10);
		//96 - 114 - gen AK15 jet eta at detector
		TH1D* recoAK15Jet_eta = new TH1D("recoAK15Jet_EtaCenter","recoAK15Jet_EtaCenter",25,-3.2,3.2);
		//97 - 115 - reco AK15 jet phi at detector
		TH1D* recoAK15Jet_phi = new TH1D("recoAK15Jet_PhiCenter","recoAK15Jet_PhiCenter",25,0.,8*atan(1));
		//98 - 116 - reco AK15 jet time at detector
		TH1D* recoAK15Jet_time = new TH1D("recoAK15Jet_TimeCenter","recoAK15Jet_TimeCenter",25,-1,1);
		//99 - 117 - reco AK15 jet pt		
		TH1D* recoAK15Jet_pt = new TH1D("recoAK15Jet_pt","recoAK15Jet_pt",25,0,1000);
		//100 - 118 - reco AK15 jet mass
		TH1D* recoAK15Jet_mass = new TH1D("recoAK15Jet_mass","recoAK15Jet_mass",50,0,250);
		//101 - 119 - reco AK15 jet energy
		TH1D* recoAK15Jet_energy = new TH1D("recoAK15Jet_energy","recoAK15Jet_energy",25,0,500);
		//102 - 120 - # constituents per reco AK15 jet
		TH1D* recoAK15Jet_nConstituents = new TH1D("recoAK15Jet_nSubclustersJet","recoAK15Jet_nSubclustersJet",30,0,30);
		//103 - 121 - dR bw reco AK15 jet and gen top its exclusively matched to
		TH1D* recoAK15JetTop_dR = new TH1D("recoAK15Jet_genTop_dR","recoAK15Jet_genTop_dR",25,0,1.5);
		//104 - 122 - E ratio bw reco AK15 jet and gen top its exclusively matched to - reco AK15 jet energy/gen top energy
		TH1D* recoAK15JetTop_Eratio = new TH1D("recoAK15Jet_genTop_Eratio","recoAK15Jet_genTop_Eratio",25,0,2);
		//105 - 123 - reco AK15 jet size
		TH1D* recoAK15Jet_jetSize = new TH1D("recoAK15Jet_jetSize","recoAK15Jet_jetSize",50,0,2.);
		//106 - 124 - n gen AK8 jets
		TH1D* nGenAK8Jets = new TH1D("genAK8Jet_nJets","genAK8Jet_nJets",10,0,10);
		//107 - 125 - gen AK15 jet eta at detector
		TH1D* genAK8Jet_eta = new TH1D("genAK8Jet_EtaCenter","genAK8Jet_EtaCenter",25,-3.2,3.2);
		//108 - 126 - gen AK8 jet phi at detector
		TH1D* genAK8Jet_phi = new TH1D("genAK8Jet_PhiCenter","genAK8Jet_PhiCenter",25,0.,8*atan(1));
		//109 - 127 - gen AK8 jet time at detector
		TH1D* genAK8Jet_time = new TH1D("genAK8Jet_TimeCenter","genAK8Jet_TimeCenter",25,-1,1);
		//110 - 128 - gen AK8 jet pt		
		TH1D* genAK8Jet_pt = new TH1D("genAK8Jet_pt","genAK8Jet_pt",25,0,1000);
		//111 - 129 - gen AK8 jet mass
		TH1D* genAK8Jet_mass = new TH1D("genAK8Jet_mass","genAK8Jet_mass",50,0,250);
		//112 - 130 - gen AK8 jet energy
		TH1D* genAK8Jet_energy = new TH1D("genAK8Jet_energy","genAK8Jet_energy",25,0,500);
		//113 - 131 - # constituents per gen AK8 jet
		TH1D* genAK8Jet_nConstituents = new TH1D("genAK8Jet_nConstituents","genAK8Jet_nConstituents",50,0,50);
		//114 - 132 - dR bw gen AK8 jet and gen top its exclusively matched to
		TH1D* genAK8JetTop_dR = new TH1D("genAK8Jet_genTop_dR","genAK8Jet_genTop_dR",25,0,1.5);
		//115 - 133 - E ratio bw gen AK8 jet and gen top its exclusively matched to - gen AK8 jet energy/gen top energy
		TH1D* genAK8JetTop_Eratio = new TH1D("genAK8Jet_genTop_Eratio","genAK8Jet_genTop_Eratio",25,0,2);
		//116 - 134 - dR bw gen AK4 jet and gen top its exclusively matched to
		TH1D* genAK4JetW_dR = new TH1D("genAK4Jet_genW_dR","genAK4Jet_genW_dR",25,0,1.5);
		//117 - 135 - E ratio bw gen AK4 jet and gen top its exclusively matched to - gen AK4 jet energy/gen top energy
		TH1D* genAK4JetW_Eratio = new TH1D("genAK4Jet_genW_Eratio","genAK4Jet_genW_Eratio",25,0,2.);
		//118 - 136 - dR bw gen AK15 jet and gen top its exclusively matched to
		TH1D* genAK15JetW_dR = new TH1D("genAK15Jet_genW_dR","genAK15Jet_genW_dR",25,0,1.5);
		//119 - 137 - E ratio bw gen AK15 jet and gen top its exclusively matched to - gen AK15 jet energy/gen top energy
		TH1D* genAK15JetW_Eratio = new TH1D("genAK15Jet_genW_Eratio","genAK15Jet_genW_Eratio",25,0,2.);
		//120 - 138 - dR bw BHC jet and gen W its exclusively matched to
		TH1D* BHCJetW_dR = new TH1D("BHCJet_genW_dR","BHCJet_genW_dR",25,0,1.5);
		//121 - 139 - E ratio bw BHC jet and gen W its exclusively matched to - BHC jet energy/gen W energy
		TH1D* BHCJetW_Eratio = new TH1D("BHCJet_genW_Eratio","BHCJet_genW_Eratio",25,0,2.);
		//122 - 140 - # subclusters in BHC jets matched to Ws
		TH1D* BHCJetW_nSubclusters = new TH1D("BHCJetW_nSubclusters","BHCJetW_nSubclusters",10,0,10);
		//123 - 141 - subcluster energy in BHC jets matched to Ws
		TH1D* BHCJetW_subClusterEnergy = new TH1D("BHCJetW_subclusterEnergy","BHCJetW_subclusterEnergy",25,0,500);
		//124 - 142 - BHC jet subcluster mass
		TH1D* BHCJet_subClusterMass = new TH1D("BHCJet_subclusterMass","BHCJet_subclusterMass",25,0,200);
		//125 - 143 - BHC jet subcluster # effective rechits
		TH1D* BHCJet_subClusterEffnRhs = new TH1D("BHCJet_subclusterEffnRhs","BHCJet_subclusterEffnRhs",25,0,150);
		//126 - 144 - reco AK4 jet subcluster mass
		TH1D* recoAK4Jet_subClusterMass = new TH1D("recoAK4Jet_subclusterMass","recoAK4Jet_subclusterMass",25,0,200);
		//127 - 145 - reco AK4 jet subcluster # effective rechits
		TH1D* recoAK4Jet_subClusterEffnRhs = new TH1D("recoAK4Jet_subclusterEffnRhs","recoAK4Jet_subclusterEffnRhs",25,0,150);
		//128 - 146 - subcluster mass in BHC jets matched to Ws
		TH1D* BHCJetW_subClusterMass = new TH1D("BHCJetW_subclusterMass","BHCJetW_subclusterMass",25,0,200);
		//129 - 147 - BHC jets - invariant mass of lead two subclusters (for jets with at least 2 subclusters)
		TH1D* BHCJetW_subClusterLeadInvMass = new TH1D("BHCJetW_subclusterLeadInvMass","BHCJetW_subclusterLeadInvMass",50,0,250);
		//130 - 148 - BHC jets - # ghost leftover after clustering
		TH1D* BHCJet_nGhosts = new TH1D("BHCJet_nGhosts","BHCJet_nGhosts",10,0,10);		
		//131 - 149 - BHC jets - ghost subcl energy
		TH1D* BHCJet_ghostSubClusterEnergy = new TH1D("BHCJet_ghostSubclusterEnergy","BHCJet_ghostSubclusterEnergy",25,0,500);
		//132 - 150 - BHC jets - ghost subcl eff # rhs
		TH1D* BHCJet_ghostSubClusterEffnRhs = new TH1D("BHCJet_ghostSubclusterEffnRhs","BHCJet_ghostSubclusterEffnRhs",25,0,200);
		//133 - 151 - BHC jets - gen-matched W - Eratio (reco/gen) of gen partons in W decay and 2 lead subclusters in BHC jet
		TH1D* BHCJetW_subclParton_dR = new TH1D("BHCJetW_subclParton_dR","BHCJetW_subclParton_dR",25,0,2);
		//134 - 152 - BHC jets - gen-matched W - Eratio (reco/gen) of gen partons in W decay and 2 lead subclusters in BHC jet
		TH1D* BHCJetW_subclParton_Eratio = new TH1D("BHCJetW_subclParton_Eratio","BHCJetW_subclParton_Eratio",25,0,2);
		//135 - 153 - # subclusters in BHC jets matched to Tops
		TH1D* BHCJetTop_nSubclusters = new TH1D("BHCJetTop_nSubclusters","BHCJetTop_nSubclusters",10,0,10);
		//136 - 154 - subcluster mass in BHC jets matched to Tops
		TH1D* BHCJetTop_subClusterMass = new TH1D("BHCJetTop_subclusterMass","BHCJetTop_subclusterMass",25,0,200);
		//137 - 155 - BHC jets - invariant mass of lead two subclusters (for jets with at least 2 subclusters)
		TH1D* BHCJetTop_subClusterLeadInvMass = new TH1D("BHCJetTop_subclusterLeadInvMass","BHCJetTop_subclusterLeadInvMass",50,0,250);
		//138 - 156 - dR bw BHC jet and gen q its exclusively matched to
		TH1D* BHCJetq_dR = new TH1D("BHCJet_genq_dR","BHCJet_genq_dR",25,0,1.5);
		//139 - 157 - E ratio bw BHC jet and gen q its exclusively matched to - BHC jet energy/gen q energy
		TH1D* BHCJetq_Eratio = new TH1D("BHCJet_genq_Eratio","BHCJet_genq_Eratio",25,0,2.5);
		//140 - 158 - # subclusters in BHC jets matched to qs
		TH1D* BHCJetq_nSubclusters = new TH1D("BHCJetq_nSubclusters","BHCJetq_nSubclusters",10,0,10);
		//141 - 159 - subcluster mass in BHC jets matched to qs
		TH1D* BHCJetq_subClusterMass = new TH1D("BHCJetq_subclusterMass","BHCJetq_subclusterMass",25,0,200);
		//142 - 160 - BHC jets - gen-matched to W - subcluster eta center		
		TH1D* BHCJetW_subclEtaCenter = new TH1D("BHCJetW_subclEtaCenter","BHCJetW_subclEtaCenter",25,-3.2,3.2);
		//143 - 161 - BHC jets - gen-matched to W - subcluster phi center		
		TH1D* BHCJetW_subclPhiCenter = new TH1D("BHCJetW_subclPhiCenter","BHCJetW_subclPhiCenter",25,-3.2,3.2);
		//144 - 162 - BHC jets - gen-matched to W - subcluster time center		
		TH1D* BHCJetW_subclTimeCenter = new TH1D("BHCJetW_subclTimeCenter","BHCJetW_subclTimeCenter",25,-1.,1);
		//145 - 163 - BHC jets - gen-matched to W - eta sigma of GMM cluster 
		TH1D* BHCJetW_subClusterEtaSig = new TH1D("BHCJetW_subclusterEtaSig","BHCJetW_subclusterEtaSig",50,0.,0.5);
		//146 - 164 - BHC jets - gen-matched to W - phi sigma of GMM cluster 
		TH1D* BHCJetW_subClusterPhiSig = new TH1D("BHCJetW_subclusterPhiSig","BHCJetW_subclusterPhiSig",50,0.,0.5);
		//147 - 165 - BHC jets - gen-matched to W - time sigma of GMM cluster 
		TH1D* BHCJetW_subClusterTimeSig = new TH1D("BHCJetW_subclusterTimeSig","BHCJetW_subclusterTimeSig",50,0.,5.);
		//148 - 166 - BHC jets - gen-matched to W - eta-phi covariance of GMM cluster 
		TH1D* BHCJetW_subClusteretaPhiCov = new TH1D("BHCJetW_subclusterEtaPhiCov","BHCJetW_subclusterEtaPhiCov",50,-0.0005,0.0005);
		//149 - 167 - BHC jets - gen-matched to W - time-eta covariance of GMM cluster 
		TH1D* BHCJetW_subClustertimeEtaCov = new TH1D("BHCJetW_subclusterTimeEtaCov","BHCJetW_subclusterTimeEtaCov",50,-0.05,0.05);
		//150 - 168 - BHC jets - gen-matched to W - time-phi covariance of GMM cluster 
		TH1D* BHCJetW_subClustertimePhiCov = new TH1D("BHCJetW_subclusterTimePhiCov","BHCJetW_subclusterTimePhiCov",50,-0.05,0.05);
		//151 - 169 - high mass + W-matched BHC jets - pt of subclusters gen-matched to W partons
		TH1D* BHCJetW_highMass_partonMatchSubclPt = new TH1D("BHCJetW_highMass_partonMatchSubclPt","BHCJetW_highMass_partonMatchSubclPt;SubclPt",25,0,500);
		//152 - 170 - high mass + W-matched BHC jets - pt of subclusters NOT gen-matched W partons
		TH1D* BHCJetW_highMass_partonNoMatchSubclPt = new TH1D("BHCJetW_highMass_partonNoMatchSubclPt","BHCJetW_highMass_partonNoMatchSubclPt;SubclPt",25,0,500);
		//153 - 171 - high mass + W-matched BHC jets - subclSize of subclusters gen-matched to W partons
		TH1D* BHCJetW_highMass_partonMatchSubclSize = new TH1D("BHCJetW_highMass_partonMatchSubclSize","BHCJetW_highMass_partonMatchSubclSize;SubclSize",25,0,1.5);
		//154 - 172 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched W partons
		TH1D* BHCJetW_highMass_partonNoMatchSubclSize = new TH1D("BHCJetW_highMass_partonNoMatchSubclSize","BHCJetW_highMass_partonNoMatchSubclSize;SubclSize",25,0,1.5);
		//155 - 173 - dR bw reco AK8 jet and gen W its exclusively matched to
		TH1D* recoAK8JetW_dR = new TH1D("recoAK8Jet_genW_dR","recoAK8Jet_genW_dR",25,0,1.5);
		//156 - 174 - E ratio bw reco AK8 jet and gen W its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK8JetW_Eratio = new TH1D("recoAK8Jet_genW_Eratio","recoAK8Jet_genW_Eratio",25,0,2.);
		//157 - 175 - dR bw reco AK8 jet and gen gluon its exclusively matched to
		TH1D* recoAK8JetGluon_dR = new TH1D("recoAK8Jet_genGluon_dR","recoAK8Jet_genGluon_dR",25,0,1.5);
		//158 - 176 - E ratio bw reco AK8 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK8JetGluon_Eratio = new TH1D("recoAK8Jet_genGluon_Eratio","recoAK8Jet_genGluon_Eratio",25,0,2);
		//159 - 177 - dR bw bhc jet and gen gluon its exclusively matched to
		TH1D* BHCJetGluon_dR = new TH1D("BHCJet_genGluon_dR","BHCJet_genGluon_dR",25,0,1.5);
		//160 - 178 - E ratio bw bhc jet and gen gluon its exclusively matched to - bhc jet energy/gen top energy
		TH1D* BHCJetGluon_Eratio = new TH1D("BHCJet_genGluon_Eratio","BHCJet_genGluon_Eratio",25,0,2);
		//161 - 179 - high mass + W-matched BHC jets - # subclusters
		TH1D* BHCJetW_highMass_nSubclustersJet = new TH1D("BHCJetW_highMass_nSubclustersJet","BHCJetW_highMass_nSubclustersJet",10,0,10);
		//162 - 180 - dR bw reco AK8 jet and gen q its exclusively matched to
		TH1D* recoAK8Jetq_dR = new TH1D("recoAK8Jet_genq_dR","recoAK8Jet_genq_dR",25,0,1.5);
		//163 - 181 - E ratio bw reco AK8 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK8Jetq_Eratio = new TH1D("recoAK8Jet_genq_Eratio","recoAK8Jet_genq_Eratio",25,0,2.5);
		//164 - 182 - dR bw reco AK15 jet and gen W its exclusively matched to
		TH1D* recoAK15JetW_dR = new TH1D("recoAK15Jet_genW_dR","recoAK15Jet_genW_dR",25,0,1.5);
		//165 - 183 - E ratio bw reco AK15 jet and gen W its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK15JetW_Eratio = new TH1D("recoAK15Jet_genW_Eratio","recoAK15Jet_genW_Eratio",25,0,2.);
		//166 - 184 - dR bw reco AK15 jet and gen gluon its exclusively matched to
		TH1D* recoAK15JetGluon_dR = new TH1D("recoAK15Jet_genGluon_dR","recoAK15Jet_genGluon_dR",25,0,1.5);
		//167 - 185 - E ratio bw reco AK15 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK15JetGluon_Eratio = new TH1D("recoAK15Jet_genGluon_Eratio","recoAK15Jet_genGluon_Eratio",25,0,2);
		//168 - 186 - dR bw reco AK15 jet and gen q its exclusively matched to
		TH1D* recoAK15Jetq_dR = new TH1D("recoAK15Jet_genq_dR","recoAK15Jet_genq_dR",25,0,1.5);
		//169 - 187 - E ratio bw reco AK15 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK15Jetq_Eratio = new TH1D("recoAK15Jet_genq_Eratio","recoAK15Jet_genq_Eratio",25,0,2.5);
		//170 - 188 - dR bw reco AK4 jet and gen gluon its exclusively matched to
		TH1D* recoAK4JetGluon_dR = new TH1D("recoAK4Jet_genGluon_dR","recoAK4Jet_genGluon_dR",25,0,1.5);
		//171 - 189 - E ratio bw reco AK4 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK4JetGluon_Eratio = new TH1D("recoAK4Jet_genGluon_Eratio","recoAK4Jet_genGluon_Eratio",25,0,2);
		//172 - 190 - dR bw reco AK4 jet and gen q its exclusively matched to
		TH1D* recoAK4Jetq_dR = new TH1D("recoAK4Jet_genq_dR","recoAK4Jet_genq_dR",25,0,1.5);
		//173 - 191 - E ratio bw reco AK4 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK4Jetq_Eratio = new TH1D("recoAK4Jet_genq_Eratio","recoAK4Jet_genq_Eratio",25,0,2.5);
		//174 - 192 - dR bw reco AK4 jet and gen W its exclusively matched to
		TH1D* recoAK4JetW_dR = new TH1D("recoAK4Jet_genW_dR","recoAK4Jet_genW_dR",25,0,1.5);
		//175 - 193 - E ratio bw reco AK4 jet and gen gluon its exclusively matched to - reco jet energy/gen top energy
		TH1D* recoAK4JetW_Eratio = new TH1D("recoAK4Jet_genW_Eratio","recoAK4Jet_genW_Eratio",25,0,2.);
		//176 - 194 - high mass + q-matched BHC jets - pt of subclusters gen-matched to W partons
		TH1D* BHCJetq_highMass_partonMatchSubclPt = new TH1D("BHCJetq_highMass_partonMatchSubclPt","BHCJetq_highMass_partonMatchSubclPt;SubclPt",25,0,500);
		//177 - 195 - high mass + q-matched BHC jets - pt of subclusters NOT gen-matched W partons
		TH1D* BHCJetq_highMass_partonNoMatchSubclPt = new TH1D("BHCJetq_highMass_partonNoMatchSubclPt","BHCJetq_highMass_partonNoMatchSubclPt;SubclPt",25,0,500);
		//178 - 196 - high mass + q-matched BHC jets - subclSize of subclusters gen-matched to W partons
		TH1D* BHCJetq_highMass_partonMatchSubclSize = new TH1D("BHCJetq_highMass_partonMatchSubclSize","BHCJetq_highMass_partonMatchSubclSize;SubclSize",25,0,1.5);
		//179 - 197 - high mass + q-matched BHC jets - subclSize of subclusters NOT gen-matched W partons
		TH1D* BHCJetq_highMass_partonNoMatchSubclSize = new TH1D("BHCJetq_highMass_partonNoMatchSubclSize","BHCJetq_highMass_partonNoMatchSubclSize;SubclSize",25,0,1.5);
		//180 - 198 - low mass + W-matched BHC jets - # subclusters
		TH1D* BHCJetW_lowMass_nSubclustersJet = new TH1D("BHCJetW_lowMass_nSubclustersJet","BHCJetW_lowMass_nSubclustersJet",10,0,10);
		//181 - 199 - high mass + W-matched BHC jets - # subclusters
		TH1D* BHCJetW_Wmass_nSubclustersJet = new TH1D("BHCJetW_Wmass_nSubclustersJet","BHCJetW_Wmass_nSubclustersJet",10,0,10);
		//182 - 200 - # subclusters in BHC jets matched to Gluons
		TH1D* BHCJetGluon_nSubclusters = new TH1D("BHCJetGluon_nSubclusters","BHCJetGluon_nSubclusters",10,0,10);
		//183 - 201 - BHC jets - gen-matched q - Eratio (reco/gen) of gen q to lead subcluster in BHC jet
		TH1D* BHCJetq_subclParton_dR = new TH1D("BHCJetq_subclParton_dR","BHCJetq_subclParton_dR",25,0,2);
		//184 - 202 - BHC jets - gen-matched q - Eratio (reco/gen) of gen q to lead subcluster in BHC jet
		TH1D* BHCJetq_subclParton_Eratio = new TH1D("BHCJetq_subclParton_Eratio","BHCJetq_subclParton_Eratio",25,0,2);
		//185 - 203 - BHC jets - gen-matched gluon - Eratio (reco/gen) of gen gluon to lead subcluster in BHC jet
		TH1D* BHCJetGluon_subclParton_dR = new TH1D("BHCJetGluon_subclParton_dR","BHCJetGluon_subclParton_dR",25,0,2);
		//186 - 204 - BHC jets - gen-matched gluon - Eratio (reco/gen) of gen gluon to lead subcluster in BHC jet
		TH1D* BHCJetGluon_subclParton_Eratio = new TH1D("BHCJetGluon_subclParton_Eratio","BHCJetGluon_subclParton_Eratio",25,0,2);
		//187 - 205 - high mass + W-matched BHC jets - E of subclusters gen-matched to W partons / E jet
		TH1D* BHCJetW_highMass_partonMatchSubclEOvJetE = new TH1D("BHCJetW_highMass_partonMatchSubclEOvJetE","BHCJetW_highMass_partonMatchSubclEOvJetE;SubclEOvJetE",25,0,1.2);
		//188 - 206 - high mass + W-matched BHC jets - E of subclusters NOT gen-matched W partons / E jet
		TH1D* BHCJetW_highMass_partonNoMatchSubclEOvJetE = new TH1D("BHCJetW_highMass_partonNoMatchSubclEOvJetE","BHCJetW_highMass_partonNoMatchSubclEOvJetE;SubclEOvJetE",25,0,1.2);
		//189 - 207 - high mass + W-matched BHC jets - subclSize of subclusters gen-matched to W partons / size jet
		TH1D* BHCJetW_highMass_partonMatchSubclSizeOvJetSize = new TH1D("BHCJetW_highMass_partonMatchSubclSizeOvJetSize","BHCJetW_highMass_partonMatchSubclSizeOvJetSize;SubclSizeOvJetSize",25,0,5.);
		//190 - 208 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched W partons / size jet
		TH1D* BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize = new TH1D("BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize","BHCJetW_highMass_partonNoMatchSubclSizeOvJetSize;SubclSizeOvJetSize",25,0,5.);
		//191 - 209 - PU-downweighted bhc jet mass for jets within dR = 1 of gen W
		TH1D* BHCJet_PUdownweighted_genFrameW_mass = new TH1D("BHCJet_PUdownweighted_genFrameW_mass","BHCJet_PUdownweighted_genFrameW_mass",50,0,250);
		//192 - 210 - PU-removed bhc jet mass for jets within dR = 1 of gen W
		TH1D* BHCJet_PUremoved_mass = new TH1D("BHCJet_PUremoved_mass","BHCJet_PUremoved_mass",50,0,250);
		//193 - 211 - dR bw parton and PU-cleaned subclusters for bhc jets matched to Ws
		TH1D* BHCJetW_subclParton_PUcleaned_dR = new TH1D("BHCJetW_subclParton_PUcleaned_dR","BHCJetW_subclParton_PUcleaned_dR",25,0,2);
		//194 - 212 - # bhc jets with ==2 cleaned subclusters matched to Ws / # bhc jets matched to Ws
		TH1D* nBHCJetsW_eq2cleanedSubcls = new TH1D("BHCJetW_nJets_eq2cleanedSubcls","BHCJetW_nJets_eq2cleanedSubcls",30,0,30);
		//195 - 213 - time of PU-like subclusters from bhc jets
		TH1D* BHCJet_subclTimeCenter_PUlike = new TH1D("BHCJet_subclTimeCenter_PUlike","BHCJet_subclTimeCenter_PUlike",25,-1.,1.);
		//196 - 214 - time of PU-cleaned suclusters from bhc jets
		TH1D* BHCJet_subclTimeCenter_PUcleaned = new TH1D("BHCJet_subclTimeCenter_PUcleaned","BHCJet_subclTimeCenter_PUcleaned",25,-1.,1);
		//197 - 215 - time var of PU-like subclusters from bhc jets
		TH1D* BHCJet_subclTimeSig_PUlike = new TH1D("BHCJet_subclTimeSig_PUlike","BHCJet_subclusterTimeSig_PUlike",50,0.,5.);
		//198 - 216 - time var of PU-cleaned suclusters from bhc jets
		TH1D* BHCJet_subclTimeSig_PUcleaned = new TH1D("BHCJet_subclTimeSig_PUcleaned","BHCJet_subclusterTimeSig_PUcleaned",50,0.,5.);
		//199 - PU-downweighted bhc jet mass for jets within dR = 1 of gen W for # subclusters = 2
		TH1D* BHCJet_PUdownweighted_nSubcleq2_mass = new TH1D("BHCJet_PUdownweighted_nSubcleq2_mass","BHCJet_PUdownweighted_nSubcleq2_mass",50,0,250);
		//200 - PU-downweighted bhc jet mass for jets within dR = 1 of gen W for # subclusters < 2
		TH1D* BHCJet_PUdownweighted_nSubcllt2_mass = new TH1D("BHCJet_PUdownweighted_nSubcllt2_mass","BHCJet_PUdownweighted_nSubcllt2_mass",50,0,250);
		//201 - PU-downweighted bhc jet mass for jets within dR = 1 of gen W for # subclusters > 2
		TH1D* BHCJet_PUdownweighted_nSubclgt2_mass = new TH1D("BHCJet_PUdownweighted_nSubclgt2_mass","BHCJet_PUdownweighted_nSubclgt2_mass",50,0,250);
		//202 - high mass + W-matched BHC jets - time sig of subclusters gen-matched to W partons / time sig jet
		TH1D* BHCJetW_highMass_partonMatchSubclTimeSigOvJetTimeSig = new TH1D("BHCJetW_highMass_partonMatchSubclTimeSigOvJetTimeSig","BHCJetW_highMass_partonMatchSubclTimeSigOvJetTimeSig;SubclTimeSigOvJetTimeSig",25,0,5.);
		//203 - high mass + W-matched BHC jets - time sig of subclusters NOT gen-matched W partons / time sig jet
		TH1D* BHCJetW_highMass_partonNoMatchSubclTimeSigOvJetTimeSig = new TH1D("BHCJetW_highMass_partonNoMatchSubclTimeSigOvJetTimeSig","BHCJetW_highMass_partonNoMatchSubclTimeSigOvJetTimeSig;SubclTimeSigOvJetTimeSig",25,0,5.);
		//204 - highMass + q-matched BHC jets - time sig of subclusters gen-matched to qs / time sig jet
		TH1D* BHCJetq_highMass_partonMatchSubclTimeSigOvJetTimeSig = new TH1D("BHCJetq_highMass_partonMatchSubclTimeSigOvJetTimeSig","BHCJetq_highMass_partonMatchSubclTimeSigOvJetTimeSig;SubclTimeSigOvJetTimeSig",25,0,5.);
		//205 - highMass + q-matched BHC jets - time sig of subclusters NOT gen-matched qs / time sig jet
		TH1D* BHCJetq_highMass_partonNoMatchSubclTimeSigOvJetTimeSig = new TH1D("BHCJetq_highMass_partonNoMatchSubclTimeSigOvJetTimeSig","BHCJetq_highMass_partonNoMatchSubclTimeSigOvJetTimeSig;SubclTimeSigOvJetTimeSig",25,0,5.);
		//206 - highMass + q-matched BHC jets - E of subclusters gen-matched to qs / E jet
		TH1D* BHCJetq_highMass_partonMatchSubclEOvJetE = new TH1D("BHCJetq_highMass_partonMatchSubclEOvJetE","BHCJetq_highMass_partonMatchSubclEOvJetE;SubclEOvJetE",25,0,1.2);
		//207 - highMass + q-matched BHC jets - E of subclusters NOT gen-matched qs / E jet
		TH1D* BHCJetq_highMass_partonNoMatchSubclEOvJetE = new TH1D("BHCJetq_highMass_partonNoMatchSubclEOvJetE","BHCJetq_highMass_partonNoMatchSubclEOvJetE;SubclEOvJetE",25,0,1.2);
		//208 - highMass + top-matched BHC jets - time sig of subclusters gen-matched to tops / time sig jet
		TH1D* BHCJetTop_highMass_partonMatchSubclTimeSigOvJetTimeSig = new TH1D("BHCJetTop_highMass_partonMatchSubclTimeSigOvJetTimeSig","BHCJetTop_highMass_partonMatchSubclTimeSigOvJetTimeSig;SubclTimeSigOvJetTimeSig",25,0,5.);
		//209 - highMass + top-matched BHC jets - time sig of subclusters NOT gen-matched tops / time sig jet
		TH1D* BHCJetTop_highMass_partonNoMatchSubclTimeSigOvJetTimeSig = new TH1D("BHCJetTop_highMass_partonNoMatchSubclTimeSigOvJetTimeSig","BHCJetTop_highMass_partonNoMatchSubclTimeSigOvJetTimeSig;SubclTimeSigOvJetTimeSig",25,0,5.);
		//210 - highMass + top-matched BHC jets - E of subclusters gen-matched to tops / E jet
		TH1D* BHCJetTop_highMass_partonMatchSubclEOvJetE = new TH1D("BHCJetTop_highMass_partonMatchSubclEOvJetE","BHCJetTop_highMass_partonMatchSubclEOvJetE;SubclEOvJetE",25,0,1.2);
		//211 - highMass + top-matched BHC jets - E of subclusters NOT gen-matched tops / E jet
		TH1D* BHCJetTop_highMass_partonNoMatchSubclEOvJetE = new TH1D("BHCJetTop_highMass_partonNoMatchSubclEOvJetE","BHCJetTop_highMass_partonNoMatchSubclEOvJetE;SubclEOvJetE",25,0,1.2);
		//212 - PU-downweighted bhc jet mass for jets within dR = 1 of gen W
		TH1D* BHCJet_PUdownweighted_mass = new TH1D("BHCJet_PUdownweighted_mass","BHCJet_PUdownweighted_mass",50,0,250);
		//213 - # jets per event with mass < 30
		TH1D* BHCJet_massle30_nJets = new TH1D("BHCJet_massle30_nJets","BHCJet_massle30_nJets",10,0,10);
		//214 - invariant mass bw jets with mass < 30 and PU cleaned
		TH1D* BHCJet_massle30_PUremoved_SubclInvMass = new TH1D("BHCJet_massle30_PUremoved_SubclInvMass","BHCJet_massle30_PUremoved_SubclInvMass",50,0,250);
		//215 - # subclusters for jets with mass < 30 and PU cleaned	
		TH1D* BHCJet_massle30_PUremoved_nSubclusters = new TH1D("BHCJet_massle30_PUremoved_nSubclusters","BHCJet_massle30_PUremoved_nSubclusters",10,0,10);


		////////////////////////////////////////////////////////////////////////
		////////////////////////////////2D plots////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		
		//0 - 2D histogram for recoGen pT resolution as a function of gen jet energy 
		TH2D* jetGenE_diffDeltaPt_recoGen = new TH2D("jetGenE_diffDeltaPt_recoGen","jetGenE_diffDeltaPt_recoGen;jet_{gen} E (GeV);#Delta p_{T}_{reco, gen} (GeV)",4,&xbins_recoGenPt[0],50,-50,50);
		//BREAK!!!! -17
		//1 - 18 - geo energy avg vs difference in time for adjacent crystals in same obj w/in 10% energy 
		TH2D* geoEavg_diffDeltaTime_adjRhs = new TH2D("geoEavg_diffDeltaTime_adjRhs","geoEavg_diffDeltaTime_adjRhs;geoEavg;diffDeltaTime;a.u.",xbins.size()-1,&xbins[0],25,-5,5);
		//2 - 22 - BHC jet multiplicity vs jet size
		TH2D* BHCJet_nJets_jetSize = new TH2D("BHCJet_nJets_jetSize","BHCJet_nJets_jetSize;nJets;jetSize",30,0,30,50,0,2);
		//3 - 27 - # subclusters vs jet mass for BHC jets
		TH2D* BHCJet_nSubclustersJet_mass = new TH2D("BHCJet_nSubclustersJet_mass","BHCJet_nSubclustersJet_mass;nSubclustersJet;mass",30,0,30,50,0,250);
		//4 - 28 - # subclusters vs jet energy for BHC jets
		TH2D* BHCJet_nSubclustersJet_energy = new TH2D("BHCJet_nSubclustersJet_energy","BHCJet_nSubclustersJet_energy;nSubclusters;energy",30,0,30,50,0,2000);
		//5 - 29 - # subclusters/evt vs # subclusters/jet for BHC jets
		TH2D* BHCJet_nSubclustersEvt_nJet = new TH2D("BHCJet_nSubclustersEvt_nJet","BHCJet_nSubclustersEvt_nJet;nSubclustersEvt;nJet",30,0,30,10,0,10);
		//6 - 31 - reco AK4 jet energy vs jet mass
		TH2D* recoAK4Jet_jetEnergy_jetMass = new TH2D("recoAK4Jet_jetEnergy_jetMass","recoAK4Jet_jetEnergy_jetMass;jetEnergy;jetMass",50,0,2000,50,0,250);
		//7 - 32 - BHC jet energy vs mass
		TH2D* BHCJet_jetEnergy_jetMass = new TH2D("BHCJet_jetEnergy_jetMass","BHCJet_jetEnergy_jetMass;jetEnergy;jetMass",50,0,2000,50,0,250);
		//8 - 33 - reco AK4 jets # subclusters vs jet size
		TH2D* recoAK4Jet_nSubclusters_jetSize = new TH2D("recoAK4Jet_nSubclustersJet_jetSize","recoAK4Jet_nSubclustersJet_jetSize;nSubclustersJet;jetsize",30,0,30,50,0,2);
		//9 - 34 - reco AK4 jet energy vs jet size
		TH2D* recoAK4Jet_jetEnergy_jetSize = new TH2D("recoAK4Jet_jetEnergy_jetSize","recoAK4Jet_jetEnergy_jetSize;jetEnergy;jetSize",50,0,500,50,0,2);
		//10 - 35 - BHC jet energy vs jet size
		TH2D* BHCJet_jetEnergy_jetSize = new TH2D("BHCJet_jetEnergy_jetSize","BHCJet_jetEnergy_jetSize;jetEnergy;jetSize",50,0,500,50,0,2);
		//11 - 56 - BHC jet pt vs gen top pt 
		TH2D* BHCJetPt_genTopPt = new TH2D("BHCJetPt_genTopPt","BHCJetPt_genTopPt;BHCJetPt;genTopPt",25,0,2000,25,0,2000);
		//12 - 57- BHC jet E vs gen top pt
		TH2D* BHCJetE_genTopE = new TH2D("BHCJetE_genTopE","BHCJetE_genTopE;BHCJetE;genTopE",25,0,2000,25,0,2000);
		//13 - 58 - BHC jet mass vs gen top jet mass 
		TH2D* BHCJetMass_genTopMass = new TH2D("BHCJetMass_genTopMass","BHCJetMass_genTopMass;BHCJetMass;genTopMass",25,0,250,25,0,250);
		//14 - 59 - BHC jet eta vs gen top eta 
		TH2D* BHCJetEtaCenter_genTopEtaCenter = new TH2D("BHCJetEtaCenter_genTopEtaCenter","BHCJetEtaCenter_genTopEtaCenter;BHCJetEtaCenter;genTopEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//15 - 60 - BHC jet phi vs gen top phi 
		TH2D* BHCJetPhiCenter_genTopPhiCenter = new TH2D("BHCJetPhiCenter_genTopPhiCenter","BHCJetPhiCenter_genTopPhiCenter;BHCJetPhiCenter;genTopPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//16 - 81 - BHC jet pt vs gen W pt 
		TH2D* BHCJetPt_genWPt = new TH2D("BHCJetPt_genWPt","BHCJetPt_genWPt;BHCJetPt;genWPt",25,0,500,25,0,500);
		//17 - 82 - BHC jet E vs gen W pt
		TH2D* BHCJetE_genWE = new TH2D("BHCJetE_genWE","BHCJetE_genWE;BHCJetE;genWE",25,0,2000,25,0,2000);
		//18 - 83 - BHC jet mass vs gen W jet mass 
		TH2D* BHCJetMass_genWMass = new TH2D("BHCJetMass_genWMass","BHCJetMass_genWMass;BHCJetMass;genWMass",25,0,250,25,0,250);
		//19 - 84 - BHC jet eta vs gen W eta 
		TH2D* BHCJetEtaCenter_genWEtaCenter = new TH2D("BHCJetEtaCenter_genWEtaCenter","BHCJetEtaCenter_genWEtaCenter;BHCJetEtaCenter;genWEtaCenter",25,-3.2,3.2,25,-3.2,3.2);
		//20 - 85 - BHC jet phi vs gen W phi 
		TH2D* BHCJetPhiCenter_genWPhiCenter = new TH2D("BHCJetPhiCenter_genWPhiCenter","BHCJetPhiCenter_genWPhiCenter;BHCJetPhiCenter;genWPhiCenter",25,0.,8*atan(1),25,0.,8*atan(1));
		//21 - 106 - BHC jet subcluster energy vs subcluster mass
		TH2D* BHCJet_subclusterEnergy_subclusterMass = new TH2D("BHCJet_subclusterEnergy_subclusterMass","BHCJet_subclusterEnergy_subclusterMass;subclusterEnergy;subclusterMass",25,0,500,25,0,100);
		//22 - 107 - BHC jet subcluster energy vs # effective rec hits
		TH2D* BHCJet_subclusterEnergy_subclusterEffnRhs = new TH2D("BHCJet_subclusterEnergy_subclusterEffnRhs","BHCJet_subclusterEnergy_subclusterEffnRhs;subclusterEnergy;subclusterEffnRhs",25,0,500,25,0,500);
		//23 - 108 - BHC jet subcluster mass vs # effective rechits 	
		TH2D* BHCJet_subclusterMass_subclusterEffnRhs = new TH2D("BHCJet_subclusterMass_subclusterEffnRhs","BHCJet_subclusterMass_subclusterEffnRhs;subclusterMass;subclusterEffnRhs",25,0,100,25,0,200);
		//24 - 109 - reco AK4 jet subcluster energy vs subcluster mass
		TH2D* recoAK4Jet_subclusterEnergy_subclusterMass = new TH2D("recoAK4Jet_subclusterEnergy_subclusterMass","recoAK4Jet_subclusterEnergy_subclusterMass;subclusterEnergy;subclusterMass",25,0,500,25,0,100);
		//25 - 110 - recoAK4 jet subcluster energy vs # effective rec hits
		TH2D* recoAK4Jet_subclusterEnergy_subclusterEffnRhs = new TH2D("recoAK4Jet_subclusterEnergy_subclusterEffnRhs","recoAK4Jet_subclusterEnergy_subclusterEffnRhs;subclusterEnergy;subclusterEffnRhs",25,0,500,25,0,200);
		//26 - 111 - recoAK4 jet subcluster mass vs # effective rechits 	
		TH2D* recoAK4Jet_subclusterMass_subclusterEffnRhs = new TH2D("recoAK4Jet_subclusterMass_subclusterEffnRhs","recoAK4Jet_subclusterMass_subclusterEffnRhs;subclusterMass;subclusterEffnRhs",25,0,100,25,0,200);
		//27 - 112 - subcluster energy vs dR of subcluster to jet center
		TH2D* BHCJet_subclusterEnergy_subclusterdRToJet = new TH2D("BHCJet_subclusterEnergy_subclusterdRToJet","BHCJet_subclusterEnergy_subclusterdRToJet;subclusterEnergy;subclusterdRToJet",25,0,500,25,0,0.5);
		//28 - 113 - subcluster eff # rhs vs dR of subcluster to jet center
		TH2D* BHCJet_subclusterEffnRhs_subclusterdRToJet = new TH2D("BHCJet_subclusterEffnRhs_subclusterdRToJet","BHCJet_subclusterEffnRhs_subclusterdRToJet;subclusterEffnRhs;subclusterdRToJet",25,0,100,25,0,0.5);
		//29 - 114 - BHC jets gen-matched to Ws - subcluster mass vs # subclusters/jet
		TH2D* BHCJetW_nSubclustersJet_mass = new TH2D("BHCJetW_nSubclustersJet_mass","BHCJetW_nSubclustersJet_mass;nSubclustersJet;mass",10,0,10,50,0,250);
		//30 - 115 - BHC jets gen-matched to Ws - dR bw gen partons of W vs # subclusters/jet
		TH2D* BHCJetW_dRGenPartons_nSubclustersJet = new TH2D("BHCJetW_dRGenPartons_nSubclustersJet","BHCJetW_dRGenPartons_nSubclustersJet;dRGenPartons;nSubclustersJet",50,0,2.,10,0,10);
		//31 - 116 - BHC jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* BHCJetW_dRGenPartons_jetSize = new TH2D("BHCJetW_dRGenPartons_jetSize","BHCJetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//32 - 117 - BHC jets gen-matched to Ws - Eratio of subcl E/gen parton E vs # subclusters/jet
		TH2D* BHCJetW_EratioSubclGenPart_nSubclustersJet = new TH2D("BHCJetW_EratioSubclGenPart_nSubclustersJet","BHCJetW_EratioSubclGenPart_nSubclustersJet;EratioSubclGenPart;nSubclustersJet",50,0,2.,10,0,10);
		//33 - 118 - BHC jets gen-matched to Tops - subcluster mass vs # subclusters/jet
		TH2D* BHCJetTop_nSubclustersJet_mass = new TH2D("BHCJetTop_nSubclustersJet_mass","BHCJetTop_nSubclustersJet_mass;nSubclustersJet;subclMass",10,0,10,50,0,250);
		//34 - 119 - # bhc jets as a function of opening angle of partons from W
		TH2D* BHCJetW_openAng_nJets = new TH2D("BHCJetW_openAng_nJets","BHCJetW_openAng_nJets;openAng;nJets",25,0,3.3,30,0,30);
		//35 - 120 - # bhc subclusters/jet as a function of opening angle of partons from W
		TH2D* BHCJetW_openAng_nSubclustersJet = new TH2D("BHCJetW_openAng_nSubclustersJet","BHCJetW_openAng_nSubclustersJet;openAng;nSubclustersJet",25,0,3.3,10,0,10);
		//36 - 121 - bhc subcluster mass as a function of opening angle of partons from W
		TH2D* BHCJetW_openAng_subclMass = new TH2D("BHCJetW_openAng_subclMass","BHCJetW_openAng_subclMass;openAng;subclMass",25,0,3.3,50,0,200);
		//37 - 122 - BHC jets gen-matched to Ws with exactly 1 subcluster - dR bw gen partons of W vs jet size
		TH2D* BHCJetW_1subcl_dRGenPartons_jetSize = new TH2D("BHCJetW_1subcl_dRGenPartons_jetSize","BHCJetW_1subcl_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//38 - 123 - BHC jets gen-matched to Ws with 2+ subclusters - dR bw gen partons of W vs jet size
		TH2D* BHCJetW_ge2subcl_dRGenPartons_jetSize = new TH2D("BHCJetW_ge2subcl_dRGenPartons_jetSize","BHCJetW_ge2subcl_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//39 - 124 - BHC jets gen-matched to Ws with exactly 1 subcluster - dR bw gen partons of W vs avg. parton energy
		TH2D* BHCJetW_1subcl_dRGenPartons_avgPartonEnergy = new TH2D("BHCJetW_1subcl_dRGenPartons_avgPartonEnergy","BHCJetW_1subcl_dRGenPartons_avgPartonEnergy;dRGenPartons;avgPartonEnergy",50,0,2.,50,0,1000.);
		//40 - 125 - BHC jets gen-matched to Ws with 2+ subclusters - dR bw gen partons of W vs avg. parton energy
		TH2D* BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy = new TH2D("BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy","BHCJetW_ge2subcl_dRGenPartons_avgPartonEnergy;dRGenPartons;avgPartonEnergy",50,0,2.,50,0,1000.);
		//41 - 126 - BHC jets gen-matched to Ws with 2+ subclusters - subcluster energy vs lead index of subcluster (ie 0 = lead, 1 = sublead, etc)
		TH2D* BHCJetW_ge2subcl_subclEnergy_subclLeadIdx = new TH2D("BHCJetW_ge2subcl_subclEnergy_subclLeadIdx","BHCJetW_ge2subcl_subclEnergy_subclLeadIdx;subclEnergy;subclLeadIdx",5,0,5,25,0,500);
		//42 - 127 - BHC jets gen-matched to Ws - Eratio of jet E/gen W E vs # subclusters/jet
		TH2D* BHCJetW_EratioJetGenW_nSubclustersJet = new TH2D("BHCJetW_ge2subcl_EratioJetGenW_nSubclustersJet","BHCJetW_EratioJetGenW_nSubclustersJet;EratioJetGenW;nSubclustersJet",50,0,2.,10,0,10);
		//43 - 128 - BHC jets - jet mass vs jet size
		TH2D* BHCJet_jetMass_jetSize = new TH2D("BHCJet_jetMass_jetSize","BHCJet_jetMass_jetSize;jetMass;jetSize",50,0,250.,50,0,2.);
		//44 - 129 - eta-phi event display of rechits for specified _evt2disp with cell energy on the z axis (overall event)
		TH2D* EvtDisplay_etaCell_phiCell = new TH2D("EvtDisplay_etaCell_phiCell","EvtDisplay_etaCell_phiCell;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,0,8*atan(1));
		//45 - 130 - reco AK8 jet mass vs reco jet pt
		TH2D* recoAK8JetMass_recoAK8JetSize = new TH2D("recoAK8Jet_jetMass_jetSize","recoAK8Jet_jetMass_jetSize;recoAK8JetMass;recoAK8JetSize",50,0,250,50,0,2.);
		//46 - 131 - reco AK15 jet mass vs reco jet pt
		TH2D* recoAK15JetMass_recoAK15JetSize = new TH2D("recoAK15Jet_jetMass_jetSize","recoAK15Jet_jetMass_jetSize;recoAK15JetMass;recoAK15JetSize",50,0,250,50,0,2.);
		//47 - 132 - matched reco AK8 jet mass vs BHC jet mass
		TH2D* recoAK8JetMass_BHCJetMass_matched = new TH2D("recoAK8JetMass_BHCJetMass_BHCtoAK8matched","recoAK8JetMass_BHCJetMass_BHCtoAK8matched;recoAK8 Jet Mass;BHC Jet Mass",25,0,250,25,0,250);
		//48 - 133 - recoAK8 jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* recoAK8JetW_dRGenPartons_jetSize = new TH2D("recoAK8JetW_dRGenPartons_jetSize","recoAK8JetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//49 - 136 - recoAK15 jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* recoAK15JetW_dRGenPartons_jetSize = new TH2D("recoAK15JetW_dRGenPartons_jetSize","recoAK15JetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//50 - 137 - recoAK4 jets gen-matched to Ws - dR bw gen partons of W vs jet size
		TH2D* recoAK4JetW_dRGenPartons_jetSize = new TH2D("recoAK4JetW_dRGenPartons_jetSize","recoAK4JetW_dRGenPartons_jetSize;dRGenPartons;jetSize",50,0,2.,50,0,2.);
		//51 - 55 - high mass + W-matched BHC jets - subclSize of subclsuters gen-matched to W partons / size jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt = new TH2D("BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt","BHCJetW_highMass_partonMatchRelSubclSize_RelSubclPt;RelSubclSize;RelSubclPt",25,0,5,25,0,1.2);
		//52 - 56 - high mass + W-matched BHC jets - subclSize of subclusters NOT gen-matched to W partons / size jet vs pt of subclusters NOT gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt","BHCJetW_highMass_partonNoMatchRelSubclSize_RelSubclPt;RelSubclSize;RelSubclPt",25,0,5,25,0,1.2);
		//BREAK!!!!! -6 
		//53 - 59 - 146 - high mass + W-matched BHC jets - subcl time sig of subclusters gen-matched to W partons / time sig jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclE = new TH2D("BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclE","BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclE;RelSubclTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//54 - 60 - 147 - high mass + W-matched BHC jets - subcl time sig of subclusters NOT gen-matched to W partons / time sig jet vs E of subclusters NOT gen-matched to W partons / E jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclE = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclE","BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclE;RelSubclTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//55 - 61 - 148 - high mass + W-matched BHC jets - subcl time sig of subclusters gen-matched to W partons / time sig jet vs size of subclusters gen-matched to W partons / size jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclSize = new TH2D("BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclSize","BHCJetW_highMass_partonMatchRelSubclTimeSig_RelSubclSize;RelSubclTimeSig;RelSubclSize",25,0,5,25,0,5);
		//56 - 62 - 149 - high mass + W-matched BHC jets - subcl time sig of subclusters NOT gen-matched to W partons / time sig jet vs size of subclusters NOT gen-matched to W partons / size jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclSize = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclSize","BHCJetW_highMass_partonNoMatchRelSubclTimeSig_RelSubclSize;RelSubclTimeSig;RelSubclSize",25,0,5,25,0,5);
		//BREAK!!!!! -10 
		//57 - 67 - highMass + q-matched BHC jets - subcl time sig of subclusters gen-matched to qs / time sig jet vs E of subclusters gen-matched to qs / E jet
		TH2D* BHCJetq_highMass_partonMatchRelSubclTimeSig_RelSubclE = new TH2D("BHCJetq_highMass_partonMatchRelSubclTimeSig_RelSubclE","BHCJetq_highMass_partonMatchRelSubclTimeSig_RelSubclE;RelSubclTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//58 - 68 - highMass + q-matched BHC jets - subcl time sig of subclusters NOT gen-matched to qs / time sig jet vs E of subclusters NOT gen-matched to qs / E jet
		TH2D* BHCJetq_highMass_partonNoMatchRelSubclTimeSig_RelSubclE = new TH2D("BHCJetq_highMass_partonNoMatchRelSubclTimeSig_RelSubclE","BHCJetq_highMass_partonNoMatchRelSubclTimeSig_RelSubclE;RelSubclTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//59 - 69 - highMass + top-matched BHC jets - subcl time sig of subclusters gen-matched to tops / time sig jet vs E of subclusters gen-matched to tops / E jet
		TH2D* BHCJetTop_highMass_partonMatchRelSubclTimeSig_RelSubclE = new TH2D("BHCJetTop_highMass_partonMatchRelSubclTimeSig_RelSubclE","BHCJetTop_highMass_partonMatchRelSubclTimeSig_RelSubclE;RelSubclTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//60 - 70 - highMass + top-matched BHC jets - subcl time sig of subclusters NOT gen-matched to tops / time sig jet vs E of subclusters NOT gen-matched to tops / E jet
		TH2D* BHCJetTop_highMass_partonNoMatchRelSubclTimeSig_RelSubclE = new TH2D("BHCJetTop_highMass_partonNoMatchRelSubclTimeSig_RelSubclE","BHCJetTop_highMass_partonNoMatchRelSubclTimeSig_RelSubclE;RelSubclTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//BREAK!!!!! -22 
		//61 - 83 - high mass + W-matched BHC jets - rel sucbl eta sig vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclEtaSig_RelSubclE = new TH2D("BHCJetW_highMass_partonMatchRelSubclEtaSig_RelSubclE","BHCJetW_highMass_partonMatchRelSubclEtaSig_RelSubclE;RelSubclEtaSig;RelSubclE",25,0,5,25,0,1.2);
		//62 - 84 - high mass + W-matched BHC jets - rel sucbl eta sig vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclEtaSig_RelSubclE = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclEtaSig_RelSubclE","BHCJetW_highMass_partonNoMatchRelSubclEtaSig_RelSubclE;RelSubclEtaSig;RelSubclE",25,0,5,25,0,1.2);
		//63 - 85 - high mass + W-matched BHC jets - rel sucbl phi sig vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclPhiSig_RelSubclE = new TH2D("BHCJetW_highMass_partonMatchRelSubclPhiSig_RelSubclE","BHCJetW_highMass_partonMatchRelSubclPhiSig_RelSubclE;RelSubclPhiSig;RelSubclE",25,0,5,25,0,1.2);
		//64 - 86 - high mass + W-matched BHC jets - rel sucbl phi sig vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclPhiSig_RelSubclE = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclPhiSig_RelSubclE","BHCJetW_highMass_partonNoMatchRelSubclPhiSig_RelSubclE;RelSubclPhiSig;RelSubclE",25,0,5,25,0,1.2);
		//65 - 87 - high mass + W-matched BHC jets - subcl time var of subclusters gen-matched to W partons / time var jet vs pt of subclusters gen-matched to W partons / pt jet
		TH2D* BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclE = new TH2D("BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclE","BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclE;RelSubclTimeVar;RelSubclE",25,0,5,25,0,1.2);
		//66 - 88 - high mass + W-matched BHC jets - subcl time var of subclusters NOT gen-matched to W partons / time var jet vs E of subclusters NOT gen-matched to W partons / E jet
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclE = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclE","BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclE;RelSubclTimeVar;RelSubclE",25,0,5,25,0,1.2);
		//67 - 89 - high mass + W-matched BHC jets - rel sucbl geo avg of spatial size and time var vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE = new TH2D("BHCJetW_highMass_partonMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE","BHCJetW_highMass_partonMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE;RelSubclGeoAvgSpatialSizeTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//68 - 90 - high mass + W-matched BHC jets - rel sucbl geo avg of spatial size and time var vs rel subcl energy for NOT matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE","BHCJetW_highMass_partonNoMatchRelSubclGeoAvgSpatialSizeTimeSig_RelSubclE;RelSubclGeoAvgSpatialSizeTimeSig;RelSubclE",25,0,5,25,0,1.2);
		//69 - 91 - high mass + W-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE = new TH2D("BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE","BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE;RelSubclGeoAvgAllVar;RelSubclE",50,0,1.2,50,0,1.2);
		//70 - 92 - high mass + W-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for NOT matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE","BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE;RelSubclGeoAvgAllVar;RelSubclE",50,0,1.2,50,0,1.2);
		//BREAK!!!!! -24 
		//71 - 95 - high mass + Top-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for matched subcls
		TH2D* BHCJetTop_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE = new TH2D("BHCJetTop_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE","BHCJetTop_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE;RelSubclGeoAvgAllVar;RelSubclE",50,0,1.2,50,0,1.2);
		//72 - 96 - high mass + Top-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for NOT matched subcls
		TH2D* BHCJetTop_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE = new TH2D("BHCJetTop_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE","BHCJetTop_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE;RelSubclGeoAvgAllVar;RelSubclE",50,0,1.2,50,0,1.2);
		//73 - 97 - high mass + q-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for matched subcls
		TH2D* BHCJetq_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE = new TH2D("BHCJetq_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE","BHCJetq_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE;RelSubclGeoAvgAllVar;RelSubclE",50,0,1.2,50,0,1.2);
		//74 - 98 - high mass + q-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for NOT matched subcls
		TH2D* BHCJetq_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE = new TH2D("BHCJetq_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE","BHCJetq_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE;RelSubclGeoAvgAllVar;RelSubclE",50,0,1.2,50,0,1.2);
		//75 - 99 - high mass + W-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE_zoom = new TH2D("BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE_zoom","BHCJetW_highMass_partonMatchRelSubclGeoAvgAllVar_RelSubclE_zoom;RelSubclGeoAvgAllVar;RelSubclE_zoom",50,0,0.5,50,0,0.5);
		//76 - 100 - high mass + W-matched BHC jets - sucbl geo avg of rel time var, rel eta var, rel phi var vs rel subcl energy for NOT matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE_zoom = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE_zoom","BHCJetW_highMass_partonNoMatchRelSubclGeoAvgAllVar_RelSubclE_zoom;RelSubclGeoAvgAllVar;RelSubclE_zoom",50,0,0.5,50,0,0.5);
		//77 - 101 - high mass + W-matched BHC jets - rel subcl eta var vs rel subcl phi for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclPhiVar = new TH2D("BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclPhiVar","BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclPhiVar;RelSubclEtaVar;RelSubclPhiVar",50,0,1.2,50,0,1.2);
		//78 - 102 - high mass + W-matched BHC jets - rel subcl eta var vs rel subcl phi for NOT matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclPhiVar = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclPhiVar","BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclPhiVar;RelSubclEtaVar;RelSubclPhiVar",50,0,1.2,50,0,1.2);
		//79 - 103 - high mass + W-matched BHC jets - rel subcl eta var vs rel subcl phi for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclTimeVar = new TH2D("BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclTimeVar","BHCJetW_highMass_partonMatchRelSubclEtaVar_RelSubclTimeVar;RelSubclEtaVar;RelSubclTimeVar",50,0,1.2,50,0,1.2);
		//80 - 104 - high mass + W-matched BHC jets - rel subcl eta var vs rel subcl phi for NOT matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclTimeVar = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclTimeVar","BHCJetW_highMass_partonNoMatchRelSubclEtaVar_RelSubclTimeVar;RelSubclEtaVar;RelSubclTimeVar",50,0,1.2,50,0,1.2);
		//81 - 105 - high mass + W-matched BHC jets - rel subcl eta var vs rel subcl phi for matched subcls
		TH2D* BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclPhiVar = new TH2D("BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclPhiVar","BHCJetW_highMass_partonMatchRelSubclTimeVar_RelSubclPhiVar;RelSubclTimeVar;RelSubclPhiVar",50,0,1.2,50,0,1.2);
		//82 - 106 - high mass + W-matched BHC jets - rel subcl eta var vs rel subcl phi for NOT matched subcls
		TH2D* BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclPhiVar = new TH2D("BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclPhiVar","BHCJetW_highMass_partonNoMatchRelSubclTimeVar_RelSubclPhiVar;RelSubclTimeVar;RelSubclPhiVar",50,0,1.2,50,0,1.2);



		//event display histograms (per gen object)
		vector<TH2D*> _evtdisps_obj;
		//these should be centered on their respective gen particles
		//0 - single W, W+gluon, first W in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_W = new TH2D("EvtDisplay_etaCell_phiCell_W","EvtDisplay_etaCell_phiCell_W;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//1 - second W in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_W2 = new TH2D("EvtDisplay_etaCell_phiCell_W2","EvtDisplay_etaCell_phiCell_W2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//2 - gluon in W+gluon
		TH2D* EvtDisplay_etaCell_phiCell_gluon = new TH2D("EvtDisplay_etaCell_phiCell_gluon","EvtDisplay_etaCell_phiCell_gluon;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//3 - q1 in QCD dijets
		TH2D* EvtDisplay_etaCell_phiCell_q1 = new TH2D("EvtDisplay_etaCell_phiCell_q1","EvtDisplay_etaCell_phiCell_q1;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//4 - q2 in QCD dijets
		TH2D* EvtDisplay_etaCell_phiCell_q2 = new TH2D("EvtDisplay_etaCell_phiCell_q2","EvtDisplay_etaCell_phiCell_q2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//5 - b1 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_b1 = new TH2D("EvtDisplay_etaCell_phiCell_b1","EvtDisplay_etaCell_phiCell_b1;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//6 - b2 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_b2 = new TH2D("EvtDisplay_etaCell_phiCell_b2","EvtDisplay_etaCell_phiCell_b2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//7 - top1 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_top1 = new TH2D("EvtDisplay_etaCell_phiCell_top1","EvtDisplay_etaCell_phiCell_top1;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//8 - top2 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_top2 = new TH2D("EvtDisplay_etaCell_phiCell_top2","EvtDisplay_etaCell_phiCell_top2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));


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
	void SetReducePU(bool t){ _zoom_window = t; if(_zoom_window) cout << "Reducing PU with zoom window." << endl; }//draws rectangle around hard scattering to reduce # of rechits to cluster

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
		bool _data, _debug, _zoom_window;
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

		double CalcSubleadSize(Matrix& cov, bool time = false){
			if(cov.GetDims()[0] != 3 || cov.GetDims()[1] != 3){
				cout << "Error: can't calculate size for matrix of size " << cov.GetDims()[0] << " x " << cov.GetDims()[1] << endl;
				return -1;
			}
			vector<double> eigvals;
			vector<Matrix> eigvecs;
			if(time){
				cov.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[1]);
			}
			else{
				Matrix cov2D(2,2);
				Get2DMat(cov,cov2D);	
				cov2D.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[0]);
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




		//PU cleaning stuff
		double _minRelPt, _maxRelTimeVar;
		//function that takes in jet and returns new, PU-cleaned jet with bool option for subcl downweight or full removal

		




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
