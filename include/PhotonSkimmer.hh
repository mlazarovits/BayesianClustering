#ifndef PHOTONSKIMMER_HH
#define PHOTONSKIMMER_HH

#include "JetPoint.hh"
#include <TFile.h>
#include "BaseSkimmer.hh"
#include "PhotonProducer.hh"
#include "JetProducer.hh"
#include "TSystem.h"
#include <math.h>
#include <fstream>

using procCat = BaseSkimmer::procCat;
class PhotonSkimmer : public BaseSkimmer{
	public:
		PhotonSkimmer(){
			_evti = 0;
			_evtj = 0;
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_gev = 1/30.;
			_applyFrac = false;
			_jetprod = nullptr;

			//reqs on iso bkg sample
			_isoBkgSel = false;
			_minPhoPt_isoBkg = 70;
			_minHt_isoBkg = 50;
			_minJetPt_isoBkg = 50;
			_maxMet_isoBkg = 150;

			_weight = 1;
			_cell = 0;
			_tresCte = 0;
			_tresNoise = 0;
			_tresStoch = 0;

		};
		virtual ~PhotonSkimmer(){ };

		//get rechits from file to cluster
		PhotonSkimmer(TFile* file) : BaseSkimmer(file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time
			_prod = new PhotonProducer(file);

			//set producer to get jets with different kin reqs - can't use same file pointer ig?
			string fname = file->GetName();
			TFile* f2 = TFile::Open(fname.c_str());
			_jetprod = new JetProducer(f2);
			

			//set histogram weights for HT slices, etc
			_weight = 1;
			if(_data || fname.find("QCD") != string::npos){ _weight = 1.; }
			else{
				cout << "Getting weights from info/EventWeights_AL1IsoPho.txt for GJets" << endl;
			        ifstream weights("info/EventWeights_AL1IsoPho.txt", std::ios::in);
			        string filein;
			        string filename = file->GetName();
			        double jet_weight, pho_weight;
			        while( weights >> filein >> jet_weight >> pho_weight){
			                if(filename.find(filein) == string::npos) continue;
			                else{
			                        _weight = pho_weight;
			                        break;
			                }
			        }
			}		


			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/photon_skims_"+_cms_label+".root";
			
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_gev = 1/30.;
			_applyFrac = false;
			
			//reqs on iso bkg sample
			_isoBkgSel = false;
			_minPhoPt_isoBkg = 70;
			_minHt_isoBkg = 50;
			_minJetPt_isoBkg = 50;
			_maxMet_isoBkg = 150;
			
			objE->SetTitle("totphoE");
			objE->SetName("totphoE");
			
			objE_clusterE->SetTitle("phoE_clusterE");
			objE_clusterE->SetName("phoE_clusterE");

			SetupDetIDsEB( _detIDmap, _ietaiphiID );
			InitHists();
		}
		
		//get rechits from file to cluster
		PhotonSkimmer(string filelist) : BaseSkimmer(filelist){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time
			//this is just the type of producer, there is a GetSuperCluster fcn in the base producer class
                        TChain* ch = MakeTChain(filelist);
                        if(ch == nullptr) return;
			_prod = new PhotonProducer(ch);

			//set producer to get jets with different kin reqs - can't use same file pointer ig?
                        TChain* ch2 = MakeTChain(filelist);
			_jetprod = new JetProducer(ch2);
			

			//set histogram weights for HT slices, etc
			_weight = 1;
			if(_data || filelist.find("QCD") != string::npos){ _weight = 1.; }
			else{
				cout << "Getting weights from info/EventWeights_AL1IsoPho.txt for GJets" << endl;
			        ifstream weights("info/EventWeights_AL1IsoPho.txt", std::ios::in);
			        string filein;
			        double jet_weight, pho_weight;
			        while( weights >> filein >> jet_weight >> pho_weight){
			                if(filelist.find(filein) == string::npos) continue;
			                else{
			                        _weight = pho_weight;
			                        break;
			                }
			        }
			}		


			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/photon_skims_"+_cms_label+".root";
			
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_gev = 1/30.;
			_applyFrac = false;
			
			//reqs on iso bkg sample
			_isoBkgSel = false;
			_minPhoPt_isoBkg = 70;
			_minHt_isoBkg = 50;
			_minJetPt_isoBkg = 50;
			_maxMet_isoBkg = 150;
			
			objE->SetTitle("totphoE");
			objE->SetName("totphoE");
			
			objE_clusterE->SetTitle("phoE_clusterE");
			objE_clusterE->SetName("phoE_clusterE");

			SetupDetIDsEB( _detIDmap, _ietaiphiID );
			InitHists();
		}
		void InitHists(){
			//add photon specific histograms
                        _hists1D.push_back(slope_space);
                        _hists1D.push_back(slope_etaT);
                        _hists1D.push_back(slope_phiT);
                        _hists1D.push_back(polar_ang);
                        _hists1D.push_back(azimuth_ang);
                        _hists1D.push_back(e_subcl);
                        _hists1D.push_back(rotundity_3D);
                        _hists1D.push_back(rotundity_2D);
                        _hists1D.push_back(velocity);
                        _hists1D.push_back(eigen2D_ratio);
                        _hists1D.push_back(etaSig);
                        _hists1D.push_back(phiSig);
                        _hists1D.push_back(timeSig);
                        _hists1D.push_back(fracE);
                        _hists1D.push_back(phiEll);
                        _hists1D.push_back(etaSig_pos);
                        _hists1D.push_back(etaSig_neg);
                        _hists1D.push_back(etaphi_cov);
                        _hists1D.push_back(timeeta_cov);
                        _hists1D.push_back(timephi_cov);
			_hists1D.push_back(timemaj_cov);
			_hists1D.push_back(timemin_cov);
                	_hists1D.push_back(timeSig_gamPtLow);
                	_hists1D.push_back(timeSig_gamPtMed);
                	_hists1D.push_back(timeSig_gamPtHi);
			_hists1D.push_back(timeeta_covUnnorm);
			_hists1D.push_back(timephi_covUnnorm);
			_hists1D.push_back(timeMaj_covUnnorm);
			_hists1D.push_back(timeMin_covUnnorm);
			_hists1D.push_back(logE_timeeta_cov);
			_hists1D.push_back(logE_timephi_cov);
			_hists1D.push_back(logE_etaphi_cov);
                	_hists1D.push_back(logE_etaSig);
                	_hists1D.push_back(logE_phiSig);
                	_hists1D.push_back(logE_timeSig);
			_hists1D.push_back(noE_timeeta_cov);
			_hists1D.push_back(noE_timephi_cov);
			_hists1D.push_back(noE_etaphi_cov);
                	_hists1D.push_back(noE_etaSig);
                	_hists1D.push_back(noE_phiSig);
                	_hists1D.push_back(noE_timeSig);
			_hists1D.push_back(noE_timeeta_covUnnorm);
			_hists1D.push_back(noE_timephi_covUnnorm);
			_hists1D.push_back(noE_etaphi_covUnnorm);
			_hists1D.push_back(logE_timeeta_covUnnorm);
			_hists1D.push_back(logE_timephi_covUnnorm);
			_hists1D.push_back(logE_etaphi_covUnnorm);
			_hists1D.push_back(etaphi_covUnnorm);
			_hists1D.push_back(noE_smaj);		
			_hists1D.push_back(noE_smin);		
			_hists1D.push_back(logE_smaj);		
			_hists1D.push_back(logE_smin);		
			_hists1D.push_back(noE_time_center);
			_hists1D.push_back(noE_eta_center);
			_hists1D.push_back(noE_phi_center);
			_hists1D.push_back(noE_phiEll);		
			_hists1D.push_back(noE_timeSmaj_cov);
			_hists1D.push_back(noE_timeSmin_cov);
			_hists1D.push_back(noE_timeSmaj_covUnnorm);
			_hists1D.push_back(noE_timeSmin_covUnnorm);
			_hists1D.push_back(noE_rotundity_2D);
			_hists1D.push_back(logE_time_center);
			_hists1D.push_back(logE_eta_center);
			_hists1D.push_back(logE_phi_center);
			_hists1D.push_back(logE_phiEll);	
			_hists1D.push_back(logE_timeSmaj_cov);
			_hists1D.push_back(logE_timeSmin_cov);
			_hists1D.push_back(logE_timeSmaj_covUnnorm);
			_hists1D.push_back(logE_timeSmin_covUnnorm);
			_hists1D.push_back(logE_rotundity_2D);
			_hists1D.push_back(rotundity_2D_Elo);
			_hists1D.push_back(rotundity_2D_Emed);
			_hists1D.push_back(rotundity_2D_Ehi);
			_hists1D.push_back(rotundity_2D_1Q);
			_hists1D.push_back(rotundity_2D_2Q);
			_hists1D.push_back(rotundity_2D_3Q);
			_hists1D.push_back(rotundity_2D_4Q);
			_hists1D.push_back(phiE_2D_1Q);
			_hists1D.push_back(phiE_2D_2Q);
			_hists1D.push_back(phiE_2D_3Q);
			_hists1D.push_back(phiE_2D_4Q);
			_hists1D.push_back(etaPhiCov_1Q);
			_hists1D.push_back(etaPhiCov_2Q);
			_hists1D.push_back(etaPhiCov_3Q);
			_hists1D.push_back(etaPhiCov_4Q);
			_hists1D.push_back(timeEtaCov_1Q);
			_hists1D.push_back(timeEtaCov_2Q);
			_hists1D.push_back(timeEtaCov_3Q);
			_hists1D.push_back(timeEtaCov_4Q);
                	_hists1D.push_back(etaCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(etaCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(phiCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(phiCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(timeCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(timeCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(phoE_phiE2Deq0PiOv2);
                	_hists1D.push_back(phoE_phiE2Dneq0PiOv2);
                	_hists1D.push_back(nSubClusters_phiE2Deq0PiOv2);
                	_hists1D.push_back(nSubClusters_phiE2Dneq0PiOv2);
                	_hists1D.push_back(etaSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(etaSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(phiSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(phiSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(timeSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(timeSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(etaPhiCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(etaPhiCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(timeEtaCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(timeEtaCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(timePhiCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(timePhiCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(timeMajCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(timeMajCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(timeMinCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(timeMinCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(phiE2D_phiE2Deq0PiOv2);
                	_hists1D.push_back(phiE2D_phiE2Dneq0PiOv2);
                	_hists1D.push_back(rot2D_phiE2Deq0PiOv2);
                	_hists1D.push_back(rot2D_phiE2Dneq0PiOv2);
			_hists1D.push_back(logErotundity_2D_Elo);
			_hists1D.push_back(logErotundity_2D_Emed);
			_hists1D.push_back(logErotundity_2D_Ehi);
			_hists1D.push_back(logErotundity_2D_1Q);
			_hists1D.push_back(logErotundity_2D_2Q);
			_hists1D.push_back(logErotundity_2D_3Q);
			_hists1D.push_back(logErotundity_2D_4Q);
			_hists1D.push_back(logEphiE_2D_1Q);
			_hists1D.push_back(logEphiE_2D_2Q);
			_hists1D.push_back(logEphiE_2D_3Q);
			_hists1D.push_back(logEphiE_2D_4Q);
			_hists1D.push_back(logEetaPhiCov_1Q);
			_hists1D.push_back(logEetaPhiCov_2Q);
			_hists1D.push_back(logEetaPhiCov_3Q);
			_hists1D.push_back(logEetaPhiCov_4Q);
			_hists1D.push_back(logEtimeEtaCov_1Q);
			_hists1D.push_back(logEtimeEtaCov_2Q);
			_hists1D.push_back(logEtimeEtaCov_3Q);
			_hists1D.push_back(logEtimeEtaCov_4Q);
                	_hists1D.push_back(logEetaCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEetaCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEphiCenter_phiE2Deq0PiOv2); 
                	_hists1D.push_back(logEphiCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEtimeCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEtimeCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEphoE_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEphoE_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEnSubClusters_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEnSubClusters_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEetaSig_phiE2Deq0PiOv2); 
                	_hists1D.push_back(logEetaSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEphiSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEphiSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEtimeSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEtimeSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEetaPhiCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEetaPhiCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEtimeEtaCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEtimeEtaCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEtimePhiCov_phiE2Deq0PiOv2); 
                	_hists1D.push_back(logEtimePhiCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEtimeMajCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEtimeMajCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEtimeMinCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEtimeMinCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logEphiE2D_phiE2Deq0PiOv2);
                	_hists1D.push_back(logEphiE2D_phiE2Dneq0PiOv2);
                	_hists1D.push_back(logErot2D_phiE2Deq0PiOv2);
                	_hists1D.push_back(logErot2D_phiE2Dneq0PiOv2);
			_hists1D.push_back(noErotundity_2D_Elo);
			_hists1D.push_back(noErotundity_2D_Emed);
			_hists1D.push_back(noErotundity_2D_Ehi);
			_hists1D.push_back(noErotundity_2D_1Q);
			_hists1D.push_back(noErotundity_2D_2Q);
			_hists1D.push_back(noErotundity_2D_3Q);
			_hists1D.push_back(noErotundity_2D_4Q);
			_hists1D.push_back(noEphiE_2D_1Q);
			_hists1D.push_back(noEphiE_2D_2Q);
			_hists1D.push_back(noEphiE_2D_3Q);
			_hists1D.push_back(noEphiE_2D_4Q);
			_hists1D.push_back(noEetaPhiCov_1Q);
			_hists1D.push_back(noEetaPhiCov_2Q);
			_hists1D.push_back(noEetaPhiCov_3Q);
			_hists1D.push_back(noEetaPhiCov_4Q);
			_hists1D.push_back(noEtimeEtaCov_1Q);
			_hists1D.push_back(noEtimeEtaCov_2Q);
			_hists1D.push_back(noEtimeEtaCov_3Q);
			_hists1D.push_back(noEtimeEtaCov_4Q);
                	_hists1D.push_back(noEetaCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEetaCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEphiCenter_phiE2Deq0PiOv2); 
                	_hists1D.push_back(noEphiCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEtimeCenter_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEtimeCenter_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEphoE_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEphoE_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEnSubClusters_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEnSubClusters_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEetaSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEetaSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEphiSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEphiSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEtimeSig_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEtimeSig_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEetaPhiCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEetaPhiCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEtimeEtaCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEtimeEtaCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEtimePhiCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEtimePhiCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEtimeMajCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEtimeMajCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEtimeMinCov_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEtimeMinCov_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noEphiE2D_phiE2Deq0PiOv2);
                	_hists1D.push_back(noEphiE2D_phiE2Dneq0PiOv2);
                	_hists1D.push_back(noErot2D_phiE2Deq0PiOv2);
                	_hists1D.push_back(noErot2D_phiE2Dneq0PiOv2);
			_hists1D.push_back(swCross); 
			_hists1D.push_back(trueSmaj);
			_hists1D.push_back(trueSmin);
			_hists1D.push_back(trueSiEtaiEta);
			_hists1D.push_back(trueSiPhiiPhi);
			_hists1D.push_back(phoNrhs);		
			_hists1D.push_back(phoNrhs_Ebin1);		
			_hists1D.push_back(phoNrhs_Ebin2);		
			_hists1D.push_back(phoNrhs_Ebin3);		
			_hists1D.push_back(phoNrhs_Ebin4);		
			_hists1D.push_back(swCrossPrime);
			_hists1D.push_back(etaDiff);
			_hists1D.push_back(dPhi_phiCenterMet);
			_hists1D.push_back(dPhi_phiCenterMet_etaSigge0p3ANDphiSigle0p3);
			_hists1D.push_back(dR_trackSubcl);	
			_hists1D.push_back(dE_trackSubcl);	
			_hists1D.push_back(E_IsoBkgSel);
                	_hists1D.push_back(etaSig_IsoBkgSel);
                	_hists1D.push_back(phiSig_IsoBkgSel);
			_hists1D.push_back(etaPhiCov_IsoBkgSel);
			_hists1D.push_back(timeEtaCov_IsoBkgSel);
			_hists1D.push_back(timePhiCov_IsoBkgSel);
                	_hists1D.push_back(majLength_IsoBkgSel);
                	_hists1D.push_back(minLength_IsoBkgSel);
			_hists1D.push_back(phi2D_IsoBkgSel);		
			_hists1D.push_back(rot2D_IsoBkgSel);
			_hists1D.push_back(rot3D_IsoBkgSel);
			_hists1D.push_back(nRhs_rnkThresh);
			_hists1D.push_back(etaAngle3D);
			_hists1D.push_back(phiAngle3D);
			_hists1D.push_back(etaAngle2D);
			_hists1D.push_back(phiAngle2D);
			_hists1D.push_back(majLength3D);
			_hists1D.push_back(majLength2D);
			_hists1D.push_back(timesSigSq_measErr);
			_hists1D.push_back(rhTime);
			_hists1D.push_back(photonPt);
			_hists1D.push_back(etaAngle2D_absTimeMajCovge0p1);
			_hists1D.push_back(sinDiffEtaAngle2D_3D);
			_hists1D.push_back(timePhiCovOvtimeMaj2DCov);
			_hists1D.push_back(timeEtaCovOvtimeMaj2DCov);
			_hists1D.push_back(timeMaj2DCovOvtimeMaj3DCov);	
			
			_hists2D.push_back(time_E);
                        _hists2D.push_back(az_E);
                        _hists2D.push_back(rot2D_E);
                        _hists2D.push_back(eta_phi);
                        _hists2D.push_back(t_eta);
                        _hists2D.push_back(t_phi);
                        _hists2D.push_back(t_mixcoeff);
                        _hists2D.push_back(E_mixcoeff);
                        _hists2D.push_back(rot3D_E);
                        _hists2D.push_back(npts_E);
                        _hists2D.push_back(nsubcl_nrhs);
                        _hists2D.push_back(nsubcl_mmcoeff);
                        _hists2D.push_back(etaSig_phiSig);
                        _hists2D.push_back(timeSig_etaSig);
                        _hists2D.push_back(timeSig_phiSig);
                        _hists2D.push_back(fracE_mmcoeff);
                        _hists2D.push_back(nsubcl_fracE);
                        _hists2D.push_back(timeCenter_timeSig);
                        _hists2D.push_back(rot2D_az2D);
                        _hists2D.push_back(timeSig_fracE);
                        _hists2D.push_back(timeSig_totE);
                        _hists2D.push_back(timeEtaCov_totE);
                        _hists2D.push_back(timePhiCov_totE);
                	_hists2D.push_back(timeSig_timeMajCov);
                	_hists2D.push_back(timeSig_timeMinCov);
                	_hists2D.push_back(phiEll2D_timeMajCov);
                	_hists2D.push_back(phiEll2D_timeMinCov);
                	_hists2D.push_back(rot2D_timeMajCov);
                	_hists2D.push_back(rot2D_timeMinCov);
                	_hists2D.push_back(rot3D_timeMajCov);
                	_hists2D.push_back(rot3D_timeMinCov);
			_hists2D.push_back(timeCenter_timeMajCov);
			_hists2D.push_back(timeCenter_timeMinCov);
			_hists2D.push_back(etaPhiCov_phiEll2D);
			_hists2D.push_back(timeMajCov_phiEll2D);
			_hists2D.push_back(timeMinCov_phiEll2D);
                	_hists2D.push_back(rot2D_timeMajCovUnnorm);
                	_hists2D.push_back(rot2D_timeMinCovUnnorm);
                	_hists2D.push_back(timeSig_timeEtaCov);
			_hists2D.push_back(timeCenter_etaPhiCov);
			_hists2D.push_back(timeCenter_timeEtaCov);
			_hists2D.push_back(timeCenter_timePhiCov);
			_hists2D.push_back(cmsSmaj_cmsTimeSig);
			_hists2D.push_back(cmsSmin_cmsTimeSig);
			_hists2D.push_back(logESmaj_logETimeSig);
			_hists2D.push_back(logESmin_logETimeSig);
			_hists2D.push_back(etaPhiCov_timeEtaCov);
			_hists2D.push_back(timeEtaCov_timePhiCov);
			_hists2D.push_back(etaPhiCov_timePhiCov);
                	_hists2D.push_back(etaPhiCov_timeMajCov);
                	_hists2D.push_back(etaPhiCov_timeMinCov);
                	_hists2D.push_back(timePhiCov_timeMajCov);
                	_hists2D.push_back(timePhiCov_timeMinCov);
                	_hists2D.push_back(timeEtaCov_timeMajCov);
                	_hists2D.push_back(timeEtaCov_timeMinCov);
			_hists2D.push_back(etaPhiCovUnnorm_timeEtaCovUnnorm); 
			_hists2D.push_back(timeEtaCovUnnorm_timePhiCovUnnorm);
			_hists2D.push_back(etaPhiCovUnnorm_timePhiCovUnnorm); 
                	_hists2D.push_back(etaPhiCovUnnorm_timeMajCovUnnorm); 
                	_hists2D.push_back(etaPhiCovUnnorm_timeMinCovUnnorm); 
                	_hists2D.push_back(timePhiCovUnnorm_timeMajCovUnnorm);
                	_hists2D.push_back(timePhiCovUnnorm_timeMinCovUnnorm);
                	_hists2D.push_back(timeEtaCovUnnorm_timeMajCovUnnorm);
                	_hists2D.push_back(timeEtaCovUnnorm_timeMinCovUnnorm);
                	_hists2D.push_back(rot2D_etaPhiCov);
                	_hists2D.push_back(rot2D_etaPhiCovUnnorm);
			_hists2D.push_back(E_phi);
                	_hists2D.push_back(rot2D_phiEll2D);
			_hists2D.push_back(etaPhiCov_timeEtaCovCounts);
			_hists2D.push_back(timeCenter_phiE2D);	
                	_hists2D.push_back(rot2D_etaPhiCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(rot2D_etaPhiCov_phiE2Dneq0PiOv2);
			_hists2D.push_back(phoE_phiE2D);	
			_hists2D.push_back(timeEtaCov_rot2D);
                	_hists2D.push_back(rot2D_timeEtaCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(rot2D_timeEtaCov_phiE2Dneq0PiOv2);
                	_hists2D.push_back(logErot2D_etaPhiCov);
                	_hists2D.push_back(logErot2D_etaPhiCovUnnorm);
			_hists2D.push_back(logEE_phi);
                	_hists2D.push_back(logErot2D_phiEll2D);
			_hists2D.push_back(logEetaPhiCov_timeEtaCovCounts);
			_hists2D.push_back(logEtimeCenter_phiE2D);	
                	_hists2D.push_back(logErot2D_etaPhiCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(logErot2D_etaPhiCov_phiE2Dneq0PiOv2);
			_hists2D.push_back(logEphoE_phiE2D);
			_hists2D.push_back(logEtimeEtaCov_rot2D);
                	_hists2D.push_back(logErot2D_timeEtaCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(logErot2D_timeEtaCov_phiE2Dneq0PiOv2);
                	_hists2D.push_back(noErot2D_etaPhiCov);
                	_hists2D.push_back(noErot2D_etaPhiCovUnnorm);
			_hists2D.push_back(noEE_phi);
                	_hists2D.push_back(noErot2D_phiEll2D);
			_hists2D.push_back(noEetaPhiCov_timeEtaCovCounts);
			_hists2D.push_back(noEtimeCenter_phiE2D);	
                	_hists2D.push_back(noErot2D_etaPhiCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(noErot2D_etaPhiCov_phiE2Dneq0PiOv2);
			_hists2D.push_back(noEphoE_phiE2D);	
			_hists2D.push_back(noEtimeEtaCov_rot2D);
                	_hists2D.push_back(noErot2D_timeEtaCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(noErot2D_timeEtaCov_phiE2Dneq0PiOv2);
                	_hists2D.push_back(etaPhiCov_timeEtaCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(etaPhiCov_timeEtaCov_phiE2Dneq0PiOv2);
			_hists2D.push_back(phoE_phiE2D_timeNeg10ToNeg2);	
			_hists2D.push_back(phoE_phiE2D_timeNeg2To5);
			_hists2D.push_back(phoE_phiE2D_time5To10);
			_hists2D.push_back(phoE_phiE2D_time10To15);
			_hists2D.push_back(timeCenter_etaCenter);
			_hists2D.push_back(timeCenter_etaCenter_phiE2Deq0PiOv2);
			_hists2D.push_back(timeCenter_etaCenter_phiE2Dneq0PiOv2);
                	_hists2D.push_back(logEetaPhiCov_timeEtaCov);
                	_hists2D.push_back(logEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(logEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2);
                	_hists2D.push_back(noEetaPhiCov_timeEtaCov);
                	_hists2D.push_back(noEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2);
                	_hists2D.push_back(noEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2);
			_hists2D.push_back(timeCenter_rot2D);
			_hists2D.push_back(timeCenter_swCross);
			_hists2D.push_back(timeCenter_swCross_phiE2Deq0PiOv2);
			_hists2D.push_back(timeCenter_swCross_phiE2Dneq0PiOv2);
                	_hists2D.push_back(rot2D_swCross);
                	_hists2D.push_back(rot2D_swCross_phiE2Deq0PiOv2);
                	_hists2D.push_back(rot2D_swCross_phiE2Dneq0PiOv2);
                	_hists2D.push_back(phiE2D_swCross);
			_hists2D.push_back(phoE_swCross);	
                	_hists2D.push_back(rot2D_phiEll2D_timeNeg10toNeg2);
                	_hists2D.push_back(rot2D_phiEll2D_timeNeg2to5);
                	_hists2D.push_back(rot2D_phiEll2D_time5to15);
			_hists2D.push_back(phoE_logEsmaj);
			_hists2D.push_back(phoE_logEsmin);
			_hists2D.push_back(phoE_logEetaSig);
			_hists2D.push_back(phoE_logEphiSig);
			_hists2D.push_back(phoE_logEetaPhiCov);
			_hists2D.push_back(phoE_logEtimeEtaCov);
			_hists2D.push_back(nRhs_logEsmaj);
			_hists2D.push_back(nRhs_logEsmin);
			_hists2D.push_back(nRhs_logEetaSig);
			_hists2D.push_back(nRhs_logEphiSig);
			_hists2D.push_back(nRhs_logEetaPhiCov);
			_hists2D.push_back(nRhs_logEtimeEtaCov);
			_hists2D.push_back(phoE_nRhs);		
			_hists2D.push_back(timeCenter_etaCenter_rot2Dge0p6le0p8);
			_hists2D.push_back(timeCenter_etaCenter_rot2Dle0p6ge0p8);
			_hists2D.push_back(rot2D_E_phiE2Deq0);
			_hists2D.push_back(rot2D_E_phiE2Dneq0);
                	_hists2D.push_back(phiCenter_rot2D_phiE2Deq0);
                	_hists2D.push_back(phiCenter_rot2D_phiE2Dneq0);
			_hists2D.push_back(phiCenter_etaCenter_rot2Dge0p6le0p8);
			_hists2D.push_back(phiCenter_etaCenter_rot2Dle0p6ge0p8);
			_hists2D.push_back(timeCenter_swCrossPrime);
			_hists2D.push_back(timeCenter_swCrossPrime_phiE2Deq0PiOv2);
			_hists2D.push_back(timeCenter_swCrossPrime_phiE2Dneq0PiOv2);
                	_hists2D.push_back(rot2D_swCrossPrime);
                	_hists2D.push_back(rot2D_swCrossPrime_phiE2Deq0PiOv2);
                	_hists2D.push_back(rot2D_swCrossPrime_phiE2Dneq0PiOv2);
                	_hists2D.push_back(phiE2D_swCrossPrime);
			_hists2D.push_back(phoE_swCrossPrime);	
			_hists2D.push_back(etaSig_phiSig_phiE2Deq0);
			_hists2D.push_back(etaSig_phiSig_phiE2Dneq0);
			_hists2D.push_back(etaCenter_phiCenter_Ege100le200_rot2Dge0p7le0p8);
			_hists2D.push_back(etaCenter_phiCenter_phiE2Deq0);
			_hists2D.push_back(timeCenter_etaCenter_Ege100le200_rot2Dge0p7le0p8);
			_hists2D.push_back(timeCenter_etaCenter_phiE2Deq0);
			_hists2D.push_back(timeCenter_etaCenter_phiE2Dneq0);
			_hists2D.push_back(rot2D_E_phiCentereq0Pi);
			_hists2D.push_back(timeCenter_etaCenter_Ele100ge200_rot2Dle0p7ge0p8);
			_hists2D.push_back(etaCenter_phiCenter_Ele100ge200_rot2Dle0p7ge0p8);
			_hists2D.push_back(etaCenter_phiCenter_phiE2Dneq0);
			_hists2D.push_back(timeCenter_etaCenter_phiE2Deq0_Ege100le200_rot2Dge0p7le0p8);
			_hists2D.push_back(timeCenter_etaCenter_phiE2Dneq0_Ege100le200_rot2Dge0p7le0p8);
			_hists2D.push_back(phiSig_phiCenter);
			_hists2D.push_back(phiSig_phiCenter_phiE2Deq0);
			_hists2D.push_back(phiSig_phiCenter_phiE2Dneq0);
			_hists2D.push_back(phiSig_etaCenter);
			_hists2D.push_back(phiSig_etaCenter_phiE2Deq0);
			_hists2D.push_back(phiSig_etaCenter_phiE2Dneq0);
			_hists2D.push_back(etaSig_phiCenter);
			_hists2D.push_back(etaSig_phiCenter_phiE2Deq0);
			_hists2D.push_back(etaSig_phiCenter_phiE2Dneq0);
			_hists2D.push_back(etaCenter_phiE2D);
			_hists2D.push_back(phoE_etaSig);
			_hists2D.push_back(phoE_etaSig_phiE2Deq0);
			_hists2D.push_back(phoE_etaSig_phiE2Dneq0);
			_hists2D.push_back(phoE_phiSig);
			_hists2D.push_back(phoE_phiSig_phiE2Deq0);
			_hists2D.push_back(phoE_phiSig_phiE2Dneq0);
			_hists2D.push_back(timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3);
			_hists2D.push_back(timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3);
			_hists2D.push_back(timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Deq0);
			_hists2D.push_back(timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Dneq0);
			_hists2D.push_back(timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Deq0);
			_hists2D.push_back(timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Dneq0);
			_hists2D.push_back(etaCenter_phiCenter_timeNeg15toNeg1);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg1to3);	
			_hists2D.push_back(etaCenter_phiCenter_time3to15);	
			_hists2D.push_back(etaCenter_phiCenter_phiSigle0p3ANDetaSigge0p3);	
			_hists2D.push_back(etaCenter_phiCenter_phiSigge0p3ORetaSigle0p3);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg15toNeg1_phiE2Deq0);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg15toNeg1_phiE2Dneq0);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg1to3_phiE2Deq0);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg1to3_phiE2Dneq0);	
			_hists2D.push_back(etaCenter_phiCenter_time3to15_phiE2Deq0);	
			_hists2D.push_back(etaCenter_phiCenter_time3to15_phiE2Dneq0);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3);	
			_hists2D.push_back(etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg15toNeg1_phiSigge0p3ORetaSigle0p3);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg1to3_phiSigge0p3ORetaSigle0p3);	
			_hists2D.push_back(etaCenter_phiCenter_time3to15_phiSigge0p3ORetaSigle0p3);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3_phiE2Deq0);	
			_hists2D.push_back(etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3_phiE2Deq0);	
			_hists2D.push_back(etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3_phiE2Deq0);	
			_hists2D.push_back(rot2D_E_phiCenterneq0Pi);
                	_hists2D.push_back(phiE2D_swCrossPrime_timeNeg15toNeg1);
                	_hists2D.push_back(phiE2D_swCrossPrime_timeNeg1to3);
                	_hists2D.push_back(phiE2D_swCrossPrime_time3to15);
			_hists2D.push_back(phiE2D_recoMet);
			_hists2D.push_back(etaSig_recoMet);
			_hists2D.push_back(phiSig_recoMet);
			_hists2D.push_back(etaCenter_recoMet);
			_hists2D.push_back(phiCenter_recoMet);
			_hists2D.push_back(timeCenter_recoMet);
			_hists2D.push_back(etaSig_timeEtaCov_timeNeg15toNeg1);
			_hists2D.push_back(etaSig_timeEtaCov_timeNeg15toNeg1_phiCenter0pi);
			_hists2D.push_back(etaSig_timeEtaCov_time3to15);
			_hists2D.push_back(etaSig_etaPhiCov_timeNeg15toNeg1);
			_hists2D.push_back(etaSig_etaPhiCov_timeNeg15toNeg1_phiCenter0pi);
			_hists2D.push_back(etaSig_etaPhiCov_time3to15);
			_hists2D.push_back(metPhi_phiCenter);
			_hists2D.push_back(metPhi_phiCenter_etaSigge0p3ANDphiSigle0p3);
			_hists2D.push_back(metPhi_phiCenter_etaSigge0p3ORphiSigle0p3);
			_hists2D.push_back(dRtrack_timeSubcl);	
			_hists2D.push_back(dEtrack_timeSubcl);	
			_hists2D.push_back(dRtrack_dEtrack);	
			_hists2D.push_back(dRtrack_dEtrack_early);	
			_hists2D.push_back(dRtrack_dEtrack_prompt);
			_hists2D.push_back(dRtrack_dEtrack_late);	
			_hists2D.push_back(rhEnergy_timesSigSqMeasErr);
			_hists2D.push_back(ENeighbors);	
			_hists2D.push_back(etaPhi_overlaidsubcl);
			_hists2D.push_back(phoE_timeSig);
			_hists2D.push_back(timeCenter_etaSig);
			_hists2D.push_back(timeCenter_phiSig);
			_hists2D.push_back(subclE_sqrtSubclEmultEtaSig);
			_hists2D.push_back(subclE_sqrtSubclEmultPhiSig);
			_hists2D.push_back(subclE_sqrtSubclEmultTimeSig);
			_hists2D.push_back(subclEmultTimeMajCov_subclEmultTimeEtaCov);
			_hists2D.push_back(timeEtaCov_timePhiCov_absTimeMajCovge0p1);
			_hists2D.push_back(timeMajCov_timeMinCov);
			_hists2D.push_back(tPlaneDist_planeEtaSig);
			_hists2D.push_back(tPlaneDist_planePhiSig);
			_hists2D.push_back(tPlaneDist_planeEtaPhiCov);
			


		};
	

		//6 - space slope
		TH1D* slope_space = new TH1D("slope_space","slope_space",50,-30,30);
		//7 - eta-time slope
		TH1D* slope_etaT = new TH1D("slope_etaT","slope_etaT",50,-2,2);
		//8 - phi-time slope
		TH1D* slope_phiT = new TH1D("slope_phiT","slope_phiT",50,-4,4);
		//9 - polar angle
		TH1D* polar_ang = new TH1D("polar_ang","polar_ang",50,-0.5,3.5);		
		//10 - azimuth angle
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",25,-0.1,6.3);		
		//11 - subcluster energy - total
		TH1D* e_subcl = new TH1D("e_subcl","e_subcl",50,0.,2000.);
		//12 - ellipsoid rotundity
		TH1D* rotundity_3D = new TH1D("rotundity_3D","rotundity_3D",10,0.74,1.01);
		//13 - spatial rotundity
		TH1D* rotundity_2D = new TH1D("rotundity_2D","rotundity_2D",50,0.45,1.1);
		//14 - velocity = z/r*k for k transfer factor to velocity units
		TH1D* velocity = new TH1D("velocity","velocity",31,-1.,30.);
		//15 - ratio of 2D eigenvals
		TH1D* eigen2D_ratio = new TH1D("eigen2D_ratio","eigen2D_ratio",50,0.,1.);
		//16 - eta sigma	
                TH1D* etaSig = new TH1D("etaSig","etaSig",50,0., 0.1);
		//17 - phi sigma	
                TH1D* phiSig = new TH1D("phiSig","phiSig",50,0.,0.1);
		//18 - time sigma	
                TH1D* timeSig = new TH1D("timeSig","timeSig",50,0,5.);
		//19 - fraction of energy in cluster
		TH1D* fracE = new TH1D("fracE","fracE",50,0.,1.1);
		//20 - azimuth angle in 2D
		TH1D* phiEll = new TH1D("phiE2D","phiE2D",50,-3.1,3.1);		
		//21 - eta sigma for positive eta clusters
                TH1D* etaSig_pos = new TH1D("etaSig_pos","etaSig_pos",25,0.01, 0.09);
		//22 - eta sigma for negative eta clusters
                TH1D* etaSig_neg = new TH1D("etaSig_neg","etaSig_neg",25,0.01, 0.09);
		//23 - normalized covariance - eta/phi
		TH1D* etaphi_cov = new TH1D("etaphi_cov","etaphi_cov",25,-0.1,0.1);
		//24 - normalized covariance - time/eta
		TH1D* timeeta_cov = new TH1D("timeeta_cov","timeeta_cov",25,-0.1,0.1);
		//25 - normalized covariance - time/phi
		TH1D* timephi_cov = new TH1D("timephi_cov","timephi_cov",25,-0.1,0.1);
		//26 - normalized covariance - time/major axis
		TH1D* timemaj_cov = new TH1D("timemaj_cov","timemaj_cov",25,-0.3,0.3);
		//27 - normalized covariance - time/minor axis
		TH1D* timemin_cov = new TH1D("timemin_cov","timemin_cov",25,-0.3,0.3);
		//28 - time sigma - low photon pt	
                TH1D* timeSig_gamPtLow = new TH1D("timeSig_Elow","timeSig_Elow",25,0,3.);
		//29 - time sigma - med photon pt
                TH1D* timeSig_gamPtMed = new TH1D("timeSig_Emed","timeSig_Emed",25,0,3.);
		//30 - time sigma - high photon pt
                TH1D* timeSig_gamPtHi = new TH1D("timeSig_Ehi","timeSig_Ehi",25,0,2.);
		//31 - unnormalized covariance - time/eta
		TH1D* timeeta_covUnnorm = new TH1D("timeeta_covUnnorm","timeeta_covUnnorm",25,-0.2,0.2);
		//32 - unnormalized covariance - time/phi
		TH1D* timephi_covUnnorm = new TH1D("timephi_covUnnorm","timephi_covUnnorm",25,-0.2,0.2);
		//33 - unnormalized covariance - time/maj
		TH1D* timeMaj_covUnnorm = new TH1D("timeMaj_covUnnorm","timeMaj_covUnnorm",25,-0.5,0.5);
		//34 - unnormalized covariance - time/min
		TH1D* timeMin_covUnnorm = new TH1D("timeMin_covUnnorm","timeMin_covUnnorm",25,-0.5,0.5);
		//35 - logE CMS normalized covariance - time/eta
		TH1D* logE_timeeta_cov = new TH1D("logE_timeeta_cov","logE_timeeta_cov",25,-1.,1.);
		//36 - logE CMS normalized covariance - time/phi
		TH1D* logE_timephi_cov = new TH1D("logE_timephi_cov","logE_timephi_cov",25,-1.,1.);
		//37 - logE CMS normalized covariance - eta/phi
		TH1D* logE_etaphi_cov = new TH1D("logE_etaphi_cov","logE_etaphi_cov",25,-1.,1);
		//38 - logE CMS eta sigma	
                TH1D* logE_etaSig = new TH1D("logE_etaSig","logE_etaSig",25,0.01, 0.09);
		//39 - logE CMS phi sigma	
                TH1D* logE_phiSig = new TH1D("logE_phiSig","logE_phiSig",25,0.01,0.09);
		//40 - logE CMS time sigma	
                TH1D* logE_timeSig = new TH1D("logE_timeSig","logE_timeSig",25,0,10.);
		//41 - CMS normalized covariance - time/eta
		TH1D* noE_timeeta_cov = new TH1D("noE_timeeta_cov","noE_timeeta_cov",25,-1.,1.);
		//42 - CMS normalized covariance - time/phi
		TH1D* noE_timephi_cov = new TH1D("noE_timephi_cov","noE_timephi_cov",25,-1.,1.);		
		//43 - CMS normalized covariance - eta/phi
		TH1D* noE_etaphi_cov = new TH1D("noE_etaphi_cov","noE_etaphi_cov",25,-1.,1.);
		//44 - CMS eta sigma	
                TH1D* noE_etaSig = new TH1D("noE_etaSig","noE_etaSig",25,0.01, 0.09);
		//45 - CMS phi sigma	
                TH1D* noE_phiSig = new TH1D("noE_phiSig","noE_phiSig",25,0.01,0.09);
		//46 - CMS time sigma	
                TH1D* noE_timeSig = new TH1D("noE_timeSig","noE_timeSig",25,0,10.);
		//47 - CMS unnormalized covariance - time/eta
		TH1D* noE_timeeta_covUnnorm = new TH1D("noE_timeeta_covUnnorm","noE_timeeta_covUnnorm",25,-1.,1.);
		//48 - CMS unnormalized covariance - time/phi
		TH1D* noE_timephi_covUnnorm = new TH1D("noE_timephi_covUnnorm","noE_timephi_covUnnorm",25,-1.,1.);
		//49 - CMS unnormalized covariance - eta/phi
		TH1D* noE_etaphi_covUnnorm = new TH1D("noE_etaphi_covUnnorm","noE_etaphi_covUnnorm",25,-1.,1.);
		//50 - logE CMS unnormalized covariance - time/eta
		TH1D* logE_timeeta_covUnnorm = new TH1D("logE_timeeta_covUnnorm","logE_timeeta_covUnnorm",25,-1.,1.);
		//51 - logE CMS unnormalized covariance - time/phi
		TH1D* logE_timephi_covUnnorm = new TH1D("logE_timephi_covUnnorm","logE_timephi_covUnnorm",25,-1.,1.);
		//52 - logE CMS normalized covariance - eta/phi
		TH1D* logE_etaphi_covUnnorm = new TH1D("logE_etaphi_covUnnorm","logE_etaphi_covUnnorm",25,-1.,1.);
		//53 - unnormalized covariance - eta/phi
		TH1D* etaphi_covUnnorm = new TH1D("etaphi_covUnnorm","etaphi_covUnnorm",25,-0.1,0.1);
		//54 - CMS smaj
		TH1D* noE_smaj = new TH1D("noE_smaj","noE_smaj",25,0.,0.004);		
		//55 - CMS smin
		TH1D* noE_smin = new TH1D("noE_smin","noE_smin",25,0.,0.004);		
		//56 - CMS logE smaj
		TH1D* logE_smaj = new TH1D("logE_smaj","logE_smaj",25,0.,0.004);		
		//57 - CMS logE smin
		TH1D* logE_smin = new TH1D("logE_smin","logE_smin",25,0.,0.004);		
		//58 - CMS time center
		TH1D* noE_time_center = new TH1D("noE_time_center","noE_time_center",50,-20,20);
		//59 - CMS eta center
		TH1D* noE_eta_center = new TH1D("noE_eta_center","noE_eta_center",50,-3.5,3.5);
		//60 - CMS phi center
		TH1D* noE_phi_center = new TH1D("noE_phi_center","noE_phi_center",50,-0.1,6.3);
		//61 - CMS ellipsoid angle (phi2D)
		TH1D* noE_phiEll = new TH1D("noE_phiE2D","noE_phiE2D",50,-3,3.);		
		//62 - CMS time-smaj covariance
		TH1D* noE_timeSmaj_cov = new TH1D("noE_timeSmaj_cov","noE_timeSmaj_cov",25,-1.,1.);
		//63 - CMS time-smin covariance
		TH1D* noE_timeSmin_cov = new TH1D("noE_timeSmin_cov","noE_timeSmin_cov",25,-1.,1.);
		//64 - CMS time-smaj covariance unnormalized 
		TH1D* noE_timeSmaj_covUnnorm = new TH1D("noE_timeSmaj_covUnnorm","noE_timeSmaj_covUnnorm",25,-0.5,0.5);
		//65 - CMS time-smin covariance unnormalized 
		TH1D* noE_timeSmin_covUnnorm = new TH1D("noE_timeSmin_covUnnorm","noE_timeSmin_covUnnorm",25,-0.5,0.5);
		//66 - CMS rotundity 2D
		TH1D* noE_rotundity_2D = new TH1D("noE_rotundity_2D","noE_rotundity_2D",20,0.4,1.1);
		//67 - CMS logE time center
		TH1D* logE_time_center = new TH1D("logE_time_center","logE_time_center",50,-20,20);
		//68 - CMS logE eta center
		TH1D* logE_eta_center = new TH1D("logE_eta_center","logE_eta_center",50,-3.5,3.5);
		//69 - CMS logE phi center
		TH1D* logE_phi_center = new TH1D("logE_phi_center","logE_phi_center",50,-0.1,6.3);
		//70 - CMS logE ellipsoid angle (phi2D)
		TH1D* logE_phiEll = new TH1D("logE_phiE2D","logE_phiE2D",50,-3,3.);	
		//71 - CMS logE time-smaj covariance
		TH1D* logE_timeSmaj_cov = new TH1D("logE_timeSmaj_cov","logE_timeSmaj_cov",25,-1.,1.);
		//72 - CMS logE time-smin covariance
		TH1D* logE_timeSmin_cov = new TH1D("logE_timeSmin_cov","logE_timeSmin_cov",25,-1.,1.);
		//73 - CMS logE time-smaj covariance unnormalized 
		TH1D* logE_timeSmaj_covUnnorm = new TH1D("logE_timeSmaj_covUnnorm","logE_timeSmaj_covUnnorm",25,-0.5,0.5);
		//74 - CMS logE time-smin covariance unnormalized 
		TH1D* logE_timeSmin_covUnnorm = new TH1D("logE_timeSmin_covUnnorm","logE_timeSmin_covUnnorm",25,-0.5,0.5);
		//75 - CMS logE rotundity 2D
		TH1D* logE_rotundity_2D = new TH1D("logE_rotundity_2D","logE_rotundity_2D",20,0.4,1.1);
		//76 - spatial rotundity, E low
		TH1D* rotundity_2D_Elo = new TH1D("rot2D_Elo","rot2D_Elo",20,0.4,1.1);
		//77 - spatial rotundity, E med
		TH1D* rotundity_2D_Emed = new TH1D("rot2D_Emed","rot2D_Emed",20,0.4,1.1);
		//78 - spatial rotundity, E high
		TH1D* rotundity_2D_Ehi = new TH1D("rot2D_Ehi","rot2D_Ehi",20,0.4,1.1);
		//Quadrants in etaPhiCov-timeEtaCov space
		//Q1 = ep_cov < 0 and te_cov < 0 
		//Q2 = ep_cov < 0 and te_cov > 0 
		//Q3 = ep_cov > 0 and te_cov < 0 
		//Q4 = ep_cov > 0 and te_cov > 0  
		//79 - rotundity in quadrant 1
		TH1D* rotundity_2D_1Q = new TH1D("rot2D_1Q","rot2D_1Q",20,0.4,1.1);
		//80 - rotundity in quadrant 2
		TH1D* rotundity_2D_2Q = new TH1D("rot2D_2Q","rot2D_2Q",20,0.4,1.1);
		//81 - rotundity in quadrant 3
		TH1D* rotundity_2D_3Q = new TH1D("rot2D_3Q","rot2D_3Q",20,0.4,1.1);
		//82 - rotundity in quadrant 4
		TH1D* rotundity_2D_4Q = new TH1D("rot2D_4Q","rot2D_4Q",20,0.4,1.1);
		//83 - phiE2D in quadrant 1
		TH1D* phiE_2D_1Q = new TH1D("phiE2D_1Q","phiE2D_1Q",20,-3.1,1);
		//84 - phiE2D in quadrant 2
		TH1D* phiE_2D_2Q = new TH1D("phiE2D_2Q","phiE2D_2Q",20,-3.1,1);
		//85 - phiE2D in quadrant 3
		TH1D* phiE_2D_3Q = new TH1D("phiE2D_3Q","phiE2D_3Q",20,-3.1,1);
		//86 - phiE2D in quadrant 4
		TH1D* phiE_2D_4Q = new TH1D("phiE2D_4Q","phiE2D_4Q",20,-3.1,1);
		//87 - etaphicov in quadrant 1
		TH1D* etaPhiCov_1Q = new TH1D("etaPhiCov_1Q","etaPhiCov_1Q",20,-1,0.1);
		//88 - etaphicov in quadrant 2
		TH1D* etaPhiCov_2Q = new TH1D("etaPhiCov_2Q","etaPhiCov_2Q",20,-1,0.1);
		//89 - etaphicov in quadrant 3
		TH1D* etaPhiCov_3Q = new TH1D("etaPhiCov_3Q","etaPhiCov_3Q",20,-0.1,1.1);
		//90 - etaphicov in quadrant 4
		TH1D* etaPhiCov_4Q = new TH1D("etaPhiCov_4Q","etaPhiCov_4Q",20,-0.1,1.1);
		//91 - timeetacov in quadrant 1
		TH1D* timeEtaCov_1Q = new TH1D("timeEtaCov_1Q","timeEtaCov_1Q",20,-1,0.1);
		//92 - timeetacov in quadrant 2
		TH1D* timeEtaCov_2Q = new TH1D("timeEtaCov_2Q","timeEtaCov_2Q",20,-0.1,1);
		//93 - timeetacov in quadrant 3
		TH1D* timeEtaCov_3Q = new TH1D("timeEtaCov_3Q","timeEtaCov_3Q",20,-1,0.1);
		//94 - timeetacov in quadrant 4
		TH1D* timeEtaCov_4Q = new TH1D("timeEtaCov_4Q","timeEtaCov_4Q",20,-0.1,1);
                //95 - subcluster eta center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* etaCenter_phiE2Deq0PiOv2 = new TH1D("etaCenter_phiE2Deq0PiOv2","etaCenter_phiE2Deq0PiOv2 ",25,-3.5,3.5);
                //96 - subcluster eta center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* etaCenter_phiE2Dneq0PiOv2 = new TH1D("etaCenter_phiE2Dneq0PiOv2","etaCenter_phiE2Dneq0PiOv2",25,-3.5,3.5);
                //97 - subcluster phi center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* phiCenter_phiE2Deq0PiOv2 = new TH1D("phiCenter_phiE2Deq0PiOv2","phiCenter_phiE2Deq0PiOv2",25,-0.6,6.6);
                //98 - subcluster phi center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* phiCenter_phiE2Dneq0PiOv2 = new TH1D("phiCenter_phiE2Dneq0PiOv2","phiCenter_phiE2Dneq0PiOv2",25,-0.6,6.6);
                //99 - subcluster time center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* timeCenter_phiE2Deq0PiOv2 = new TH1D("timeCenter_phiE2Deq0PiOv2","timeCenter_phiE2Deq0PiOv2",25,-20,20);
                //100 - subcluster time center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* timeCenter_phiE2Dneq0PiOv2 = new TH1D("timeCenter_phiE2Dneq0PiOv2","timeCenter_phiE2Dneq0PiOv2",25,-20,20);
                //101 - pho E (rh sum), phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* phoE_phiE2Deq0PiOv2 = new TH1D("phoE_phiE2Deq0PiOv2","phoE_phiE2Deq0PiOv2",25,0,1000);
                //102 - pho E (rh sum), phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* phoE_phiE2Dneq0PiOv2 = new TH1D("phoE_phiE2Dneq0PiOv2","phoE_phiE2Dneq0PiOv2",25,0,1000);
                //103 - nsubclusters, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* nSubClusters_phiE2Deq0PiOv2 = new TH1D("nSubClusters_phiE2Deq0PiOv2","nSubClusters_phiE2Deq0PiOv2",10,0,10);
                //104 - nsubclusters, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* nSubClusters_phiE2Dneq0PiOv2 = new TH1D("nSubClusters_phiE2Dneq0PiOv2","nSubClusters_phiE2Dneq0PiOv2",10,0,10);
                //105 - eta sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* etaSig_phiE2Deq0PiOv2 = new TH1D("etaSig_phiE2Deq0PiOv2","etaSig_phiE2Deq0PiOv2",25,0.01,0.05);
                //106 - eta sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* etaSig_phiE2Dneq0PiOv2 = new TH1D("etaSig_phiE2Dneq0PiOv2","etaSig_phiE2Dneq0PiOv2",25,0.01,0.05);
                //107 - phi sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* phiSig_phiE2Deq0PiOv2 = new TH1D("phiSig_phiE2Deq0PiOv2","phiSig_phiE2Deq0PiOv2",25,0.01,0.09);
                //108 - phi sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* phiSig_phiE2Dneq0PiOv2 = new TH1D("phiSig_phiE2Dneq0PiOv2","phiSig_phiE2Dneq0PiOv2",25,0.01,0.09);
                //109 - time sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* timeSig_phiE2Deq0PiOv2 = new TH1D("timeSig_phiE2Deq0PiOv2","timeSig_phiE2Deq0PiOv2",25,0.,4);
                //110 - time sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* timeSig_phiE2Dneq0PiOv2 = new TH1D("timeSig_phiE2Dneq0PiOv2","timeSig_phiE2Dneq0PiOv2",25,0.,4.);
                //111 - etaPhi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* etaPhiCov_phiE2Deq0PiOv2 = new TH1D("etaPhiCov_phiE2Deq0PiOv2","etaPhiCov_phiE2Deq0PiOv2",25,-1,1);
                //112 - etaPhi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* etaPhiCov_phiE2Dneq0PiOv2 = new TH1D("etaPhiCov_phiE2Dneq0PiOv2","etaPhiCov_phiE2Dneq0PiOv2",25,-1,1);
                //113 - timeEta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* timeEtaCov_phiE2Deq0PiOv2 = new TH1D("timeEtaCov_phiE2Deq0PiOv2","timeEtaCov_phiE2Deq0PiOv2",25,-1,1);
                //114 - timeEta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* timeEtaCov_phiE2Dneq0PiOv2 = new TH1D("timeEtaCov_phiE2Dneq0PiOv2","timeEtaCov_phiE2Dneq0PiOv2",25,-1,1);
                //115 - timePhi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* timePhiCov_phiE2Deq0PiOv2 = new TH1D("timePhiCov_phiE2Deq0PiOv2","timePhiCov_phiE2Deq0PiOv2",25,-1,1);
                //116 - timePhi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* timePhiCov_phiE2Dneq0PiOv2 = new TH1D("timePhiCov_phiE2Dneq0PiOv2","timePhiCov_phiE2Dneq0PiOv2",25,-1,1);
                //117 - timeMaj cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* timeMajCov_phiE2Deq0PiOv2 = new TH1D("timeMajCov_phiE2Deq0PiOv2","timeMajCov_phiE2Deq0PiOv2",25,-1,1);
                //118 - timeMaj cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* timeMajCov_phiE2Dneq0PiOv2 = new TH1D("timeMajCov_phiE2Dneq0PiOv2","timeMajCov_phiE2Dneq0PiOv2",25,-1,1);
                //119 - timeMin cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* timeMinCov_phiE2Deq0PiOv2 = new TH1D("timeMinCov_phiE2Deq0PiOv2","timeMinCov_phiE2Deq0PiOv2",25,-1,1);
                //120 - timeMin cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* timeMinCov_phiE2Dneq0PiOv2 = new TH1D("timeMinCov_phiE2Dneq0PiOv2","timeMinCov_phiE2Dneq0PiOv2",25,-1,1);
                //121 - phiE2D, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* phiE2D_phiE2Deq0PiOv2 = new TH1D("phiE2D_phiE2Deq0PiOv2","phiE2D_phiE2Deq0PiOv2",25,-3.1,3.1);
                //122 - phiE2D, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* phiE2D_phiE2Dneq0PiOv2 = new TH1D("phiE2D_phiE2Dneq0PiOv2","phiE2D_phiE2Dneq0PiOv2",25,-3.1,3.1);
                //123 - rot2D cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* rot2D_phiE2Deq0PiOv2 = new TH1D("rot2D_phiE2Deq0PiOv2","rot2D_phiE2Deq0PiOv2",25,0.4,1.1);
                //124 - rot2D cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* rot2D_phiE2Dneq0PiOv2 = new TH1D("rot2D_phiE2Dneq0PiOv2","rot2D_phiE2Dneq0PiOv2",25,0.4,1.1);
		
		
		//////////logE weighted//////////
		//125 - spatial rotundity, E low
		TH1D* logErotundity_2D_Elo = new TH1D("logErot2D_Elo","logErot2D_Elo",20,0.4,1.1);
		//126 - spatial rotundity, E med
		TH1D* logErotundity_2D_Emed = new TH1D("logErot2D_Emed","logErot2D_Emed",20,0.4,1.1);
		//127 - spatial rotundity, E high
		TH1D* logErotundity_2D_Ehi = new TH1D("logErot2D_Ehi","logErot2D_Ehi",20,0.4,1.1);
		//Quadrants in etaPhiCov-timeEtaCov space
		//Q1 = ep_cov < 0 and te_cov < 0 
		//Q2 = ep_cov < 0 and te_cov > 0 
		//Q3 = ep_cov > 0 and te_cov < 0 
		//Q4 = ep_cov > 0 and te_cov > 0  
		//128 - rotundity in quadrant 1
		TH1D* logErotundity_2D_1Q = new TH1D("logErot2D_1Q","logErot2D_1Q",20,0.4,1.1);
		//129 - rotundity in quadrant 2
		TH1D* logErotundity_2D_2Q = new TH1D("logErot2D_2Q","logErot2D_2Q",20,0.4,1.1);
		//130 - rotundity in quadrant 3
		TH1D* logErotundity_2D_3Q = new TH1D("logErot2D_3Q","logErot2D_3Q",20,0.4,1.1);
		//131 - rotundity in quadrant 4
		TH1D* logErotundity_2D_4Q = new TH1D("logErot2D_4Q","logErot2D_4Q",20,0.4,1.1);
		//132 - phiE2D in quadrant 1
		TH1D* logEphiE_2D_1Q = new TH1D("logEphiE2D_1Q","logEphiE2D_1Q",20,-3.1,1);
		//133 - phiE2D in quadrant 2
		TH1D* logEphiE_2D_2Q = new TH1D("logEphiE2D_2Q","logEphiE2D_2Q",20,-3.1,1);
		//134 - phiE2D in quadrant 3
		TH1D* logEphiE_2D_3Q = new TH1D("logEphiE2D_3Q","logEphiE2D_3Q",20,-3.1,1);
		//135 - phiE2D in quadrant 4
		TH1D* logEphiE_2D_4Q = new TH1D("logEphiE2D_4Q","logEphiE2D_4Q",20,-3.1,1);
		//136 - etaphicov in quadrant 1
		TH1D* logEetaPhiCov_1Q = new TH1D("logEetaPhiCov_1Q","logEetaPhiCov_1Q",20,-1,0.1);
		//137 - etaphicov in quadrant 2
		TH1D* logEetaPhiCov_2Q = new TH1D("logEetaPhiCov_2Q","logEetaPhiCov_2Q",20,-1,0.1);
		//138 - etaphicov in quadrant 3
		TH1D* logEetaPhiCov_3Q = new TH1D("logEetaPhiCov_3Q","logEetaPhiCov_3Q",20,-0.1,1.1);
		//139 - etaphicov in quadrant 4
		TH1D* logEetaPhiCov_4Q = new TH1D("logEetaPhiCov_4Q","logEetaPhiCov_4Q",20,-0.1,1.1);
		//140 - timeetacov in quadrant 1
		TH1D* logEtimeEtaCov_1Q = new TH1D("logEtimeEtaCov_1Q","logEtimeEtaCov_1Q",20,-1,0.1);
		//141 - timeetacov in quadrant 2
		TH1D* logEtimeEtaCov_2Q = new TH1D("logEtimeEtaCov_2Q","logEtimeEtaCov_2Q",20,-0.1,1);
		//142 - timeetacov in quadrant 3
		TH1D* logEtimeEtaCov_3Q = new TH1D("logEtimeEtaCov_3Q","logEtimeEtaCov_3Q",20,-1,0.1);
		//143 - timeetacov in quadrant 4
		TH1D* logEtimeEtaCov_4Q = new TH1D("logEtimeEtaCov_4Q","logEtimeEtaCov_4Q",20,-0.1,1);
                //144 - subcluster eta center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEetaCenter_phiE2Deq0PiOv2 = new TH1D("logEetaCenter_phiE2Deq0PiOv2","logEetaCenter_phiE2Deq0PiOv2 ",25,-3.5,3.5);
                //145 - subcluster eta center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEetaCenter_phiE2Dneq0PiOv2 = new TH1D("logEetaCenter_phiE2Dneq0PiOv2","logEetaCenter_phiE2Dneq0PiOv2",25,-3.5,3.5);
                //146 - subcluster phi center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEphiCenter_phiE2Deq0PiOv2 = new TH1D("logEphiCenter_phiE2Deq0PiOv2","logEphiCenter_phiE2Deq0PiOv2",25,-3.5,3.5);
                //147 - subcluster phi center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEphiCenter_phiE2Dneq0PiOv2 = new TH1D("logEphiCenter_phiE2Dneq0PiOv2","logEphiCenter_phiE2Dneq0PiOv2",25,-3.5,3.5);
                //148 - subcluster time center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEtimeCenter_phiE2Deq0PiOv2 = new TH1D("logEtimeCenter_phiE2Deq0PiOv2","logEtimeCenter_phiE2Deq0PiOv2",25,-20,20);
                //149 - subcluster time center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEtimeCenter_phiE2Dneq0PiOv2 = new TH1D("logEtimeCenter_phiE2Dneq0PiOv2","logEtimeCenter_phiE2Dneq0PiOv2",25,-20,20);
                //150 - pho E (rh sum), phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEphoE_phiE2Deq0PiOv2 = new TH1D("logEphoE_phiE2Deq0PiOv2","logEphoE_phiE2Deq0PiOv2",25,0,1000);
                //151 - pho E (rh sum), phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEphoE_phiE2Dneq0PiOv2 = new TH1D("logEphoE_phiE2Dneq0PiOv2","logEphoE_phiE2Dneq0PiOv2",25,0,1000);
                //152 - nsubclusters, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEnSubClusters_phiE2Deq0PiOv2 = new TH1D("logEnSubClusters_phiE2Deq0PiOv2","logEnSubClusters_phiE2Deq0PiOv2",10,0,10);
                //153 - nsubclusters, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEnSubClusters_phiE2Dneq0PiOv2 = new TH1D("logEnSubClusters_phiE2Dneq0PiOv2","logEnSubClusters_phiE2Dneq0PiOv2",10,0,10);
                //154 - eta sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEetaSig_phiE2Deq0PiOv2 = new TH1D("logEetaSig_phiE2Deq0PiOv2","logEetaSig_phiE2Deq0PiOv2",25,0.01,0.09);
                //155 - eta sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEetaSig_phiE2Dneq0PiOv2 = new TH1D("logEetaSig_phiE2Dneq0PiOv2","logEetaSig_phiE2Dneq0PiOv2",25,0.01,0.09);
                //156 - phi sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEphiSig_phiE2Deq0PiOv2 = new TH1D("logEphiSig_phiE2Deq0PiOv2","logEphiSig_phiE2Deq0PiOv2",25,0.01,0.09);
                //157 - phi sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEphiSig_phiE2Dneq0PiOv2 = new TH1D("logEphiSig_phiE2Dneq0PiOv2","logEphiSig_phiE2Dneq0PiOv2",25,0.,0.09);
                //158 - time sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEtimeSig_phiE2Deq0PiOv2 = new TH1D("logEtimeSig_phiE2Deq0PiOv2","logEtimeSig_phiE2Deq0PiOv2",25,0.,10.);
                //159 - time sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEtimeSig_phiE2Dneq0PiOv2 = new TH1D("logEtimeSig_phiE2Dneq0PiOv2","logEtimeSig_phiE2Dneq0PiOv2",25,0.,10.);
                //160 - etaPhi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEetaPhiCov_phiE2Deq0PiOv2 = new TH1D("logEetaPhiCov_phiE2Deq0PiOv2","logEetaPhiCov_phiE2Deq0PiOv2",25,-1,1);
                //161 - etaPhi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEetaPhiCov_phiE2Dneq0PiOv2 = new TH1D("logEetaPhiCov_phiE2Dneq0PiOv2","logEetaPhiCov_phiE2Dneq0PiOv2",25,-1,1);
                //162 - timeEta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEtimeEtaCov_phiE2Deq0PiOv2 = new TH1D("logEtimeEtaCov_phiE2Deq0PiOv2","logEtimeEtaCov_phiE2Deq0PiOv2",25,-1,1);
                //163 - timeEta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEtimeEtaCov_phiE2Dneq0PiOv2 = new TH1D("logEtimeEtaCov_phiE2Dneq0PiOv2","logEtimeEtaCov_phiE2Dneq0PiOv2",25,-1,1);
                //164 - timePhi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEtimePhiCov_phiE2Deq0PiOv2 = new TH1D("logEtimePhiCov_phiE2Deq0PiOv2","logEtimePhiCov_phiE2Deq0PiOv2",25,-1,1);
                //165 - timePhi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEtimePhiCov_phiE2Dneq0PiOv2 = new TH1D("logEtimePhiCov_phiE2Dneq0PiOv2","logEtimePhiCov_phiE2Dneq0PiOv2",25,-1,1);
                //166 - timeMaj cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEtimeMajCov_phiE2Deq0PiOv2 = new TH1D("logEtimeMajCov_phiE2Deq0PiOv2","logEtimeMajCov_phiE2Deq0PiOv2",25,-1,1);
                //167 - timeMaj cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEtimeMajCov_phiE2Dneq0PiOv2 = new TH1D("logEtimeMajCov_phiE2Dneq0PiOv2","logEtimeMajCov_phiE2Dneq0PiOv2",25,-1,1);
		//168 - timeMin cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEtimeMinCov_phiE2Deq0PiOv2 = new TH1D("logEtimeMinCov_phiE2Deq0PiOv2","logEtimeMinCov_phiE2Deq0PiOv2",25,-1,1);
                //169 - timeMin cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEtimeMinCov_phiE2Dneq0PiOv2 = new TH1D("logEtimeMinCov_phiE2Dneq0PiOv2","logEtimeMinCov_phiE2Dneq0PiOv2",25,-1,1);
                //170 - phiE2D, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logEphiE2D_phiE2Deq0PiOv2 = new TH1D("logEphiE2D_phiE2Deq0PiOv2","logEphiE2D_phiE2Deq0PiOv2",25,-3.1,3.1);
                //171 - phiE2D, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logEphiE2D_phiE2Dneq0PiOv2 = new TH1D("logEphiE2D_phiE2Dneq0PiOv2","logEphiE2D_phiE2Dneq0PiOv2",25,-3.1,3.1);
                //172 - rot2D cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* logErot2D_phiE2Deq0PiOv2 = new TH1D("logErot2D_phiE2Deq0PiOv2","logErot2D_phiE2Deq0PiOv2",25,0.1,1.1);
                //173 - rot2D cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* logErot2D_phiE2Dneq0PiOv2 = new TH1D("logErot2D_phiE2Dneq0PiOv2","logErot2D_phiE2Dneq0PiOv2",25,0.4,1.1);
		
		//////////noE weighted//////////
		//174 - spatial rotundity, E low
		TH1D* noErotundity_2D_Elo = new TH1D("noErot2D_Elo","noErot2D_Elo",20,0.4,1.1);
		//175 - spatial rotundity, E med
		TH1D* noErotundity_2D_Emed = new TH1D("noErot2D_Emed","noErot2D_Emed",20,0.4,1.1);
		//176 - spatial rotundity, E high
		TH1D* noErotundity_2D_Ehi = new TH1D("noErot2D_Ehi","noErot2D_Ehi",20,0.4,1.1);
		//Quadrants in etaPhiCov-timeEtaCov space
		//Q1 = ep_cov < 0 and te_cov < 0 
		//Q2 = ep_cov < 0 and te_cov > 0 
		//Q3 = ep_cov > 0 and te_cov < 0 
		//Q4 = ep_cov > 0 and te_cov > 0  
		//177 - rotundity in quadrant 1
		TH1D* noErotundity_2D_1Q = new TH1D("noErot2D_1Q","noErot2D_1Q",20,0.4,1.1);
		//178 - rotundity in quadrant 2
		TH1D* noErotundity_2D_2Q = new TH1D("noErot2D_2Q","noErot2D_2Q",20,0.4,1.1);
		//179 - rotundity in quadrant 3
		TH1D* noErotundity_2D_3Q = new TH1D("noErot2D_3Q","noErot2D_3Q",20,0.4,1.1);
		//180 - rotundity in quadrant 4
		TH1D* noErotundity_2D_4Q = new TH1D("noErot2D_4Q","noErot2D_4Q",20,0.4,1.1);
		//181 - phiE2D in quadrant 1
		TH1D* noEphiE_2D_1Q = new TH1D("noEphiE2D_1Q","noEphiE2D_1Q",20,-3.1,1);
		//182 - phiE2D in quadrant 2
		TH1D* noEphiE_2D_2Q = new TH1D("noEphiE2D_2Q","noEphiE2D_2Q",20,-3.1,1);
		//183 - phiE2D in quadrant 3
		TH1D* noEphiE_2D_3Q = new TH1D("noEphiE2D_3Q","noEphiE2D_3Q",20,-3.1,1);
		//184 - phiE2D in quadrant 4
		TH1D* noEphiE_2D_4Q = new TH1D("noEphiE2D_4Q","noEphiE2D_4Q",20,-3.1,1);
		//185 - etaphicov in quadrant 1
		TH1D* noEetaPhiCov_1Q = new TH1D("noEetaPhiCov_1Q","noEetaPhiCov_1Q",20,-1,0.1);
		//186 - etaphicov in quadrant 2
		TH1D* noEetaPhiCov_2Q = new TH1D("noEetaPhiCov_2Q","noEetaPhiCov_2Q",20,-1,0.1);
		//187 - etaphicov in quadrant 3
		TH1D* noEetaPhiCov_3Q = new TH1D("noEetaPhiCov_3Q","noEetaPhiCov_3Q",20,-0.1,1.1);
		//188 - etaphicov in quadrant 4
		TH1D* noEetaPhiCov_4Q = new TH1D("noEetaPhiCov_4Q","noEetaPhiCov_4Q",20,-0.1,1.1);
		//189 - timeetacov in quadrant 1
		TH1D* noEtimeEtaCov_1Q = new TH1D("noEtimeEtaCov_1Q","noEtimeEtaCov_1Q",20,-1,0.1);
		//190 - timeetacov in quadrant 2
		TH1D* noEtimeEtaCov_2Q = new TH1D("noEtimeEtaCov_2Q","noEtimeEtaCov_2Q",20,-0.1,1);
		//191 - timeetacov in quadrant 3
		TH1D* noEtimeEtaCov_3Q = new TH1D("noEtimeEtaCov_3Q","noEtimeEtaCov_3Q",20,-1,0.1);
		//192 - timeetacov in quadrant 4
		TH1D* noEtimeEtaCov_4Q = new TH1D("noEtimeEtaCov_4Q","noEtimeEtaCov_4Q",20,-0.1,1);
                //193 - subcluster eta center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEetaCenter_phiE2Deq0PiOv2 = new TH1D("noEetaCenter_phiE2Deq0PiOv2","noEetaCenter_phiE2Deq0PiOv2 ",25,-3.5,3.5);
                //194 - subcluster eta center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEetaCenter_phiE2Dneq0PiOv2 = new TH1D("noEetaCenter_phiE2Dneq0PiOv2","noEetaCenter_phiE2Dneq0PiOv2",25,-3.5,3.5);
                //195 - subcluster phi center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEphiCenter_phiE2Deq0PiOv2 = new TH1D("noEphiCenter_phiE2Deq0PiOv2","noEphiCenter_phiE2Deq0PiOv2",25,-3.5,3.5);
                //196 - subcluster phi center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEphiCenter_phiE2Dneq0PiOv2 = new TH1D("noEphiCenter_phiE2Dneq0PiOv2","noEphiCenter_phiE2Dneq0PiOv2",25,-3.5,3.5);
                //197 - subcluster time center, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEtimeCenter_phiE2Deq0PiOv2 = new TH1D("noEtimeCenter_phiE2Deq0PiOv2","noEtimeCenter_phiE2Deq0PiOv2",25,-20,20);
                //198 - subcluster time center, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEtimeCenter_phiE2Dneq0PiOv2 = new TH1D("noEtimeCenter_phiE2Dneq0PiOv2","noEtimeCenter_phiE2Dneq0PiOv2",25,-20,20);
                //199 - pho E (rh sum), phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEphoE_phiE2Deq0PiOv2 = new TH1D("noEphoE_phiE2Deq0PiOv2","noEphoE_phiE2Deq0PiOv2",25,0,1000);
                //200 - pho E (rh sum), phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEphoE_phiE2Dneq0PiOv2 = new TH1D("noEphoE_phiE2Dneq0PiOv2","noEphoE_phiE2Dneq0PiOv2",25,0,1000);
                //201 - nsubclusters, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEnSubClusters_phiE2Deq0PiOv2 = new TH1D("noEnSubClusters_phiE2Deq0PiOv2","noEnSubClusters_phiE2Deq0PiOv2",10,0,10);
                //202 - nsubclusters, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEnSubClusters_phiE2Dneq0PiOv2 = new TH1D("noEnSubClusters_phiE2Dneq0PiOv2","noEnSubClusters_phiE2Dneq0PiOv2",10,0,10);
                //203 - eta sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEetaSig_phiE2Deq0PiOv2 = new TH1D("noEetaSig_phiE2Deq0PiOv2","noEetaSig_phiE2Deq0PiOv2",25,0.01,0.09);
                //204 - eta sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEetaSig_phiE2Dneq0PiOv2 = new TH1D("noEetaSig_phiE2Dneq0PiOv2","noEetaSig_phiE2Dneq0PiOv2",25,0.01,0.09);
                //205 - phi sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEphiSig_phiE2Deq0PiOv2 = new TH1D("noEphiSig_phiE2Deq0PiOv2","noEphiSig_phiE2Deq0PiOv2",25,0.01,0.09);
                //206 - phi sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEphiSig_phiE2Dneq0PiOv2 = new TH1D("noEphiSig_phiE2Dneq0PiOv2","noEphiSig_phiE2Dneq0PiOv2",25,0.01,0.09);
                //207 - time sigma, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEtimeSig_phiE2Deq0PiOv2 = new TH1D("noEtimeSig_phiE2Deq0PiOv2","noEtimeSig_phiE2Deq0PiOv2",25,0.,10.);
                //208 - time sigma, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEtimeSig_phiE2Dneq0PiOv2 = new TH1D("noEtimeSig_phiE2Dneq0PiOv2","noEtimeSig_phiE2Dneq0PiOv2",25,0.,10.);
                //209 - etaPhi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEetaPhiCov_phiE2Deq0PiOv2 = new TH1D("noEetaPhiCov_phiE2Deq0PiOv2","noEetaPhiCov_phiE2Deq0PiOv2",25,-1,1);
                //210 - etaPhi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEetaPhiCov_phiE2Dneq0PiOv2 = new TH1D("noEetaPhiCov_phiE2Dneq0PiOv2","noEetaPhiCov_phiE2Dneq0PiOv2",25,-1,1);
                //211 - timeEta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEtimeEtaCov_phiE2Deq0PiOv2 = new TH1D("noEtimeEtaCov_phiE2Deq0PiOv2","noEtimeEtaCov_phiE2Deq0PiOv2",25,-1,1);
                //212 - timeEta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEtimeEtaCov_phiE2Dneq0PiOv2 = new TH1D("noEtimeEtaCov_phiE2Dneq0PiOv2","noEtimeEtaCov_phiE2Dneq0PiOv2",25,-1,1);
                //213 - timePhi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEtimePhiCov_phiE2Deq0PiOv2 = new TH1D("noEtimePhiCov_phiE2Deq0PiOv2","noEtimePhiCov_phiE2Deq0PiOv2",25,-1,1);
                //214 - timePhi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEtimePhiCov_phiE2Dneq0PiOv2 = new TH1D("noEtimePhiCov_phiE2Dneq0PiOv2","noEtimePhiCov_phiE2Dneq0PiOv2",25,-1,1);
                //215 - timeMaj cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEtimeMajCov_phiE2Deq0PiOv2 = new TH1D("noEtimeMajCov_phiE2Deq0PiOv2","noEtimeMajCov_phiE2Deq0PiOv2",25,-1,1);
                //216 - timeMaj cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEtimeMajCov_phiE2Dneq0PiOv2 = new TH1D("noEtimeMajCov_phiE2Dneq0PiOv2","noEtimeMajCov_phiE2Dneq0PiOv2",25,-1,1);
                //217 - timeMin cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEtimeMinCov_phiE2Deq0PiOv2 = new TH1D("noEtimeMinCov_phiE2Deq0PiOv2","noEtimeMinCov_phiE2Deq0PiOv2",25,-1,1);
                //218 - timeMin cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEtimeMinCov_phiE2Dneq0PiOv2 = new TH1D("noEtimeMinCov_phiE2Dneq0PiOv2","noEtimeMinCov_phiE2Dneq0PiOv2",25,-1,1);
                //219 - phiE2D, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noEphiE2D_phiE2Deq0PiOv2 = new TH1D("noEphiE2D_phiE2Deq0PiOv2","noEphiE2D_phiE2Deq0PiOv2",25,-3.1,3.1);
                //220 - phiE2D, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noEphiE2D_phiE2Dneq0PiOv2 = new TH1D("noEphiE2D_phiE2Dneq0PiOv2","noEphiE2D_phiE2Dneq0PiOv2",25,-3.1,3.1);
                //221 - rot2D cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH1D* noErot2D_phiE2Deq0PiOv2 = new TH1D("noErot2D_phiE2Deq0PiOv2","noErot2D_phiE2Deq0PiOv2",25,0.4,1.1);
                //222 - rot2D cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH1D* noErot2D_phiE2Dneq0PiOv2 = new TH1D("noErot2D_phiE2Dneq0PiOv2","noErot2D_phiE2Dneq0PiOv2",25,0.4,1.1);
		//223 - swissCross
		TH1D* swCross = new TH1D("swCross","swCross",25,0.9,1.1); 
		//224 - true smaj
		TH1D* trueSmaj = new TH1D("trueSmaj","trueSmaj",25,0,5);
		//225 - true smaj
		TH1D* trueSmin = new TH1D("trueSmin","trueSmin",25,0,2);
		//226 - true sietaieta
		TH1D* trueSiEtaiEta = new TH1D("trueSiEtaiEta","trueSiEtaiEta",25,0.01,0.09);
		//227 - true siphiiphi
		TH1D* trueSiPhiiPhi = new TH1D("trueSiPhiiPhi","trueSiPhiiPhi",25,0.01,0.09);
		//228 - nrhs
		TH1D* phoNrhs = new TH1D("phoNrhs","phoNrhs",25,0,100);		
		//229 - nrhs, 0 <= E < 200
		TH1D* phoNrhs_Ebin1 = new TH1D("phoNrhs_Ebin1","phoNrhs_Ebin1",25,0,100);		
		//230 - nrhs, 200 <= E < 400
		TH1D* phoNrhs_Ebin2 = new TH1D("phoNrhs_Ebin2","phoNrhs_Ebin2",25,0,100);		
		//231 - nrhs, 400 <= E < 600
		TH1D* phoNrhs_Ebin3 = new TH1D("phoNrhs_Ebin3","phoNrhs_Ebin3",25,0,100);		
		//232 - nrhs, 600 <= E < 1000
		TH1D* phoNrhs_Ebin4 = new TH1D("phoNrhs_Ebin4","phoNrhs_Ebin4",25,0,100);		
		//233 - swiss cross prime (swiss cross recreation from subcluster information)
		TH1D* swCrossPrime = new TH1D("swCrossPrime","swCrossPrime",25,-0.05,1.5);
		//234 - difference between eta mean from points and eta center from GMM
		TH1D* etaDiff = new TH1D("etaDiff","etaDiff",25,-0.02,0.02);
		//235 - dPhi bw phiCenter and MET	
		TH1D* dPhi_phiCenterMet = new TH1D("dPhi_phiCenterMet","dPhi_phiCenterMet",25,-3.5,3.5);
		//236 - dPhi bw phiCenter and MET, etaSig + phiSig cuts
		TH1D* dPhi_phiCenterMet_etaSigge0p3ANDphiSigle0p3 = new TH1D("dPhi_phiCenterMet_etaSigge0p3ANDphiSigle0p3","dPhi_phiCenterMet_etaSigge0p3ANDphiSigle0p3",25,-3.5,3.5);
		//237 - dR bw subcluster and closest matching track
		TH1D* dR_trackSubcl = new TH1D("dR_trackSubcl","dR_trackSubcl",50,0,5);	
		//238 - normalized dE bw subcluster and closest matching track	
		TH1D* dE_trackSubcl = new TH1D("dE_trackSubcl","dE_trackSubcl",25,-2,2);	
		//iso bkg == bkg for sig vs bkg MVA
		//239 - energy - iso bkg selection
		TH1D* E_IsoBkgSel = new TH1D("E_IsoBkgSel","E_IsoBkgSel",50,0,500);
		//240 - etaSig - iso bkg selection	
                TH1D* etaSig_IsoBkgSel = new TH1D("etaSig_IsoBkgSel","etaSig_IsoBkgSel",50,0.01,0.05);
		//241 - phiSig - iso bkg selection	
                TH1D* phiSig_IsoBkgSel = new TH1D("phiSig_IsoBkgSel","phiSig_IsoBkgSel",50,0.01,0.09);
		//242 - etaPhiCov - iso bkg selection	
		TH1D* etaPhiCov_IsoBkgSel = new TH1D("etaPhiCov_IsoBkgSel","etaPhiCov_IsoBkgSel",50,-0.01,0.01);
		//243 - timeEtaCov - iso bkg selection	
		TH1D* timeEtaCov_IsoBkgSel = new TH1D("timeEtaCov_IsoBkgSel","timeEtaCov_IsoBkgSel",50,-0.2,0.2);
		//244 - timePhiCov - iso bkg selection	
		TH1D* timePhiCov_IsoBkgSel = new TH1D("timePhiCov_IsoBkgSel","timePhiCov_IsoBkgSel",50,-0.2,0.1);
		//245 - majLength - iso bkg selection	
                TH1D* majLength_IsoBkgSel = new TH1D("majLength_IsoBkgSel","majLength_IsoBkgSel",50,0.01,1.);
		//246 - minLength - iso bkg selection	
                TH1D* minLength_IsoBkgSel = new TH1D("minLength_IsoBkgSel","minLength_IsoBkgSel",50,0.01,0.09);
		//247 - phi2D - iso bkg selection	
		TH1D* phi2D_IsoBkgSel = new TH1D("phiE2D_IsoBkgSel","phiE2D_IsoBkgSel",50,-2.5,1.);		
		//248 - rot2D - iso bkg selection	
		TH1D* rot2D_IsoBkgSel = new TH1D("rot2D_IsoBkgSel","rot2D_IsoBkgSel",50,0.4,1.1);
		//249 - rot3D - iso bkg selection 
		TH1D* rot3D_IsoBkgSel = new TH1D("rot3D_IsoBkgSel","rot3D_IsoBkgSel",50,0.95,1.01);
		//250 - # rhs above r_nk = 0.8
		TH1D* nRhs_rnkThresh = new TH1D("nRhs_rnkThresh","nRhs_rnkThresh",30,0,30);
		//251 - eta angle (angle bw maj axis + eta axis in 3D)
		TH1D* etaAngle3D = new TH1D("etaAngle3D","etaAngle3D",25,-0.4,3.4);
		//252 - phi angle (angle bw maj axis + phi axis in 3D)
		TH1D* phiAngle3D = new TH1D("phiAngle3D","phiAngle3D",25,-0.4,3.4);
		//253 - eta angle (angle bw maj axis + eta axis in 2D)
		TH1D* etaAngle2D = new TH1D("etaAngle2D","etaAngle2D",25,-0.4,3.4);
		//254 - phi angle (angle bw maj axis + phi axis in 2D)
		TH1D* phiAngle2D = new TH1D("phiAngle2D","phiAngle2D",25,-0.4,3.4);
		//255 - major axis length (in 3D)
		TH1D* majLength3D = new TH1D("majLength3D","majLength3D",50,0,5.);
		//256 - major axis length (in 2D)
		TH1D* majLength2D = new TH1D("majLength2D","majLength2D",50,0,0.2);
		//257 - sigma^2_t from meas err
		TH1D* timesSigSq_measErr = new TH1D("timesSigSq_measErr","timesSigSq_measErr",25,0,5);
		//258 - input rh times
		TH1D* rhTime = new TH1D("rhTime","rhTime",50,-3,3);
		//259 - photon ps
		TH1D* photonPt = new TH1D("photonPt","photonPt",25,30,1000);
		//260 - eta angle (angle bw maj axis + eta axis in 2D) with |time maj cov| > 0.1
		TH1D* etaAngle2D_absTimeMajCovge0p1 = new TH1D("etaAngle2D_absTimeMajCovge0p1","etaAngle2D_absTimeMajCovge0p1",25,-0.4,3.4);
		//261 - sin(etaAngle_2D - etaAngle_3D)
		TH1D* sinDiffEtaAngle2D_3D = new TH1D("sinDiffEtaAngle2D_3D","sinDiffEtaAngle2D_3D",25,-1.2,1.2);
		//262 - timePhi cov over timeMaj3D cov
		TH1D* timePhiCovOvtimeMaj2DCov = new TH1D("timePhiCovOvtimeMaj2DCov","timePhiCovOvtimeMaj2DCov",25,-2,2);
		//263 - timeEta cov over timeMaj3D cov
		TH1D* timeEtaCovOvtimeMaj2DCov = new TH1D("timeEtaCovOvtimeMaj2DCov","timeEtaCovOvtimeMaj2DCov",25,-2,2);
		//264 - time maj 2d / time maj 3d
		TH1D* timeMaj2DCovOvtimeMaj3DCov = new TH1D("timeMaj2DCovOvtimeMaj3DCov","timeMaj2DCovOvtimeMaj3DCov",25,-3,3);	
		


	
		//0 - time v subcl subcluster energy
		TH2D* time_E = new TH2D("time_subclE","time_subclE;time_center;E;a.u.", 50,-30,30,10,0,1000);
		//1 - azimuthal angle v subcl energy
		TH2D* az_E = new TH2D("az_subclE","az_subclE;azimuthal_angle;E;a.u.",50,-3.5,3.5,10,0,1000);
		//2 - rotundity (2D) v subcl energy
		TH2D* rot2D_E = new TH2D("rot2D_subclE","rot2D_subclE;rotundity2D;E;a.u.",25,0.4,1.1,25,0,1000);
		//3 - eta v phi
		TH2D* eta_phi = new TH2D("eta_phi","eta_phi;eta_center;phi_center",50,-2,2,50,-0.2,6.4);
		//4 - t v eta
		TH2D* t_eta = new TH2D("t_eta","t_eta;time_center;eta_center",50,-30,30,50,-3.5,3.5);
		//5 - t v phi
		TH2D* t_phi = new TH2D("t_phi","t_phi;time_center;phi_center;a.u.",50,-30,30,50,-0.1,6.3);
		//6 - time to mixing coeff
		TH2D* t_mixcoeff = new TH2D("t_mixcoeff","t_mixcoeff;time_center;mixing_coeff;a.u.",50,-20,20,50,0,1.1);
		//7 - subcl E vs mixing coeff
		TH2D* E_mixcoeff = new TH2D("subclE_mixcoeff","subclE_mixcoeff;E;mixing_coeff;a.u.",20,0,2000,20,0,1.1);	
		//8 - rotundity (3D) v subcluster energy
		TH2D* rot3D_E = new TH2D("rot3D_subclE","rot3D_subclE;rotundity3D;E;a.u.",50,0.4,1.1,50,0,100);
		//9 - npts effective (ie weighted) v subcl energy
		TH2D* npts_E = new TH2D("npts_E","npts_subclE;npts;E;a.u.",50,0.,100,10,0,1000);
		//10 - nsubclusters vs nrhs in cluster
		TH2D* nsubcl_nrhs = new TH2D("nsubcl_nrhs","nsubcl_nrhs;nsubclusters;nrhs",10,0,10,10,0,50);
		//11 - nsubclusters vs mm coeff
		TH2D* nsubcl_mmcoeff = new TH2D("nsubcl_mmcoeff","nsubcl_mmcoeff;nsubclusters;mmcoeff",10,0,10,20,0.,1.1);
                //12 - eta sigma v phi sigma
		TH2D* etaSig_phiSig = new TH2D("etaSig_phiSig","etaSig_phiSig;etaSig;phiSig",25,0.,0.1,25,0.,0.1);
                //13 - time sigma v phi sigma
                TH2D* timeSig_etaSig = new TH2D("timeSig_etaSig","timeSig_etaSig;timeSig;etaSig",25,0,5.,25,0.,0.1);
                //14 - time sigma v phi sigma
                TH2D* timeSig_phiSig = new TH2D("timeSig_phiSig","timeSig_phiSig;timeSig;phiSig",25,0,5.,25,0.,0.1);
		//15 - fraction of energy in subcluster vs mm coeff of subcluster
		TH2D* fracE_mmcoeff = new TH2D("fracE_mmcoeff","fracE_mmcoeff;fracE;mmcoeff",20,0.,1.,20,0,1.1);
		//16 - number of subclusters vs fraction of energy in particular subcluster (really only applicable to lead subcluster)
		TH2D* nsubcl_fracE = new TH2D("nsubcl_fracE","nsubcl_fracE;nSubClusters;fracE",10,0,10.,20,0,1.1);
		//17 - time sigma vs time center
		TH2D* timeCenter_timeSig = new TH2D("timeCenter_timeSig","timeCenter_timeSig;timeCenter;time sig",25,-3,3,25,0.,5);
		//18 - 2d rotundity vs 2d az angle
		TH2D* rot2D_az2D = new TH2D("rot2D_az2D","rot2D_az2D;rotundity2D;azangle2D",50,0.4,1.1,50,0.,3.5);
		//19 - time variation vs frac E
		TH2D* timeSig_fracE = new TH2D("timeSig_fracE","timeSig_fracE;timeSig;fracE",25,0,25,20,0.,1.1);
		//20 - time sigma vs total E
		TH2D* timeSig_totE = new TH2D("timeSig_totE","timeSig_totE;timeSig;totE",30,0,5,30,0,1400);
		//21 - time-eta covariance vs total E
		TH2D* timeEtaCov_totE = new TH2D("timeEtaCov_totE","timeEtaCov_totE;timeEtaCov;totE",25,-1,1,20,0,2000);
		//22 - time-phi covariance vs total E
		TH2D* timePhiCov_totE = new TH2D("timePhiCov_totE","timePhiCov_totE;timePhiCov;totE",25,-1,1,20,0,2000);
                //23 - time sigma vs. TimeMajCov
                TH2D* timeSig_timeMajCov = new TH2D("timeSig_timeMajCov","timeSig_timeMajCov;timeSig;timeMajCov",25,0,1.,25,-0.5,0.5);
                //24 - time sigma vs. TimeMinCov
                TH2D* timeSig_timeMinCov = new TH2D("timeSig_timeMinCov","timeSig_timeMinCov;timeSig;timeMinCov",25,0,1.,25,-0.5,0.5);
                //25 - az angle 2D vs. TimeMajCov
                TH2D* phiEll2D_timeMajCov = new TH2D("phiE2D_timeMajCov","phiE2D_timeMajCov;phiE2D;timeMajCov",25,0,3.5,25,-0.5,0.5);
                //26 - az 2D angle vs. TimeMinCov
                TH2D* phiEll2D_timeMinCov = new TH2D("phiE2D_timeMinCov","phiE2D_timeMinCov;phiE2D;timeMinCov",25,0,3.5,25,-0.5,0.5);
                //27 - rot 2D vs. TimeMajCov
                TH2D* rot2D_timeMajCov = new TH2D("rot2D_timeMajCov","rot2D_timeMajCov;rot2D;timeMajCov",25,0.45,1.1,25,-0.5,0.5);
                //28 - rot 2D angle vs. TimeMinCov
                TH2D* rot2D_timeMinCov = new TH2D("rot2D_timeMinCov","rot2D_timeMinCov;rot2D;timeMinCov",25,0.45,1.1,25,-0.5,0.5);
                //29 - rot 3D vs. TimeMajCov
                TH2D* rot3D_timeMajCov = new TH2D("rot3D_timeMajCov","rot3D_timeMajCov;rot3D;timeMajCov",25,0.8,1.1,25,-0.5,0.5);
                //30 - rot 3D angle vs. TimeMinCov
                TH2D* rot3D_timeMinCov = new TH2D("rot3D_timeMinCov","rot3D_timeMinCov;rot3D;timeMinCov",25,0.8,1.1,25,-0.5,0.5);
		//31 - timemaj cov vs time center
		TH2D* timeCenter_timeMajCov = new TH2D("timeCenter_timeMajCov","timeCenter_timeMajCov;timeCenter;timemaj cov",25,-3,3,25,-0.3,0.3);
		//32 - timemin cov vs time center
		TH2D* timeCenter_timeMinCov = new TH2D("timeCenter_timeMinCov","timeCenter_timeMinCov;timeCenter;timemin cov",25,-3,3,25,-0.3,0.3);
		//33 - eta-phi covariance vs azimuth ellipsoid angle
		TH2D* etaPhiCov_phiEll2D = new TH2D("etaPhiCov_phiE2D","etaPhiCov_phiE2D;etaPhiCov;ellipsoid phi 2D",25,-0.2,0.2,25,0,3.5);
		//34 - time-maj covariance vs azimuth ellipsoid angle
		TH2D* timeMajCov_phiEll2D = new TH2D("timeMajCov_phiE2D","timeMajCov_phiE2D;timeMajCov;ellipsoid phi 2D",25,-0.5,0.5,25,0,3.5);
		//35 - time-min covariance vs azimuth ellipsoid angle
		TH2D* timeMinCov_phiEll2D = new TH2D("timeMinCov_phiE2D","timeMinCov_phiE2D;timeMinCov;ellipsoid phi 2D",25,-0.5,0.5,25,0,3.5);
                //36 - rot 2D vs. TimeMajCov unnorm
                TH2D* rot2D_timeMajCovUnnorm = new TH2D("rot2D_timeMajCovUnnorm","rot2D_timeMajCovUnnorm;rot2D;timeMajCovUnnorm",25,0.4,1.1,25,-0.5,0.5);
                //37 - rot 2D angle vs. TimeMinCov unnorm
                TH2D* rot2D_timeMinCovUnnorm = new TH2D("rot2D_timeMinCovUnnorm","rot2D_timeMinCovUnnorm;rot2D;timeMinCovUnnorm",25,0.4,1.1,25,-0.5,0.5);
                //38 - time sigma vs. timeEtaCov
                TH2D* timeSig_timeEtaCov = new TH2D("timeSig_timeEtaCov","timeSig_timeEtaCov;timeSig;timeEtaCov",25,0,1.,25,-1,1);
		//39 - etaphi cov vs time center
		TH2D* timeCenter_etaPhiCov = new TH2D("timeCenter_etaPhiCov","timeCenter_etaPhiCov;timeCenter;etaphi cov",25,-3,3,25,-0.1,0.1);
		//40 - timeeta cov vs time center
		TH2D* timeCenter_timeEtaCov = new TH2D("timeCenter_timeEtaCov","timeCenter_timeEtaCov;timeCenter;timeeta cov",25,-3,3,25,-0.1,0.1);
		//41 - timeeta cov vs time center
		TH2D* timeCenter_timePhiCov = new TH2D("timeCenter_timePhiCov","timeCenter_timePhiCov;timeCenter;time center",25,-3,3,25,-0.1,0.1);
		//42 - CMS smaj vs cms sigma_t
		TH2D* cmsSmaj_cmsTimeSig = new TH2D("noESmaj_noETimeSig","noESmaj_noETimeSig;noE_Smaj;noE_timeSig",25,0.,0.01,25,0,10);
		//43 - CMS smin vs cms sigma_t
		TH2D* cmsSmin_cmsTimeSig = new TH2D("noESmin_noETimeSig","noESmin_noETimeSig;noE_Smin;noE_timeSig",25,0.,0.01,25,0,10);
		//44 - CMS logE smaj vs cms sigma_t
		TH2D* logESmaj_logETimeSig = new TH2D("logESmaj_logETimeSig","logESmaj_logETimeSig;logE_Smaj;logE_timeSig",25,0.,0.01,25,0,10);
		//45 - CMS logE smin vs cms sigma_t
		TH2D* logESmin_logETimeSig = new TH2D("logESmin_logETimeSig","logESmin_logETimeSig;logE_Smin;logE_timeSig",25,0.,0.01,25,0,10);
		//46 - etaphi cov vs timeeta cov
		TH2D* etaPhiCov_timeEtaCov = new TH2D("etaPhiCov_timeEtaCov","etaPhiCov_timeEtaCov;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//47 - timeeta cov vs timephi cov
		TH2D* timeEtaCov_timePhiCov = new TH2D("timeEtaCov_timePhiCov","timeEtaCov_timePhiCov;timeEtaCov;timePhiCov",25,-0.1,0.1,25,-0.1,0.1);
		//48 - etaphi cov vs timephi cov
		TH2D* etaPhiCov_timePhiCov = new TH2D("etaPhiCov_timePhiCov","etaPhiCov_timePhiCov;etaPhiCov;timePhiCov",25,-1,1,25,-1,1);
                //49 - etaphi cov vs. TimeMajCov
                TH2D* etaPhiCov_timeMajCov = new TH2D("etaPhiCov_timeMajCov","etaPhiCov_timeMajCov;etaPhiCov;timeMajCov",25,-0.1,0.1,25,-0.3,0.3);
                //50 - etaphi cov vs. TimeMinCov
                TH2D* etaPhiCov_timeMinCov = new TH2D("etaPhiCov_timeMinCov","etaPhiCov_timeMinCov;etaPhiCov;timeMinCov",25,-0.1,0.1,25,-0.3,0.3);
                //51 - timephi cov vs. TimeMajCov
                TH2D* timePhiCov_timeMajCov = new TH2D("timePhiCov_timeMajCov","timePhiCov_timeMajCov;timePhiCov;timeMajCov",25,-0.1,0.1,25,-0.3,0.3);
                //52 - timephi cov vs. TimeMinCov
                TH2D* timePhiCov_timeMinCov = new TH2D("timePhiCov_timeMinCov","timePhiCov_timeMinCov;timePhiCov;timeMinCov",25,-0.1,0.1,25,-0.3,0.3);
                //53 - timeeta cov vs. TimeMajCov
                TH2D* timeEtaCov_timeMajCov = new TH2D("timeEtaCov_timeMajCov","timeEtaCov_timeMajCov;timeEtaCov;timeMajCov",25,-0.02,0.02,25,-0.3,0.3);
                //54 - timeeta cov vs. TimeMinCov
                TH2D* timeEtaCov_timeMinCov = new TH2D("timeEtaCov_timeMinCov","timeEtaCov_timeMinCov;timeEtaCov;timeMinCov",25,-0.1,0.1,25,-0.3,0.3);
		//55 - etaphi cov unnorm vs timeeta cov unnorm
		TH2D* etaPhiCovUnnorm_timeEtaCovUnnorm = new TH2D("etaPhiCovUnnorm_timeEtaCovUnnorm","etaPhiCovUnnorm_timeEtaCovUnnorm;etaPhiCovUnnorm;timeEtaCovUnnorm",25,-0.2,0.2,25,-0.5,0.5);
		//56 - timeeta cov unnorm vs timephi cov unnorm
		TH2D* timeEtaCovUnnorm_timePhiCovUnnorm = new TH2D("timeEtaCovUnnorm_timePhiCovUnnorm","timeEtaCovUnnorm_timePhiCovUnnorm;timeEtaCovUnnorm;timePhiCovUnnorm",25,-0.5,0.5,25,-0.5,0.5);
		//57 - etaphi cov unnorm vs timephi cov unnorm
		TH2D* etaPhiCovUnnorm_timePhiCovUnnorm = new TH2D("etaPhiCovUnnorm_timePhiCovUnnorm","etaPhiCovUnnorm_timePhiCovUnnorm;etaPhiCovUnnorm;timePhiCovUnnorm",25,-0.2,0.2,25,-0.5,0.5);
                //58 - etaphi cov unnorm vs. TimeMajCovUnnorm 
                TH2D* etaPhiCovUnnorm_timeMajCovUnnorm = new TH2D("etaPhiCovUnnorm_timeMajCovUnnorm","etaPhiCovUnnorm_timeMajCovUnnorm;etaPhiCovUnnorm;timeMajCovUnnorm",25,-0.2,0.2,25,-0.5,0.5);
                //59 - etaphi cov unnorm vs. TimeMinCovUnnorm
                TH2D* etaPhiCovUnnorm_timeMinCovUnnorm = new TH2D("etaPhiCovUnnorm_timeMinCovUnnorm","etaPhiCovUnnorm_timeMinCovUnnorm;etaPhiCovUnnorm;timeMinCovUnnorm",25,-0.2,0.2,25,-0.5,0.5);
                //60 - timephi cov unnorm vs. TimeMajCovUnnorm
                TH2D* timePhiCovUnnorm_timeMajCovUnnorm = new TH2D("timePhiCovUnnorm_timeMajCovUnnorm","timePhiCovUnnorm_timeMajCovUnnorm;timePhiCovUnnorm;timeMajCovUnnorm",25,-0.5,0.5,25,-0.5,0.5);
                //61 - timephi cov unnorm vs. TimeMinCovUnnorm
                TH2D* timePhiCovUnnorm_timeMinCovUnnorm = new TH2D("timePhiCovUnnorm_timeMinCovUnnorm","timePhiCovUnnorm_timeMinCovUnnorm;timePhiCovUnnorm;timeMinCovUnnorm",25,-0.5,0.5,25,-0.5,0.5);
                //62 - timeeta cov unnorm vs. TimeMajCovUnnorm
                TH2D* timeEtaCovUnnorm_timeMajCovUnnorm = new TH2D("timeEtaCovUnnorm_timeMajCovUnnorm","timeEtaCovUnnorm_timeMajCovUnnorm;timeEtaCovUnnorm;timeMajCovUnnorm",25,-0.5,0.5,25,-0.5,0.5);
                //63 - timeeta cov unnorm vs. TimeMinCovUnnorm
                TH2D* timeEtaCovUnnorm_timeMinCovUnnorm = new TH2D("timeEtaCovUnnorm_timeMinCovUnnorm","timeEtaCovUnnorm_timeMinCovUnnorm;timeEtaCovUnnorm;timeMinCovUnnorm",25,-0.5,0.5,25,-0.5,0.5);
                //64 - rot 2D vs. etaphi cov
                TH2D* rot2D_etaPhiCov = new TH2D("rot2D_etaPhiCov","rot2D_etaPhiCov;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
                //65 - rot 2D vs. etaphi cov unnorm
                TH2D* rot2D_etaPhiCovUnnorm = new TH2D("rot2D_etaPhiCovUnnorm","rot2D_etaPhiCovUnnorm;rot2D;etaPhiCovUnnorm",25,0.4,1.1,25,-0.2,0.2);
		//66 - phi center vs. energy
		TH2D* E_phi = new TH2D("E_phi","E_phi;E_center;phi_center",50,0.,1000,50,-0.2,6.4);
		//67 - rot2D vs phiE2D
                TH2D* rot2D_phiEll2D = new TH2D("rot2D_phiE2D","rot2D_phiE2D;rot2D;phiE2D",25,0.4,1.1,25,-3.1,1);
		//68 - counts of etaphi cov vs timeeta cov
		TH2D* etaPhiCov_timeEtaCovCounts = new TH2D("etaPhiCov_timeEtaCovCounts","etaPhiCov_timeEtaCovCounts;etaPhiCovCounts;timeEtaCovCounts",2,-1,1,2,-1,1);
		//69 (nice) - time center vs phiE2D
		TH2D* timeCenter_phiE2D = new TH2D("timeCenter_phiE2D","timeCenter_phiE2D;timeCenter;phiE2D",25,-3,3,25,-3.1,1);	
                //70 - rot 2D vs. etaphi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* rot2D_etaPhiCov_phiE2Deq0PiOv2 = new TH2D("rot2D_etaPhiCov_phiE2Deq0PiOv2","rot2D_etaPhiCov_phiE2Deq0PiOv2;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
                //71 - rot 2D vs. etaphi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* rot2D_etaPhiCov_phiE2Dneq0PiOv2 = new TH2D("rot2D_etaPhiCov_phiE2Dneq0PiOv2","rot2D_etaPhiCov;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
		//72 - energy vs phiE2D
		TH2D* phoE_phiE2D = new TH2D("phoEnergy_phiE2D","phoEnergy_phiE2D;E_{pho};phiE2D",25,0,1000,25,-3.1,1);	
		//73 - timeEtaCov vs rotundity 2D
		TH2D* timeEtaCov_rot2D = new TH2D("timeEtaCov_rot2D","timeEtaCov_rot2D;timeEtaCov;rot2D",25,-1.,1.,25,0.4,1.1);
                //74 - rot 2D vs. timeeta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* rot2D_timeEtaCov_phiE2Deq0PiOv2 = new TH2D("rot2D_timeEtaCov_phiE2Deq0PiOv2","rot2D_timeEtaCov_phiE2Deq0PiOv2;rot2D;timeEtaCov",25,0.4,1.1,25,-0.2,0.2);
                //75 - rot 2D vs. timeeta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* rot2D_timeEtaCov_phiE2Dneq0PiOv2 = new TH2D("rot2D_timeEtaCov_phiE2Dneq0PiOv2","rot2D_timeEtaCov;rot2D;timeEtaCov",25,0.4,1.1,25,-0.2,0.2);
		

		//////////logE weighted//////////
                //76 - logE rot 2D vs. etaphi cov
                TH2D* logErot2D_etaPhiCov = new TH2D("logErot2D_etaPhiCov","logErot2D_etaPhiCov;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
                //77 - logE rot 2D vs. etaphi cov unnorm
                TH2D* logErot2D_etaPhiCovUnnorm = new TH2D("logErot2D_etaPhiCovUnnorm","logErot2D_etaPhiCovUnnorm;rot2D;etaPhiCovUnnorm",25,0.4,1.1,25,-0.2,0.2);
		//78 - logE phi center vs. energy
		TH2D* logEE_phi = new TH2D("logEE_phi","logEE_phi;E_center;phi_center",50,0.,1000,50,-0.2,6.4);
		//79 - logE rot2D vs phiE2D
                TH2D* logErot2D_phiEll2D = new TH2D("logErot2D_phiE2D","logErot2D_phiE2D;rot2D;phiE2D",25,0.4,1.1,25,-3.1,1);
		//80 - logE counts of etaphi cov vs timeeta cov
		TH2D* logEetaPhiCov_timeEtaCovCounts = new TH2D("logEetaPhiCov_timeEtaCovCounts","logEetaPhiCov_timeEtaCovCounts;etaPhiCovCounts;timeEtaCovCounts",2,-1,1,2,-1,1);
		//81 - logE time center vs phiE2D
		TH2D* logEtimeCenter_phiE2D = new TH2D("logEtimeCenter_phiE2D","logEtimeCenter_phiE2D;timeCenter;phiE2D",25,-15,15,25,-3.1,1);	
                //82 - logE rot 2D vs. etaphi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* logErot2D_etaPhiCov_phiE2Deq0PiOv2 = new TH2D("logErot2D_etaPhiCov_phiE2Deq0PiOv2","logErot2D_etaPhiCov_phiE2Deq0PiOv2;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
                //83 - logE rot 2D vs. etaphi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* logErot2D_etaPhiCov_phiE2Dneq0PiOv2 = new TH2D("logErot2D_etaPhiCov_phiE2Dneq0PiOv2","logErot2D_etaPhiCov;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
		//84 - logE energy vs phiE2D
		TH2D* logEphoE_phiE2D = new TH2D("logEphoEnergy_phiE2D","logEphoEnergy_phiE2D;E_{pho};phiE2D",25,0,1000,25,-3.1,1);	
		//85 - logE timeEtaCov vs rotundity 2D
		TH2D* logEtimeEtaCov_rot2D = new TH2D("logEtimeEtaCov_rot2D","logEtimeEtaCov_rot2D;timeEtaCov;rot2D",25,-1.,1.,25,0.4,1.1);
                //86 - logE rot 2D vs. timeeta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* logErot2D_timeEtaCov_phiE2Deq0PiOv2 = new TH2D("logErot2D_timeEtaCov_phiE2Deq0PiOv2","logErot2D_timeEtaCov_phiE2Deq0PiOv2;rot2D;timeEtaCov",25,0.4,1.1,25,-0.2,0.2);
                //87 - logE rot 2D vs. timeeta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* logErot2D_timeEtaCov_phiE2Dneq0PiOv2 = new TH2D("logErot2D_timeEtaCov_phiE2Dneq0PiOv2","logErot2D_timeEtaCov;rot2D;timeEtaCov",25,0.4,1.1,25,-0.2,0.2);

		

		//////////noE weighted//////////
                //88 - noE rot 2D vs. etaphi cov
                TH2D* noErot2D_etaPhiCov = new TH2D("noErot2D_etaPhiCov","noErot2D_etaPhiCov;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
                //89 - noE rot 2D vs. etaphi cov unnorm
                TH2D* noErot2D_etaPhiCovUnnorm = new TH2D("noErot2D_etaPhiCovUnnorm","noErot2D_etaPhiCovUnnorm;rot2D;etaPhiCovUnnorm",25,0.4,1.1,25,-0.2,0.2);
		//90 - noE phi center vs. energy
		TH2D* noEE_phi = new TH2D("noEE_phi","noEE_phi;E_center;phi_center",50,0.,1000,50,-0.2,6.4);
		//91 - noE rot2D vs phiE2D
                TH2D* noErot2D_phiEll2D = new TH2D("noErot2D_phiE2D","noErot2D_phiE2D;rot2D;phiE2D",25,0.4,1.1,25,-3.1,1);
		//92 - noE counts of etaphi cov vs timeeta cov
		TH2D* noEetaPhiCov_timeEtaCovCounts = new TH2D("noEetaPhiCov_timeEtaCovCounts","noEetaPhiCov_timeEtaCovCounts;etaPhiCovCounts;timeEtaCovCounts",2,-1,1,2,-1,1);
		//93 - noE time center vs phiE2D
		TH2D* noEtimeCenter_phiE2D = new TH2D("noEtimeCenter_phiE2D","noEtimeCenter_phiE2D;timeCenter;phiE2D",25,-15,15,25,-3.1,1);	
                //94 - noE rot 2D vs. etaphi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* noErot2D_etaPhiCov_phiE2Deq0PiOv2 = new TH2D("noErot2D_etaPhiCov_phiE2Deq0PiOv2","noErot2D_etaPhiCov_phiE2Deq0PiOv2;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
                //95 - noE rot 2D vs. etaphi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* noErot2D_etaPhiCov_phiE2Dneq0PiOv2 = new TH2D("noErot2D_etaPhiCov_phiE2Dneq0PiOv2","noErot2D_etaPhiCov;rot2D;etaPhiCov",25,0.4,1.1,25,-0.2,0.2);
		//96 - noE energy vs phiE2D
		TH2D* noEphoE_phiE2D = new TH2D("noEphoEnergy_phiE2D","noEphoEnergy_phiE2D;E_{pho};phiE2D",25,0,1000,25,-3.1,1);	
		//97 - noE timeEtaCov vs rotundity 2D
		TH2D* noEtimeEtaCov_rot2D = new TH2D("noEtimeEtaCov_rot2D","noEtimeEtaCov_rot2D;timeEtaCov;rot2D",25,-1.,1.,25,0.4,1.1);
                //98 - noE rot 2D vs. timeeta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* noErot2D_timeEtaCov_phiE2Deq0PiOv2 = new TH2D("noErot2D_timeEtaCov_phiE2Deq0PiOv2","noErot2D_timeEtaCov_phiE2Deq0PiOv2;rot2D;timeEtaCov",25,0.4,1.1,25,-0.2,0.2);
                //99 - noE rot 2D vs. timeeta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* noErot2D_timeEtaCov_phiE2Dneq0PiOv2 = new TH2D("noErot2D_timeEtaCov_phiE2Dneq0PiOv2","noErot2D_timeEtaCov;rot2D;timeEtaCov",25,0.4,1.1,25,-0.2,0.2);

		//linear E-weighting
		//100 - etaphi cov vs timeeta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* etaPhiCov_timeEtaCov_phiE2Deq0PiOv2 = new TH2D("etaPhiCov_timeEtaCov_phiE2Deq0PiOv2","etaPhiCov_timeEtaCov_phiE2Deq0PiOv2;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
                //101 - etaphi cov vs. timeeta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* etaPhiCov_timeEtaCov_phiE2Dneq0PiOv2 = new TH2D("etaPhiCov_timeEtaCov_phiE2Dneq0PiOv2","etaPhiCov_timeEtaCov_phiE2Dneq0PiOv2;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//102 - phoE vs phiE2D, -10 < t < -2
		TH2D* phoE_phiE2D_timeNeg10ToNeg2 = new TH2D("phoEnergy_phiE2D_timeNeg10ToNeg2","phoEnergy_phiE2D_timeNeg10ToNeg2;E_{pho};phiE2D",25,0,1000,25,-3.1,1);	
		//103 - phoE vs phiE2D, -2 < t < 5
		TH2D* phoE_phiE2D_timeNeg2To5 = new TH2D("phoEnergy_phiE2D_timeNeg2To5","phoEnergy_phiE2D_timeNeg2To5;E_{pho};phiE2D",25,0,1000,25,-3.1,1);
		//104 - phoE vs phiE2D, 5 < t < 10
		TH2D* phoE_phiE2D_time5To10 = new TH2D("phoEnergy_phiE2D_time5To10","phoEnergy_phiE2D_time5To10;E_{pho};phiE2D",25,0,1000,25,-3.1,1);
		//105 - phoE vs phiE2D, 10 < t < 15
		TH2D* phoE_phiE2D_time10To15 = new TH2D("phoEnergy_phiE2D_time10To15","phoEnergy_phiE2D_time10To15;E_{pho};phiE2D",25,0,1000,25,-3.1,1);
		//106 - eta center vs time center
		TH2D* timeCenter_etaCenter = new TH2D("timeCenter_etaCenter","timeCenter_etaCenter;timeCenter;etaCenter",25,-15,15,25,-1.6,1.6);
		//107 - eta center vs time center, phiE2D ~ 0 && phiE2D ~ pi/2 
		TH2D* timeCenter_etaCenter_phiE2Deq0PiOv2 = new TH2D("timeCenter_etaCenter_phiE2Deq0PiOv2","timeCenter_etaCenter_phiE2Deq0PiOv2;timeCenter_phiE2Deq0PiOv2;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//108 - eta center vs time center, phiE2D !~ 0 && phiE2D !~ pi/2
		TH2D* timeCenter_etaCenter_phiE2Dneq0PiOv2 = new TH2D("timeCenter_etaCenter_phiE2Dneq0PiOv2","timeCenter_etaCenter_phiE2Dneq0PiOv2;timeCenter_phiE2Dneq0PiOv2;etaCenter;a.u.",25,-15,15,25,-1.6, 1.6);

		//109 - no E etaPhiCov vs timeEtaCov
                TH2D* logEetaPhiCov_timeEtaCov = new TH2D("logEetaPhiCov_timeEtaCov","logEetaPhiCov_timeEtaCov;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//110 - log E etaPhiCov vs timeEtaCov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* logEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2 = new TH2D("logEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2","logEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//111 - log E etaPhiCov vs timeEtaCov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* logEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2 = new TH2D("logEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2","logEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);

		//112 - no E etaPhiCov vs timeEtaCov
                TH2D* noEetaPhiCov_timeEtaCov = new TH2D("noEetaPhiCov_timeEtaCov","noEetaPhiCov_timeEtaCov;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//113 - no E etaPhiCov vs timeEtaCov, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* noEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2 = new TH2D("noEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2","noEetaPhiCov_timeEtaCov_phiE2Deq0PiOv2;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//114 - no E etaPhiCov vs timeEtaCov, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* noEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2 = new TH2D("noEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2","noEetaPhiCov_timeEtaCov_phiE2Dneq0PiOv2;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);

		//115 - time center vs rot2D
		TH2D* timeCenter_rot2D = new TH2D("timeCenter_rot2D","timeCenter_rot2D;timeCenter;rot2D",25,-3,3,25,0.4,1.1);
		//116 - time center vs swCross
		TH2D* timeCenter_swCross = new TH2D("timeCenter_swCross","timeCenter_swCross;timeCenter;swCross",25,-15,15,25,0.9,1.1);
		//117 - time center vs swCross, phiE2D ~ 0 && phiE2D ~ pi/2
		TH2D* timeCenter_swCross_phiE2Deq0PiOv2 = new TH2D("timeCenter_swCross_phiE2Deq0PiOv2","timeCenter_swCross;timeCenter;swCross",25,-15,15,25,0.9,1.1);
		//118 - time center vs swCross, phiE2D !~ 0 && phiE2D !~ pi/2
		TH2D* timeCenter_swCross_phiE2Dneq0PiOv2 = new TH2D("timeCenter_swCross_phiE2Dneq0PiOv2","timeCenter_swCross;timeCenter;swCross",25,-15,15,25,0.9,1.1);
		//119 - rot2D vs swCross
                TH2D* rot2D_swCross = new TH2D("rot2D_swCross","rot2D_swCross;rot2D;swCross",25,0.4,1.1,25,0.9,1.1);
		//120 - rot2D vs swCross, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* rot2D_swCross_phiE2Deq0PiOv2 = new TH2D("rot2D_swCross_phiE2Deq0PiOv2","rot2D_swCross_phiE2Deq0PiOv2;rot2D;swCross",25,0.4,1.1,25,0.9,1.1);
		//121 - rot2D vs swCross, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* rot2D_swCross_phiE2Dneq0PiOv2 = new TH2D("rot2D_swCross_phiE2Dneq0PiOv2","rot2D_swCross_phiE2Dneq0PiOv2;rot2D;swCross",25,0.4,1.1,25,0.9,1.1);
		//122 - phiE2D vs swCross
                TH2D* phiE2D_swCross = new TH2D("phiE2D_swCross","phiE2D_swCross;phiE2D;swCross",25,-3.,3.,25,0.9,1.1);
		//123 - phoE vs swCross
		TH2D* phoE_swCross = new TH2D("phoEnergy_swCross","phoEnergy_swCross;phoE;swCross",25,0,1000,25,0.9,1.1);	
		//124 - rot2D vs phiE2D, -10 < t < -2
                TH2D* rot2D_phiEll2D_timeNeg10toNeg2 = new TH2D("rot2D_phiE2D_timeNeg10toNeg2","rot2D_phiE2D_timeNeg10toNeg2;rot2D;phiE2D_timeNeg10toNeg2",25,0.4,1.1,25,-3.1,1);
		//125 - rot2D vs phiE2D, -2 < t < 5
                TH2D* rot2D_phiEll2D_timeNeg2to5 = new TH2D("rot2D_phiE2D_timeNeg2to5","rot2D_phiE2D_timeNeg2to5;rot2D;phiE2D_timeNeg2to5",25,0.4,1.1,25,-3.1,1);
		//126 - rot2D vs phiE2D, 5 < t < 15
                TH2D* rot2D_phiEll2D_time5to15 = new TH2D("rot2D_phiE2D_time5to15","rot2D_phiE2D_time5to15;rot2D;phiE2D_time5to15",25,0.4,1.1,25,-3.1,1);


		//127 - phoE vs logE smaj
		TH2D* phoE_logEsmaj = new TH2D("phoE_logEsmaj","phoE_logEsmaj;phoE;logEsmaj;a.u.",25,0,1000,25,0,0.004);
		//128 - phoE vs logE smin
		TH2D* phoE_logEsmin = new TH2D("phoE_logEsmin","phoE_logEsmin;phoE;logEsmin;a.u.",25,0,1000,25,0,0.004);
		//129 - phoE vs logE sieie
		TH2D* phoE_logEetaSig = new TH2D("phoE_logEetaSig","phoE_logEetaSig;phoE;logEetaSig",25,0,1000,25,0.01,0.09);
		//130 - phoE vs logE sipip
		TH2D* phoE_logEphiSig = new TH2D("phoE_logEphiSig","phoE_logEphiSig;phoE;logEphiSig",25,0,1000,25,0.01,0.09);
		//131 - phoE vs logE etaphicov
		TH2D* phoE_logEetaPhiCov = new TH2D("phoE_logEetaPhiCov","phoE_logEetaPhiCov",25,0,1000,25,-1,1);
		//132 - phoE vs logE timeetacov
		TH2D* phoE_logEtimeEtaCov = new TH2D("phoE_logEtimeEtaCov","phoE_logEtimeEtaCov",25,0,1000,25,-1,1);
		//133 - nRhs vs logE smaj
		TH2D* nRhs_logEsmaj = new TH2D("nRhs_logEsmaj","nRhs_logEsmaj;nRhs;logEsmaj;a.u.",25,0,100,25,0,0.004);
		//134 - nRhs vs logE smin
		TH2D* nRhs_logEsmin = new TH2D("nRhs_logEsmin","nRhs_logEsmin;nRhs;logEsmin;a.u.",25,0,100,25,0,0.004);
		//135 - nRhs vs logE sieie
		TH2D* nRhs_logEetaSig = new TH2D("nRhs_logEetaSig","nRhs_logEetaSig;nRhs;logEetaSig;a.u.",25,0,100,25,0.01,0.09);
		//136 - nRhs vs logE sipip
		TH2D* nRhs_logEphiSig = new TH2D("nRhs_logEphiSig","nRhs_logEphiSig;nRhs;logEphiSig;a.u.",25,0,100,25,0.01,0.09);
		//137 - phoE vs logE etaphicov
		TH2D* nRhs_logEetaPhiCov = new TH2D("nRhs_logEetaPhiCov","nRhs_logEetaPhiCov;nRhs;etaPhiCov;a.u.",25,0,100,25,-1,1);
		//138 - phoE vs logE timeetacov
		TH2D* nRhs_logEtimeEtaCov = new TH2D("nRhs_logEtimeEtaCov","nRhs_logEtimeEtaCov;nRhs;timeEtaCov;a.u.",25,0,100,25,-1,1);
		//139 - nrhs vs phoE
		TH2D* phoE_nRhs = new TH2D("phoE_nRhs","phoE_nRhs;phoEnergy;nRhs;a.u.",25,0,1000,25,0,100);		
		
		//BEAM HALO CR PLOTS
		//140 - time center vs eta center, 0.6 < rot2D < 0.8 
		TH2D* timeCenter_etaCenter_rot2Dge0p6le0p8 = new TH2D("timeCenter_etaCenter_rot2Dge0p6le0p8","timeCenter_etaCenter_rot2Dge0p6le0p8;timeCenter_rot2Dge0p6le0p8;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//141 - time center vs eta center, 0.6 > rotE2D && rot2D > 0.8 
		TH2D* timeCenter_etaCenter_rot2Dle0p6ge0p8 = new TH2D("timeCenter_etaCenter_rot2Dle0p6ge0p8","timeCenter_etaCenter_rot2Dle0p6ge0p8;timeCenter_rot2Dle0p6ge0p8;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//142 - rotundity (2D) v subcl energy, phiE2D ~ 0 
		TH2D* rot2D_E_phiE2Deq0 = new TH2D("rot2D_subclE_phiE2Deq0","rot2D_subclE_phiE2Deq0;rotundity2D_phiE2Deq0;E;a.u.",25,0.4,1.1,25,0,1000);
		//143 - rot2D v subcl energy, phiE2D !~ 0 
		TH2D* rot2D_E_phiE2Dneq0 = new TH2D("rot2D_subclE_phiE2Dneq0","rot2D_subclE_phiE2Dneq0;rotundity2D_phiE2Dneq0;E;a.u.",25,0.4,1.1,25,0,1000);
                //144 - subcluster phi center, phiE2D ~ 0 
                TH2D* phiCenter_rot2D_phiE2Deq0 = new TH2D("phiCenter_rot2D_phiE2Deq0","phiCenter_rot2D_phiE2Deq0;phiCenter_phiE2Deq0;rot2D;a.u.",25,-0.6,6.6,25,0.4,1.1);
                //145 - subcluster phi center, phiE2D !~ 0 
                TH2D* phiCenter_rot2D_phiE2Dneq0 = new TH2D("phiCenter_rot2D_phiE2Dneq0","phiCenter_rot2D_phiE2Dneq0;phiCenter_phiE2Dneq0;rot2D;a.u.",25,-0.6,6.6,25,0.4,1.1);
		//146 - phi center vs eta center, 0.6 < rot2D < 0.8
		TH2D* phiCenter_etaCenter_rot2Dge0p6le0p8 = new TH2D("phiCenter_etaCenter_rot2Dge0p6le0p8","phiCenter_etaCenter_rot2Dge0p6le0p8;phiCenter_rot2Dge0p6le0p8;etaCenter;a.u.",25,-0.6,6.6,25,-1.6,1.6);
		//147 - eta center vs time center, 0.6 > rotE2D && rot2D > 0.8 
		TH2D* phiCenter_etaCenter_rot2Dle0p6ge0p8 = new TH2D("phiCenter_etaCenter_rot2Dle0p6ge0p8","phiCenter_etaCenter_rot2Dle0p6ge0p8;phiCenter_rot2Dle0p6ge0p8;etaCenter;a.u.",25,-0.6,6.6,25,-1.6,1.6);
		
		//SPIKE CR PLOTS
		//148 - time center vs swCrossPrime
		TH2D* timeCenter_swCrossPrime = new TH2D("timeCenter_swCrossPrime","timeCenter_swCrossPrime;timeCenter;swCrossPrime",25,-15,15,25,-0.05,1.5);
		//150 - time center vs swCrossPrime, phiE2D ~ 0 && phiE2D ~ pi/2
		TH2D* timeCenter_swCrossPrime_phiE2Deq0PiOv2 = new TH2D("timeCenter_swCrossPrime_phiE2Deq0PiOv2","timeCenter_swCrossPrime;timeCenter;swCrossPrime",25,-15,15,25,-0.05,1.5);
		//151 - time center vs swCrossPrime, phiE2D !~ 0 && phiE2D !~ pi/2
		TH2D* timeCenter_swCrossPrime_phiE2Dneq0PiOv2 = new TH2D("timeCenter_swCrossPrime_phiE2Dneq0PiOv2","timeCenter_swCrossPrime;timeCenter;swCrossPrime",25,-15,15,25,-0.05,1.5);
		//151 - rot2D vs swCrossPrime
                TH2D* rot2D_swCrossPrime = new TH2D("rot2D_swCrossPrime","rot2D_swCrossPrime;rot2D;swCrossPrime",25,0.4,1.1,25,-0.05,1.5);
		//152 - rot2D vs swCrossPrime, phiE2D ~ 0 && phiE2D ~ pi/2
                TH2D* rot2D_swCrossPrime_phiE2Deq0PiOv2 = new TH2D("rot2D_swCrossPrime_phiE2Deq0PiOv2","rot2D_swCrossPrime_phiE2Deq0PiOv2;rot2D;swCrossPrime",25,0.4,1.1,25,-0.05,1.5);
		//153 - rot2D vs swCrossPrime, phiE2D !~ 0 && phiE2D !~ pi/2
                TH2D* rot2D_swCrossPrime_phiE2Dneq0PiOv2 = new TH2D("rot2D_swCrossPrime_phiE2Dneq0PiOv2","rot2D_swCrossPrime_phiE2Dneq0PiOv2;rot2D;swCrossPrime",25,0.4,1.1,25,-0.05,1.5);
		//154 - phiE2D vs swCrossPrime
                TH2D* phiE2D_swCrossPrime = new TH2D("phiE2D_swCrossPrime","phiE2D_swCrossPrime;phiE2D;swCrossPrime",25,-3.,3.,25,-0.05,1.5);
		//155 - phoE vs swCrossPrime
		TH2D* phoE_swCrossPrime = new TH2D("phoEnergy_swCrossPrime","phoEnergy_swCrossPrime;phoE;swCrossPrime",25,0,1000,25,-0.05,1.5);	

		//BEAM HALO CR PLOTS
                //156 - eta sigma v phi sigma, phiE2D ~ 0
		TH2D* etaSig_phiSig_phiE2Deq0 = new TH2D("etaSig_phiSig_phiE2Deq0","etaSig_phiSig_phiE2Deq0;etaSig;phiSig_phiE2Deq0",25,0.01,0.09,25,0.01,0.09);
                //157 - eta sigma v phi sigma, phiE2D !~ 0
		TH2D* etaSig_phiSig_phiE2Dneq0 = new TH2D("etaSig_phiSig_phiE2Dneq0","etaSig_phiSig_phiE2Dneq0;etaSig;phiSig_phiE2Dneq0",25,0.01,0.09,25,0.01,0.09);
		//158 - eta center vs phi center, E > 100 && E < 200 && rot2D > 0.7 && rot2D < 0.8 (isolate population in rot2D vs E plot)
		TH2D* etaCenter_phiCenter_Ege100le200_rot2Dge0p7le0p8 = new TH2D("etaCenter_phiCenter_Ege100le200_rot2Dge0p7le0p8","etaCenter_phiCenter_Ege100le200_rot2Dge0p7le0p8;etaCenter;phiCenter_Ege100le200_rot2Dge0p7le0p8",25,-1.6,1.6,25,-0.2,6.4);	
		//159 - eta center vs phi center, phiE2D ~ 0 
		TH2D* etaCenter_phiCenter_phiE2Deq0 = new TH2D("etaCenter_phiCenter_phiE2Deq0","etaCenter_phiCenter_phiE2Deq0;etaCenter;phiCenter_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);
		//160 - time center vs eta center, E > 100 && E < 200 && rot2D > 0.7 && rot2D < 0.8 (isolate population in rot2D vs E plot)
		TH2D* timeCenter_etaCenter_Ege100le200_rot2Dge0p7le0p8 = new TH2D("timeCenter_etaCenter_Ege100le200_rot2Dge0p7le0p8","timeCenter_etaCenter_Ege100le200_rot2Dge0p7le0p8;timeCenter;etaCenter_Ege100le200_rot2Dge0p7le0p8",25,-15,15,25,-1.6,1.6);
		//161 - eta center vs time center, phiE2D ~ 0 
		TH2D* timeCenter_etaCenter_phiE2Deq0 = new TH2D("timeCenter_etaCenter_phiE2Deq0","timeCenter_etaCenter_phiE2Deq0;timeCenter_phiE2Deq0;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//162 - eta center vs time center, phiE2D !~ 0
		TH2D* timeCenter_etaCenter_phiE2Dneq0 = new TH2D("timeCenter_etaCenter_phiE2Dneq0","timeCenter_etaCenter_phiE2Dneq0;timeCenter_phiE2Dneq0;etaCenter;a.u.",25,-15,15,25,-1.6, 1.6);
		//163 - rot2D vs subcl E, phi center ~ 0 && phi center ~ pi 
		TH2D* rot2D_E_phiCentereq0Pi = new TH2D("rot2D_subclE_phiCentereq0Pi","rot2D_subclE_phiCentereq0Pi;rotundity2D_phiCentereq0Pi;E;a.u.",25,0.4,1.1,25,0,1000);
		//164 - eta center vs time center, E > 100 && E < 200, rot2D > 0.7 && rot2D < 0.8 
		TH2D* timeCenter_etaCenter_Ele100ge200_rot2Dle0p7ge0p8 = new TH2D("timeCenter_etaCenter_Ele100ge200_rot2Dle0p7ge0p8","timeCenter_etaCenter_Ele100ge200_rot2Dle0p7ge0p8;timeCenter_Ele100ge200_rot2Dle0p7ge0p8;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//165 - eta center vs phi center, E < 100 && E > 200 && rot2D < 0.7 && rot2D > 0.8 
		TH2D* etaCenter_phiCenter_Ele100ge200_rot2Dle0p7ge0p8 = new TH2D("etaCenter_phiCenter_Ele100ge200_rot2Dle0p7ge0p8","etaCenter_phiCenter_Ele100ge200_rot2Dle0p7ge0p8;etaCenter;phiCenter_Ele100ge200_rot2Dle0p7ge0p8",25,-1.6,1.6,25,-0.2,6.4);	
		//166 - eta center vs phi center, phiE2D !~ 0 
		TH2D* etaCenter_phiCenter_phiE2Dneq0 = new TH2D("etaCenter_phiCenter_phiE2Dneq0","etaCenter_phiCenter_phiE2Dneq0;etaCenter;phiCenter_phiE2Dneq0",25,-1.6,1.6,25,-0.2,6.4);
		//167 - eta center vs time center, phiE2D ~ 0, E > 100 && E < 200, rot2D > 0.7 && rot2D < 0.8 
		TH2D* timeCenter_etaCenter_phiE2Deq0_Ege100le200_rot2Dge0p7le0p8 = new TH2D("timeCenter_etaCenter_phiE2Deq0_Ege100le200_rot2Dge0p7le0p8","timeCenter_etaCenter_phiE2Deq0_Ege100le200_rot2Dge0p7le0p8;timeCenter_phiE2Deq0_Ege100le200_rot2Dge0p7le0p8;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//168 - eta center vs time center, phiE2D !~ 0, E > 100 && E < 200, rot2D > 0.7 && rot2D < 0.8 
		TH2D* timeCenter_etaCenter_phiE2Dneq0_Ege100le200_rot2Dge0p7le0p8 = new TH2D("timeCenter_etaCenter_phiE2Dneq0_Ege100le200_rot2Dge0p7le0p8","timeCenter_etaCenter_phiE2Dneq0_Ege100le200_rot2Dge0p7le0p8;timeCenter_phiE2Dneq0_Ege100le200_rot2Dge0p7le0p8;etaCenter;a.u.",25,-15,15,25,-1.6,1.6);
		//169 - phi sigma vs phi center
		TH2D* phiSig_phiCenter = new TH2D("phiSig_phiCenter","phiSig_phiCenter;phiSig;phiCenter",25,0.01,0.09,25,-0.2,6.4);
		//170 - phi sigma vs phi center, phiE2D ~ 0
		TH2D* phiSig_phiCenter_phiE2Deq0 = new TH2D("phiSig_phiCenter_phiE2Deq0","phiSig_phiCenter;phiSig_phiE2Deq0;phiCenter",25,0.01,0.09,25,-0.2,6.4);
		//171 - phi sigma vs phi center, phiE2D !~ 0
		TH2D* phiSig_phiCenter_phiE2Dneq0 = new TH2D("phiSig_phiCenter_phiE2Dneq0","phiSig_phiCenter;phiSig_phiE2Dneq0;phiCenter",25,0.01,0.09,25,-0.2,6.4);
		//172 - phi sig vs eta center
		TH2D* phiSig_etaCenter = new TH2D("phiSig_etaCenter","phiSig_etaCenter;phiSig;etaCenter",25,0.01,0.09,25,-1.6,1.6);
		//173 - phi sig vs eta center, phiE2D ~ 0
		TH2D* phiSig_etaCenter_phiE2Deq0 = new TH2D("phiSig_etaCenter_phiE2Deq0","phiSig_etaCenter_phiE2Deq0;phiSig;etaCenter_phiE2Deq0",25,0.01,0.09,25,-1.6,1.6);
		//174 - phi sig vs eta center, phiE2D !~ 0
		TH2D* phiSig_etaCenter_phiE2Dneq0 = new TH2D("phiSig_etaCenter_phiE2Dneq0","phiSig_etaCenter_phiE2Dneq0;phiSig;etaCenter_phiE2Dneq0",25,0.01,0.09,25,-1.6,1.6);
		//175 - eta sig vs phi center
		TH2D* etaSig_phiCenter = new TH2D("etaSig_phiCenter","etaSig_phiCenter;etaSig;phiCenter",25,0.01,0.09,25,-0.2,6.4);
		//176 - eta sig vs phi center, phiE2D ~ 0
		TH2D* etaSig_phiCenter_phiE2Deq0 = new TH2D("etaSig_phiCenter_phiE2Deq0","etaSig_phiCenter_phiE2Deq0;etaSig;phiCenter_phiE2Deq0",25,0.01,0.09,25,-0.2,6.4);
		//177 - eta sig vs phi center, phiE2D !~ 0
		TH2D* etaSig_phiCenter_phiE2Dneq0 = new TH2D("etaSig_phiCenter_phiE2Dneq0","etaSig_phiCenter_phiE2Dneq0;etaSig;phiCenter_phiE2Dneq0",25,0.01,0.09,25,-0.2,6.4);
		//178 - eta center vs phiE2D
		TH2D* etaCenter_phiE2D = new TH2D("etaCenter_phiE2D","etaCenter_phiE2D;etaCenter;phiE2D",25,-1.6,1.6,25,-3.1,1.);
		//179 - energy vs eta sig
		TH2D* phoE_etaSig = new TH2D("subclE_etaSig","subclE_etaSig;subclEnergy;etaSig",25,0,1000,25,0.01,0.1);
		//180 - energy vs eta sig, phiE2D ~ 0
		TH2D* phoE_etaSig_phiE2Deq0 = new TH2D("subclE_etaSig_phiE2Deq0","subclE_etaSig_phiE2Deq0;subclEnergy;etaSig_phiE2Deq0",25,0,1000,25,0.01,0.09);
		//181 - energy vs eta sig, phiE2D !~ 0
		TH2D* phoE_etaSig_phiE2Dneq0 = new TH2D("subclE_etaSig_phiE2Dneq0","subclE_etaSig_phiE2Dneq0;subclEnergy;etaSig_phiE2Dneq0",25,0,1000,25,0.01,0.09);
		//182 - energy vs phi sig
		TH2D* phoE_phiSig = new TH2D("subclE_phiSig","subclE_phiSig;subclEnergy;phiSig",25,0,1000,25,0.01,0.1);
		//183 - energy vs phi sig, phiE2D ~ 0
		TH2D* phoE_phiSig_phiE2Deq0 = new TH2D("subclE_phiSig_phiE2Deq0","subclE_phiSig_phiE2Deq0;subclEnergy;phiSig_phiE2Deq0",25,0,1000,25,0.01,0.09);
		//184 - energy vs phi sig, phiE2D !~ 0
		TH2D* phoE_phiSig_phiE2Dneq0 = new TH2D("subclE_phiSig_phiE2Dneq0","subclE_phiSig_phiE2Dneq0;subclEnergy;phiSig_phiE2Dneq0",25,0,1000,25,0.01,0.09);
		//185 - time center vs eta center, phiSigle0p3ANDetaSigge0p3
		TH2D* timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3 = new TH2D("timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3","timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3;timeCenter;etaCenter_phiSigle0p3ANDetaSigge0p3;a.u.",25,-15,15,25,-1.6,1.6);
		//186 - time center vs eta center, phiSigge0p3ORetaSigle0p3
		TH2D* timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3 = new TH2D("timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3","timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3;timeCenter;etaCenter_phiSigge0p3ORetaSigle0p3;a.u.",25,-15,15,25,-1.6,1.6);
		//187 - time center vs eta center, phiE2D ~ 0, phiSigle0p3ANDetaSigge0p3
		TH2D* timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Deq0 = new TH2D("timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Deq0","timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Deq0;timeCenter;etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Deq0;a.u.",25,-15,15,25,-1.6,1.6);
		//188 - time center vs eta center, phiE2D !~ 0, phiSigle0p3ANDetaSigge0p3
		TH2D* timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Dneq0 = new TH2D("timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Dneq0","timeCenter_etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Dneq0;timeCenter;etaCenter_phiSigle0p3ANDetaSigge0p3_phiE2Dneq0;a.u.",25,-15,15,25,-1.6,1.6);
		//189 - time center vs eta center, phiE2D ~ 0, phiSigge0p3ORetaSigle0p3
		TH2D* timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Deq0 = new TH2D("timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Deq0","timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Deq0;timeCenter;etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Deq0;a.u.",25,-15,15,25,-1.6,1.6);
		//190 - time center vs eta center, phiE2D !~ 0, phiSigge0p3ORetaSigle0p3
		TH2D* timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Dneq0 = new TH2D("timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Dneq0","timeCenter_etaCenter_phiSigge0p3ORetaSigle0p3_phiE2Dneq0;timeCenter;etaCenter_phiSigge0p3ORetaSigle0p3;a.u.",25,-15,15,25,-1.6,1.6);
		//191 - eta center vs phi center, timeNeg15toNeg1
		TH2D* etaCenter_phiCenter_timeNeg15toNeg1 = new TH2D("etaCenter_phiCenter_timeNeg15toNeg1","etaCenter_phiCenter_timeNeg15toNeg1;etaCenter;phiCenter_timeNeg15toNeg1",25,-1.6,1.6,25,-0.2,6.4);	
		//192 - eta center vs phi center, timeNeg1to3
		TH2D* etaCenter_phiCenter_timeNeg1to3 = new TH2D("etaCenter_phiCenter_timeNeg1to3","etaCenter_phiCenter_timeNeg1to3;etaCenter;phiCenter_timeNeg1to3",25,-1.6,1.6,25,-0.2,6.4);	
		//193 - eta center vs phi center, time3to15
		TH2D* etaCenter_phiCenter_time3to15 = new TH2D("etaCenter_phiCenter_time3to15","etaCenter_phiCenter_time3to15;etaCenter;phiCenter_time3to15",25,-1.6,1.6,25,-0.2,6.4);	
		//194 - eta center vs phi center, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_phiSigle0p3ANDetaSigge0p3 = new TH2D("etaCenter_phiCenter_phiSigle0p3ANDetaSigge0p3","etaCenter_phiCenter_phiSigle0p3ANDetaSigge0p3;etaCenter;phiCenter_phiSigle0p3ANDetaSigge0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//195 - eta center vs phi center, phiSigge0p3ORetaSigle0p3
		TH2D* etaCenter_phiCenter_phiSigge0p3ORetaSigle0p3 = new TH2D("etaCenter_phiCenter_phiSigge0p3ORetaSigle0p3","etaCenter_phiCenter_phiSigge0p3ORetaSigle0p3;etaCenter;phiCenter_phiSigge0p3ORetaSigle0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//196 - eta center vs phi center, timeNeg15toNeg1, phiE2D ~ 0
		TH2D* etaCenter_phiCenter_timeNeg15toNeg1_phiE2Deq0 = new TH2D("etaCenter_phiCenter_timeNeg15toNeg1_phiE2Deq0","etaCenter_phiCenter_timeNeg15toNeg1_phiE2Deq0;etaCenter;phiCenter_timeNeg15toNeg1_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);	
		//197 - eta center vs phi center, timeNeg15toNeg1, phiE2D !~ 0
		TH2D* etaCenter_phiCenter_timeNeg15toNeg1_phiE2Dneq0 = new TH2D("etaCenter_phiCenter_timeNeg15toNeg1_phiE2Dneq0","etaCenter_phiCenter_timeNeg15toNeg1_phiE2Dneq0;etaCenter;phiCenter_timeNeg15toNeg1_phiE2Dneq0",25,-1.6,1.6,25,-0.2,6.4);	
		//198 - eta center vs phi center, timeNeg1to3, phiE2D ~ 0
		TH2D* etaCenter_phiCenter_timeNeg1to3_phiE2Deq0 = new TH2D("etaCenter_phiCenter_timeNeg1to3_phiE2Deq0","etaCenter_phiCenter_timeNeg1to3_phiE2Deq0;etaCenter;phiCenter_timeNeg1to3_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);	
		//199 - eta center vs phi center, timeNeg1to3, phiE2D !~ 0
		TH2D* etaCenter_phiCenter_timeNeg1to3_phiE2Dneq0 = new TH2D("etaCenter_phiCenter_timeNeg1to3_phiE2Dneq0","etaCenter_phiCenter_timeNeg1to3_phiE2Dneq0;etaCenter;phiCenter_timeNeg1to3_phiE2Dneq0",25,-1.6,1.6,25,-0.2,6.4);	
		//200 - eta center vs phi center, time3to15, phiE2D ~ 0
		TH2D* etaCenter_phiCenter_time3to15_phiE2Deq0 = new TH2D("etaCenter_phiCenter_time3to15_phiE2Deq0","etaCenter_phiCenter_time3to15_phiE2Deq0;etaCenter;phiCenter_time3to15_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);	
		//201 - eta center vs phi center, time3to15, phiE2D !~ 0
		TH2D* etaCenter_phiCenter_time3to15_phiE2Dneq0 = new TH2D("etaCenter_phiCenter_time3to15_phiE2Dneq0","etaCenter_phiCenter_time3to15_phiE2Dneq0;etaCenter;phiCenter_time3to15_phiE2Dneq0",25,-1.6,1.6,25,-0.2,6.4);	
		//202 - eta center vs phi center, timeNeg15toNeg1, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3 = new TH2D("etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3","etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3;etaCenter;phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//203 - eta center vs phi center, timeNeg1to3, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3 = new TH2D("etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3","etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3;etaCenter;phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//204 - eta center vs phi center, time3to15, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3 = new TH2D("etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3","etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3;etaCenter;phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//205 - eta center vs phi center, timeNeg15toNeg1, phiSigge0p3ORetaSigle0p3
		TH2D* etaCenter_phiCenter_timeNeg15toNeg1_phiSigge0p3ORetaSigle0p3 = new TH2D("etaCenter_phiCenter_timeNeg15toNeg1_phiSigge0p3ORetaSigle0p3","etaCenter_phiCenter_timeNeg15toNeg1_phiSigge0p3ORetaSigle0p3;etaCenter;phiCenter_timeNeg15toNeg1_phiSigge0p3ORetaSigle0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//206 - eta center vs phi center, timeNeg1to3, phiSigge0p3ORetaSigle0p3
		TH2D* etaCenter_phiCenter_timeNeg1to3_phiSigge0p3ORetaSigle0p3 = new TH2D("etaCenter_phiCenter_timeNeg1to3_phiSigge0p3ORetaSigle0p3","etaCenter_phiCenter_timeNeg1to3_phiSigge0p3ORetaSigle0p3;etaCenter;phiCenter_timeNeg1to3_phiSigge0p3ORetaSigle0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//207 - eta center vs phi center, time3to15, phiSigge0p3ORetaSigle0p3
		TH2D* etaCenter_phiCenter_time3to15_phiSigge0p3ORetaSigle0p3 = new TH2D("etaCenter_phiCenter_time3to15_phiSigge0p3ORetaSigle0p3","etaCenter_phiCenter_time3to15_phiSigge0p3ORetaSigle0p3;etaCenter;phiCenter_time3to15_phiSigge0p3ORetaSigle0p3",25,-1.6,1.6,25,-0.2,6.4);	
		//208 - eta center vs phi center, timeNeg15toNeg1, phiE2D ~ 0, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3_phiE2Deq0 = new TH2D("etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3_phiE2Deq0","etaCenter_phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3_phiE2Deq0;etaCenter;phiCenter_timeNeg15toNeg1_phiSigle0p3ANDetaSigge0p3_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);	
		//209 - eta center vs phi center, timeNeg1to3, phiE2D ~ 0, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3_phiE2Deq0 = new TH2D("etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3_phiE2Deq0","etaCenter_phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3_phiE2Deq0;etaCenter;phiCenter_timeNeg1to3_phiSigle0p3ANDetaSigge0p3_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);	
		//210 - eta center vs phi center, time3to15, phiE2D ~ 0, phiSigle0p3ANDetaSigge0p3
		TH2D* etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3_phiE2Deq0 = new TH2D("etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3_phiE2Deq0","etaCenter_phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3_phiE2Deq0;etaCenter;phiCenter_time3to15_phiSigle0p3ANDetaSigge0p3_phiE2Deq0",25,-1.6,1.6,25,-0.2,6.4);	
		//211 - rot2D vs subcl E, !(phi center ~ 0 && phi center ~ pi) 
		TH2D* rot2D_E_phiCenterneq0Pi = new TH2D("rot2D_subclE_phiCenterneq0Pi","rot2D_subclE_phiCenterneq0Pi;rotundity2D_phiCenterneq0Pi;E;a.u.",25,0.4,1.1,25,0,1000);
		//212 - phiE2D vs sw+', early times
                TH2D* phiE2D_swCrossPrime_timeNeg15toNeg1 = new TH2D("phiE2D_swCrossPrime_timeNeg15toNeg1","phiE2D_swCrossPrime_timeNeg15toNeg1;phiE2D;swCrossPrime_timeNeg15toNeg1",25,-3.,3.,25,-0.05,1.5);
		//213 - phiE2D vs sw+', prompt times
                TH2D* phiE2D_swCrossPrime_timeNeg1to3 = new TH2D("phiE2D_swCrossPrime_timeNeg1to3","phiE2D_swCrossPrime_timeNeg1to3;phiE2D;swCrossPrime_timeNeg1to3",25,-3.,3.,25,-0.05,1.5);
		//214 - phiE2D vs sw+', late times
                TH2D* phiE2D_swCrossPrime_time3to15 = new TH2D("phiE2D_swCrossPrime_time3to15","phiE2D_swCrossPrime_time3to15;phiE2D;swCrossPrime_time3to15",25,-3,3.,25,-0.05,1.5);
		//215 - phiE2D vs reco MET
		TH2D* phiE2D_recoMet = new TH2D("phiE2D_recoMet","phiE2D_recoMet;phiE2D;recoMet",25,-3,3,25,0,1000);
		//216 - etaSig vs reco MET
		TH2D* etaSig_recoMet = new TH2D("etaSig_recoMet","etaSig_recoMet;etaSig;recoMet",25,-0.01,0.09,25,0,1000);
		//217 - phiSig vs reco MET
		TH2D* phiSig_recoMet = new TH2D("phiSig_recoMet","phiSig_recoMet;phiSig;recoMet",25,-0.01,0.09,25,0,1000);
		//218 - etaCenter vs reco MET
		TH2D* etaCenter_recoMet = new TH2D("etaCenter_recoMet","etaCenter_recoMet;etaCenter;recoMet",25,-1.6,1.6,25,0,1000);
		//219 - phiCenter vs reco MET
		TH2D* phiCenter_recoMet = new TH2D("phiCenter_recoMet","phiCenter_recoMet;phiCenter;recoMet",25,-0.2,6.4,25,0,1000);
		//220 - time center vs reco MET
		TH2D* timeCenter_recoMet = new TH2D("timeCenter_recoMet","timeCenter_recoMet;timeCenter;recoMet",25,-15,15,25,0,1000);
		//221 - etaSig vs timeEtaCov, early times
		TH2D* etaSig_timeEtaCov_timeNeg15toNeg1 = new TH2D("etaSig_timeEtaCov_timeNeg15toNeg1","etaSig_timeEtaCov_timeNeg15toNeg1;etaSig;timeEtaCov_timeNeg15toNeg1",25,-0.01,0.09,25,-1,1);
		//222 - etaSig vs timeEtaCov, early times, phiCenter ~ 0, pi
		TH2D* etaSig_timeEtaCov_timeNeg15toNeg1_phiCenter0pi = new TH2D("etaSig_timeEtaCov_timeNeg15toNeg1_phiCenter0pi","etaSig_timeEtaCov_timeNeg15toNeg1_phiCenter0pi;etaSig;timeEtaCov_timeNeg15toNeg1_phiCenter0pi",25,-0.01,0.09,25,-1,1);
		//223 - etaSig vs timeEtaCov, late times
		TH2D* etaSig_timeEtaCov_time3to15 = new TH2D("etaSig_timeEtaCov_time3to15","etaSig_timeEtaCov_time3to15;etaSig;timeEtaCov_time3to15",25,-0.01,0.09,25,-1,1);
		//224 - etaSig vs etaPhiCov, early times
		TH2D* etaSig_etaPhiCov_timeNeg15toNeg1 = new TH2D("etaSig_etaPhiCov_timeNeg15toNeg1","etaSig_etaPhiCov_timeNeg15toNeg1;etaSig;etaPhiCov_timeNeg15toNeg1",25,-0.01,0.09,25,-1,1);
		//225 - etaSig vs etaPhiCov, early times, phiCenter ~ 0, pi
		TH2D* etaSig_etaPhiCov_timeNeg15toNeg1_phiCenter0pi = new TH2D("etaSig_etaPhiCov_timeNeg15toNeg1_phiCenter0pi","etaSig_etaPhiCov_timeNeg15toNeg1_phiCenter0pi;etaSig;etaPhiCov_timeNeg15toNeg1_phiCenter0pi",25,-0.01,0.09,25,-1,1);
		//226 - etaSig vs etaPhiCov, late times
		TH2D* etaSig_etaPhiCov_time3to15 = new TH2D("etaSig_etaPhiCov_time3to15","etaSig_etaPhiCov_time3to15;etaSig;etaPhiCov_time3to15",25,-0.01,0.09,25,-1,1);
		//227 - MET phi vs phiCenter
		TH2D* metPhi_phiCenter = new TH2D("metPhi_phiCenter","metPhi_phiCenter;metPhi;phiCenter",25,-3.5,3.5,25,-0.2,6.2);
		//228 - MET phi vs phiCenter, etaSig + phiSig cuts
		TH2D* metPhi_phiCenter_etaSigge0p3ANDphiSigle0p3 = new TH2D("metPhi_phiCenter_etaSigge0p3ANDphiSigle0p3","metPhi_phiCenter_etaSigge0p3ANDphiSigle0p3;metPhi;phiCenter_etaSigge0p3ANDphiSigle0p3",25,-3.5,3.5,25,-0.2,6.2);
		//229 - MET phi vs phiCenter, !(etaSig + phiSig cuts)
		TH2D* metPhi_phiCenter_etaSigge0p3ORphiSigle0p3 = new TH2D("metPhi_phiCenter_etaSigge0p3ORphiSigle0p3","metPhi_phiCenter_etaSigge0p3ORphiSigle0p3;metPhi;phiCenter_etaSigge0p3ORphiSigle0p3",25,-3.5,3.5,25,-0.2,6.2);
		//230 - dR trackSubcl vs subcl time
		TH2D* dRtrack_timeSubcl = new TH2D("dRtrack_timeSubcl","dRtrack_timeSubcl;dRtrack;timeSubcl",50,0,5,50,-10,10);	
		//231 - dE trackSubcl vs subcl time	
		TH2D* dEtrack_timeSubcl = new TH2D("dEtrack_timeSubcl","dEtrack_timeSubcl;dEtrack;timeSubcl",25,-2,2,50,-10,10);	
		//232 - dR trackSubcl vs dE trackSubcl	
		TH2D* dRtrack_dEtrack = new TH2D("dRtrack_dEtrack","dRtrack_dEtrack;dRtrack;dEtrack",25,0,5,25,-2,2);	
		//233 - dR trackSubcl vs dE trackSubck, -10 < time subclust < -2	
		TH2D* dRtrack_dEtrack_early = new TH2D("dRtrack_dEtrack_early","dRtrack_dEtrack_timeSubclNeg10toNeg2;dRtrack;dEtrack",25,0,5,25,-2,2);	
		//234 - dR trackSubcl vs dE trackSubck, -2 < time subclust < 2	
		TH2D* dRtrack_dEtrack_prompt = new TH2D("dRtrack_dEtrack_prompt","dRtrack_dEtrack_timeSubclNeg2to2;dRtrack;dEtrack",25,0,5,25,-2,2);
		//235 - dR trackSubcl vs dE trackSubck, 2 < time subclust < 10	
		TH2D* dRtrack_dEtrack_late = new TH2D("dRtrack_dEtrack_late","dRtrack_dEtrack_timeSubcl2to10;dRtrack;dEtrack",25,0,5,25,-2,2);
		//236 - sigma^2_t from meas err
		TH2D* rhEnergy_timesSigSqMeasErr = new TH2D("rhEnergy_timesSigSqMeasErr","rhEnergy_timesSigSqMeasErr;rhEnergy;#sigma^2_t",50,0,10,50,0,5);
		//237 - eta-phi view of overlaid subcl (energy = z axis) in 9x9 grid (oversized)
		TH2D* ENeighbors = new TH2D("ENeighbors","ENeighbors;local ieta;local iphi",9,-4,5,9,-4,5);	
		//238 - eta-phi view of overlaid subcl (energy = z axis) in 30x30 grid (oversized)
		//deta = dphi = pi/180 -> etamax = deta * nbins / 2 = -etamin
		TH2D* etaPhi_overlaidsubcl = new TH2D("etaPhi_overlaidsubcl","etaPhi_overlaidsubcl;eta;phi;energy",30,-0.2618,0.2618,30,-0.2618,0.2618);
		//239 - phoE vs timeSig
		TH2D* phoE_timeSig = new TH2D("subclE_timeSig","subclE_timeSig;subclE;timeSig",25,0,1000,25,0,5);
		//240 - time center vs eta sig
		TH2D* timeCenter_etaSig = new TH2D("timeCenter_etaSig","timeCenter_etaSig;timeCenter;etaSig",25,-3,3,25,0.,0.1);
		//241 - time center vs phi sig
		TH2D* timeCenter_phiSig = new TH2D("timeCenter_phiSig","timeCenter_phiSig;timeCenter;phiSig",25,-3,3,25,0.,0.1);
		//242 - sqrt(subcl E) vs subcl E*eta sig
		TH2D* subclE_sqrtSubclEmultEtaSig = new TH2D("subclE_sqrtSubclEmultEtaSig","subclE_sqrtSubclEmultEtaSig;subclE;subclEmultEtaSig",25,0,1000,25,0,0.02);
		//243 - sqrt(subcl E) vs subcl E*phi sig
		TH2D* subclE_sqrtSubclEmultPhiSig = new TH2D("subclE_sqrtSubclEmultPhiSig","subclE_sqrtSubclEmultPhiSig;subclE;subclEmultPhiSig",25,0,1000,25,0.2,0.8);
		//244 - sqrt(subcl E) vs subcl E*time sig
		TH2D* subclE_sqrtSubclEmultTimeSig = new TH2D("subclE_sqrtSubclEmultTimeSig","subclE_sqrtSubclEmultTimeSig;subclE;subclEmultTimeSig",25,0,1000,25,15,35);
		//245 - timeMajCov*subclE^2 vs timeEtaCov*subclE^2
		TH2D* subclEmultTimeMajCov_subclEmultTimeEtaCov = new TH2D("subclEmultTimeMajCov_subclEmultTimeEtaCov","subclEmultTimeMajCov_subclEmultTimeEtaCov;subclEmultTimeMajCov;subclEmultTimeEtaCov",50,-500,500,50,-500,500);
		//246 - timeeta cov vs timephi cov with |timemaj cov| > 0.1
		TH2D* timeEtaCov_timePhiCov_absTimeMajCovge0p1 = new TH2D("timeEtaCov_timePhiCov_absTimeMajCovge0p1","timeEtaCov_timePhiCov_absTimeMajCovge0p1;timeEtaCov;timePhiCov_absTimeMajCovge0p1",25,-0.1,0.1,25,-0.1,0.1);
		//247 - time maj cov vs time min cov
		TH2D* timeMajCov_timeMinCov = new TH2D("timeMajCov_timeMinCov","timeMajCov_timeMinCov;timeMajCov;timeMinCov",25,-0.1,0.1,25,-0.1,0.1);
		//248 - distance to plane in time vs plane eta rad
		TH2D* tPlaneDist_planeEtaSig = new TH2D("tPlaneDist_planeMin1Sig","tPlaneDist_planeMin1Sig;tPlaneDist;planeMin1Sig",50,-3.,3.,50,0.5,3.5);
		//249 - distance to plane in time vs plane phi rad
		TH2D* tPlaneDist_planePhiSig = new TH2D("tPlaneDist_planeMin2Sig","tPlaneDist_planeMin2Sig;tPlaneDist;planeMin2Sig",50,-3.,3.,50,0.1,0.2);
		//250 - distance to plane in time vs plane etaphi cov
		TH2D* tPlaneDist_planeEtaPhiCov = new TH2D("tPlaneDist_planeEtaPhiCov","tPlaneDist_planeEtaPhiCov;tPlaneDist;planeEtaPhiCov",50,-3.,3.,50,0.05,0.2);


		enum weightScheme{
			noWeight = 0,
			Eweight = 1,
			logEweight = 2
		};
		

		void Skim();

		void AddSample(TFile* file);

		void SetIsoCuts(){ _isocuts = true; }
		bool _isocuts;
		//set skip for outstream
		void SetSkip(int i){ _oskip = i; }
		int _oskip;
		void SetThresh(double t){ _thresh = t; }
		void SetBHCAlpha(double a){ _alpha = a; }
		void SetEMAlpha(double a){ _emAlpha = a; }
		double _thresh, _alpha, _emAlpha, _timeoffset, _swcross; 
		void ApplyFractions(bool a){ _applyFrac = a; if(_applyFrac) cout << "Applying RH fractions" << endl; }
		bool _applyFrac;




		void WritePlotCat1D(TFile* ofile, const procCat& pc){
			ofile->cd();
			string name;
			vector<vector<TH1D*>> hists1D = pc.hists1D;
			//write 1D hists
			for(int i = 0; i < (int)hists1D.size(); i++){
				for(int j = 0; j < hists1D[i].size(); j++){
					if(pc.hists1D[i][j] == nullptr) continue;
					name = pc.hists1D[i][j]->GetName();
					if(pc.hists1D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram: " << name << " not filled." << endl; continue; }
					pc.hists1D[i][j]->Write();
				}
			}
		}
		void WritePlotCat2D(TFile* ofile, const procCat& pc){
			ofile->cd();
			string name;
			vector<vector<TH2D*>> hists2D = pc.hists2D;
			

			//write 2D hists
			for(int i = 0; i < (int)hists2D.size(); i++){
				for(int j = 0; j < hists2D[i].size(); j++){
					if(pc.hists2D[i][j] == nullptr) continue;
					//write total hist to file
					name = pc.hists2D[i][j]->GetName();
					//if ends in "_" remove
					if(strcmp(&name[name.size() - 1],"_") == 0) name.pop_back();
					name += "2D";
					pc.hists2D[i][j]->SetName(name.c_str());		
					pc.hists2D[i][j]->SetTitle(pc.plotName.c_str());
					if(pc.hists2D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram: " << name << " not filled." << endl; continue; }
					pc.hists2D[i][j]->Write();
				}
			}

		}



		void WritePlotCatStack(TFile* ofile, const vector<procCat>& pcs){
			ofile->cd();
			string name;
			//number of histogram categories (ie leading, !leading, etc)
			int nhistCats = pcs[0].hists1D.size();
			//write 1D hists
			//variables
			for(int j = 0; j < (int)_hists1D.size(); j++){
				if(pcs[0].hists1D[0][j]->GetEntries() == 0 && pcs[0].hists1D[1][j]->GetEntries() == 0){ continue; }//cout << "Histogram for proc " << pcs[k].plotName << " not filled." << endl; continue; }
				//lead, not lead, etc.
				name = _hists1D[j]->GetName();
				TDirectory* dir = ofile->mkdir((name+"_stack").c_str());
				dir->cd(); 
				for(int i = 0; i < nhistCats; i++){
					if(pcs[0].hists1D[i][j] == nullptr) continue;
					//needs to be reset for each category
					//name = _hists1D[j]->GetName();
					if(!pcs[0].histcatnames[i].empty()) name += "_"+pcs[0].histcatnames[i];
					//cout << "	category: " << pcs[0].histcatnames[i] << endl; 
					//proc
					for(int k = 0; k < pcs.size(); k++){
						if(pcs[k].hists1D[i][j] == nullptr) continue;
						if(pcs[k].hists1D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram for proc " << pcs[k].plotName << " not filled." << endl; continue; }
						//cout << "		adding proc " << pcs[k].plotName << " to plot with hist " << pcs[k].hists1D[i][j]->GetName() << endl;
						pcs[k].hists1D[i][j]->SetTitle(pcs[k].plotName.c_str());
						pcs[k].hists1D[i][j]->Write();
					}
				}
			}
			
		}


		void WriteHists(TFile* ofile){
			//normalize histograms	
			for(int i = 0; i < (int)_procCats.size(); i++){
				for(int j = 0; j < _procCats[i].hists1D.size(); j++){
					//relative fraction histograms
					//nSubClusters
					//_procCats[i].hists1D[j][0]->Scale(1./_procCats[i].hists1D[j][0]->Integral());
					//ellipsoid center coordinates
					_procCats[i].hists1D[j][1]->Scale(1./_procCats[i].hists1D[j][1]->Integral());
					_procCats[i].hists1D[j][2]->Scale(1./_procCats[i].hists1D[j][2]->Integral());
					_procCats[i].hists1D[j][3]->Scale(1./_procCats[i].hists1D[j][3]->Integral());
					//theta + azimuthal angles
					_procCats[i].hists1D[j][7]->Scale(1./_procCats[i].hists1D[j][7]->Integral());
					_procCats[i].hists1D[j][8]->Scale(1./_procCats[i].hists1D[j][8]->Integral());
					//rotundity
					_procCats[i].hists1D[j][10]->Scale(1./_procCats[i].hists1D[j][10]->Integral());
					_procCats[i].hists1D[j][11]->Scale(1./_procCats[i].hists1D[j][11]->Integral());
					//velocity
					_procCats[i].hists1D[j][14]->Scale(1./_procCats[i].hists1D[j][14]->Integral());
				}
			}

			WritePlotCat1D(ofile, _procCats[0]);
			for(int i = 0; i < (int)_procCats.size(); i++)
				WritePlotCat2D(ofile, _procCats[i]);
			vector<procCat> id_cats(_procCats.begin()+1, _procCats.end());
			
			if(_procCats.size() > 1) WritePlotCatStack(ofile, id_cats);

			ofile->Close();

		}




		//k = sum_n(E_n)/N
		//void FillModelHists(BasePDFMixture* model, int id_idx, vector<double>& obs){
		void FillModelHists(BasePDFMixture* model, int id_idx, map<string,double>& obs){
			obs.clear();	
			map<string, Matrix> params;
			vector<double> eigenvals, eigenvals_space, norms;
			vector<Matrix> eigenvecs, eigenvecs_space; 
			Matrix space_mat = Matrix(2,2);
			Matrix rotmat2D = Matrix(3,3);		
			Matrix majmin3DCovMat = Matrix(3,3);
			Matrix majmin2DCovMat = Matrix(3,3);
			PointCollection majminpts3D, majminpts2D;

			double npts = (double)model->GetData()->GetNPoints();
		//	cout << "FillHists - starting subcluster loop" << endl;	
			double E_k, phi, rot2D, ec, pc, tc, pi, E_lead, phi2D;
			//double theta, r, rot3D, vel, v_x, v_y, v_z;
			double ep_cov, te_cov, tp_cov, e_var, p_var, t_var;
			double ep_cov_unnorm, te_cov_unnorm, tp_cov_unnorm;
			double ep_cov_norm, te_cov_norm, tp_cov_norm;
			double majtime_cov_3d, mintime_cov_3d, majtime_cov_unnorm, mintime_cov_unnorm;
			double majtime_cov_2d, mintime_cov_2d;
			double timespace_cov;		
	
			//for swiss cross prime - wmax/N_k	
			PointCollection* points = model->GetData();
			vector<double> spikeObs;
			SpikeObs(points, spikeObs);
			double swCP = spikeObs[0];
			//calculate weighted eta mean of points to check eta calculation
			double etaCentroid = points->Centroid(0);

	
			int nclusters = model->GetNClusters();
			//get leading cluster index
			vector<int> idxs;
			//sort by mixing coeffs in ascending order (smallest first)
			model->SortIdxs(idxs);
			int leadidx = idxs[nclusters-1];
			int k = leadidx;
			
			double E_tot = 0.;
			int nrhs_thresh = 0;
			double rnk_thresh = 0.8;
			Matrix post = model->GetPosterior();
			for(int i = 0; i < npts; i++){
				E_tot += model->GetData()->at(i).w()/_gev;
				if(post.at(i,k) > rnk_thresh) nrhs_thresh++;
			}
			_procCats[id_idx].hists1D[0][250]->Fill(nrhs_thresh);
			Matrix cov, lead_eigenvec, lead_eigenvec_space;
			
			_procCats[id_idx].hists1D[0][0]->Fill(nclusters);
			_procCats[id_idx].hists2D[0][10]->Fill((double)nclusters,npts);
			

			//fill for lead subcluster only
			//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
			model->GetNorms(norms);
			E_k = norms[k]/_gev; 
			
			
			params = model->GetLHPosteriorParameters(k);
			//params = model->GetDataStatistics(k);
			ec = params["mean"].at(0,0);
			pc = params["mean"].at(1,0);
			if(isnan(pc)) cout << "pc is nan" << endl;
			if(isinf(pc)) cout << "pc is inf" << endl;
			if(pc < 0 || pc > 2*acos(-1)) cout << "pc out of bounds " << pc << endl;
			tc = params["mean"].at(2,0);
			pi = params["pi"].at(0,0);
			cov = params["cov"];	
			//distance from xmax to mean_k
			///double dist = 0;
			///for(int d = 0; d < xmax.Dim(); d++)
			///	dist += (xmax.at(d) - params["mean"].at(d,0))*(xmax.at(d) - params["mean"].at(d,0));
			///dist = sqrt(dist);
			///
			///cout << "wmax " << wmax << " norms_k " << norms[k] << " Nk/wmax " << -(swCP-1) << " wmax/Nk " << wmax/norms[k] << " E_k " << E_k << " Emax " << wmax/_gev << " distance to mean " << dist << " wmax/Nk * (1/dist) " << wmax/(norms[k]*dist) << endl;
			//eta - time sign convention
			//define relative sign for eta and time components
			//based on where the cluster is in the detector
			if(ec < 0){
				//time sign does NOT match eta sign
				//flip sign of eta-time entry
				cov.SetEntry(-cov.at(0,2),0,2);	
				cov.SetEntry(-cov.at(2,0),2,0);
				//flip sign of etaphi cov
				cov.SetEntry(-cov.at(0,1),0,1);	
				cov.SetEntry(-cov.at(1,0),1,0);
			}
			//else time sign matches eta sign - no change
			//phi sign convention - etaphi-cov should always be positive
			if(cov.at(1,0) < 0){
				cov.SetEntry(-cov.at(0,1),0,1);
				cov.SetEntry(-cov.at(1,0),1,0);
				cov.SetEntry(-cov.at(1,2),1,2);
				cov.SetEntry(-cov.at(2,1),2,1);
			}
			
			e_var = sqrt(cov.at(0,0));
			p_var = sqrt(cov.at(1,1));
			t_var = sqrt(cov.at(2,2));

			//normalized covs
			ep_cov = CalcCov(cov, 1, 0);
			te_cov = CalcCov(cov, 2, 0);
			tp_cov = CalcCov(cov, 2, 1);
			ep_cov_unnorm = CalcCov(cov, 1, 0, false);
			te_cov_unnorm = CalcCov(cov, 2, 0, false);
			tp_cov_unnorm = CalcCov(cov, 2, 1, false);


//cout << "lead subcluster has e sig " << e_var << " p sig " << p_var << " t sig " << t_var << " ep cov " << ep_cov << " te cov " << te_cov << " tp cov " << tp_cov << endl; 
//			cout << "exp post cov/e cov" << endl; cov.Print();
			
			vector<Matrix> eigvecs;
			vector<double> eigvals;
			cov.eigenCalc(eigvals,eigvecs);
			double majLength = sqrt(eigvals[2]);
			if(eigvals[1] < 0) cout << "negative eigenvalue " << eigvals[1] << endl;
			double minLength; 
			if(eigvals[1] < 0) minLength = -sqrt(-eigvals[1]);
			else minLength = sqrt(eigvals[1]);	

			//rotate into 3D eigenvector space
			Matrix rotmat3D(3,3);
			Get3DRotationMatrix(eigvecs,rotmat3D);
			RotatePoints(model->GetData(), rotmat3D, majminpts3D);
			MakeCovMat(&majminpts3D, majmin3DCovMat, weightScheme(1));
			majtime_cov_3d = CalcCov(majmin3DCovMat,2,0, false);
			if(majtime_cov_3d < 0){
				eigvecs[2].SetEntry(-eigvecs[2].at(0,0),0,0);
				eigvecs[2].SetEntry(-eigvecs[2].at(1,0),1,0);
				majmin3DCovMat.SetEntry(-majtime_cov_3d,2,0);
				majmin3DCovMat.SetEntry(-majtime_cov_3d,0,2);
			}	
			//angle bw major axis and eta (3D) - eigenvectors normalized
			double eta_angle_3d = acos(eigvecs[2].at(0,0));
			double phi_angle_3d = acos(eigvecs[2].at(1,0));
			double time_angle_3d = acos(eigvecs[2].at(2,0));
			//rotundity - 2D
			//take upper 2x2 submatrix from covariance
			Get2DMat(cov,space_mat);
			space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
			phi2D = PhiEll(space_mat);			
			rot2D = Rotundity(space_mat);
			
			//rotate points into 2D (spatial only) maj/min axes
			Get2DRotationMatrix(eigenvecs_space,rotmat2D);
			RotatePoints(model->GetData(), rotmat2D, majminpts2D);
			MakeCovMat(&majminpts2D, majmin2DCovMat, weightScheme(1));

			
			//set time covariance from GMM for major/minor 2D covariance
			majmin2DCovMat.SetEntry(cov.at(2,2),2,2);
			majtime_cov_2d = CalcCov(majmin2DCovMat,2,0, false);
			mintime_cov_2d = CalcCov(majmin2DCovMat,2,1, false);
			//switch sign of maj axis based on sign of time-maj cov for 2d + 3d major axes
			if(majtime_cov_2d < 0){
				eigenvecs_space[1].SetEntry(-eigenvecs_space[1].at(0,0),0,0);
				eigenvecs_space[1].SetEntry(-eigenvecs_space[1].at(1,0),1,0);
				majmin2DCovMat.SetEntry(-majtime_cov_2d,2,0);
				majmin2DCovMat.SetEntry(-majtime_cov_2d,0,2);
			}
			majtime_cov_unnorm = CalcCov(majmin2DCovMat,2,0,false);
			mintime_cov_unnorm = CalcCov(majmin2DCovMat,2,1,false);
	
			//angle bw major axis and eta (2D)
			double eta_angle_2d = acos(eigenvecs_space[1].at(0,0));
			double phi_angle_2d = acos(eigenvecs_space[1].at(1,0));
			double majLength_2d = sqrt(eigenvals_space[1]);				

		
			//project 2D eta-phi components of 3D lead eigenvector into 2D spatial maj-min space
			//x2d = v_3D_maj \dot v_2D_maj
			//y2d = v_3D_maj \dot v_2D_min
			double x2d = eigvecs[2].at(0,0)*eigenvecs_space[1].at(0,0) + eigvecs[2].at(1,0)*eigenvecs_space[1].at(1,0);
			double y2d = eigvecs[2].at(0,0)*eigenvecs_space[0].at(0,0) + eigvecs[2].at(1,0)*eigenvecs_space[0].at(1,0);
	
			double angleDiff = atan2(y2d,x2d);
			double x2d_eta = 1*eigenvecs_space[1].at(0,0) + 0.*eigenvecs_space[1].at(1,0);
			double y2d_eta = 1*eigenvecs_space[0].at(0,0) + 0.*eigenvecs_space[0].at(1,0);
		
			double angleDiff_eta = atan2(y2d_eta,x2d_eta);
			//if(id_idx == 0) cout << "x2d " << x2d << " y2d " << y2d << " angle to maj of 3d lead " << angleDiff << " sin(angle) " << sin(angleDiff) << " angleDiff_eta " << angleDiff_eta << " eta_angle_2d " << eta_angle_2d << " eta_angle_3d " << eta_angle_3d << " time_angle_3d " << time_angle_3d << endl;
			//if(id_idx == 0) cout << " eigenvals 3d " << eigvals[0] << " " << eigvals[1] << " " << eigvals[2] << endl;
			//if(id_idx == 0) cout << " eigenvals 2d " << eigenvals_space[0] << " " << eigenvals_space[1] << endl;	
			_procCats[id_idx].hists1D[1][261]->Fill(sin(angleDiff));
	
			
			_procCats[id_idx].hists1D[1][262]->Fill(tp_cov/majtime_cov_2d);
			_procCats[id_idx].hists1D[1][263]->Fill(te_cov/majtime_cov_2d);
			_procCats[id_idx].hists1D[1][264]->Fill(majtime_cov_2d/majtime_cov_3d);

			//calculate time-space covariance
			//make space-time rotation matrix from space only eigenvalues + 1 for time
			Matrix spacetimeR(3,3);
			spacetimeR.SetEntry(rotmat2D.at(0,0),0,0);
			spacetimeR.SetEntry(rotmat2D.at(0,1),0,1);
			spacetimeR.SetEntry(rotmat2D.at(1,0),1,0);
			spacetimeR.SetEntry(rotmat2D.at(1,1),1,1);
			spacetimeR.SetEntry(1,2,2);
			//cout << "original covariance" << endl; cov.Print();
			//cout << "maj/min cov" << endl; majminCovMat.Print();
			//cout << "spacetime rotation matrix" << endl; spacetimeR.Print();
			Matrix spacetimeRT;
			spacetimeRT.transpose(spacetimeR);
			Matrix cov_spacetimeR(3,3);
			cov_spacetimeR.mult(spacetimeRT,cov);
			//cout << "spacetime rotated covariance" << endl; cov_spacetimeR.Print();
			cov_spacetimeR.mult(cov_spacetimeR,spacetimeR);

			
			//planar sections - on average (over many subclusters)
			int nstep_tplane = 2;
			double eplane_sig, pplane_sig, t0;
			//testing
			cout << "major axis length " << majLength << " other lengths " << sqrt(eigvals[1]) << " " << sqrt(eigvals[0]) << endl;
			cout << "major axis" << endl; eigvecs[2].Print();
			cout << "rotated major axis" << endl;
			Matrix test(3,1);
			test.mult(rotmat3D,eigvecs[2]);
			test.Print();

			//setting distance to plane along time axis based on time component of major axis
			if(id_idx == 0){
			for(int i = 0; i < 2*majLength*eigvecs[2].at(2,0)*nstep_tplane; i++){
				//create plane at t0 in eta-phi according to major axis (scaled appropriately)
				t0 = majLength*eigvecs[2].at(2,0) - (double)i/nstep_tplane;
				Matrix plane(3,1);
				Matrix planecov(2,2);
				plane.SetEntry(1,2,0); //0*^eta + 0*^phi + ^time = t0
				cout << "distance to plane is " << t0 << endl;	
				//rotate plane s.t. ellipsoid is aligned in coord system (ie is along maj-min axes) bc thats how the plane sections are defined
				RotateNormPlane(rotmat3D, plane, t0);
				cout << "rotated, norm distance," << t0 << " plane " << endl; plane.Print();
				//trad is time-proj of major
				//get eta-phi cov matrix in planar section at (rotated) distance of t0
				//use data rotated into maj-min axes
				EtaPhiPlaneCov(&majminpts3D, majLength*eigvecs[2].at(2,0), plane, t0, planecov);
				
				//EtaPhiPlaneSection(e_var, p_var, t_var, t0, eplane_sig, pplane_sig);
				//cout << "trad " << t_var << " t0 " << t0 << " eplanerad " << eplane_sig << " pplanerad " << pplane_sig << endl;
				if(!planecov.empty()){
				cout << "etaphi cov valid for this d = " << t0 << ", eta sig " << sqrt(planecov.at(0,0)) << " phi sig " << sqrt(planecov.at(1,1)) << " etaphi cov " << planecov.at(0,1) << endl;
					_procCats[id_idx].hists2D[1][248]->Fill(t0,sqrt(planecov.at(0,0)));
					_procCats[id_idx].hists2D[1][249]->Fill(t0,sqrt(planecov.at(1,1)));
					_procCats[id_idx].hists2D[1][250]->Fill(t0,planecov.at(0,1));
				}	
			cout << endl;	
			cout << endl;	
		}
cout << endl;
cout << endl;
cout << endl;
}	

			//do track matching
			int nTracks = _base->ECALTrack_nTracks;
			double bestTrackDr = 999;
			//double maxTrackDr;
			double dr, teta, tphi, de;
			unsigned int detid;
			int ieta, iphi;
			pair<double, double> tcoords;
			double pc_02pi = pc;
			if(pc_02pi < 0) pc_02pi += 2*pi;
			for(int t = 0; t < nTracks; t++){
				//use TrackDetId to see where in ECAL track was propagated to
				detid = _base->ECALTrackDetID_detId->at(t);
				//check if detid in map
				if(_detIDmap.count(detid) == 0) continue;
				iphi = _detIDmap[detid].i1;
				ieta = _detIDmap[detid].i2;
		
				tcoords = iEtaiPhi2EtaPhi(ieta, iphi);
				teta = tcoords.first;
				tphi = tcoords.second;			
	
				dr = sqrt((teta - ec)*(teta - ec) + (tphi - pc_02pi)*(tphi - pc_02pi));
				
				if(dr < bestTrackDr){
					bestTrackDr = dr;
					//E = p for photons
					de = (E_k - _base->ECALTrack_p->at(t))/E_k;
				}
			}

                        obs["eta_center"] = ec;
                        obs["phi_center"] = pc;
                        obs["time_center"] = tc;
                        obs["eta_sig"] = e_var;
                        obs["phi_sig"] = p_var;
                        obs["etaphi_cov"] = ep_cov;
                        obs["timeeta_cov"] = te_cov;
                        obs["energy"] = E_k;
                        obs["sw+"] = swCP;
                        obs["major_length"] = majLength;
                        obs["minor_length"] = minLength;	
			//get max weighted point
			points->Sort();
			double maxE = points->at(points->GetNPoints() - 1).w();
			obs["maxOvtotE"] = maxE/E_k;	
	
			//fill hists - lead only
			//centers
			_procCats[id_idx].hists1D[1][1]->Fill(tc);
			_procCats[id_idx].hists1D[1][2]->Fill(ec);
			_procCats[id_idx].hists1D[1][3]->Fill(pc);
			//4 - phoE filled in .cc
			//cluster E
			_procCats[id_idx].hists1D[1][5]->Fill(E_tot);
			//subcluster energy
			_procCats[id_idx].hists1D[1][11]->Fill(E_k);
			//rotundity measures
			//_procCats[id_idx].hists1D[1][12]->Fill(rot3D);
			_procCats[id_idx].hists1D[1][13]->Fill(rot2D);
			//get variances
			_procCats[id_idx].hists1D[1][16]->Fill(e_var);
			_procCats[id_idx].hists1D[1][17]->Fill(p_var);
			_procCats[id_idx].hists1D[1][18]->Fill(t_var);
			//fractional E
			_procCats[id_idx].hists1D[1][19]->Fill(E_k/E_tot);
			//"azimuth" angle in 2D (angle from x axis)
			_procCats[id_idx].hists1D[1][20]->Fill(phi2D);
			//pos/neg eta split sigma
			if(ec > 0) _procCats[id_idx].hists1D[1][21]->Fill(e_var);
			else _procCats[id_idx].hists1D[1][22]->Fill(e_var);
			//covariances
			_procCats[id_idx].hists1D[1][23]->Fill(ep_cov);
			_procCats[id_idx].hists1D[1][24]->Fill(te_cov);
			_procCats[id_idx].hists1D[1][25]->Fill(tp_cov);
			//major/minor covariances with time
			_procCats[id_idx].hists1D[1][26]->Fill(majtime_cov_2d);
			_procCats[id_idx].hists1D[1][27]->Fill(mintime_cov_2d);
			if(E_tot <= 100){
				_procCats[id_idx].hists1D[1][28]->Fill(t_var);
				_procCats[id_idx].hists1D[1][76]->Fill(rot2D);	
			}
			else if(E_tot > 100 && E_tot <= 200){
				_procCats[id_idx].hists1D[1][29]->Fill(t_var);
				_procCats[id_idx].hists1D[1][77]->Fill(rot2D);	
			}
			else{
				_procCats[id_idx].hists1D[1][30]->Fill(t_var);
				_procCats[id_idx].hists1D[1][78]->Fill(rot2D);	
			}
			_procCats[id_idx].hists1D[1][31]->Fill(te_cov_unnorm);
			_procCats[id_idx].hists1D[1][32]->Fill(tp_cov_unnorm);
			_procCats[id_idx].hists1D[1][33]->Fill(majtime_cov_unnorm);
			_procCats[id_idx].hists1D[1][34]->Fill(mintime_cov_unnorm);
			_procCats[id_idx].hists1D[1][53]->Fill(ep_cov_unnorm);
			_procCats[id_idx].hists2D[1][72]->Fill(E_tot,phi2D);
			if(ep_cov < 0){
				if(te_cov < 0){
					_procCats[id_idx].hists1D[1][79]->Fill(rot2D);
					_procCats[id_idx].hists1D[1][83]->Fill(phi2D);
					_procCats[id_idx].hists1D[1][87]->Fill(ep_cov);
					_procCats[id_idx].hists1D[1][91]->Fill(te_cov);
				}	
				else{
					_procCats[id_idx].hists1D[1][80]->Fill(rot2D);
					_procCats[id_idx].hists1D[1][84]->Fill(phi2D);
					_procCats[id_idx].hists1D[1][88]->Fill(ep_cov);
					_procCats[id_idx].hists1D[1][92]->Fill(te_cov);
				}
			}
			else{
				if(te_cov < 0){
					_procCats[id_idx].hists1D[1][81]->Fill(rot2D);
					_procCats[id_idx].hists1D[1][85]->Fill(phi2D);
					_procCats[id_idx].hists1D[1][89]->Fill(ep_cov);
					_procCats[id_idx].hists1D[1][93]->Fill(te_cov);
				}	
				else{
					_procCats[id_idx].hists1D[1][82]->Fill(rot2D);
					_procCats[id_idx].hists1D[1][86]->Fill(phi2D);
					_procCats[id_idx].hists1D[1][90]->Fill(ep_cov);
					_procCats[id_idx].hists1D[1][94]->Fill(te_cov);

				}
			}
			if((phi2D < 0.2 && phi2D > -0.2) || (fabs(phi2D) < 0.2+acos(-1)/2. && fabs(phi2D) > -0.2+acos(-1)/2.)){
				_procCats[id_idx].hists1D[1][95]->Fill(ec);
				_procCats[id_idx].hists1D[1][97]->Fill(pc);
				_procCats[id_idx].hists1D[1][99]->Fill(tc);
				_procCats[id_idx].hists1D[1][101]->Fill(E_tot);
				_procCats[id_idx].hists1D[1][103]->Fill(nclusters);
				_procCats[id_idx].hists1D[1][105]->Fill(e_var);
				_procCats[id_idx].hists1D[1][107]->Fill(p_var);
				_procCats[id_idx].hists1D[1][109]->Fill(t_var);
				_procCats[id_idx].hists1D[1][111]->Fill(ep_cov);
				_procCats[id_idx].hists1D[1][113]->Fill(te_cov);
				_procCats[id_idx].hists1D[1][115]->Fill(tp_cov);
				_procCats[id_idx].hists1D[1][117]->Fill(majtime_cov_2d);
				_procCats[id_idx].hists1D[1][119]->Fill(mintime_cov_2d);
				_procCats[id_idx].hists1D[1][121]->Fill(phi2D);
				_procCats[id_idx].hists1D[1][123]->Fill(rot2D);
			}
			else{
				_procCats[id_idx].hists1D[1][96]->Fill(ec);
				_procCats[id_idx].hists1D[1][98]->Fill(pc);
				_procCats[id_idx].hists1D[1][100]->Fill(tc);
				_procCats[id_idx].hists1D[1][102]->Fill(E_tot);
				_procCats[id_idx].hists1D[1][104]->Fill(nclusters);
				_procCats[id_idx].hists1D[1][106]->Fill(e_var);
				_procCats[id_idx].hists1D[1][108]->Fill(p_var);
				_procCats[id_idx].hists1D[1][110]->Fill(t_var);
				_procCats[id_idx].hists1D[1][112]->Fill(ep_cov);
				_procCats[id_idx].hists1D[1][114]->Fill(te_cov);
				_procCats[id_idx].hists1D[1][116]->Fill(tp_cov);
				_procCats[id_idx].hists1D[1][118]->Fill(majtime_cov_2d);
				_procCats[id_idx].hists1D[1][120]->Fill(mintime_cov_2d);
				_procCats[id_idx].hists1D[1][122]->Fill(phi2D);
				_procCats[id_idx].hists1D[1][124]->Fill(rot2D);
			}
			_procCats[id_idx].hists1D[1][223]->Fill(_swcross);
			_procCats[id_idx].hists1D[1][233]->Fill(swCP);
			_procCats[id_idx].hists1D[1][234]->Fill(etaCentroid - ec);
			double dPhiMetpc = _base->Met_phi - pc;
			if(dPhiMetpc > acos(-1)) dPhiMetpc -= 2*acos(-1);
			else if(dPhiMetpc < -acos(-1)) dPhiMetpc += 2*acos(-1);
			_procCats[id_idx].hists1D[1][235]->Fill(dPhiMetpc);
			if(e_var > 0.03 && p_var < 0.03) _procCats[id_idx].hists1D[1][236]->Fill(dPhiMetpc);
			_procCats[id_idx].hists1D[1][237]->Fill(bestTrackDr);
			_procCats[id_idx].hists1D[1][238]->Fill(de);

			_procCats[id_idx].hists1D[1][251]->Fill(eta_angle_3d);
			_procCats[id_idx].hists1D[1][252]->Fill(phi_angle_3d);
			_procCats[id_idx].hists1D[1][253]->Fill(eta_angle_2d);
			_procCats[id_idx].hists1D[1][254]->Fill(phi_angle_2d);
			_procCats[id_idx].hists1D[1][255]->Fill(majLength);
			_procCats[id_idx].hists1D[1][256]->Fill(majLength_2d);

			if(fabs(majtime_cov_2d) > 0.1){
				_procCats[id_idx].hists1D[1][260]->Fill(eta_angle_2d);
				_procCats[id_idx].hists2D[1][246]->Fill(te_cov,tp_cov);
			}


	
			//2D hists
			_procCats[id_idx].hists2D[1][0]->Fill(tc, E_k);
			_procCats[id_idx].hists2D[1][1]->Fill(phi,E_k);
			_procCats[id_idx].hists2D[1][2]->Fill(rot2D,E_k);
			_procCats[id_idx].hists2D[1][3]->Fill(ec,pc);
			_procCats[id_idx].hists2D[1][4]->Fill(tc,ec);
			_procCats[id_idx].hists2D[1][5]->Fill(tc,pc);
			_procCats[id_idx].hists2D[1][6]->Fill(tc,pi);
			_procCats[id_idx].hists2D[1][7]->Fill(E_k,pi);
			//_procCats[id_idx].hists2D[1][8]->Fill(rot3D,E_k);
			_procCats[id_idx].hists2D[1][9]->Fill(norms[k], E_k);
			_procCats[id_idx].hists2D[1][11]->Fill(nclusters, pi);
			_procCats[id_idx].hists2D[1][12]->Fill(e_var, p_var);
			_procCats[id_idx].hists2D[1][13]->Fill(t_var, e_var);
			_procCats[id_idx].hists2D[1][14]->Fill(t_var, p_var);
			_procCats[id_idx].hists2D[1][15]->Fill(E_k/E_tot, pi);
			_procCats[id_idx].hists2D[1][16]->Fill(nclusters, E_k/E_tot);
			_procCats[id_idx].hists2D[1][17]->Fill(tc, t_var);
			_procCats[id_idx].hists2D[1][18]->Fill(rot2D, phi2D);
			_procCats[id_idx].hists2D[1][19]->Fill(t_var, E_k/E_tot);
			_procCats[id_idx].hists2D[1][20]->Fill(t_var, E_tot);
			_procCats[id_idx].hists2D[1][21]->Fill(te_cov, E_tot);
			_procCats[id_idx].hists2D[1][22]->Fill(tp_cov, E_tot);
			_procCats[id_idx].hists2D[1][23]->Fill(t_var, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][24]->Fill(t_var, mintime_cov_2d);
			_procCats[id_idx].hists2D[1][25]->Fill(phi2D, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][26]->Fill(phi2D, mintime_cov_2d);
			_procCats[id_idx].hists2D[1][27]->Fill(rot2D, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][28]->Fill(rot2D, mintime_cov_2d);
			//_procCats[id_idx].hists2D[1][29]->Fill(rot3D, majtime_cov);
			//_procCats[id_idx].hists2D[1][30]->Fill(rot3D, mintime_cov);
			_procCats[id_idx].hists2D[1][31]->Fill(tc, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][32]->Fill(tc, mintime_cov_2d);
			_procCats[id_idx].hists2D[1][33]->Fill(ep_cov, phi2D);
			_procCats[id_idx].hists2D[1][34]->Fill(majtime_cov_2d, phi2D);
			_procCats[id_idx].hists2D[1][35]->Fill(mintime_cov_2d, phi2D);
			_procCats[id_idx].hists2D[1][36]->Fill(rot2D, majtime_cov_unnorm);
			_procCats[id_idx].hists2D[1][37]->Fill(rot2D, mintime_cov_unnorm);
			_procCats[id_idx].hists2D[1][38]->Fill(t_var, te_cov);
			_procCats[id_idx].hists2D[1][39]->Fill(tc, ep_cov);
			_procCats[id_idx].hists2D[1][40]->Fill(tc, te_cov);
			_procCats[id_idx].hists2D[1][41]->Fill(tc, tp_cov);
			_procCats[id_idx].hists2D[1][46]->Fill(ep_cov, te_cov);
			_procCats[id_idx].hists2D[1][47]->Fill(te_cov, tp_cov);
			_procCats[id_idx].hists2D[1][48]->Fill(ep_cov, tp_cov);
			_procCats[id_idx].hists2D[1][49]->Fill(ep_cov, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][50]->Fill(ep_cov, mintime_cov_2d);
			_procCats[id_idx].hists2D[1][51]->Fill(tp_cov, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][52]->Fill(tp_cov, mintime_cov_2d);
			_procCats[id_idx].hists2D[1][53]->Fill(te_cov, majtime_cov_2d);
			_procCats[id_idx].hists2D[1][54]->Fill(te_cov, mintime_cov_2d);
			_procCats[id_idx].hists2D[1][55]->Fill(ep_cov_unnorm, te_cov_unnorm);
			_procCats[id_idx].hists2D[1][56]->Fill(te_cov_unnorm, tp_cov_unnorm);
			_procCats[id_idx].hists2D[1][57]->Fill(ep_cov_unnorm, tp_cov_unnorm);
			_procCats[id_idx].hists2D[1][58]->Fill(ep_cov_unnorm, majtime_cov_unnorm);
			_procCats[id_idx].hists2D[1][59]->Fill(ep_cov_unnorm, mintime_cov_unnorm);
			_procCats[id_idx].hists2D[1][60]->Fill(tp_cov_unnorm, majtime_cov_unnorm);
			_procCats[id_idx].hists2D[1][61]->Fill(tp_cov_unnorm, mintime_cov_unnorm);
			_procCats[id_idx].hists2D[1][62]->Fill(te_cov_unnorm, majtime_cov_unnorm);
			_procCats[id_idx].hists2D[1][63]->Fill(te_cov_unnorm, mintime_cov_unnorm);
			_procCats[id_idx].hists2D[1][64]->Fill(rot2D, ep_cov);
			_procCats[id_idx].hists2D[1][65]->Fill(rot2D, ep_cov_unnorm);
			_procCats[id_idx].hists2D[1][66]->Fill(E_k, pc);
			_procCats[id_idx].hists2D[1][67]->Fill(rot2D,phi2D);
			_procCats[id_idx].hists2D[1][68]->Fill(ep_cov,te_cov);
			_procCats[id_idx].hists2D[1][69]->Fill(tc,phi2D);
			_procCats[id_idx].hists2D[1][73]->Fill(te_cov,rot2D);
			_procCats[id_idx].hists2D[1][239]->Fill(E_k,t_var);
			_procCats[id_idx].hists2D[1][240]->Fill(tc,e_var);
			_procCats[id_idx].hists2D[1][241]->Fill(tc,p_var);
			_procCats[id_idx].hists2D[1][242]->Fill(E_k,e_var);
			_procCats[id_idx].hists2D[1][243]->Fill(E_k,sqrt(E_k)*p_var);
			_procCats[id_idx].hists2D[1][244]->Fill(E_k,sqrt(E_k)*t_var);
			_procCats[id_idx].hists2D[1][245]->Fill(E_k*majtime_cov_2d,E_k*te_cov);
			_procCats[id_idx].hists2D[1][247]->Fill(majtime_cov_2d,mintime_cov_2d);
	
			


	

			if((phi2D < 0.2 && phi2D > -0.2) || (fabs(phi2D) < 0.2+acos(-1)/2. && fabs(phi2D) > -0.2+acos(-1)/2.)){
				_procCats[id_idx].hists2D[1][70]->Fill(rot2D,ep_cov);
				_procCats[id_idx].hists2D[1][74]->Fill(rot2D,te_cov);
				_procCats[id_idx].hists2D[1][100]->Fill(ep_cov,te_cov);
				_procCats[id_idx].hists2D[1][107]->Fill(tc, ec);
				_procCats[id_idx].hists2D[1][117]->Fill(tc, _swcross);
				_procCats[id_idx].hists2D[1][120]->Fill(rot2D, _swcross);
				_procCats[id_idx].hists2D[1][149]->Fill(tc,swCP);
				_procCats[id_idx].hists2D[1][152]->Fill(rot2D,swCP);
			}
			else{
				_procCats[id_idx].hists2D[1][71]->Fill(rot2D,ep_cov);
				_procCats[id_idx].hists2D[1][75]->Fill(rot2D,te_cov);
				_procCats[id_idx].hists2D[1][101]->Fill(ep_cov,te_cov);
				_procCats[id_idx].hists2D[1][108]->Fill(tc, ec);
				_procCats[id_idx].hists2D[1][118]->Fill(tc, _swcross);
				_procCats[id_idx].hists2D[1][121]->Fill(rot2D, _swcross);
				_procCats[id_idx].hists2D[1][150]->Fill(tc,swCP);
				_procCats[id_idx].hists2D[1][153]->Fill(rot2D,swCP);
			}

			
			if(-10 <= tc && tc < -2){
				_procCats[id_idx].hists2D[1][102]->Fill(E_tot, phi2D);
				_procCats[id_idx].hists2D[1][124]->Fill(rot2D, phi2D);
			}
			if(-2 <= tc && tc < 5){
				_procCats[id_idx].hists2D[1][103]->Fill(E_tot, phi2D);
				_procCats[id_idx].hists2D[1][125]->Fill(rot2D, phi2D);
			}
			if(5 <= tc && tc < 10){
				_procCats[id_idx].hists2D[1][104]->Fill(E_tot, phi2D);
				_procCats[id_idx].hists2D[1][126]->Fill(rot2D, phi2D);
			}
			if(10 <= tc && tc < 15){
				_procCats[id_idx].hists2D[1][105]->Fill(E_tot, phi2D);
				//below is not a typo - the window for this histogram is 5 < t < 15
				_procCats[id_idx].hists2D[1][126]->Fill(rot2D, phi2D);
			}
				_procCats[id_idx].hists2D[1][106]->Fill(tc, ec);
				_procCats[id_idx].hists2D[1][115]->Fill(tc,rot2D);
				_procCats[id_idx].hists2D[1][116]->Fill(tc, _swcross);
				_procCats[id_idx].hists2D[1][119]->Fill(rot2D, _swcross);
				_procCats[id_idx].hists2D[1][122]->Fill(phi2D, _swcross);
				_procCats[id_idx].hists2D[1][123]->Fill(E_tot, _swcross);
			if(rot2D > 0.6 && rot2D < 0.8){
				_procCats[id_idx].hists2D[1][140]->Fill(tc, ec);
				_procCats[id_idx].hists2D[1][146]->Fill(pc, ec);
			}
			else{
				_procCats[id_idx].hists2D[1][141]->Fill(tc, ec);
				_procCats[id_idx].hists2D[1][147]->Fill(pc, ec);
			}		
			_procCats[id_idx].hists2D[1][148]->Fill(tc,swCP);
			_procCats[id_idx].hists2D[1][151]->Fill(rot2D,swCP);
			_procCats[id_idx].hists2D[1][154]->Fill(phi2D,swCP);
			_procCats[id_idx].hists2D[1][155]->Fill(E_tot,swCP);
			//wider phiE2D window & only for phiE2D ~ 0 for beam halo
			if((phi2D < 0.5 && phi2D > -0.5)){
				_procCats[id_idx].hists2D[1][142]->Fill(rot2D, E_k);
				_procCats[id_idx].hists2D[1][144]->Fill(pc, rot2D);
				_procCats[id_idx].hists2D[1][156]->Fill(e_var,p_var);	
				_procCats[id_idx].hists2D[1][159]->Fill(ec, pc);
				_procCats[id_idx].hists2D[1][161]->Fill(tc,ec);
				_procCats[id_idx].hists2D[1][170]->Fill(p_var,pc);
				_procCats[id_idx].hists2D[1][173]->Fill(p_var,ec);
				_procCats[id_idx].hists2D[1][176]->Fill(e_var,pc);
				_procCats[id_idx].hists2D[1][180]->Fill(E_k,e_var);
				_procCats[id_idx].hists2D[1][183]->Fill(E_k,p_var);
			}
			else{
				_procCats[id_idx].hists2D[1][143]->Fill(rot2D, E_k);
				_procCats[id_idx].hists2D[1][145]->Fill(pc, rot2D);
				_procCats[id_idx].hists2D[1][157]->Fill(e_var,p_var);	
				_procCats[id_idx].hists2D[1][162]->Fill(tc,ec);
				_procCats[id_idx].hists2D[1][166]->Fill(ec,pc);
				_procCats[id_idx].hists2D[1][171]->Fill(p_var,pc);
				_procCats[id_idx].hists2D[1][174]->Fill(p_var,ec);
				_procCats[id_idx].hists2D[1][177]->Fill(e_var,pc);
				_procCats[id_idx].hists2D[1][181]->Fill(E_k,e_var);
				_procCats[id_idx].hists2D[1][184]->Fill(E_k,p_var);
			}
			//isolate population in E and rot2D
			if((E_k > 100 && E_k < 200) && (rot2D > 0.7 && rot2D < 0.8)){
				_procCats[id_idx].hists2D[1][158]->Fill(ec, pc);
				_procCats[id_idx].hists2D[1][160]->Fill(tc,ec);
				if((phi2D < 0.5 && phi2D > -0.5)){
					_procCats[id_idx].hists2D[1][167]->Fill(tc,ec);
				}
				else{
					_procCats[id_idx].hists2D[1][168]->Fill(tc,ec);
		
				}
			}
			else{
				_procCats[id_idx].hists2D[1][164]->Fill(tc,ec);
				_procCats[id_idx].hists2D[1][165]->Fill(ec,pc);
			}
			//phi center around 0, 2pi and pi
			if(pc < 0.2 || pc > 2*acos(-1)-0.2 || (pc > acos(-1)-0.2 && pc < acos(-1)+0.2)){
				_procCats[id_idx].hists2D[1][163]->Fill(rot2D,E_k);
			}
			else{
				_procCats[id_idx].hists2D[1][211]->Fill(rot2D,E_k);

			} 
			_procCats[id_idx].hists2D[1][169]->Fill(p_var,pc);
			_procCats[id_idx].hists2D[1][172]->Fill(p_var,ec);
			_procCats[id_idx].hists2D[1][175]->Fill(e_var,pc);
			_procCats[id_idx].hists2D[1][178]->Fill(ec,phi2D);
			_procCats[id_idx].hists2D[1][179]->Fill(E_k,e_var);
			_procCats[id_idx].hists2D[1][182]->Fill(E_k,p_var);
			if(p_var < 0.03 && e_var > 0.03){
				_procCats[id_idx].hists2D[1][185]->Fill(tc,ec);
				_procCats[id_idx].hists2D[1][194]->Fill(ec,pc);
				if((phi2D < 0.5 && phi2D > -0.5)){
					_procCats[id_idx].hists2D[1][187]->Fill(tc,ec);
				
				}
				else{
					_procCats[id_idx].hists2D[1][188]->Fill(tc,ec);

				}
			}
			else{
				_procCats[id_idx].hists2D[1][186]->Fill(tc,ec);
				_procCats[id_idx].hists2D[1][195]->Fill(ec,pc);
				if((phi2D < 0.5 && phi2D > -0.5)){
					_procCats[id_idx].hists2D[1][189]->Fill(tc,ec);
				
				}
				else{
					_procCats[id_idx].hists2D[1][190]->Fill(tc,ec);

				}

			}

			if(tc > -15 && tc <= -1){
				_procCats[id_idx].hists2D[1][191]->Fill(ec,pc);
				if((phi2D < 0.5 && phi2D > -0.5)){
					_procCats[id_idx].hists2D[1][196]->Fill(ec,pc);
				
				}
				else{
					_procCats[id_idx].hists2D[1][197]->Fill(ec,pc);

				}
				if(p_var < 0.03 && e_var > 0.03){
					_procCats[id_idx].hists2D[1][202]->Fill(ec,pc);
					if((phi2D < 0.5 && phi2D > -0.5)){
						_procCats[id_idx].hists2D[1][208]->Fill(ec,pc);
					
					}
				}
				else{
					_procCats[id_idx].hists2D[1][205]->Fill(ec,pc);
				
				}	
			}
			if(tc > -1 && tc <= 3){
				_procCats[id_idx].hists2D[1][192]->Fill(ec,pc);
				if((phi2D < 0.5 && phi2D > -0.5)){
					_procCats[id_idx].hists2D[1][198]->Fill(ec,pc);
				
				}
				else{
					_procCats[id_idx].hists2D[1][199]->Fill(ec,pc);

				}
				if(p_var < 0.03 && e_var > 0.03){
					_procCats[id_idx].hists2D[1][203]->Fill(ec,pc);
					if((phi2D < 0.5 && phi2D > -0.5)){
						_procCats[id_idx].hists2D[1][209]->Fill(ec,pc);
					
					}
				}
				else{
					_procCats[id_idx].hists2D[1][206]->Fill(ec,pc);
				
				}	

	
			}
			if(tc > 3 && tc <= 15){
				_procCats[id_idx].hists2D[1][193]->Fill(ec,pc);
				if((phi2D < 0.5 && phi2D > -0.5)){
					_procCats[id_idx].hists2D[1][200]->Fill(ec,pc);
				
				}
				else{
					_procCats[id_idx].hists2D[1][201]->Fill(ec,pc);

				}
				if(p_var < 0.03 && e_var > 0.03){
					_procCats[id_idx].hists2D[1][204]->Fill(ec,pc);
					if((phi2D < 0.5 && phi2D > -0.5)){
						_procCats[id_idx].hists2D[1][210]->Fill(ec,pc);
					
					}
				}
				else{
					_procCats[id_idx].hists2D[1][207]->Fill(ec,pc);
				
				}	

			}
			if(tc > -15 && tc < -1) _procCats[id_idx].hists2D[1][212]->Fill(phi2D,swCP);
			if(tc > -1 && tc < 3) _procCats[id_idx].hists2D[1][213]->Fill(phi2D,swCP);		
			if(tc > 3 && tc < 15) _procCats[id_idx].hists2D[1][214]->Fill(phi2D,swCP);		

			_procCats[id_idx].hists2D[1][215]->Fill(phi2D,_base->Met_pt);
			_procCats[id_idx].hists2D[1][216]->Fill(e_var,_base->Met_pt);
			_procCats[id_idx].hists2D[1][217]->Fill(p_var,_base->Met_pt);
			_procCats[id_idx].hists2D[1][218]->Fill(ec,_base->Met_pt);
			_procCats[id_idx].hists2D[1][219]->Fill(pc,_base->Met_pt);
			_procCats[id_idx].hists2D[1][220]->Fill(tc,_base->Met_pt);
			if(tc > -15 && tc < -1) _procCats[id_idx].hists2D[1][221]->Fill(e_var,te_cov);
			if(tc > -15 && tc < -1){
				if(pc > -0.2 && pc < 0.2) _procCats[id_idx].hists2D[1][222]->Fill(e_var,te_cov);
				if(pc > acos(-1)-0.2 && pc < acos(-1)+0.2) _procCats[id_idx].hists2D[1][222]->Fill(e_var,te_cov);
				if(pc > 2*acos(-1)-0.2 && pc < 2*acos(-1)+0.2) _procCats[id_idx].hists2D[1][222]->Fill(e_var,te_cov);
			}
			if(tc > 3 && tc < 15) _procCats[id_idx].hists2D[1][223]->Fill(e_var,te_cov);
			if(tc > -15 && tc < -1) _procCats[id_idx].hists2D[1][224]->Fill(e_var,ep_cov);
			if(tc > -15 && tc < -1){
				if(pc > -0.2 && pc < 0.2) _procCats[id_idx].hists2D[1][225]->Fill(e_var,ep_cov);
				if(pc > acos(-1)-0.2 && pc < acos(-1)+0.2) _procCats[id_idx].hists2D[1][225]->Fill(e_var,ep_cov);
				if(pc > 2*acos(-1)-0.2 && pc < 2*acos(-1)+0.2) _procCats[id_idx].hists2D[1][225]->Fill(e_var,ep_cov);
			}
			if(tc > 3 && tc < 15) _procCats[id_idx].hists2D[1][226]->Fill(e_var,ep_cov);
			_procCats[id_idx].hists2D[1][227]->Fill(_base->Met_phi,pc);	
			if(p_var < 0.03 && e_var > 0.03) _procCats[id_idx].hists2D[1][228]->Fill(_base->Met_phi,pc);	
			else _procCats[id_idx].hists2D[1][229]->Fill(_base->Met_phi,pc);	
			_procCats[id_idx].hists2D[1][230]->Fill(bestTrackDr,tc);	
			_procCats[id_idx].hists2D[1][231]->Fill(de,tc);	
			_procCats[id_idx].hists2D[1][232]->Fill(bestTrackDr,de);	
			if(tc > -10 && tc < -2) _procCats[id_idx].hists2D[1][233]->Fill(bestTrackDr,de);
			if(tc > -2 && tc < 2) _procCats[id_idx].hists2D[1][234]->Fill(bestTrackDr,de);
			if(tc > 2 && tc < 10) _procCats[id_idx].hists2D[1][235]->Fill(bestTrackDr,de);

		}


	BayesPoint GetCMSSWMean(PointCollection* pc, bool logw = false){
		int npts = pc->GetNPoints();
		int maxd = pc->Dim();
		double E_tot = 0.;
		//CMSSW value for super clusters
		double w0 = 4.7;
		//zero suppression involves a hitsAndFractions transfer factor (noZS in cmssw)
		//not user here (is 1)
		for(int i = 0; i < npts; i++){
			E_tot += pc->at(i).w()/_gev;
		}
		PointCollection pcnew;
		BayesPoint pt;
		for(int i = 0; i < npts; i ++){
			pt = pc->at(i);
			if(logw) pt.SetWeight( log( w0 + (pc->at(i).w()/_gev)/E_tot ) );
			else pt.SetWeight( 1.0 );
			//do phi wraparound
			if(pt.at(1) > acos(-1)) pt.SetValue(-pt.at(1),1);
			pcnew.AddPoint(pt);
		}
		//is weighted mean
		BayesPoint mean = BayesPoint(maxd);
		for(int d = 0; d < maxd; d++)
			mean.SetValue(pcnew.Centroid(d),d);

		return mean;
	}

	void MakeCovMat(PointCollection* pc, Matrix& outcov, const weightScheme& ws){
		if(!outcov.square()) return;
		if(outcov.GetDims()[0] != pc->Dim()) return;
		
		//set weights to logE
		//og w_i = _gev*E_i
		int npts = pc->GetNPoints();
		int maxd = pc->Dim();
		double E_tot = 0.;
		//CMSSW value for super clusters
		double w0 = 4.7;
		//zero suppression involves a hitsAndFractions transfer factor (noZS in cmssw)
		//not user here (is 1)
		for(int i = 0; i < npts; i++){
			E_tot += pc->at(i).w()/_gev;
		}
		PointCollection pcnew;
		BayesPoint pt;
		double denom = 0;
		BayesPoint mean = BayesPoint(maxd);
		double meta = 0;
		double mphi = 0;
		double mtime = 0;
		for(int i = 0; i < npts; i ++){
			pt = pc->at(i);
			if(ws == 0) pt.SetWeight( 1.0 );
			//already e-weighted
			if(ws == 2) pt.SetWeight( log( w0 + (pc->at(i).w()/_gev)/E_tot ) );
			pcnew.AddPoint(pt);
		}
		//center at 0 - [-pi, pi]
		pcnew.Center();
		//calculate weighted mean - with phi wraparound
		mean.SetValue(pcnew.Centroid(0),0);
		mean.SetValue(pcnew.Centroid(1),1);
		mean.SetValue(pcnew.Centroid(2),2);
		double ent;
		double ent_pt, dd1, dd2;
		double pi = acos(-1);
		for(int d1 = 0; d1 < maxd; d1++){
			for(int d2 = d1; d2 < maxd; d2++){ 
				ent = 0;
				for(int i = 0; i < npts; i++){
					dd1 = pcnew.at(i).at(d1) - mean.at(d1);
					//phi wraparound
					if(d1 == 1){
						//if(dd1 > pi) dd1 = 2*pi - dd1;
						dd1 = acos(cos(dd1));
						if(dd1 > pi) cout << "dd1: " << dd1 << endl;
						
					}
					dd2 = pcnew.at(i).at(d2) - mean.at(d2);
					if(d2 == 1){
						//if(dd2 > pi) dd2 = 2*pi - dd2;
						dd2 = acos(cos(dd2));
						if(dd2 > pi) cout << "dd2: " << dd2 << endl;
					}
					ent_pt = pcnew.at(i).w() * (dd1)*(dd2) / pcnew.Sumw();
					//if(d1 == 1) ent_pt *= sqrt(phiCorrectionFactor);
					//if(d2 == 1) ent_pt *= sqrt(phiCorrectionFactor);
					ent += ent_pt;
				}
				outcov.SetEntry(ent,d1,d2);
				if(d1 != d2) outcov.SetEntry(ent,d2,d1);
			}
		}
		//eta time sign convention
		//if(mean.at(0) < 0){
		//	//time sign does NOT match eta sign
		//	//flip sign of eta-time entry
		//	outcov.SetEntry(-outcov.at(0,2),0,2);	
		//	outcov.SetEntry(-outcov.at(2,0),2,0);	
		//}



	}


	void FillCMSHists(const vector<Jet>& rhs, int id_idx){
			double cmsLogE_phi2D, cmsLogE_rot2D, cmsLogE_ec, cmsLogE_pc, cmsLogE_tc;
			double cmsLogE_ep_cov, cmsLogE_te_cov, cmsLogE_tp_cov, cmsLogE_e_var, cmsLogE_p_var, cmsLogE_t_var;
			double cmsLogE_ep_cov_unnorm, cmsLogE_te_cov_unnorm, cmsLogE_tp_cov_unnorm;
			double cmsLogE_majtime_cov, cmsLogE_mintime_cov, cmsLogE_majtime_cov_unnorm, cmsLogE_mintime_cov_unnorm;
			double cmsLogE_smaj, cmsLogE_smin, cmsNoE_smaj, cmsNoE_smin;
			
			double cmsNoE_ep_cov_unnorm, cmsNoE_te_cov_unnorm, cmsNoE_tp_cov_unnorm;
			double cmsNoE_phi2D, cmsNoE_rot2D, cmsNoE_ec, cmsNoE_pc, cmsNoE_tc;
			double cmsNoE_ep_cov, cmsNoE_te_cov, cmsNoE_tp_cov, cmsNoE_e_var, cmsNoE_p_var, cmsNoE_t_var;
			double cmsNoE_majtime_cov, cmsNoE_mintime_cov, cmsNoE_majtime_cov_unnorm, cmsNoE_mintime_cov_unnorm;

			Matrix space_mat = Matrix(2,2);
			Matrix rotmat2D = Matrix(3,3);		
			Matrix majminCovMat = Matrix(3,3);
			vector<double> eigenvals, eigenvals_space;
			vector<Matrix> eigenvecs, eigenvecs_space; 
			PointCollection majminpts;

			//create PointCollection of (ieta, iphi, time)
			PointCollection* ipts = new PointCollection();
			double eta, phi, time;
			double E_tot = 0;
			for(int i = 0; i < rhs.size(); i++){
				eta = rhs[i].eta();
				phi = rhs[i].phi();
				//phi wraparound
				if(phi > acos(-1)) phi -= 2*acos(-1);	
				time = rhs[i].time();

				BayesPoint pt({eta, phi, time});
				pt.SetWeight(rhs[i].E()*_gev);
				E_tot += rhs[i].E();
				ipts->AddPoint(pt);
			}


			//////////////////////calculate CMS Log energy weighting variables//////////////////////
			Matrix logEw = Matrix(3,3);
			//calculate covariances like in CMSSW
			MakeCovMat(ipts, logEw, weightScheme(2));	
			cmsLogE_e_var = sqrt(logEw.at(0,0));
			cmsLogE_p_var = sqrt(logEw.at(1,1));
			cmsLogE_t_var = sqrt(logEw.at(2,2));

			cmsLogE_ep_cov = CalcCov(logEw, 1, 0);
			cmsLogE_te_cov = CalcCov(logEw, 2, 0);
			cmsLogE_tp_cov = CalcCov(logEw, 2, 1);
			cmsLogE_ep_cov_unnorm = CalcCov(logEw, 1, 0, false);
			cmsLogE_te_cov_unnorm = CalcCov(logEw, 2, 0, false);
			cmsLogE_tp_cov_unnorm = CalcCov(logEw, 2, 1, false);
			

			Get2DMat(logEw,space_mat);
			space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
			cmsLogE_phi2D = PhiEll(space_mat);			
			cmsLogE_rot2D = Rotundity(space_mat);
			cmsLogE_smaj = eigenvals_space[1];		
			cmsLogE_smin = eigenvals_space[0];		
			
		
			//rotate points into maj/min axes
			Get2DRotationMatrix(eigenvecs_space,rotmat2D);
			RotatePoints(ipts, rotmat2D, majminpts);
			MakeCovMat(&majminpts, majminCovMat, weightScheme(2));
			//set time covariance from GMM
			majminCovMat.SetEntry(logEw.at(2,2),2,2);
			cmsLogE_majtime_cov = CalcCov(majminCovMat,2,0);
			cmsLogE_mintime_cov = CalcCov(majminCovMat,2,1);
			cmsLogE_majtime_cov_unnorm = CalcCov(majminCovMat,2,0,false);
			cmsLogE_mintime_cov_unnorm = CalcCov(majminCovMat,2,1,false);
			//for centers
			BayesPoint logEwCenter = GetCMSSWMean(ipts,true);
			cmsLogE_ec = logEwCenter.at(0);
			cmsLogE_pc = logEwCenter.at(1);
			cmsLogE_tc = logEwCenter.at(2);


			//////////logE weighted//////////
			_procCats[id_idx].hists1D[0][35]->Fill(cmsLogE_te_cov);
			_procCats[id_idx].hists1D[0][36]->Fill(cmsLogE_tp_cov);
			_procCats[id_idx].hists1D[0][37]->Fill(cmsLogE_ep_cov);
			_procCats[id_idx].hists1D[0][38]->Fill(cmsLogE_e_var);
			_procCats[id_idx].hists1D[0][39]->Fill(cmsLogE_p_var);
			_procCats[id_idx].hists1D[0][40]->Fill(cmsLogE_t_var);
			_procCats[id_idx].hists1D[0][50]->Fill(cmsLogE_te_cov_unnorm);
			_procCats[id_idx].hists1D[0][51]->Fill(cmsLogE_tp_cov_unnorm);
			_procCats[id_idx].hists1D[0][52]->Fill(cmsLogE_ep_cov_unnorm);
			_procCats[id_idx].hists1D[0][56]->Fill(cmsLogE_smaj);
			_procCats[id_idx].hists1D[0][57]->Fill(cmsLogE_smin);
			_procCats[id_idx].hists1D[0][67]->Fill(logEwCenter.at(2));
			_procCats[id_idx].hists1D[0][68]->Fill(logEwCenter.at(0));
			_procCats[id_idx].hists1D[0][69]->Fill(logEwCenter.at(1));
			_procCats[id_idx].hists1D[0][70]->Fill(cmsLogE_phi2D);
			_procCats[id_idx].hists1D[0][71]->Fill(cmsLogE_majtime_cov);
			_procCats[id_idx].hists1D[0][72]->Fill(cmsLogE_mintime_cov);
			_procCats[id_idx].hists1D[0][73]->Fill(cmsLogE_majtime_cov_unnorm);
			_procCats[id_idx].hists1D[0][74]->Fill(cmsLogE_mintime_cov_unnorm);
			_procCats[id_idx].hists1D[0][75]->Fill(cmsLogE_rot2D);
			if(E_tot <= 100){
				_procCats[id_idx].hists1D[0][125]->Fill(cmsLogE_rot2D);	
			}
			else if(E_tot > 100 && E_tot <= 200){
				_procCats[id_idx].hists1D[0][126]->Fill(cmsLogE_rot2D);	
			}
			else{
				_procCats[id_idx].hists1D[0][127]->Fill(cmsLogE_rot2D);	
			}
			//Quadrants in etaPhiCov-timeEtaCov space
			//Q1 = ep_cov < 0 and te_cov < 0 
			//Q2 = ep_cov < 0 and te_cov > 0 
			//Q3 = ep_cov > 0 and te_cov < 0 
			//Q4 = ep_cov > 0 and te_cov > 0  
			if(cmsLogE_ep_cov < 0){
				if(cmsLogE_te_cov < 0){
					_procCats[id_idx].hists1D[0][128]->Fill(cmsLogE_rot2D);
					_procCats[id_idx].hists1D[0][132]->Fill(cmsLogE_phi2D);
					_procCats[id_idx].hists1D[0][136]->Fill(cmsLogE_ep_cov);
					_procCats[id_idx].hists1D[0][140]->Fill(cmsLogE_te_cov);
				}	
				else{
					_procCats[id_idx].hists1D[0][129]->Fill(cmsLogE_rot2D);
					_procCats[id_idx].hists1D[0][133]->Fill(cmsLogE_phi2D);
					_procCats[id_idx].hists1D[0][137]->Fill(cmsLogE_ep_cov);
					_procCats[id_idx].hists1D[0][141]->Fill(cmsLogE_te_cov);
				}
			}
			else{
				if(cmsLogE_te_cov < 0){
					_procCats[id_idx].hists1D[0][130]->Fill(cmsLogE_rot2D);
					_procCats[id_idx].hists1D[0][134]->Fill(cmsLogE_phi2D);
					_procCats[id_idx].hists1D[0][138]->Fill(cmsLogE_ep_cov);
					_procCats[id_idx].hists1D[0][142]->Fill(cmsLogE_te_cov);
				}	
				else{
					_procCats[id_idx].hists1D[0][131]->Fill(cmsLogE_rot2D);
					_procCats[id_idx].hists1D[0][135]->Fill(cmsLogE_phi2D);
					_procCats[id_idx].hists1D[0][139]->Fill(cmsLogE_ep_cov);
					_procCats[id_idx].hists1D[0][143]->Fill(cmsLogE_te_cov);

				}
			}
			if((cmsLogE_phi2D < 0.2 && cmsLogE_phi2D > -0.2) || (fabs(cmsLogE_phi2D) < 0.2+acos(-1)/2. && fabs(cmsLogE_phi2D) > -0.2+acos(-1)/2.)){
				_procCats[id_idx].hists1D[0][144]->Fill(cmsLogE_ec);
				_procCats[id_idx].hists1D[0][146]->Fill(cmsLogE_pc);
				_procCats[id_idx].hists1D[0][148]->Fill(cmsLogE_tc);
				_procCats[id_idx].hists1D[0][150]->Fill(E_tot);
				_procCats[id_idx].hists1D[0][152]->Fill(1.);
				_procCats[id_idx].hists1D[0][154]->Fill(cmsLogE_e_var);
				_procCats[id_idx].hists1D[0][156]->Fill(cmsLogE_p_var);
				_procCats[id_idx].hists1D[0][158]->Fill(cmsLogE_t_var);
				_procCats[id_idx].hists1D[0][160]->Fill(cmsLogE_ep_cov);
				_procCats[id_idx].hists1D[0][162]->Fill(cmsLogE_te_cov);
				_procCats[id_idx].hists1D[0][164]->Fill(cmsLogE_tp_cov);
				_procCats[id_idx].hists1D[0][166]->Fill(cmsLogE_majtime_cov);
				_procCats[id_idx].hists1D[0][168]->Fill(cmsLogE_mintime_cov);
				_procCats[id_idx].hists1D[0][170]->Fill(cmsLogE_phi2D);
				_procCats[id_idx].hists1D[0][172]->Fill(cmsLogE_rot2D);
			}
		
		//168 - timeMin cov, phiE2D ~ 0 && phiE2D ~ pi/2
                //169 - timeMin cov, phiE2D !~ 0 && phiE2D !~ pi/2
                //170 - phiE2D, phiE2D ~ 0 && phiE2D ~ pi/2
                //171 - phiE2D, phiE2D !~ 0 && phiE2D !~ pi/2
                //172 - rot2D cov, phiE2D ~ 0 && phiE2D ~ pi/2
                //173 - rot2D cov, phiE2D !~ 0 && phiE2D !~ pi/2
		
			else{
				_procCats[id_idx].hists1D[0][145]->Fill(cmsLogE_ec);
				_procCats[id_idx].hists1D[0][147]->Fill(cmsLogE_pc);
				_procCats[id_idx].hists1D[0][149]->Fill(cmsLogE_tc);
				_procCats[id_idx].hists1D[0][151]->Fill(E_tot);
				_procCats[id_idx].hists1D[0][153]->Fill(1.);
				_procCats[id_idx].hists1D[0][155]->Fill(cmsLogE_e_var);
				_procCats[id_idx].hists1D[0][157]->Fill(cmsLogE_p_var);
				_procCats[id_idx].hists1D[0][159]->Fill(cmsLogE_t_var);
				_procCats[id_idx].hists1D[0][161]->Fill(cmsLogE_ep_cov);
				_procCats[id_idx].hists1D[0][163]->Fill(cmsLogE_te_cov);
				_procCats[id_idx].hists1D[0][165]->Fill(cmsLogE_tp_cov);
				_procCats[id_idx].hists1D[0][167]->Fill(cmsLogE_majtime_cov);
				_procCats[id_idx].hists1D[0][169]->Fill(cmsLogE_mintime_cov);
				_procCats[id_idx].hists1D[0][171]->Fill(cmsLogE_phi2D);
				_procCats[id_idx].hists1D[0][173]->Fill(cmsLogE_rot2D);
			}

			

			//////////////////////calculate CMS No energy weighting variables//////////////////////
			Matrix noEw = Matrix(3,3);
			//calculate covariances like in CMSSW
			MakeCovMat(ipts, noEw, weightScheme(0));	
			cmsNoE_e_var = sqrt(noEw.at(0,0));
			cmsNoE_p_var = sqrt(noEw.at(1,1));
			cmsNoE_t_var = sqrt(noEw.at(2,2));

			cmsNoE_ep_cov = CalcCov(noEw, 1, 0);
			cmsNoE_te_cov = CalcCov(noEw, 2, 0);
			cmsNoE_tp_cov = CalcCov(noEw, 2, 1);
			cmsNoE_ep_cov_unnorm = CalcCov(noEw, 1, 0, false);
			cmsNoE_te_cov_unnorm = CalcCov(noEw, 2, 0, false);
			cmsNoE_tp_cov_unnorm = CalcCov(noEw, 2, 1, false);
			
			Get2DMat(noEw,space_mat);
			space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
			cmsNoE_phi2D = PhiEll(space_mat);			
			cmsNoE_rot2D = Rotundity(space_mat);
			cmsNoE_smaj = eigenvals_space[1];		
			cmsNoE_smin = eigenvals_space[0];		
		
			//rotate points into maj/min axes
			Get2DRotationMatrix(eigenvecs_space,rotmat2D);
			RotatePoints(ipts, rotmat2D, majminpts);
			MakeCovMat(&majminpts, majminCovMat, weightScheme(0));
			//set time covariance from GMM
			majminCovMat.SetEntry(noEw.at(2,2),2,2);
			cmsNoE_majtime_cov = CalcCov(majminCovMat,2,0);
			cmsNoE_mintime_cov = CalcCov(majminCovMat,2,1);
			cmsNoE_majtime_cov_unnorm = CalcCov(majminCovMat,2,0,false);
			cmsNoE_mintime_cov_unnorm = CalcCov(majminCovMat,2,1,false);
			//for centers
			BayesPoint noEwCenter = GetCMSSWMean(ipts,false);
			cmsNoE_ec = noEwCenter.at(0);
			cmsNoE_pc = noEwCenter.at(1);
			cmsNoE_tc = noEwCenter.at(2);


			_procCats[id_idx].hists1D[0][41]->Fill(cmsNoE_te_cov);	
			_procCats[id_idx].hists1D[0][42]->Fill(cmsNoE_tp_cov);	
			_procCats[id_idx].hists1D[0][43]->Fill(cmsNoE_ep_cov);	
			_procCats[id_idx].hists1D[0][44]->Fill(cmsNoE_e_var);
			_procCats[id_idx].hists1D[0][45]->Fill(cmsNoE_p_var);
			_procCats[id_idx].hists1D[0][46]->Fill(cmsNoE_t_var);
			_procCats[id_idx].hists1D[0][47]->Fill(cmsNoE_te_cov_unnorm);
			_procCats[id_idx].hists1D[0][48]->Fill(cmsNoE_tp_cov_unnorm);
			_procCats[id_idx].hists1D[0][49]->Fill(cmsNoE_ep_cov_unnorm);
			_procCats[id_idx].hists1D[0][54]->Fill(cmsNoE_smaj);
			_procCats[id_idx].hists1D[0][55]->Fill(cmsNoE_smin);
			_procCats[id_idx].hists1D[0][58]->Fill(noEwCenter.at(2));
			_procCats[id_idx].hists1D[0][59]->Fill(noEwCenter.at(0));
			_procCats[id_idx].hists1D[0][60]->Fill(noEwCenter.at(1));
			_procCats[id_idx].hists1D[0][61]->Fill(cmsNoE_phi2D);
			_procCats[id_idx].hists1D[0][62]->Fill(cmsNoE_majtime_cov);
			_procCats[id_idx].hists1D[0][63]->Fill(cmsNoE_mintime_cov);
			_procCats[id_idx].hists1D[0][64]->Fill(cmsNoE_majtime_cov_unnorm);
			_procCats[id_idx].hists1D[0][65]->Fill(cmsNoE_mintime_cov_unnorm);
			_procCats[id_idx].hists1D[0][66]->Fill(cmsNoE_rot2D);
			if(E_tot <= 100){
				_procCats[id_idx].hists1D[0][174]->Fill(cmsNoE_rot2D);	
			}
			else if(E_tot > 100 && E_tot <= 200){
				_procCats[id_idx].hists1D[0][175]->Fill(cmsNoE_rot2D);	
			}
			else{
				_procCats[id_idx].hists1D[0][176]->Fill(cmsNoE_rot2D);	
			}
			//Quadrants in etaPhiCov-timeEtaCov space
			//Q1 = ep_cov < 0 and te_cov < 0 
			//Q2 = ep_cov < 0 and te_cov > 0 
			//Q3 = ep_cov > 0 and te_cov < 0 
			//Q4 = ep_cov > 0 and te_cov > 0  
			if(cmsNoE_ep_cov < 0){
				if(cmsNoE_te_cov < 0){
					_procCats[id_idx].hists1D[0][177]->Fill(cmsNoE_rot2D);
					_procCats[id_idx].hists1D[0][181]->Fill(cmsNoE_phi2D);
					_procCats[id_idx].hists1D[0][185]->Fill(cmsNoE_ep_cov);
					_procCats[id_idx].hists1D[0][189]->Fill(cmsNoE_te_cov);
				}	
				else{
					_procCats[id_idx].hists1D[0][178]->Fill(cmsNoE_rot2D);
					_procCats[id_idx].hists1D[0][182]->Fill(cmsNoE_phi2D);
					_procCats[id_idx].hists1D[0][186]->Fill(cmsNoE_ep_cov);
					_procCats[id_idx].hists1D[0][190]->Fill(cmsNoE_te_cov);
				}
			}
			else{
				if(cmsNoE_te_cov < 0){
					_procCats[id_idx].hists1D[0][179]->Fill(cmsNoE_rot2D);
					_procCats[id_idx].hists1D[0][183]->Fill(cmsNoE_phi2D);
					_procCats[id_idx].hists1D[0][187]->Fill(cmsNoE_ep_cov);
					_procCats[id_idx].hists1D[0][191]->Fill(cmsNoE_te_cov);
				}	
				else{
					_procCats[id_idx].hists1D[0][180]->Fill(cmsNoE_rot2D);
					_procCats[id_idx].hists1D[0][184]->Fill(cmsNoE_phi2D);
					_procCats[id_idx].hists1D[0][188]->Fill(cmsNoE_ep_cov);
					_procCats[id_idx].hists1D[0][192]->Fill(cmsNoE_te_cov);

				}
			}
			if((cmsNoE_phi2D < 0.2 && cmsNoE_phi2D > -0.2) || (fabs(cmsNoE_phi2D) < 0.2+acos(-1)/2. && fabs(cmsNoE_phi2D) > -0.2+acos(-1)/2.)){
				_procCats[id_idx].hists1D[0][193]->Fill(cmsNoE_ec);
				_procCats[id_idx].hists1D[0][195]->Fill(cmsNoE_pc);
				_procCats[id_idx].hists1D[0][197]->Fill(cmsNoE_tc);
				_procCats[id_idx].hists1D[0][199]->Fill(E_tot);
				_procCats[id_idx].hists1D[0][201]->Fill(1.);
				_procCats[id_idx].hists1D[0][203]->Fill(cmsNoE_e_var);
				_procCats[id_idx].hists1D[0][205]->Fill(cmsNoE_p_var);
				_procCats[id_idx].hists1D[0][207]->Fill(cmsNoE_t_var);
				_procCats[id_idx].hists1D[0][209]->Fill(cmsNoE_ep_cov);
				_procCats[id_idx].hists1D[0][211]->Fill(cmsNoE_te_cov);
				_procCats[id_idx].hists1D[0][213]->Fill(cmsNoE_tp_cov);
				_procCats[id_idx].hists1D[0][215]->Fill(cmsNoE_majtime_cov);
				_procCats[id_idx].hists1D[0][217]->Fill(cmsNoE_mintime_cov);
				_procCats[id_idx].hists1D[0][219]->Fill(cmsNoE_phi2D);
				_procCats[id_idx].hists1D[0][221]->Fill(cmsNoE_rot2D);
			}
			else{
				_procCats[id_idx].hists1D[0][194]->Fill(cmsNoE_ec);
				_procCats[id_idx].hists1D[0][196]->Fill(cmsNoE_pc);
				_procCats[id_idx].hists1D[0][198]->Fill(cmsNoE_tc);
				_procCats[id_idx].hists1D[0][200]->Fill(E_tot);
				_procCats[id_idx].hists1D[0][202]->Fill(1.);
				_procCats[id_idx].hists1D[0][204]->Fill(cmsNoE_e_var);
				_procCats[id_idx].hists1D[0][206]->Fill(cmsNoE_p_var);
				_procCats[id_idx].hists1D[0][208]->Fill(cmsNoE_t_var);
				_procCats[id_idx].hists1D[0][210]->Fill(cmsNoE_ep_cov);
				_procCats[id_idx].hists1D[0][212]->Fill(cmsNoE_te_cov);
				_procCats[id_idx].hists1D[0][214]->Fill(cmsNoE_tp_cov);
				_procCats[id_idx].hists1D[0][216]->Fill(cmsNoE_majtime_cov);
				_procCats[id_idx].hists1D[0][218]->Fill(cmsNoE_mintime_cov);
				_procCats[id_idx].hists1D[0][220]->Fill(cmsNoE_phi2D);
				_procCats[id_idx].hists1D[0][222]->Fill(cmsNoE_rot2D);
			}
		
		//2D histograms
		//////////logE weighted//////////
                //76 - logE rot 2D vs. etaphi cov
                _procCats[id_idx].hists2D[0][76]->Fill(cmsLogE_rot2D, cmsLogE_ep_cov);
		//77 - logE rot 2D vs. etaphi cov unnorm
                _procCats[id_idx].hists2D[0][77]->Fill(cmsLogE_rot2D, cmsLogE_ep_cov_unnorm);
		//78 - logE phi center vs. energy
                _procCats[id_idx].hists2D[0][78]->Fill(cmsLogE_pc, E_tot);
		//79 - logE rot2D vs phiE2D
                _procCats[id_idx].hists2D[0][79]->Fill(cmsLogE_rot2D, cmsLogE_phi2D);
		//80 - logE counts of etaphi cov vs timeeta cov counts
                _procCats[id_idx].hists2D[0][80]->Fill(cmsLogE_ep_cov, cmsLogE_te_cov);
		//81 - logE time center vs phiE2D
                _procCats[id_idx].hists2D[0][81]->Fill(cmsLogE_tc, cmsLogE_phi2D);
		if((cmsNoE_phi2D < 0.2 && cmsNoE_phi2D > -0.2) || (fabs(cmsNoE_phi2D) < 0.2+acos(-1)/2. && fabs(cmsNoE_phi2D) > -0.2+acos(-1)/2.)){
                	//82 - logE rot 2D vs. etaphi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                	_procCats[id_idx].hists2D[0][82]->Fill(cmsLogE_rot2D, cmsLogE_ep_cov);
	                //86 - logE rot 2D vs. timeeta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                	_procCats[id_idx].hists2D[0][86]->Fill(cmsLogE_rot2D, cmsLogE_te_cov);

		}
		else{
                	//83 - logE rot 2D vs. etaphi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                	_procCats[id_idx].hists2D[0][83]->Fill(cmsLogE_rot2D, cmsLogE_ep_cov);
	                //87 - logE rot 2D vs. timeeta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                	_procCats[id_idx].hists2D[0][84]->Fill(cmsLogE_rot2D, cmsLogE_te_cov);

		}
		//84 - logE energy vs phiE2D
                _procCats[id_idx].hists2D[0][84]->Fill(E_tot, cmsLogE_phi2D);
		//85 - logE timeEtaCov vs rotundity 2D
                _procCats[id_idx].hists2D[0][85]->Fill(cmsLogE_te_cov, cmsLogE_rot2D);

		

		//////////noE weighted//////////
                //88 - noE rot 2D vs. etaphi cov
                _procCats[id_idx].hists2D[0][88]->Fill(cmsNoE_rot2D, cmsNoE_ep_cov);
                //89 - noE rot 2D vs. etaphi cov unnorm
                _procCats[id_idx].hists2D[0][89]->Fill(cmsNoE_rot2D, cmsNoE_ep_cov_unnorm);
		//90 - noE phi center vs. energy
                _procCats[id_idx].hists2D[0][90]->Fill(cmsNoE_pc, E_tot);
		//91 - noE rot2D vs phiE2D
                _procCats[id_idx].hists2D[0][91]->Fill(cmsNoE_rot2D, cmsNoE_phi2D);
		//92 - noE counts of etaphi cov vs timeeta cov
                _procCats[id_idx].hists2D[0][92]->Fill(cmsNoE_ep_cov, cmsNoE_te_cov);
		//93 - noE time center vs phiE2D
                _procCats[id_idx].hists2D[0][93]->Fill(cmsNoE_tc, cmsNoE_phi2D);
		if((cmsNoE_phi2D < 0.2 && cmsNoE_phi2D > -0.2) || (fabs(cmsNoE_phi2D) < 0.2+acos(-1)/2. && fabs(cmsNoE_phi2D) > -0.2+acos(-1)/2.)){
                	//94 - noE rot 2D vs. etaphi cov, phiE2D ~ 0 && phiE2D ~ pi/2
                	_procCats[id_idx].hists2D[0][94]->Fill(cmsNoE_rot2D, cmsNoE_ep_cov);
                	//98 - noE rot 2D vs. timeeta cov, phiE2D ~ 0 && phiE2D ~ pi/2
                	_procCats[id_idx].hists2D[0][98]->Fill(cmsNoE_rot2D, cmsNoE_te_cov);
			//110 - log E etaPhiCov vs timeEtaCov, phiE2D ~ 0 && phiE2D ~ pi/2
			_procCats[id_idx].hists2D[0][110]->Fill(cmsLogE_ep_cov, cmsLogE_te_cov);
			//113 - no E etaPhiCov vs timeEtaCov, phiE2D ~ 0 && phiE2D ~ pi/2
			_procCats[id_idx].hists2D[0][113]->Fill(cmsNoE_ep_cov, cmsNoE_te_cov);

		}
		else{
                	//95 - noE rot 2D vs. etaphi cov, phiE2D !~ 0 && phiE2D !~ pi/2
                	_procCats[id_idx].hists2D[0][95]->Fill(cmsNoE_rot2D, cmsNoE_ep_cov);
                	//99 - noE rot 2D vs. timeeta cov, phiE2D !~ 0 && phiE2D !~ pi/2
                	_procCats[id_idx].hists2D[0][96]->Fill(cmsNoE_rot2D, cmsNoE_te_cov);
			//111 - log E etaPhiCov vs timeEtaCov, phiE2D !~ 0 && phiE2D !~ pi/2
			_procCats[id_idx].hists2D[0][111]->Fill(cmsLogE_ep_cov, cmsLogE_te_cov);
			//114 - no E etaPhiCov vs timeEtaCov, phiE2D !~ 0 && phiE2D !~ pi/2
			_procCats[id_idx].hists2D[0][114]->Fill(cmsNoE_ep_cov, cmsNoE_te_cov);

		}
		//96 - noE energy vs phiE2D
                _procCats[id_idx].hists2D[0][96]->Fill(E_tot, cmsNoE_phi2D);
		//97 - noE timeEtaCov vs rotundity 2D
                _procCats[id_idx].hists2D[0][97]->Fill(cmsNoE_te_cov, cmsNoE_rot2D);
		//109 - log E etaPhiCov vs timeEtaCov
		_procCats[id_idx].hists2D[0][109]->Fill(cmsLogE_ep_cov, cmsLogE_te_cov);
		//112 - no E etaPhiCov vs timeEtaCov
		_procCats[id_idx].hists2D[0][112]->Fill(cmsNoE_ep_cov, cmsNoE_te_cov);

		//127 - phoE vs logE smaj
		_procCats[id_idx].hists2D[0][127]->Fill(E_tot, cmsLogE_smaj);
		//128 - phoE vs logE smin
		_procCats[id_idx].hists2D[0][128]->Fill(E_tot, cmsLogE_smin);
		//129 - phoE vs logE sieie
		_procCats[id_idx].hists2D[0][129]->Fill(E_tot, cmsLogE_e_var);
		//130 - phoE vs logE sipip
		_procCats[id_idx].hists2D[0][130]->Fill(E_tot, cmsLogE_p_var);
		//131 - phoE vs logE etaphicov
		_procCats[id_idx].hists2D[0][131]->Fill(E_tot, cmsLogE_ep_cov);
		//132 - phoE vs logE timeetacov
		_procCats[id_idx].hists2D[0][132]->Fill(E_tot, cmsLogE_te_cov);
		//133 - nRhs vs logE smaj
		_procCats[id_idx].hists2D[0][133]->Fill(rhs.size(), cmsLogE_smaj);
		//134 - nRhs vs logE smin
		_procCats[id_idx].hists2D[0][134]->Fill(rhs.size(), cmsLogE_smin);
		//135 - nRhs vs logE sieie
		_procCats[id_idx].hists2D[0][135]->Fill(rhs.size(), cmsLogE_e_var);
		//136 - nRhs vs logE sipip
		_procCats[id_idx].hists2D[0][136]->Fill(rhs.size(), cmsLogE_p_var);
		//137 - phoE vs logE etaphicov
		_procCats[id_idx].hists2D[0][137]->Fill(rhs.size(), cmsLogE_ep_cov);
		//138 - phoE vs logE timeetacov
		_procCats[id_idx].hists2D[0][138]->Fill(rhs.size(), cmsLogE_te_cov);
		//139 - phoE vs nRhs
		_procCats[id_idx].hists2D[0][139]->Fill(E_tot, rhs.size());


	}




	double sqrtcov(double c){
		if(c > 0)
			return sqrt(c);
		else
			return -sqrt(-c);

	}


	double CalcCov(const Matrix& cov, int i, int j, bool norm = true){
		if(i >= cov.GetDims()[0]) return -999;
		if(j >= cov.GetDims()[1]) return -999;
		if(cov.empty()) return -999;
		double denom = 1;
		if(norm){
			denom = sqrt(cov.at(i,i)*cov.at(j,j));
		}
		//return sqrtcov(cov.at(i,j)/denom);
		return cov.at(i,j)/denom;
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

	double PhiEll(Matrix& inmat){
		vector<Matrix> eigenvecs;
		vector<double> eigenvals;
		inmat.eigenCalc(eigenvals, eigenvecs);
		int maxd = inmat.GetDims()[0] - 1;
		double v_x = eigenvecs[maxd].at(0,0);	
		double v_y = eigenvecs[maxd].at(1,0);	
		//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
		double phi = atan2( v_y , v_x );
		return phi;
	}

	void SpikeObs(PointCollection* pc, vector<double>& obs){
		obs.clear();
		pc->Sort();
		int npts = pc->GetNPoints();
		double wmax = pc->at(npts-1).w();
		BayesPoint xmax = pc->at(npts-1);
		//find 4 closest neighbors to xmax
		//map<double, Point> distToMax;
		double dist, deta, dphi;
		//distance threshold to include rh in sw+'
		double distThresh = 0.02;
		int nNeighbors = 0;
		double sumNeighbors = 0;
		for(int i = 0; i < npts; i++){
			//skip xmax
			if(xmax == pc->at(i)) continue;
			deta = pc->at(i).at(0) - xmax.at(0);
			dphi = pc->at(i).at(1) - xmax.at(1);
			dist = sqrt(deta*deta + dphi*dphi);
			if(dist < distThresh){
				sumNeighbors += pc->at(i).w();
				nNeighbors++;
			}
			//distToMax[dist] = points->at(i);
		}
		//swCP
		obs.push_back(sumNeighbors/wmax);

	}


	//potential BH variables
	void BeamHaloObs(PointCollection* pc, vector<double>& obs){
		obs.clear();	
		//find seed crystal -> largest weight
		pc->Sort();
		//seed crystal is last one
		BayesPoint seed = pc->at(pc->GetNPoints()-1);
		//ie looking at neighbor eta energy ratio + neighbor phi energy ratio
		
			// double ratio
		//ratio of center crystal to 2 neighbors in eta (smaller), phi (larger)
		//ratio of eta strips to surrounding eta strips in phi
		

	}
	
	void Get3DRotationMatrix(vector<Matrix> eigenvecs, Matrix& rotmat){
		if(rotmat.GetDims()[0] != 3) return;
		if(eigenvecs.size() != 3) return;
		rotmat.reset();
		rotmat.SetEntry(eigenvecs[2].at(0,0),0,0);
		rotmat.SetEntry(eigenvecs[2].at(1,0),1,0);
		rotmat.SetEntry(eigenvecs[2].at(2,0),2,0);

		rotmat.SetEntry(eigenvecs[1].at(0,0),0,1);
		rotmat.SetEntry(eigenvecs[1].at(1,0),1,1);
		rotmat.SetEntry(eigenvecs[1].at(2,0),2,1);

		rotmat.SetEntry(eigenvecs[0].at(0,0),0,2);
		rotmat.SetEntry(eigenvecs[0].at(1,0),1,2);
		rotmat.SetEntry(eigenvecs[0].at(2,0),2,2);
	}

	//is a 3D rotation matrix but only does 2D rotation
	void Get2DRotationMatrix(vector<Matrix> eigenvecs, Matrix& rotmat){
		if(rotmat.GetDims()[0] != 3) return;
		if(eigenvecs.size() != 2) return;
		rotmat.reset();
		rotmat.SetEntry(eigenvecs[1].at(0,0),0,0);
		rotmat.SetEntry(eigenvecs[1].at(1,0),1,0);
		rotmat.SetEntry(eigenvecs[0].at(0,0),0,1);
		rotmat.SetEntry(eigenvecs[0].at(1,0),1,1);
		rotmat.SetEntry(1.,2,2);
	}


	void RotatePoints(PointCollection* pc, const Matrix& rotmat, PointCollection& newpts){
		Matrix oldpt = Matrix(3,1);	
		Matrix newpt = Matrix(3,1);
		newpts.Clear();	
		//calculate points rotated into major/minor axes
		for(int i = 0; i < pc->GetNPoints(); i++){
			//put eta/phi coordinates into major/minor coordinates
			BayesPoint pt = pc->at(i);
			oldpt.SetEntry(pt.at(0),0,0);
			oldpt.SetEntry(pt.at(1),1,0);
			oldpt.SetEntry(pt.at(2),2,0);	

			//rotate by above transformation
			newpt.mult(rotmat,oldpt);
			//may need to preserve weight
			PointCollection tmp_newpts = newpt.MatToPoints();
			BayesPoint newpoint = tmp_newpts.at(0);
			newpoint.SetWeight(pt.w());
			newpts.AddPoint(newpoint);	
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


	double swissCross(const vector<Jet>& jets){
		//find seed crystal (highest weight, E = w*_gev)
		PointCollection sc(1); //save cmsswIds with associated weights for rec hit id
		JetPoint rh;
		for(int j = 0; j < jets.size(); j++){
			//should only have 1 rh per jet
			if(jets[j].GetNRecHits() > 1) continue;
			rh = jets[j].GetJetPoints()[0];
			BayesPoint pt(1);
			pt.SetValue(int(rh.rhId()),0);
			pt.SetWeight(rh.GetWeight());
			sc += pt;	
		}
		sc.Sort();
		double e1, e4;
		e1 = sc.at(sc.GetNPoints()-1).w()/_gev;
		int e1id = sc.at(sc.GetNPoints()-1).at(0);
		int e1_iphi = _detIDmap[e1id].i1;
		int e1_ieta = _detIDmap[e1id].i2;
		//find up, left, right, down crystals (e4)
		//eta on x-axis - shift in eta = left/right
		//phi on y-axis - shift in phi = up/down
		int e4upid = offsetBy(e1_ieta, e1_iphi, 0, 1);
		int e4downid = offsetBy(e1_ieta, e1_iphi, 0, -1);
		int e4rightid = offsetBy(e1_ieta, e1_iphi, 1, 0);
		int e4leftid = offsetBy(e1_ieta, e1_iphi, -1, 0);
		e4 = 0;
		for(int i = 0; i < sc.GetNPoints(); i++){
			if(sc.at(i).at(0) == e4upid && e4upid != -1)
				e4 += sc.at(i).w()/_gev;
			if(sc.at(i).at(0) == e4downid && e4downid != -1)
				e4 += sc.at(i).w()/_gev;
			if(sc.at(i).at(0) == e4leftid && e4leftid != -1)
				e4 += sc.at(i).w()/_gev;
			if(sc.at(i).at(0) == e4rightid && e4rightid != -1)
				e4 += sc.at(i).w()/_gev;
		}
		return 1 - e4/e1;

	}

	//from CMSSW: https://cmssdt.cern.ch/lxr/source/DataFormats/EcalDetId/src/EBDetId.cc
	UInt_t offsetBy(int ieta, int iphi, int nrStepsEta, int nrStepsPhi) const {
		int newEta = ieta + nrStepsEta;
		if (newEta * ieta <= 0) {
		  if (ieta < 0) {
		    newEta++;
		  } else if (ieta > 0) {
		    newEta--;
		  }
		}
		int newPhi = iphi + nrStepsPhi;
		while (newPhi > 360)
		  newPhi -= 360;
		while (newPhi <= 0)
		  newPhi += 360;
		pair<int, int> newpair = make_pair(newEta, newPhi);	
		if (validDetId(newEta, newPhi)) {
		  UInt_t id = _ietaiphiID.at(newpair);
		  return id;
		} else {
		  return -1;
		}
	}

	//from CMSSW: https://cmssdt.cern.ch/lxr/source/DataFormats/EcalDetId/interface/EBDetId.h
	/// check if a valid index combination
	static bool validDetId(int i, int j) {
		int max_ieta = 85;
		int max_iphi = 360;
		int min_ieta = 1;
		int min_iphi = 1;
		return i != 0 && (std::abs(i) <= max_ieta) && (j >= min_iphi) && (j <= max_iphi);
	}



	pair<double, double> iEtaiPhi2EtaPhi(int ieta, int iphi){
		//offset by 10 (+1 for starting at 1)?
		if(iphi > 11) iphi -= 11;
		else iphi += 349; 	
		double phinew = iphi*acos(-1)/180;
		double etanew = ieta*0.017453292519943295;
		return make_pair(etanew, phinew);
	}


	//assume centered at zero, t0 = center of plane section
	void EtaPhiPlaneSection(double erad, double prad, double trad, double t0, double& etaplanerad, double& phiplanerad){
		//for a horizontal planar section where z = d = t0, the eta and phi radii are orthogonal to the plane along their respective directions
		//the eta/phi plane radii depends on the distance to the plane (t0 = d), the radius in time (trad), scaled by the respective overall radius (erad, prad)
		//see https://en.wikipedia.org/wiki/Ellipsoid#Determining_the_ellipse_of_a_plane_section
		if(fabs(t0) > trad){
			etaplanerad = -1;
			phiplanerad = -1;
			return;
		}
		double rho = sqrt(1 - (t0/trad)*(t0/trad));
		etaplanerad = erad*rho;
		phiplanerad = prad*rho;
	}

	//rotations preserve length so plane = d -> plane' = d
	//but d gets scaled into Hessian normal form
	void RotateNormPlane(Matrix& rotcov, Matrix& plane, double& d){
		cout << "original plane" << endl; plane.Print();
		plane.mult(rotcov, plane);
		cout << "rotated plane" << endl; plane.Print();
		//put into Hessian normal form 
		//https://mathworld.wolfram.com/HessianNormalForm.html
		double dist = sqrt(plane.at(0,0)*plane.at(0,0) + plane.at(1,0)*plane.at(1,0) + plane.at(2,0)*plane.at(2,0));
		cout << "dist " << dist << endl;
		plane.mult(plane,1/dist);
		cout << "normal plane" << endl; plane.Print();
		d = d/dist;
	}


	//plane is plane_x*x + plane_y*y + plane_z*z = d
	//trad should be time-component of major axis
	void EtaPhiPlaneCov(PointCollection* pc, double trad, Matrix& plane, double d, Matrix& planecov){
		planecov.clear();
		//plane cannot be above/below time-component of major axis
		if(fabs(d) > trad){
			return;
		}
		planecov = Matrix(2,2);
		//smoothly weight pts according to distance to time slice
		double w = -1;
		double norm = 0;
		Matrix planemean = Matrix(pc->mean());
		planemean.SetEntry(pc->CircularMean(1),1,0);
		double etadist, phidist, planedist;	
		double plnorm = sqrt(plane.at(0,0)*plane.at(0,0) + plane.at(1,0)*plane.at(1,0) + plane.at(2,0)*plane.at(2,0));
		Matrix normplane = plane;
		normplane.mult(normplane,1/plnorm);
		for(int i = 0; i < pc->GetNPoints(); i++){
			//weight by inverse of distance to t0 plane
			//min distance of point to plane - distance vector bw point and plane is perpendicular to plane
			//use dot product to find projection of point onto plane
			//planedist = a*pt_x + b*pt_y + c*pt_z - d in Hesse normal
			//plane-point distance = ^n \dot ^plane + d 
			//https://mathworld.wolfram.com/Point-PlaneDistance.html
			planedist = fabs(normplane.at(0,0)*pc->at(i).at(0) + normplane.at(1,0)*pc->at(i).at(1) + normplane.at(2,0)*pc->at(i).at(2) - d);
			//cout << " planedist " << planedist << " for pt " << endl; pc->at(i).Print(); 
			w = pc->at(i).w();//1/planedist; 
			norm += w;

			etadist = pc->at(i).at(0) - planemean.at(0,0);
			phidist = pc->at(i).at(1) - planemean.at(1,0);
			planecov.SetEntry(w*etadist*etadist + planecov.at(0,0), 0, 0);
			planecov.SetEntry(w*etadist*phidist + planecov.at(0,1), 0, 1);
			planecov.SetEntry(w*phidist*etadist + planecov.at(1,0), 1, 0);
			planecov.SetEntry(w*phidist*phidist + planecov.at(1,1), 1, 1);

		}
		planecov.mult(planecov,1/norm);
	cout << "d " << d << " plane" << endl; plane.Print(); cout << "cov at plane" << endl; planecov.Print();
	} 



	//used for sig/bkg MVA, iso/!iso MVA
	int GetTrainingLabel(int nobj, int ncl, BasePDFMixture* gmm){
		//labels
		//unmatched = -1

		//sig vs bkg
		//signal = 0
		//iso bkg = 4

		//iso vs !iso
		//iso sig = 5
		//!iso bkg = 6

		//phys bkg vs det bkg - these are done in SuperClusterSkimmer
		//phys bkg = 1 
		//BH = 2
		//spike = 3
		
		double ec, pc, tc;
		auto params = gmm->GetLHPosteriorParameters(ncl);
		ec = params["mean"].at(0,0);
		pc = params["mean"].at(1,0);
		tc = params["mean"].at(2,0);

		vector<double> norms;
		gmm->GetNorms(norms);
		double Ek = norms[ncl]/_gev;

		Matrix cov = params["cov"];
		vector<Matrix> eigvecs;
		vector<double> eigvals;
		cov.eigenCalc(eigvals,eigvecs);
		double majLength = sqrt(eigvals[2]);
		if(eigvals[1] < 0) cout << "negative eigenvalue " << eigvals[1] << endl;
		double minLength; 
		if(eigvals[1] < 0) minLength = -sqrt(-eigvals[1]);
		else minLength = sqrt(eigvals[1]);	
		
		double rot3D = Rotundity(cov);
		Matrix space_mat(2,2);
		Get2DMat(cov,space_mat);
		double phi2D = PhiEll(space_mat);			
		double rot2D = Rotundity(space_mat);
		
		bool trksum, ecalrhsum, htowoverem, iso;	
		int label = -1;
		//signal
		if(!_data){
			//find photon associated with subcluster
			int phoidx = _base->SuperCluster_PhotonIndx->at(nobj);
			//matched to photon
			if(phoidx != -1){
                		trksum = _base->Photon_trkSumPtSolidConeDR04->at(phoidx) < 6.0;
                		ecalrhsum = _base->Photon_ecalRHSumEtConeDR04->at(phoidx) < 10.0;
                		htowoverem = _base->Photon_hadTowOverEM->at(phoidx) < 0.02;
                		iso = trksum && ecalrhsum && htowoverem;
				if(_base->Photon_genIdx->at(phoidx) != -1){
					int genidx = _base->Photon_genIdx->at(phoidx);
                			//needs to be isolated
					if(_isocuts){
						if(!iso) label = -1;
						else{
							//sig vs bkg
							//iso bkg = 0
							//obj selection for iso bkg
							if(_base->Gen_susId->at(genidx) == 22)
								label = 0;
							else{
								//sig vs bkg
								//iso bkg = 4
								//obj selection for iso bkg
								if(_isoBkgSel){
									if(_base->Photon_pt->at(nobj) < _minPhoPt_isoBkg) label = -1; //failed pho pt req for iso bkg
									//pt asymmetry bw photon + jet - put in if modelling bw data/MC is not good enough
									else{
										label = 4; //selection for iso bkg is on (event sel)
									}
								}
								else label = -1; //removal of GMSB !sig photons is done in data processing for NN
							}	
						}
					}
					//not applying isolation - use for iso network
					else{
						//photon from C2
						if(_base->Gen_susId->at(genidx) == 22)
							label = 0;
						//photon from hard subprocess - isolated sig 
						else if(_base->Gen_status->at(genidx) == 23)
							label = 5;
						else if(_base->Gen_motherIdx->at(genidx) != -1 && _base->Gen_status->at(_base->Gen_motherIdx->at(genidx)) == 23)
							label = 5;
						//photons from QCD are all nonisolated - !iso bkg
						else if(_oname.find("QCD") != string::npos)
							label = 6;
						else
							label = 1; //removal of GMSB !sig photons is done in data processing for NN	
					}
				}
				else //no gen match
					label = -1;

			}
			else
				label = -1;
		}
		//else in data - could be spikes or BH
		else{
			if(_isocuts){
				int phoidx = _base->SuperCluster_PhotonIndx->at(nobj);
				if(phoidx == -1){
					//not matched to a photon so cant be bkg
					label = -1;
				}
				else{
                			trksum = _base->Photon_trkSumPtSolidConeDR04->at(phoidx) < 6.0;
                			ecalrhsum = _base->Photon_ecalRHSumEtConeDR04->at(phoidx) < 10.0;
                			htowoverem = _base->Photon_hadTowOverEM->at(phoidx) < 0.02;
                			iso = trksum && ecalrhsum && htowoverem;
                			if(!iso) label = -1; //not isolated photon - won't make it into analysis anyway
					//sig vs bkg
					//iso bkg = 4
					//obj selection for iso bkg
					if(_isoBkgSel){
						if(_base->Photon_pt->at(nobj) < _minPhoPt_isoBkg) label = -1; //failed pho pt req for iso bkg
						//pt asymmetry bw photon + jet - put in if modelling bw data/MC is not good enough
						else label = 4; //selection for iso bkg is on (event sel)
					}
					else label = -1; 
					
		
				}
			}
			else{
				//no iso cuts -> photons for iso network -> can only be evaluated in MC gen-matching above (n/a data)
			}
		
		}
		//iso bkg CR distributions
		if(label == 4){
			//cout << "prompt Ek " << Ek << " norm " << norms[ncl] << " for subcl " << ncl << " of " << norms.size() << endl;
			for(int i = 0; i < (int)_procCats.size(); i++){
				_procCats[i].hists1D[0][239]->Fill(Ek, _weight);
				_procCats[i].hists1D[0][240]->Fill(sqrt(cov.at(0,0)), _weight);
				_procCats[i].hists1D[0][241]->Fill(sqrt(cov.at(1,1)), _weight);
				_procCats[i].hists1D[0][242]->Fill(cov.at(0,1), _weight);
				_procCats[i].hists1D[0][243]->Fill(cov.at(0,2), _weight);
				_procCats[i].hists1D[0][244]->Fill(cov.at(1,2), _weight);
				_procCats[i].hists1D[0][245]->Fill(majLength, _weight);
				_procCats[i].hists1D[0][246]->Fill(minLength, _weight);
				_procCats[i].hists1D[0][247]->Fill(phi2D, _weight);
				_procCats[i].hists1D[0][248]->Fill(rot2D, _weight);
				_procCats[i].hists1D[0][249]->Fill(rot3D, _weight);


			}
		}
		return label;
	}


	void SetMinPt_IsoBkg(double p){ _minPhoPt_isoBkg = p; _prod->SetMinPt(p); }
	void SetMinHt_IsoBkg(double p){ _minHt_isoBkg = p; }
	void SetMinJetPt_IsoBkg(double p){ _minJetPt_isoBkg = p; _jetprod->SetMinPt(p); }
	void SetMaxMet_IsoBkg(double p){ _maxMet_isoBkg = p; }
	void SetIsoBkgSel(bool b){ _isoBkgSel = b;}

	private:
		JetProducer* _jetprod;
		double _minPhoPt_isoBkg, _minHt_isoBkg, _minJetPt_isoBkg, _maxMet_isoBkg;
		bool _isoBkgSel; 

};
#endif
