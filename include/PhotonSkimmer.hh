#ifndef PHOTONSKIMMER_HH
#define PHOTONSKIMMER_HH

#include "JetPoint.hh"
#include <TFile.h>
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include "PhotonProducer.hh"
#include "TSystem.h"
#include <math.h>

using plotCat = BaseSkimmer::plotCat;
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
		};
		virtual ~PhotonSkimmer(){ };

		//get rechits from file to cluster
		PhotonSkimmer(TFile* file) : BaseSkimmer(file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time
			_prod = new PhotonProducer(file);
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
			objE->SetTitle("phoE");
			objE->SetName("phoE");
			
			objE_clusterE->SetTitle("phoE_clusterE");
			objE_clusterE->SetName("phoE_clusterE");

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
			_hists1D.push_back(cmsLogE_timeeta_cov);
			_hists1D.push_back(cmsLogE_timephi_cov);
			_hists1D.push_back(cmsLogE_etaphi_cov);
                	_hists1D.push_back(cmsLogE_etaSig);
                	_hists1D.push_back(cmsLogE_phiSig);
                	_hists1D.push_back(cmsLogE_timeSig);
			_hists1D.push_back(cmsNoE_timeeta_cov);
			_hists1D.push_back(cmsNoE_timephi_cov);
			_hists1D.push_back(cmsNoE_etaphi_cov);
                	_hists1D.push_back(cmsNoE_etaSig);
                	_hists1D.push_back(cmsNoE_phiSig);
                	_hists1D.push_back(cmsNoE_timeSig);
			_hists1D.push_back(cmsNoE_timeeta_covUnnorm);
			_hists1D.push_back(cmsNoE_timephi_covUnnorm);
			_hists1D.push_back(cmsNoE_etaphi_covUnnorm);
			_hists1D.push_back(cmsLogE_timeeta_covUnnorm);
			_hists1D.push_back(cmsLogE_timephi_covUnnorm);
			_hists1D.push_back(cmsLogE_etaphi_covUnnorm);
			_hists1D.push_back(etaphi_covUnnorm);
			_hists1D.push_back(cmsNoE_smaj);		
			_hists1D.push_back(cmsNoE_smin);		
			_hists1D.push_back(cmsLogE_smaj);		
			_hists1D.push_back(cmsLogE_smin);		
			_hists1D.push_back(cmsNoE_time_center);
			_hists1D.push_back(cmsNoE_eta_center);
			_hists1D.push_back(cmsNoE_phi_center);
			_hists1D.push_back(cmsNoE_phiEll);		
			_hists1D.push_back(cmsNoE_timeSmaj_cov);
			_hists1D.push_back(cmsNoE_timeSmin_cov);
			_hists1D.push_back(cmsNoE_timeSmaj_covUnnorm);
			_hists1D.push_back(cmsNoE_timeSmin_covUnnorm);
			_hists1D.push_back(cmsNoE_rotundity_2D);
			_hists1D.push_back(cmsLogE_time_center);
			_hists1D.push_back(cmsLogE_eta_center);
			_hists1D.push_back(cmsLogE_phi_center);
			_hists1D.push_back(cmsLogE_phiEll);	
			_hists1D.push_back(cmsLogE_timeSmaj_cov);
			_hists1D.push_back(cmsLogE_timeSmin_cov);
			_hists1D.push_back(cmsLogE_timeSmaj_covUnnorm);
			_hists1D.push_back(cmsLogE_timeSmin_covUnnorm);
			_hists1D.push_back(cmsLogE_rotundity_2D);
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
                        _hists2D.push_back(timeSig_timeCenter);
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
			_hists2D.push_back(timeMajCov_timeCenter);
			_hists2D.push_back(timeMinCov_timeCenter);
			_hists2D.push_back(etaPhiCov_phiEll2D);
			_hists2D.push_back(timeMajCov_phiEll2D);
			_hists2D.push_back(timeMinCov_phiEll2D);
                	_hists2D.push_back(rot2D_timeMajCovUnnorm);
                	_hists2D.push_back(rot2D_timeMinCovUnnorm);
                	_hists2D.push_back(timeSig_timeEtaCov);
			_hists2D.push_back(etaPhiCov_timeCenter);
			_hists2D.push_back(timeEtaCov_timeCenter);
			_hists2D.push_back(timePhiCov_timeCenter);
			_hists2D.push_back(cmsSmaj_cmsTimeSig);
			_hists2D.push_back(cmsSmin_cmsTimeSig);
			_hists2D.push_back(cmsLogESmaj_cmsLogETimeSig);
			_hists2D.push_back(cmsLogESmin_cmsLogETimeSig);
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
		TH1D* rotundity_2D = new TH1D("rotundity_2D","rotundity_2D",20,0.4,1.1);
		//14 - velocity = z/r*k for k transfer factor to velocity units
		TH1D* velocity = new TH1D("velocity","velocity",31,-1.,30.);
		//15 - ratio of 2D eigenvals
		TH1D* eigen2D_ratio = new TH1D("eigen2D_ratio","eigen2D_ratio",50,0.,1.);
		//16 - eta sigma	
                TH1D* etaSig = new TH1D("etaSig","etaSig",25,0.01, 0.09);
		//17 - phi sigma	
                TH1D* phiSig = new TH1D("phiSig","phiSig",25,0.01,0.09);
		//18 - time sigma	
                TH1D* timeSig = new TH1D("timeSig","timeSig",25,0,4.);
		//19 - fraction of energy in cluster
		TH1D* fracE = new TH1D("fracE","fracE",50,0.,1.1);
		//20 - azimuth angle in 2D
		TH1D* phiEll = new TH1D("phi_ell","phi_ell",50,-3.1,3.1);		
		//21 - eta sigma for positive eta clusters
                TH1D* etaSig_pos = new TH1D("etaSig_pos","etaSig_pos",25,0.01, 0.09);
		//22 - eta sigma for negative eta clusters
                TH1D* etaSig_neg = new TH1D("etaSig_neg","etaSig_neg",25,0.01, 0.09);
		//23 - normalized covariance - eta/phi
		TH1D* etaphi_cov = new TH1D("etaphi_cov","etaphi_cov",25,-1.,1.);
		//24 - normalized covariance - time/eta
		TH1D* timeeta_cov = new TH1D("timeeta_cov","timeeta_cov",25,-1.,1.);
		//25 - normalized covariance - time/phi
		TH1D* timephi_cov = new TH1D("timephi_cov","timephi_cov",25,-1.,1.);
		//26 - normalized covariance - time/major axis
		TH1D* timemaj_cov = new TH1D("timemaj_cov","timemaj_cov",25,-5.,5.);
		//27 - normalized covariance - time/minor axis
		TH1D* timemin_cov = new TH1D("timemin_cov","timemin_cov",25,-5.,5.);
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
		TH1D* cmsLogE_timeeta_cov = new TH1D("cmsLogE_timeeta_cov","cmsLogE_timeeta_cov",25,-1.,1.);
		//36 - logE CMS normalized covariance - time/phi
		TH1D* cmsLogE_timephi_cov = new TH1D("cmsLogE_timephi_cov","cmsLogE_timephi_cov",25,-1.,1.);
		//37 - logE CMS normalized covariance - eta/phi
		TH1D* cmsLogE_etaphi_cov = new TH1D("cmsLogE_etaphi_cov","cmsLogE_etaphi_cov",25,-50.,50.);
		//38 - logE CMS eta sigma	
                TH1D* cmsLogE_etaSig = new TH1D("cmsLogE_etaSig","cmsLogE_etaSig",25,0.01, 0.09);
		//39 - logE CMS phi sigma	
                TH1D* cmsLogE_phiSig = new TH1D("cmsLogE_phiSig","cmsLogE_phiSig",25,0.01,0.09);
		//40 - logE CMS time sigma	
                TH1D* cmsLogE_timeSig = new TH1D("cmsLogE_timeSig","cmsLogE_timeSig",25,0,10.);
		//41 - CMS normalized covariance - time/eta
		TH1D* cmsNoE_timeeta_cov = new TH1D("cmsNoE_timeeta_cov","cmsNoE_timeeta_cov",25,-1.,1.);
		//42 - CMS normalized covariance - time/phi
		TH1D* cmsNoE_timephi_cov = new TH1D("cmsNoE_timephi_cov","cmsNoE_timephi_cov",25,-1.,1.);		
		//43 - CMS normalized covariance - eta/phi
		TH1D* cmsNoE_etaphi_cov = new TH1D("cmsNoE_etaphi_cov","cmsNoE_etaphi_cov",25,-1.,1.);
		//44 - CMS eta sigma	
                TH1D* cmsNoE_etaSig = new TH1D("cmsNoE_etaSig","cmsNoE_etaSig",25,0.01, 0.09);
		//45 - CMS phi sigma	
                TH1D* cmsNoE_phiSig = new TH1D("cmsNoE_phiSig","cmsNoE_phiSig",25,0.01,0.09);
		//46 - CMS time sigma	
                TH1D* cmsNoE_timeSig = new TH1D("cmsNoE_timeSig","cmsNoE_timeSig",25,0,10.);
		//47 - CMS unnormalized covariance - time/eta
		TH1D* cmsNoE_timeeta_covUnnorm = new TH1D("cmsNoE_timeeta_covUnnorm","cmsNoE_timeeta_covUnnorm",25,-1.,1.);
		//48 - CMS unnormalized covariance - time/phi
		TH1D* cmsNoE_timephi_covUnnorm = new TH1D("cmsNoE_timephi_covUnnorm","cmsNoE_timephi_covUnnorm",25,-1.,1.);
		//49 - CMS unnormalized covariance - eta/phi
		TH1D* cmsNoE_etaphi_covUnnorm = new TH1D("cmsNoE_etaphi_covUnnorm","cmsNoE_etaphi_covUnnorm",25,-1.,1.);
		//50 - logE CMS unnormalized covariance - time/eta
		TH1D* cmsLogE_timeeta_covUnnorm = new TH1D("cmsLogE_timeeta_covUnnorm","cmsLogE_timeeta_covUnnorm",25,-1.,1.);
		//51 - logE CMS unnormalized covariance - time/phi
		TH1D* cmsLogE_timephi_covUnnorm = new TH1D("cmsLogE_timephi_covUnnorm","cmsLogE_timephi_covUnnorm",25,-1.,1.);
		//52 - logE CMS normalized covariance - eta/phi
		TH1D* cmsLogE_etaphi_covUnnorm = new TH1D("cmsLogE_etaphi_covUnnorm","cmsLogE_etaphi_covUnnorm",25,-1.,1.);
		//53 - unnormalized covariance - eta/phi
		TH1D* etaphi_covUnnorm = new TH1D("etaphi_covUnnorm","etaphi_covUnnorm",25,-0.1,0.1);
		//54 - CMS smaj
		TH1D* cmsNoE_smaj = new TH1D("cmsNoE_smaj","cmsNoE_smaj",25,0.,0.004);		
		//55 - CMS smin
		TH1D* cmsNoE_smin = new TH1D("cmsNoE_smin","cmsNoE_smin",25,0.,0.004);		
		//56 - CMS logE smaj
		TH1D* cmsLogE_smaj = new TH1D("cmsLogE_smaj","cmsLogE_smaj",25,0.,0.004);		
		//57 - CMS logE smin
		TH1D* cmsLogE_smin = new TH1D("cmsLogE_smin","cmsLogE_smin",25,0.,0.004);		
		//58 - CMS time center
		TH1D* cmsNoE_time_center = new TH1D("cmsNoE_time_center","cmsNoE_time_center",50,-20,20);
		//59 - CMS eta center
		TH1D* cmsNoE_eta_center = new TH1D("cmsNoE_eta_center","cmsNoE_eta_center",50,-3.5,3.5);
		//60 - CMS phi center
		TH1D* cmsNoE_phi_center = new TH1D("cmsNoE_phi_center","cmsNoE_phi_center",50,-0.1,6.3);
		//61 - CMS ellipsoid angle (phi2D)
		TH1D* cmsNoE_phiEll = new TH1D("cmsNoE_phi_ell","cmsNoE_phi_ell",50,-0.1,3.);		
		//62 - CMS time-smaj covariance
		TH1D* cmsNoE_timeSmaj_cov = new TH1D("cmsNoE_timeSmaj_cov","cmsNoE_timeSmaj_cov",25,-3.,3.);
		//63 - CMS time-smin covariance
		TH1D* cmsNoE_timeSmin_cov = new TH1D("cmsNoE_timeSmin_cov","cmsNoE_timeSmin_cov",25,-3.,3.);
		//64 - CMS time-smaj covariance unnormalized 
		TH1D* cmsNoE_timeSmaj_covUnnorm = new TH1D("cmsNoE_timeSmaj_covUnnorm","cmsNoE_timeSmaj_covUnnorm",25,-0.5,0.5);
		//65 - CMS time-smin covariance unnormalized 
		TH1D* cmsNoE_timeSmin_covUnnorm = new TH1D("cmsNoE_timeSmin_covUnnorm","cmsNoE_timeSmin_covUnnorm",25,-0.5,0.5);
		//66 - CMS rotundity 2D
		TH1D* cmsNoE_rotundity_2D = new TH1D("cmsNoE_rotundity_2D","cmsNoE_rotundity_2D",20,0.4,1.1);
		//67 - CMS logE time center
		TH1D* cmsLogE_time_center = new TH1D("cmsLogE_time_center","cmsLogE_time_center",50,-20,20);
		//68 - CMS logE eta center
		TH1D* cmsLogE_eta_center = new TH1D("cmsLogE_eta_center","cmsLogE_eta_center",50,-3.5,3.5);
		//69 - CMS logE phi center
		TH1D* cmsLogE_phi_center = new TH1D("cmsLogE_phi_center","cmsLogE_phi_center",50,-0.1,6.3);
		//70 - CMS logE ellipsoid angle (phi2D)
		TH1D* cmsLogE_phiEll = new TH1D("cmsLogE_phi_ell","cmsLogE_phi_ell",50,-0.1,3.);	
		//71 - CMS logE time-smaj covariance
		TH1D* cmsLogE_timeSmaj_cov = new TH1D("cmsLogE_timeSmaj_cov","cmsLogE_timeSmaj_cov",25,-3.,3.);
		//72 - CMS logE time-smin covariance
		TH1D* cmsLogE_timeSmin_cov = new TH1D("cmsLogE_timeSmin_cov","cmsLogE_timeSmin_cov",25,-3.,3.);
		//73 - CMS logE time-smaj covariance unnormalized 
		TH1D* cmsLogE_timeSmaj_covUnnorm = new TH1D("cmsLogE_timeSmaj_covUnnorm","cmsLogE_timeSmaj_covUnnorm",25,-0.5,0.5);
		//74 - CMS logE time-smin covariance unnormalized 
		TH1D* cmsLogE_timeSmin_covUnnorm = new TH1D("cmsLogE_timeSmin_covUnnorm","cmsLogE_timeSmin_covUnnorm",25,-0.5,0.5);
		//75 - CMS logE rotundity 2D
		TH1D* cmsLogE_rotundity_2D = new TH1D("cmsLogE_rotundity_2D","cmsLogE_rotundity_2D",20,0.4,1.1);
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
		TH2D* etaSig_phiSig = new TH2D("etaSig_phiSig","etaSig_phiSig;etaSig;phiSig",25,0.1,0.09,25,0.01,0.09);
                //13 - time sigma v phi sigma
                TH2D* timeSig_etaSig = new TH2D("timeSig_etaSig","timeSig_etaSig;timeSig;etaSig",25,0,5.,25,0.01,0.09);
                //14 - time sigma v phi sigma
                TH2D* timeSig_phiSig = new TH2D("timeSig_phiSig","timeSig_phiSig;timeSig;phiSig",25,0,5.,25,0.01,0.09);
		//15 - fraction of energy in subcluster vs mm coeff of subcluster
		TH2D* fracE_mmcoeff = new TH2D("fracE_mmcoeff","fracE_mmcoeff;fracE;mmcoeff",20,0.,1.,20,0,1.1);
		//16 - number of subclusters vs fraction of energy in particular subcluster (really only applicable to lead subcluster)
		TH2D* nsubcl_fracE = new TH2D("nsubcl_fracE","nsubcl_fracE;nSubClusters;fracE",10,0,10.,20,0,1.1);
		//17 - time sigma vs time center
		TH2D* timeSig_timeCenter = new TH2D("timeSig_timeCenter","timeSig_timeCenter;timeSig;time center",25,0,1,25,-5,15);
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
                TH2D* phiEll2D_timeMajCov = new TH2D("azAngle2D_timeMajCov","azAngle2D_timeMajCov;azAngle2D;timeMajCov",25,0,3.5,25,-0.5,0.5);
                //26 - az 2D angle vs. TimeMinCov
                TH2D* phiEll2D_timeMinCov = new TH2D("azAngle2D_timeMinCov","azAngle2D_timeMinCov;azAngle2D;timeMinCov",25,0,3.5,25,-0.5,0.5);
                //27 - rot 2D vs. TimeMajCov
                TH2D* rot2D_timeMajCov = new TH2D("rot2D_timeMajCov","rot2D_timeMajCov;rot2D;timeMajCov",25,0.4,1.1,25,-0.5,0.5);
                //28 - rot 2D angle vs. TimeMinCov
                TH2D* rot2D_timeMinCov = new TH2D("rot2D_timeMinCov","rot2D_timeMinCov;rot2D;timeMinCov",25,0.4,1.1,25,-0.5,0.5);
                //29 - rot 3D vs. TimeMajCov
                TH2D* rot3D_timeMajCov = new TH2D("rot3D_timeMajCov","rot3D_timeMajCov;rot3D;timeMajCov",25,0.8,1.1,25,-0.5,0.5);
                //30 - rot 3D angle vs. TimeMinCov
                TH2D* rot3D_timeMinCov = new TH2D("rot3D_timeMinCov","rot3D_timeMinCov;rot3D;timeMinCov",25,0.8,1.1,25,-0.5,0.5);
		//31 - timemaj cov vs time center
		TH2D* timeMajCov_timeCenter = new TH2D("timeMajCov_timeCenter","timeMajCov_timeCenter;timeMajCov;time center",25,-0.5,0.5,25,-0.1,0.1);
		//32 - timemin cov vs time center
		TH2D* timeMinCov_timeCenter = new TH2D("timeMinCov_timeCenter","timeMinCov_timeCenter;timeMinCov;time center",25,-0.5,0.5,25,-0.1,0.1);
		//33 - eta-phi covariance vs azimuth ellipsoid angle
		TH2D* etaPhiCov_phiEll2D = new TH2D("etaPhiCov_phiE2D","etaPhiCov_phiE2D;etaPhiCov;ellipsoid phi 2D",25,-10.,10.,25,0,3.5);
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
		TH2D* etaPhiCov_timeCenter = new TH2D("etaPhiCov_timeCenter","etaPhiCov_timeCenter;etaPhiCov;time center",25,-1,1,25,-5,15);
		//40 - timeeta cov vs time center
		TH2D* timeEtaCov_timeCenter = new TH2D("timeEtaCov_timeCenter","timeEtaCov_timeCenter;timeEtaCov;time center",25,-1,1,25,-5,15);
		//41 - timeeta cov vs time center
		TH2D* timePhiCov_timeCenter = new TH2D("timePhiCov_timeCenter","timePhiCov_timeCenter;timePhiCov;time center",25,-1,1,25,-5,15);
		//42 - CMS smaj vs cms sigma_t
		TH2D* cmsSmaj_cmsTimeSig = new TH2D("cmsNoESmaj_cmsNoETimeSig","cmsNoESmaj_cmsNoETimeSig;cmsNoE_Smaj;cmsNoE_timeSig",25,0.,0.01,25,0,10);
		//43 - CMS smin vs cms sigma_t
		TH2D* cmsSmin_cmsTimeSig = new TH2D("cmsNoESmin_cmsNoETimeSig","cmsNoESmin_cmsNoETimeSig;cmsNoE_Smin;cmsNoE_timeSig",25,0.,0.01,25,0,10);
		//44 - CMS logE smaj vs cms sigma_t
		TH2D* cmsLogESmaj_cmsLogETimeSig = new TH2D("cmsLogESmaj_cmsLogETimeSig","cmsLogESmaj_cmsLogETimeSig;cmsLogE_Smaj;cmsLogE_timeSig",25,0.,0.01,25,0,10);
		//45 - CMS logE smin vs cms sigma_t
		TH2D* cmsLogESmin_cmsLogETimeSig = new TH2D("cmsLogESmin_cmsLogETimeSig","cmsLogESmin_cmsLogETimeSig;cmsLogE_Smin;cmsLogE_timeSig",25,0.,0.01,25,0,10);
		//46 - etaphi cov vs timeeta cov
		TH2D* etaPhiCov_timeEtaCov = new TH2D("etaPhiCov_timeEtaCov","etaPhiCov_timeEtaCov;etaPhiCov;timeEtaCov",25,-1,1,25,-1,1);
		//47 - timeeta cov vs timephi cov
		TH2D* timeEtaCov_timePhiCov = new TH2D("timeEtaCov_timePhiCov","timeEtaCov_timePhiCov;timeEtaCov;timePhiCov",25,-1,1,25,-1,1);
		//48 - etaphi cov vs timephi cov
		TH2D* etaPhiCov_timePhiCov = new TH2D("etaPhiCov_timePhiCov","etaPhiCov_timePhiCov;etaPhiCov;timePhiCov",25,-1,1,25,-1,1);
                //49 - etaphi cov vs. TimeMajCov
                TH2D* etaPhiCov_timeMajCov = new TH2D("etaPhiCov_timeMajCov","etaPhiCov_timeMajCov;etaPhiCov;timeMajCov",25,-1,1.,25,-0.5,0.5);
                //50 - etaphi cov vs. TimeMinCov
                TH2D* etaPhiCov_timeMinCov = new TH2D("etaPhiCov_timeMinCov","etaPhiCov_timeMinCov;etaPhiCov;timeMinCov",25,-1,1.,25,-0.5,0.5);
                //51 - timephi cov vs. TimeMajCov
                TH2D* timePhiCov_timeMajCov = new TH2D("timePhiCov_timeMajCov","timePhiCov_timeMajCov;timePhiCov;timeMajCov",25,-1,1.,25,-0.5,0.5);
                //52 - timephi cov vs. TimeMinCov
                TH2D* timePhiCov_timeMinCov = new TH2D("timePhiCov_timeMinCov","timePhiCov_timeMinCov;timePhiCov;timeMinCov",25,-1,1.,25,-0.5,0.5);
                //53 - timeeta cov vs. TimeMajCov
                TH2D* timeEtaCov_timeMajCov = new TH2D("timeEtaCov_timeMajCov","timeEtaCov_timeMajCov;timeEtaCov;timeMajCov",25,-1,1.,25,-0.5,0.5);
                //54 - timeeta cov vs. TimeMinCov
                TH2D* timeEtaCov_timeMinCov = new TH2D("timeEtaCov_timeMinCov","timeEtaCov_timeMinCov;timeEtaCov;timeMinCov",25,-1,1.,25,-0.5,0.5);
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
		double _thresh, _alpha, _emAlpha;

		vector<plotCat> plotCats;


		enum weightScheme{
			noWeight = 0,
			Eweight = 1,
			logEweight = 2
		};
		

		void MakeIDHists(string sample){
			//total
			plotCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			plotCats.push_back(tot);	
			
			if(sample.find("GMSB") != string::npos){
				//notSunm
				plotCat notSunm(_hists1D, _hists2D, "notSunm","notSunm");
				//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
				notSunm.ids = {29, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
				plotCats.push_back(notSunm);
				
				//signal
				plotCat sig(_hists1D, _hists2D, "chiGam","#Chi^{0} #rightarrow #gamma");
				sig.ids = {22};
				plotCats.push_back(sig);
			}
			else if(sample.find("JetHT") != string::npos){
				//data
				plotCat jetht(_hists1D, _hists2D, "JetHT", "JetHT");
				jetht.ids = {-999};
				plotCats.push_back(jetht);
			}
			else if(sample.find("GJets") != string::npos){
				//data
				plotCat gjets(_hists1D, _hists2D, "GJets", "GJets");
				gjets.ids = {-999};
				plotCats.push_back(gjets);
			}
			else return;

		}
	


		void WritePlotCat1D(TFile* ofile, const plotCat& pc){
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
		void WritePlotCat2D(TFile* ofile, const plotCat& pc){
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



		void WritePlotCatStack(TFile* ofile, const vector<plotCat>& pcs){
			ofile->cd();
			string name;
			//number of histogram categories (ie leading, !leading, etc)
			int nhistCats = pcs[0].hists1D.size();
			//write 1D hists
			//variables
			for(int j = 0; j < (int)_hists1D.size(); j++){
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
			//			cout << "		adding proc " << pcs[k].plotName << " to plot with hist " << pcs[k].hists1D[i][j]->GetName() << endl;
						pcs[k].hists1D[i][j]->SetTitle(pcs[k].plotName.c_str());
						pcs[k].hists1D[i][j]->Write();
					}
				}
			}
			
		}


		void WriteHists(TFile* ofile){
			//normalize histograms	
			for(int i = 0; i < (int)plotCats.size(); i++){
				for(int j = 0; j < plotCats[i].hists1D.size(); j++){
					//relative fraction histograms
					//nSubClusters
					//plotCats[i].hists1D[j][0]->Scale(1./plotCats[i].hists1D[j][0]->Integral());
					//ellipsoid center coordinates
					plotCats[i].hists1D[j][1]->Scale(1./plotCats[i].hists1D[j][1]->Integral());
					plotCats[i].hists1D[j][2]->Scale(1./plotCats[i].hists1D[j][2]->Integral());
					plotCats[i].hists1D[j][3]->Scale(1./plotCats[i].hists1D[j][3]->Integral());
					//theta + azimuthal angles
					plotCats[i].hists1D[j][7]->Scale(1./plotCats[i].hists1D[j][7]->Integral());
					plotCats[i].hists1D[j][8]->Scale(1./plotCats[i].hists1D[j][8]->Integral());
					//rotundity
					plotCats[i].hists1D[j][10]->Scale(1./plotCats[i].hists1D[j][10]->Integral());
					plotCats[i].hists1D[j][11]->Scale(1./plotCats[i].hists1D[j][11]->Integral());
					//velocity
					plotCats[i].hists1D[j][14]->Scale(1./plotCats[i].hists1D[j][14]->Integral());
				}
			}

			WritePlotCat1D(ofile, plotCats[0]);
			for(int i = 0; i < (int)plotCats.size(); i++)
				WritePlotCat2D(ofile, plotCats[i]);
			vector<plotCat> id_cats(plotCats.begin()+1, plotCats.end());
			WritePlotCatStack(ofile, id_cats);

			ofile->Close();

		}




		//k = sum_n(E_n)/N
		void FillModelHists(BasePDFMixture* model, int id_idx){
			map<string, Matrix> params;
			vector<double> eigenvals, eigenvals_space, norms, eigenvals_noEw, eigenvals_logEw;
			vector<Matrix> eigenvecs, eigenvecs_space, eigenvecs_noEw, eigenvecs_logEw;
			Matrix space_mat = Matrix(2,2);

			double npts = (double)model->GetData()->GetNPoints();
		//	cout << "FillHists - starting subcluster loop" << endl;	
			double E_k, phi, rot2D, ec, pc, tc, pi, E_lead, phi2D;
			//double theta, r, rot3D, vel, v_x, v_y, v_z;
			double ep_cov, te_cov, tp_cov, e_var, p_var, t_var;
			double ep_cov_unnorm, te_cov_unnorm, tp_cov_unnorm;
			double majtime_cov, mintime_cov, majtime_cov_unnorm, mintime_cov_unnorm;
			
			double cmsLogE_phi2D, cmsLogE_rot2D, cmsLogE_ec, cmsLogE_pc, cmsLogE_tc;
			double cmsLogE_ep_cov, cmsLogE_te_cov, cmsLogE_tp_cov, cmsLogE_e_var, cmsLogE_p_var, cmsLogE_t_var;
			double cmsLogE_ep_cov_unnorm, cmsLogE_te_cov_unnorm, cmsLogE_tp_cov_unnorm;
			double cmsLogE_majtime_cov, cmsLogE_mintime_cov, cmsLogE_majtime_cov_unnorm, cmsLogE_mintime_cov_unnorm;
			double cmsLogE_smaj, cmsLogE_smin, cmsNoE_smaj, cmsNoE_smin;
			
			double cmsNoE_ep_cov_unnorm, cmsNoE_te_cov_unnorm, cmsNoE_tp_cov_unnorm;
			double cmsNoE_phi2D, cmsNoE_rot2D, cmsNoE_ec, cmsNoE_pc, cmsNoE_tc;
			double cmsNoE_ep_cov, cmsNoE_te_cov, cmsNoE_tp_cov, cmsNoE_e_var, cmsNoE_p_var, cmsNoE_t_var;
			double cmsNoE_majtime_cov, cmsNoE_mintime_cov, cmsNoE_majtime_cov_unnorm, cmsNoE_mintime_cov_unnorm;

			double E_tot = 0.;
			for(int i = 0; i < npts; i++){
				E_tot += model->GetData()->at(i).w()/_gev;
			}
			Matrix cov, lead_eigenvec, lead_eigenvec_space;
			
			int nclusters = model->GetNClusters();
			plotCats[id_idx].hists1D[0][0]->Fill(nclusters);
			plotCats[id_idx].hists2D[0][10]->Fill((double)nclusters,npts);
			model->GetNorms(norms);
			
			//get leading cluster index
			vector<int> idxs;
			//sort by mixing coeffs in ascending order (smallest first)
			model->SortIdxs(idxs);
			int leadidx = idxs[nclusters-1];

			for(int k = 0; k < nclusters; k++){
				//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
				E_k = norms[k]/_gev; 
				
				params = model->GetPriorParameters(k);
				ec = params["mean"].at(0,0);
				pc = params["mean"].at(1,0);
				tc = params["mean"].at(2,0);
				pi = params["pi"].at(0,0);
				cov = params["cov"];	
				
				//eta - time sign convention
				//define relative sign for eta and time components
				//based on where the cluster is in the detector
				for(int i = 0; i < 3; i++){
					if(ec < 0){
						//time sign does NOT match eta sign
						//flip sign of eta-time entry
						cov.SetEntry(-cov.at(0,2),0,2);	
						cov.SetEntry(-cov.at(2,0),2,0);	
					}
					//else time sign matches eta sign - no change
				}


				e_var = sqrt(cov.at(0,0));
				p_var = sqrt(cov.at(1,1));
				t_var = sqrt(cov.at(2,2));

				ep_cov = CalcCov(cov, 1, 0);
				te_cov = CalcCov(cov, 2, 0);
				tp_cov = CalcCov(cov, 2, 1);
				ep_cov_unnorm = CalcCov(cov, 1, 0, false);
				te_cov_unnorm = CalcCov(cov, 2, 0, false);
				tp_cov_unnorm = CalcCov(cov, 2, 1, false);
	

				//calculate slopes from eigenvectors
				//cov.eigenCalc(eigenvals, eigenvecs);
				//lead_eigenvec = eigenvecs[2];			
				//v_x = lead_eigenvec.at(0,0);	
				//v_y = lead_eigenvec.at(1,0);	
				//v_z = lead_eigenvec.at(2,0);	
				//r = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
				//polar angle with lead eigenvector
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				//theta = acos( v_z / r );
				//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
				//phi = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				//phi = atan2( v_y , v_x  );
				//if(signbit(v_y)) phi *= -1;
				////rotundity - 3D
				//rot3D = 0;
				//for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				//rot3D = eigenvals[2]/rot3D;
				////velocity = z/r * rad/deg * deg/cm => ns/cm
				//vel = (lead_eigenvec.at(2,0)/sqrt(v_x*v_x + v_y*v_y)) * (acos(-1)/180.) * (1./2.2);
				//vel = fabs(1./vel);
				//if(isnan(vel) || isinf(vel)) vel = -999;
				
				//rotundity - 2D
				//take upper 2x2 submatrix from covariance
				Get2DMat(cov,space_mat);
				space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
				phi2D = PhiEll(space_mat);			
				rot2D = Rotundity(space_mat);
	
				//rotate points into maj/min axes
				Matrix rotmat2D = Matrix(3,3);		
				Get2DRotationMatrix(eigenvecs_space,rotmat2D);
				PointCollection majminpts;
				RotatePoints(model->GetData(), rotmat2D, majminpts);
				Matrix majminCovMat = Matrix(3,3);
				MakeCovMat(&majminpts, majminCovMat, weightScheme(1));
				//set time covariance from GMM
				majminCovMat.SetEntry(cov.at(2,2),2,2);
				majtime_cov = CalcCov(majminCovMat,2,0);
				mintime_cov = CalcCov(majminCovMat,2,1);
				majtime_cov_unnorm = CalcCov(majminCovMat,2,0,false);
				mintime_cov_unnorm = CalcCov(majminCovMat,2,1,false);

				//////////////////////calculate CMS No energy weighting variables//////////////////////
				Matrix noEw = Matrix(3,3);
				//calculate covariances like in CMSSW
				MakeCovMat(model->GetData(), noEw, weightScheme(0));	
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
				RotatePoints(model->GetData(), rotmat2D, majminpts);
				MakeCovMat(&majminpts, majminCovMat, weightScheme(2));
				//set time covariance from GMM
				majminCovMat.SetEntry(cov.at(2,2),2,2);
				cmsNoE_majtime_cov = CalcCov(majminCovMat,2,0);
				cmsNoE_mintime_cov = CalcCov(majminCovMat,2,1);
				cmsNoE_majtime_cov_unnorm = CalcCov(majminCovMat,2,0,false);
				cmsNoE_mintime_cov_unnorm = CalcCov(majminCovMat,2,1,false);
				//for centers
				Point noEwCenter = GetCMSSWMean(model->GetData(),false);
				cmsNoE_ec = noEwCenter.at(0);
				cmsNoE_pc = noEwCenter.at(1);
				cmsNoE_tc = noEwCenter.at(2);

				//////////////////////calculate CMS Log energy weighting variables//////////////////////
				Matrix logEw = Matrix(3,3);
				//calculate covariances like in CMSSW
				MakeCovMat(model->GetData(), logEw, weightScheme(2));	
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
				RotatePoints(model->GetData(), rotmat2D, majminpts);
				MakeCovMat(&majminpts, majminCovMat, weightScheme(1));
				//set time covariance from GMM
				majminCovMat.SetEntry(cov.at(2,2),2,2);
				cmsLogE_majtime_cov = CalcCov(majminCovMat,2,0);
				cmsLogE_mintime_cov = CalcCov(majminCovMat,2,1);
				cmsLogE_majtime_cov_unnorm = CalcCov(majminCovMat,2,0,false);
				cmsLogE_mintime_cov_unnorm = CalcCov(majminCovMat,2,1,false);
				//for centers
				Point logEwCenter = GetCMSSWMean(model->GetData(),true);
				cmsLogE_ec = logEwCenter.at(0);
				cmsLogE_pc = logEwCenter.at(1);
				cmsLogE_tc = logEwCenter.at(2);

			
				//fill hists
				//centers
				plotCats[id_idx].hists1D[0][1]->Fill(tc);
				plotCats[id_idx].hists1D[0][2]->Fill(ec);
				plotCats[id_idx].hists1D[0][3]->Fill(pc);
				//4 - phoE filled in .cc
				//cluster E
				plotCats[id_idx].hists1D[0][5]->Fill(E_tot);
				//slope space - phi/eta
				//plotCats[id_idx].hists1D[0][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        			////slope - eta/time
				//plotCats[id_idx].hists1D[0][7]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
				////slope - phi/time
				//plotCats[id_idx].hists1D[0][8]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
				//polar angle in 3D space
				//plotCats[id_idx].hists1D[0][9]->Fill(theta);
				//azimuthal angle in 2D space
				//plotCats[id_idx].hists1D[0][10]->Fill(phi);
				//subcluster energy
				plotCats[id_idx].hists1D[0][11]->Fill(E_k);
				//rotundity measures
				//plotCats[id_idx].hists1D[0][12]->Fill(rot3D);
				plotCats[id_idx].hists1D[0][13]->Fill(rot2D);
				//velocity	
				//plotCats[id_idx].hists1D[0][14]->Fill(vel);
				//2D eigenval ratio
				//plotCats[id_idx].hists1D[0][15]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				//get variances
				plotCats[id_idx].hists1D[0][16]->Fill(e_var);
				plotCats[id_idx].hists1D[0][17]->Fill(p_var);
				plotCats[id_idx].hists1D[0][18]->Fill(t_var);
				//fractional E
				plotCats[id_idx].hists1D[0][19]->Fill(E_k/E_tot);
				//"azimuth" angle in 2D (angle from x axis)
				plotCats[id_idx].hists1D[0][20]->Fill(phi2D);
				//pos/neg eta split sigma
				if(ec > 0) plotCats[id_idx].hists1D[0][21]->Fill(e_var);
				else plotCats[id_idx].hists1D[0][22]->Fill(e_var);
				//covariances
				plotCats[id_idx].hists1D[0][23]->Fill(ep_cov);
				plotCats[id_idx].hists1D[0][24]->Fill(te_cov);
				plotCats[id_idx].hists1D[0][25]->Fill(tp_cov);
				//major/minor covariances with time
				plotCats[id_idx].hists1D[0][26]->Fill(majtime_cov);
				plotCats[id_idx].hists1D[0][27]->Fill(mintime_cov);
				if(E_tot <= 100){
					plotCats[id_idx].hists1D[0][28]->Fill(t_var);
					plotCats[id_idx].hists1D[0][76]->Fill(rot2D);	
				}
				else if(E_tot > 100 && E_tot <= 200){
					plotCats[id_idx].hists1D[0][29]->Fill(t_var);
					plotCats[id_idx].hists1D[0][77]->Fill(rot2D);	
				}
				else{
					plotCats[id_idx].hists1D[0][30]->Fill(t_var);
					plotCats[id_idx].hists1D[0][78]->Fill(rot2D);	
				}
				plotCats[id_idx].hists1D[0][31]->Fill(te_cov_unnorm);
				plotCats[id_idx].hists1D[0][32]->Fill(tp_cov_unnorm);
				plotCats[id_idx].hists1D[0][33]->Fill(majtime_cov_unnorm);
				plotCats[id_idx].hists1D[0][34]->Fill(mintime_cov_unnorm);
				plotCats[id_idx].hists1D[0][35]->Fill(cmsLogE_te_cov);
				plotCats[id_idx].hists1D[0][36]->Fill(cmsLogE_tp_cov);
				plotCats[id_idx].hists1D[0][37]->Fill(cmsLogE_ep_cov);
				plotCats[id_idx].hists1D[0][38]->Fill(cmsLogE_e_var);
				plotCats[id_idx].hists1D[0][39]->Fill(cmsLogE_p_var);
				plotCats[id_idx].hists1D[0][40]->Fill(cmsLogE_t_var);
				plotCats[id_idx].hists1D[0][41]->Fill(cmsNoE_te_cov);	
				plotCats[id_idx].hists1D[0][42]->Fill(cmsNoE_tp_cov);	
				plotCats[id_idx].hists1D[0][43]->Fill(cmsNoE_ep_cov);	
				plotCats[id_idx].hists1D[0][44]->Fill(cmsNoE_e_var);
				plotCats[id_idx].hists1D[0][45]->Fill(cmsNoE_p_var);
				plotCats[id_idx].hists1D[0][46]->Fill(cmsNoE_t_var);
				plotCats[id_idx].hists1D[0][47]->Fill(cmsNoE_te_cov_unnorm);
				plotCats[id_idx].hists1D[0][48]->Fill(cmsNoE_tp_cov_unnorm);
				plotCats[id_idx].hists1D[0][49]->Fill(cmsNoE_ep_cov_unnorm);
				plotCats[id_idx].hists1D[0][50]->Fill(cmsLogE_te_cov_unnorm);
				plotCats[id_idx].hists1D[0][51]->Fill(cmsLogE_tp_cov_unnorm);
				plotCats[id_idx].hists1D[0][52]->Fill(cmsLogE_ep_cov_unnorm);
				plotCats[id_idx].hists1D[0][53]->Fill(ep_cov_unnorm);
				plotCats[id_idx].hists1D[0][54]->Fill(cmsNoE_smaj);
				plotCats[id_idx].hists1D[0][55]->Fill(cmsNoE_smin);
				plotCats[id_idx].hists1D[0][56]->Fill(cmsLogE_smaj);
				plotCats[id_idx].hists1D[0][57]->Fill(cmsLogE_smin);
				plotCats[id_idx].hists1D[0][58]->Fill(noEwCenter.at(2));
				plotCats[id_idx].hists1D[0][59]->Fill(noEwCenter.at(0));
				plotCats[id_idx].hists1D[0][60]->Fill(noEwCenter.at(1));
				plotCats[id_idx].hists1D[0][61]->Fill(cmsNoE_phi2D);
				plotCats[id_idx].hists1D[0][62]->Fill(cmsNoE_majtime_cov);
				plotCats[id_idx].hists1D[0][63]->Fill(cmsNoE_mintime_cov);
				plotCats[id_idx].hists1D[0][64]->Fill(cmsNoE_majtime_cov_unnorm);
				plotCats[id_idx].hists1D[0][65]->Fill(cmsNoE_mintime_cov_unnorm);
				plotCats[id_idx].hists1D[0][66]->Fill(cmsNoE_rot2D);
				plotCats[id_idx].hists1D[0][67]->Fill(logEwCenter.at(2));
				plotCats[id_idx].hists1D[0][68]->Fill(logEwCenter.at(0));
				plotCats[id_idx].hists1D[0][69]->Fill(logEwCenter.at(1));
				plotCats[id_idx].hists1D[0][70]->Fill(cmsLogE_phi2D);
				plotCats[id_idx].hists1D[0][71]->Fill(cmsLogE_majtime_cov);
				plotCats[id_idx].hists1D[0][72]->Fill(cmsLogE_mintime_cov);
				plotCats[id_idx].hists1D[0][73]->Fill(cmsLogE_majtime_cov_unnorm);
				plotCats[id_idx].hists1D[0][74]->Fill(cmsLogE_mintime_cov_unnorm);
				plotCats[id_idx].hists1D[0][75]->Fill(cmsLogE_rot2D);
				if(ep_cov < 0){
					if(te_cov < 0){
						plotCats[id_idx].hists1D[0][79]->Fill(rot2D);
						plotCats[id_idx].hists1D[0][83]->Fill(phi2D);
						plotCats[id_idx].hists1D[0][87]->Fill(ep_cov);
						plotCats[id_idx].hists1D[0][91]->Fill(te_cov);
					}	
					else{
						plotCats[id_idx].hists1D[0][80]->Fill(rot2D);
						plotCats[id_idx].hists1D[0][84]->Fill(phi2D);
						plotCats[id_idx].hists1D[0][88]->Fill(ep_cov);
						plotCats[id_idx].hists1D[0][92]->Fill(te_cov);
					}
				}
				else{
					if(te_cov < 0){
						plotCats[id_idx].hists1D[0][81]->Fill(rot2D);
						plotCats[id_idx].hists1D[0][85]->Fill(phi2D);
						plotCats[id_idx].hists1D[0][89]->Fill(ep_cov);
						plotCats[id_idx].hists1D[0][93]->Fill(te_cov);
					}	
					else{
						plotCats[id_idx].hists1D[0][82]->Fill(rot2D);
						plotCats[id_idx].hists1D[0][86]->Fill(phi2D);
						plotCats[id_idx].hists1D[0][90]->Fill(ep_cov);
						plotCats[id_idx].hists1D[0][94]->Fill(te_cov);

					}
				}


	
				//2D hists
				plotCats[id_idx].hists2D[0][0]->Fill(tc, E_k);
				plotCats[id_idx].hists2D[0][1]->Fill(phi,E_k);
				plotCats[id_idx].hists2D[0][2]->Fill(rot2D,E_k);
				plotCats[id_idx].hists2D[0][3]->Fill(ec,pc);
				plotCats[id_idx].hists2D[0][4]->Fill(tc,ec);
				plotCats[id_idx].hists2D[0][5]->Fill(tc,pc);
				plotCats[id_idx].hists2D[0][6]->Fill(tc,pi);
				plotCats[id_idx].hists2D[0][7]->Fill(E_k,pi);
				//plotCats[id_idx].hists2D[0][8]->Fill(rot3D,E_k);
				plotCats[id_idx].hists2D[0][9]->Fill(norms[k], E_k);
				plotCats[id_idx].hists2D[0][11]->Fill(nclusters, pi);
				plotCats[id_idx].hists2D[0][12]->Fill(e_var, p_var);
				plotCats[id_idx].hists2D[0][13]->Fill(t_var, e_var);
				plotCats[id_idx].hists2D[0][14]->Fill(t_var, p_var);
				plotCats[id_idx].hists2D[0][15]->Fill(E_k/E_tot, pi);
				plotCats[id_idx].hists2D[0][16]->Fill(nclusters, E_k/E_tot);
				plotCats[id_idx].hists2D[0][17]->Fill(t_var, tc);
				plotCats[id_idx].hists2D[0][18]->Fill(rot2D, phi2D);
				plotCats[id_idx].hists2D[0][19]->Fill(t_var, E_k/E_tot);
				plotCats[id_idx].hists2D[0][20]->Fill(t_var, E_tot);
				plotCats[id_idx].hists2D[0][21]->Fill(te_cov, E_tot);
				plotCats[id_idx].hists2D[0][22]->Fill(tp_cov, E_tot);
				plotCats[id_idx].hists2D[0][23]->Fill(t_var, majtime_cov);
				plotCats[id_idx].hists2D[0][24]->Fill(t_var, mintime_cov);
				plotCats[id_idx].hists2D[0][25]->Fill(phi2D, majtime_cov);
				plotCats[id_idx].hists2D[0][26]->Fill(phi2D, mintime_cov);
				plotCats[id_idx].hists2D[0][27]->Fill(rot2D, majtime_cov);
				plotCats[id_idx].hists2D[0][28]->Fill(rot2D, mintime_cov);
				//plotCats[id_idx].hists2D[0][29]->Fill(rot3D, majtime_cov);
				//plotCats[id_idx].hists2D[0][30]->Fill(rot3D, mintime_cov);
				plotCats[id_idx].hists2D[0][31]->Fill(majtime_cov, tc);
				plotCats[id_idx].hists2D[0][32]->Fill(mintime_cov, tc);
				plotCats[id_idx].hists2D[0][33]->Fill(ep_cov, phi2D);
				plotCats[id_idx].hists2D[0][34]->Fill(majtime_cov, phi2D);
				plotCats[id_idx].hists2D[0][35]->Fill(mintime_cov, phi2D);
				plotCats[id_idx].hists2D[0][36]->Fill(rot2D, majtime_cov_unnorm);
				plotCats[id_idx].hists2D[0][37]->Fill(rot2D, mintime_cov_unnorm);
				plotCats[id_idx].hists2D[0][38]->Fill(t_var, te_cov);
				plotCats[id_idx].hists2D[0][39]->Fill(ep_cov, tc);
				plotCats[id_idx].hists2D[0][40]->Fill(te_cov, tc);
				plotCats[id_idx].hists2D[0][41]->Fill(tp_cov, tc);
				plotCats[id_idx].hists2D[0][42]->Fill(cmsNoE_smaj, cmsNoE_t_var);
				plotCats[id_idx].hists2D[0][43]->Fill(cmsNoE_smin, cmsNoE_t_var);
				plotCats[id_idx].hists2D[0][44]->Fill(cmsLogE_smaj, cmsLogE_t_var);
				plotCats[id_idx].hists2D[0][45]->Fill(cmsLogE_smin, cmsLogE_t_var);
				plotCats[id_idx].hists2D[0][46]->Fill(ep_cov, te_cov);
				plotCats[id_idx].hists2D[0][47]->Fill(te_cov, tp_cov);
				plotCats[id_idx].hists2D[0][48]->Fill(ep_cov, tp_cov);
				plotCats[id_idx].hists2D[0][49]->Fill(ep_cov, majtime_cov);
				plotCats[id_idx].hists2D[0][50]->Fill(ep_cov, mintime_cov);
				plotCats[id_idx].hists2D[0][51]->Fill(tp_cov, majtime_cov);
				plotCats[id_idx].hists2D[0][52]->Fill(tp_cov, mintime_cov);
				plotCats[id_idx].hists2D[0][53]->Fill(te_cov, majtime_cov);
				plotCats[id_idx].hists2D[0][54]->Fill(te_cov, mintime_cov);
				plotCats[id_idx].hists2D[0][55]->Fill(ep_cov_unnorm, te_cov_unnorm);
				plotCats[id_idx].hists2D[0][56]->Fill(te_cov_unnorm, tp_cov_unnorm);
				plotCats[id_idx].hists2D[0][57]->Fill(ep_cov_unnorm, tp_cov_unnorm);
				plotCats[id_idx].hists2D[0][58]->Fill(ep_cov_unnorm, majtime_cov_unnorm);
				plotCats[id_idx].hists2D[0][59]->Fill(ep_cov_unnorm, mintime_cov_unnorm);
				plotCats[id_idx].hists2D[0][60]->Fill(tp_cov_unnorm, majtime_cov_unnorm);
				plotCats[id_idx].hists2D[0][61]->Fill(tp_cov_unnorm, mintime_cov_unnorm);
				plotCats[id_idx].hists2D[0][62]->Fill(te_cov_unnorm, majtime_cov_unnorm);
				plotCats[id_idx].hists2D[0][63]->Fill(te_cov_unnorm, mintime_cov_unnorm);
				plotCats[id_idx].hists2D[0][64]->Fill(rot2D, ep_cov);
				plotCats[id_idx].hists2D[0][65]->Fill(rot2D, ep_cov_unnorm);
				plotCats[id_idx].hists2D[0][66]->Fill(E_k, pc);
				plotCats[id_idx].hists2D[0][67]->Fill(rot2D,phi2D);


				//histograms for leading/subleading clusters
				if(k == leadidx){
					//centers
					plotCats[id_idx].hists1D[1][1]->Fill(tc);
					plotCats[id_idx].hists1D[1][2]->Fill(ec);
					plotCats[id_idx].hists1D[1][3]->Fill(pc);
					//4 - phoE filled in .cc
					//cluster E
					plotCats[id_idx].hists1D[1][5]->Fill(E_tot);
					//slope space - phi/eta
					//plotCats[id_idx].hists1D[1][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				////slope - eta/time
					//plotCats[id_idx].hists1D[1][7]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					////slope - phi/time
					//plotCats[id_idx].hists1D[1][8]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//polar angle in 3D space
					//plotCats[id_idx].hists1D[1][9]->Fill(theta);
					//azimuthal angle in 2D space
					//plotCats[id_idx].hists1D[1][10]->Fill(phi);
					//subcluster energy
					plotCats[id_idx].hists1D[1][11]->Fill(E_k);
					//rotundity measures
					//plotCats[id_idx].hists1D[1][12]->Fill(rot3D);
					plotCats[id_idx].hists1D[1][13]->Fill(rot2D);
					//velocity	
					//plotCats[id_idx].hists1D[1][14]->Fill(vel);
					//2D eigenval ratio
					//plotCats[id_idx].hists1D[1][15]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					//get variances
					plotCats[id_idx].hists1D[1][16]->Fill(e_var);
					plotCats[id_idx].hists1D[1][17]->Fill(p_var);
					plotCats[id_idx].hists1D[1][18]->Fill(t_var);
					//fractional E
					plotCats[id_idx].hists1D[1][19]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (angle from x axis)
					plotCats[id_idx].hists1D[1][20]->Fill(phi2D);
					//pos/neg eta split sigma
					if(ec > 0) plotCats[id_idx].hists1D[1][21]->Fill(e_var);
					else plotCats[id_idx].hists1D[1][22]->Fill(e_var);
					//covariances
					plotCats[id_idx].hists1D[1][23]->Fill(ep_cov);
					plotCats[id_idx].hists1D[1][24]->Fill(te_cov);
					plotCats[id_idx].hists1D[1][25]->Fill(tp_cov);
					//major/minor covariances with time
					plotCats[id_idx].hists1D[1][26]->Fill(majtime_cov);
					plotCats[id_idx].hists1D[1][27]->Fill(mintime_cov);
					if(E_tot <= 100){
						plotCats[id_idx].hists1D[1][28]->Fill(t_var);
						plotCats[id_idx].hists1D[1][76]->Fill(rot2D);	
					}
					else if(E_tot > 100 && E_tot <= 200){
						plotCats[id_idx].hists1D[1][29]->Fill(t_var);
						plotCats[id_idx].hists1D[1][77]->Fill(rot2D);	
					}
					else{
						plotCats[id_idx].hists1D[1][30]->Fill(t_var);
						plotCats[id_idx].hists1D[1][78]->Fill(rot2D);	
					}
					plotCats[id_idx].hists1D[1][31]->Fill(te_cov_unnorm);
					plotCats[id_idx].hists1D[1][32]->Fill(tp_cov_unnorm);
					plotCats[id_idx].hists1D[1][33]->Fill(majtime_cov_unnorm);
					plotCats[id_idx].hists1D[1][34]->Fill(mintime_cov_unnorm);
					plotCats[id_idx].hists1D[1][35]->Fill(cmsLogE_te_cov);
					plotCats[id_idx].hists1D[1][36]->Fill(cmsLogE_tp_cov);
					plotCats[id_idx].hists1D[1][37]->Fill(cmsLogE_ep_cov);
					plotCats[id_idx].hists1D[1][38]->Fill(cmsLogE_e_var);
					plotCats[id_idx].hists1D[1][39]->Fill(cmsLogE_p_var);
					plotCats[id_idx].hists1D[1][40]->Fill(cmsLogE_t_var);
					plotCats[id_idx].hists1D[1][41]->Fill(cmsNoE_te_cov);	
					plotCats[id_idx].hists1D[1][42]->Fill(cmsNoE_tp_cov);	
					plotCats[id_idx].hists1D[1][43]->Fill(cmsNoE_ep_cov);	
					plotCats[id_idx].hists1D[1][44]->Fill(cmsNoE_e_var);
					plotCats[id_idx].hists1D[1][45]->Fill(cmsNoE_p_var);
					plotCats[id_idx].hists1D[1][46]->Fill(cmsNoE_t_var);
					plotCats[id_idx].hists1D[1][47]->Fill(cmsNoE_te_cov_unnorm);
					plotCats[id_idx].hists1D[1][48]->Fill(cmsNoE_tp_cov_unnorm);
					plotCats[id_idx].hists1D[1][49]->Fill(cmsNoE_ep_cov_unnorm);
					plotCats[id_idx].hists1D[1][50]->Fill(cmsLogE_te_cov_unnorm);
					plotCats[id_idx].hists1D[1][51]->Fill(cmsLogE_tp_cov_unnorm);
					plotCats[id_idx].hists1D[1][52]->Fill(cmsLogE_ep_cov_unnorm);
					plotCats[id_idx].hists1D[1][53]->Fill(ep_cov_unnorm);
					plotCats[id_idx].hists1D[1][54]->Fill(cmsNoE_smaj);
					plotCats[id_idx].hists1D[1][55]->Fill(cmsNoE_smin);
					plotCats[id_idx].hists1D[1][56]->Fill(cmsLogE_smaj);
					plotCats[id_idx].hists1D[1][57]->Fill(cmsLogE_smin);
					plotCats[id_idx].hists1D[1][58]->Fill(noEwCenter.at(2));
					plotCats[id_idx].hists1D[1][59]->Fill(noEwCenter.at(0));
					plotCats[id_idx].hists1D[1][60]->Fill(noEwCenter.at(1));
					plotCats[id_idx].hists1D[1][61]->Fill(cmsNoE_phi2D);
					plotCats[id_idx].hists1D[1][62]->Fill(cmsNoE_majtime_cov);
					plotCats[id_idx].hists1D[1][63]->Fill(cmsNoE_mintime_cov);
					plotCats[id_idx].hists1D[1][64]->Fill(cmsNoE_majtime_cov_unnorm);
					plotCats[id_idx].hists1D[1][65]->Fill(cmsNoE_mintime_cov_unnorm);
					plotCats[id_idx].hists1D[1][66]->Fill(cmsNoE_rot2D);
					plotCats[id_idx].hists1D[1][67]->Fill(logEwCenter.at(2));
					plotCats[id_idx].hists1D[1][68]->Fill(logEwCenter.at(0));
					plotCats[id_idx].hists1D[1][69]->Fill(logEwCenter.at(1));
					plotCats[id_idx].hists1D[1][70]->Fill(cmsLogE_phi2D);
					plotCats[id_idx].hists1D[1][71]->Fill(cmsLogE_majtime_cov);
					plotCats[id_idx].hists1D[1][72]->Fill(cmsLogE_mintime_cov);
					plotCats[id_idx].hists1D[1][73]->Fill(cmsLogE_majtime_cov_unnorm);
					plotCats[id_idx].hists1D[1][74]->Fill(cmsLogE_mintime_cov_unnorm);
					plotCats[id_idx].hists1D[1][75]->Fill(cmsLogE_rot2D);
				if(ep_cov < 0){
					if(te_cov < 0){
						plotCats[id_idx].hists1D[1][79]->Fill(rot2D);
						plotCats[id_idx].hists1D[1][83]->Fill(phi2D);
						plotCats[id_idx].hists1D[1][87]->Fill(ep_cov);
						plotCats[id_idx].hists1D[1][91]->Fill(te_cov);
					}	
					else{
						plotCats[id_idx].hists1D[1][80]->Fill(rot2D);
						plotCats[id_idx].hists1D[1][84]->Fill(phi2D);
						plotCats[id_idx].hists1D[1][88]->Fill(ep_cov);
						plotCats[id_idx].hists1D[1][92]->Fill(te_cov);
					}
				}
				else{
					if(te_cov < 0){
						plotCats[id_idx].hists1D[1][81]->Fill(rot2D);
						plotCats[id_idx].hists1D[1][85]->Fill(phi2D);
						plotCats[id_idx].hists1D[1][89]->Fill(ep_cov);
						plotCats[id_idx].hists1D[1][93]->Fill(te_cov);
					}	
					else{
						plotCats[id_idx].hists1D[1][82]->Fill(rot2D);
						plotCats[id_idx].hists1D[1][86]->Fill(phi2D);
						plotCats[id_idx].hists1D[1][90]->Fill(ep_cov);
						plotCats[id_idx].hists1D[1][94]->Fill(te_cov);

					}
				}



	
					//2D hists
					plotCats[id_idx].hists2D[1][0]->Fill(tc, E_k);
					plotCats[id_idx].hists2D[1][1]->Fill(phi,E_k);
					plotCats[id_idx].hists2D[1][2]->Fill(rot2D,E_k);
					plotCats[id_idx].hists2D[1][3]->Fill(ec,pc);
					plotCats[id_idx].hists2D[1][4]->Fill(tc,ec);
					plotCats[id_idx].hists2D[1][5]->Fill(tc,pc);
					plotCats[id_idx].hists2D[1][6]->Fill(tc,pi);
					plotCats[id_idx].hists2D[1][7]->Fill(E_k,pi);
					//plotCats[id_idx].hists2D[1][8]->Fill(rot3D,E_k);
					plotCats[id_idx].hists2D[1][9]->Fill(norms[k], E_k);
					plotCats[id_idx].hists2D[1][11]->Fill(nclusters, pi);
					plotCats[id_idx].hists2D[1][12]->Fill(e_var, p_var);
					plotCats[id_idx].hists2D[1][13]->Fill(t_var, e_var);
					plotCats[id_idx].hists2D[1][14]->Fill(t_var, p_var);
					plotCats[id_idx].hists2D[1][15]->Fill(E_k/E_tot, pi);
					plotCats[id_idx].hists2D[1][16]->Fill(nclusters, E_k/E_tot);
					plotCats[id_idx].hists2D[1][17]->Fill(t_var, tc);
					plotCats[id_idx].hists2D[1][18]->Fill(rot2D, phi2D);
					plotCats[id_idx].hists2D[1][19]->Fill(t_var, E_k/E_tot);
					plotCats[id_idx].hists2D[1][20]->Fill(t_var, E_tot);
					plotCats[id_idx].hists2D[1][21]->Fill(te_cov, E_tot);
					plotCats[id_idx].hists2D[1][22]->Fill(tp_cov, E_tot);
					plotCats[id_idx].hists2D[1][23]->Fill(t_var, majtime_cov);
					plotCats[id_idx].hists2D[1][24]->Fill(t_var, mintime_cov);
					plotCats[id_idx].hists2D[1][25]->Fill(phi2D, majtime_cov);
					plotCats[id_idx].hists2D[1][26]->Fill(phi2D, mintime_cov);
					plotCats[id_idx].hists2D[1][27]->Fill(rot2D, majtime_cov);
					plotCats[id_idx].hists2D[1][28]->Fill(rot2D, mintime_cov);
					//plotCats[id_idx].hists2D[1][29]->Fill(rot3D, majtime_cov);
					//plotCats[id_idx].hists2D[1][30]->Fill(rot3D, mintime_cov);
					plotCats[id_idx].hists2D[1][31]->Fill(majtime_cov, tc);
					plotCats[id_idx].hists2D[1][32]->Fill(mintime_cov, tc);
					plotCats[id_idx].hists2D[1][33]->Fill(ep_cov, phi2D);
					plotCats[id_idx].hists2D[1][34]->Fill(majtime_cov, phi2D);
					plotCats[id_idx].hists2D[1][35]->Fill(mintime_cov, phi2D);
					plotCats[id_idx].hists2D[1][36]->Fill(rot2D, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[1][37]->Fill(rot2D, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[1][38]->Fill(t_var, te_cov);
					plotCats[id_idx].hists2D[1][39]->Fill(ep_cov, tc);
					plotCats[id_idx].hists2D[1][40]->Fill(te_cov, tc);
					plotCats[id_idx].hists2D[1][41]->Fill(tp_cov, tc);
					plotCats[id_idx].hists2D[1][42]->Fill(cmsNoE_smaj, cmsNoE_t_var);
					plotCats[id_idx].hists2D[1][43]->Fill(cmsNoE_smin, cmsNoE_t_var);
					plotCats[id_idx].hists2D[1][44]->Fill(cmsLogE_smaj, cmsLogE_t_var);
					plotCats[id_idx].hists2D[1][45]->Fill(cmsLogE_smin, cmsLogE_t_var);
					plotCats[id_idx].hists2D[1][46]->Fill(ep_cov, te_cov);
					plotCats[id_idx].hists2D[1][47]->Fill(te_cov, tp_cov);
					plotCats[id_idx].hists2D[1][48]->Fill(ep_cov, tp_cov);
					plotCats[id_idx].hists2D[1][49]->Fill(ep_cov, majtime_cov);
					plotCats[id_idx].hists2D[1][50]->Fill(ep_cov, mintime_cov);
					plotCats[id_idx].hists2D[1][51]->Fill(tp_cov, majtime_cov);
					plotCats[id_idx].hists2D[1][52]->Fill(tp_cov, mintime_cov);
					plotCats[id_idx].hists2D[1][53]->Fill(te_cov, majtime_cov);
					plotCats[id_idx].hists2D[1][54]->Fill(te_cov, mintime_cov);
					plotCats[id_idx].hists2D[1][55]->Fill(ep_cov_unnorm, te_cov_unnorm);
					plotCats[id_idx].hists2D[1][56]->Fill(te_cov_unnorm, tp_cov_unnorm);
					plotCats[id_idx].hists2D[1][57]->Fill(ep_cov_unnorm, tp_cov_unnorm);
					plotCats[id_idx].hists2D[1][58]->Fill(ep_cov_unnorm, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[1][59]->Fill(ep_cov_unnorm, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[1][60]->Fill(tp_cov_unnorm, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[1][61]->Fill(tp_cov_unnorm, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[1][62]->Fill(te_cov_unnorm, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[1][63]->Fill(te_cov_unnorm, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[1][64]->Fill(rot2D, ep_cov);
					plotCats[id_idx].hists2D[1][65]->Fill(rot2D, ep_cov_unnorm);
					plotCats[id_idx].hists2D[1][66]->Fill(E_k, pc);
					plotCats[id_idx].hists2D[1][67]->Fill(rot2D,phi2D);
				

				}

				else if(k != leadidx){
					//centers
					plotCats[id_idx].hists1D[2][1]->Fill(tc);
					plotCats[id_idx].hists1D[2][2]->Fill(ec);
					plotCats[id_idx].hists1D[2][3]->Fill(pc);
					//4 - phoE filled in .cc
					//cluster E
					plotCats[id_idx].hists1D[2][5]->Fill(E_tot);
					//slope space - phi/eta
					//plotCats[id_idx].hists1D[2][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				////slope - eta/time
					//plotCats[id_idx].hists1D[2][7]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					////slope - phi/time
					//plotCats[id_idx].hists1D[2][8]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//polar angle in 3D space
					//plotCats[id_idx].hists1D[2][9]->Fill(theta);
					//azimuthal angle in 2D space
					//plotCats[id_idx].hists1D[2][10]->Fill(phi);
					//subcluster energy
					plotCats[id_idx].hists1D[2][11]->Fill(E_k);
					//rotundity measures
					//plotCats[id_idx].hists1D[2][12]->Fill(rot3D);
					plotCats[id_idx].hists1D[2][13]->Fill(rot2D);
					//velocity	
					//plotCats[id_idx].hists1D[2][14]->Fill(vel);
					//2D eigenval ratio
					//plotCats[id_idx].hists1D[2][15]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					//get variances
					plotCats[id_idx].hists1D[2][16]->Fill(e_var);
					plotCats[id_idx].hists1D[2][17]->Fill(p_var);
					plotCats[id_idx].hists1D[2][18]->Fill(t_var);
					//fractional E
					plotCats[id_idx].hists1D[2][19]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (angle from x axis)
					plotCats[id_idx].hists1D[2][20]->Fill(phi2D);
					//pos/neg eta split sigma
					if(ec > 0) plotCats[id_idx].hists1D[2][21]->Fill(e_var);
					else plotCats[id_idx].hists1D[2][22]->Fill(e_var);
					//covariances
					plotCats[id_idx].hists1D[2][23]->Fill(ep_cov);
					plotCats[id_idx].hists1D[2][24]->Fill(te_cov);
					plotCats[id_idx].hists1D[2][25]->Fill(tp_cov);
					//major/minor covariances with time
					plotCats[id_idx].hists1D[2][26]->Fill(majtime_cov);
					plotCats[id_idx].hists1D[2][27]->Fill(mintime_cov);
					if(E_tot <= 100) plotCats[id_idx].hists1D[2][28]->Fill(t_var);
					else if(E_tot > 100 && E_tot <= 200) plotCats[id_idx].hists1D[2][29]->Fill(t_var);
					else plotCats[id_idx].hists1D[2][30]->Fill(t_var);
					plotCats[id_idx].hists1D[2][31]->Fill(te_cov_unnorm);
					plotCats[id_idx].hists1D[2][32]->Fill(tp_cov_unnorm);
					plotCats[id_idx].hists1D[2][33]->Fill(majtime_cov_unnorm);
					plotCats[id_idx].hists1D[2][34]->Fill(mintime_cov_unnorm);
					plotCats[id_idx].hists1D[2][35]->Fill(cmsLogE_te_cov);
					plotCats[id_idx].hists1D[2][36]->Fill(cmsLogE_tp_cov);
					plotCats[id_idx].hists1D[2][37]->Fill(cmsLogE_ep_cov);
					plotCats[id_idx].hists1D[2][38]->Fill(cmsLogE_e_var);
					plotCats[id_idx].hists1D[2][39]->Fill(cmsLogE_p_var);
					plotCats[id_idx].hists1D[2][40]->Fill(cmsLogE_t_var);
					plotCats[id_idx].hists1D[2][41]->Fill(cmsNoE_te_cov);	
					plotCats[id_idx].hists1D[2][42]->Fill(cmsNoE_tp_cov);	
					plotCats[id_idx].hists1D[2][43]->Fill(cmsNoE_ep_cov);	
					plotCats[id_idx].hists1D[2][44]->Fill(cmsNoE_e_var);
					plotCats[id_idx].hists1D[2][45]->Fill(cmsNoE_p_var);
					plotCats[id_idx].hists1D[2][46]->Fill(cmsNoE_t_var);
					plotCats[id_idx].hists1D[2][47]->Fill(cmsNoE_te_cov_unnorm);
					plotCats[id_idx].hists1D[2][48]->Fill(cmsNoE_tp_cov_unnorm);
					plotCats[id_idx].hists1D[2][49]->Fill(cmsNoE_ep_cov_unnorm);
					plotCats[id_idx].hists1D[2][50]->Fill(cmsLogE_te_cov_unnorm);
					plotCats[id_idx].hists1D[2][51]->Fill(cmsLogE_tp_cov_unnorm);
					plotCats[id_idx].hists1D[2][52]->Fill(cmsLogE_ep_cov_unnorm);
					plotCats[id_idx].hists1D[2][53]->Fill(ep_cov_unnorm);
					plotCats[id_idx].hists1D[2][54]->Fill(cmsNoE_smaj);
					plotCats[id_idx].hists1D[2][55]->Fill(cmsNoE_smin);
					plotCats[id_idx].hists1D[2][56]->Fill(cmsLogE_smaj);
					plotCats[id_idx].hists1D[2][57]->Fill(cmsLogE_smin);
					plotCats[id_idx].hists1D[2][58]->Fill(noEwCenter.at(2));
					plotCats[id_idx].hists1D[2][59]->Fill(noEwCenter.at(0));
					plotCats[id_idx].hists1D[2][60]->Fill(noEwCenter.at(1));
					plotCats[id_idx].hists1D[2][61]->Fill(cmsNoE_phi2D);
					plotCats[id_idx].hists1D[2][62]->Fill(cmsNoE_majtime_cov);
					plotCats[id_idx].hists1D[2][63]->Fill(cmsNoE_mintime_cov);
					plotCats[id_idx].hists1D[2][64]->Fill(cmsNoE_majtime_cov_unnorm);
					plotCats[id_idx].hists1D[2][65]->Fill(cmsNoE_mintime_cov_unnorm);
					plotCats[id_idx].hists1D[2][66]->Fill(cmsNoE_rot2D);
					plotCats[id_idx].hists1D[2][67]->Fill(logEwCenter.at(2));
					plotCats[id_idx].hists1D[2][68]->Fill(logEwCenter.at(0));
					plotCats[id_idx].hists1D[2][69]->Fill(logEwCenter.at(1));
					plotCats[id_idx].hists1D[2][70]->Fill(cmsLogE_phi2D);
					plotCats[id_idx].hists1D[2][71]->Fill(cmsLogE_majtime_cov);
					plotCats[id_idx].hists1D[2][72]->Fill(cmsLogE_mintime_cov);
					plotCats[id_idx].hists1D[2][73]->Fill(cmsLogE_majtime_cov_unnorm);
					plotCats[id_idx].hists1D[2][74]->Fill(cmsLogE_mintime_cov_unnorm);
					plotCats[id_idx].hists1D[2][75]->Fill(cmsLogE_rot2D);
				if(ep_cov < 0){
					if(te_cov < 0){
						plotCats[id_idx].hists1D[2][79]->Fill(rot2D);
						plotCats[id_idx].hists1D[2][83]->Fill(phi2D);
						plotCats[id_idx].hists1D[2][87]->Fill(ep_cov);
						plotCats[id_idx].hists1D[2][91]->Fill(te_cov);
					}	
					else{
						plotCats[id_idx].hists1D[2][80]->Fill(rot2D);
						plotCats[id_idx].hists1D[2][84]->Fill(phi2D);
						plotCats[id_idx].hists1D[2][88]->Fill(ep_cov);
						plotCats[id_idx].hists1D[2][92]->Fill(te_cov);
					}
				}
				else{
					if(te_cov < 0){
						plotCats[id_idx].hists1D[2][81]->Fill(rot2D);
						plotCats[id_idx].hists1D[2][85]->Fill(phi2D);
						plotCats[id_idx].hists1D[2][89]->Fill(ep_cov);
						plotCats[id_idx].hists1D[2][93]->Fill(te_cov);
					}	
					else{
						plotCats[id_idx].hists1D[0][82]->Fill(rot2D);
						plotCats[id_idx].hists1D[0][86]->Fill(phi2D);
						plotCats[id_idx].hists1D[0][90]->Fill(ep_cov);
						plotCats[id_idx].hists1D[0][94]->Fill(te_cov);

					}
				}

	
					//2D hists
					plotCats[id_idx].hists2D[2][0]->Fill(tc, E_k);
					plotCats[id_idx].hists2D[2][1]->Fill(phi,E_k);
					plotCats[id_idx].hists2D[2][2]->Fill(rot2D,E_k);
					plotCats[id_idx].hists2D[2][3]->Fill(ec,pc);
					plotCats[id_idx].hists2D[2][4]->Fill(tc,ec);
					plotCats[id_idx].hists2D[2][5]->Fill(tc,pc);
					plotCats[id_idx].hists2D[2][6]->Fill(tc,pi);
					plotCats[id_idx].hists2D[2][7]->Fill(E_k,pi);
					//plotCats[id_idx].hists2D[2][8]->Fill(rot3D,E_k);
					plotCats[id_idx].hists2D[2][9]->Fill(norms[k], E_k);
					plotCats[id_idx].hists2D[2][11]->Fill(nclusters, pi);
					plotCats[id_idx].hists2D[2][12]->Fill(e_var, p_var);
					plotCats[id_idx].hists2D[2][13]->Fill(t_var, e_var);
					plotCats[id_idx].hists2D[2][14]->Fill(t_var, p_var);
					plotCats[id_idx].hists2D[2][15]->Fill(E_k/E_tot, pi);
					plotCats[id_idx].hists2D[2][16]->Fill(nclusters, E_k/E_tot);
					plotCats[id_idx].hists2D[2][17]->Fill(t_var, tc);
					plotCats[id_idx].hists2D[2][18]->Fill(rot2D, phi2D);
					plotCats[id_idx].hists2D[2][19]->Fill(t_var, E_k/E_tot);
					plotCats[id_idx].hists2D[2][20]->Fill(t_var, E_tot);
					plotCats[id_idx].hists2D[2][21]->Fill(te_cov, E_tot);
					plotCats[id_idx].hists2D[2][22]->Fill(tp_cov, E_tot);
					plotCats[id_idx].hists2D[2][23]->Fill(t_var, majtime_cov);
					plotCats[id_idx].hists2D[2][24]->Fill(t_var, mintime_cov);
					plotCats[id_idx].hists2D[2][25]->Fill(phi2D, majtime_cov);
					plotCats[id_idx].hists2D[2][26]->Fill(phi2D, mintime_cov);
					plotCats[id_idx].hists2D[2][27]->Fill(rot2D, majtime_cov);
					plotCats[id_idx].hists2D[2][28]->Fill(rot2D, mintime_cov);
					//plotCats[id_idx].hists2D[2][29]->Fill(rot3D, majtime_cov);
					//plotCats[id_idx].hists2D[2][30]->Fill(rot3D, mintime_cov);
					plotCats[id_idx].hists2D[2][31]->Fill(majtime_cov, tc);
					plotCats[id_idx].hists2D[2][32]->Fill(mintime_cov, tc);
					plotCats[id_idx].hists2D[2][33]->Fill(ep_cov, phi2D);
					plotCats[id_idx].hists2D[2][34]->Fill(majtime_cov, phi2D);
					plotCats[id_idx].hists2D[2][35]->Fill(mintime_cov, phi2D);
					plotCats[id_idx].hists2D[2][36]->Fill(rot2D, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[2][37]->Fill(rot2D, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[2][38]->Fill(t_var, te_cov);
					plotCats[id_idx].hists2D[2][39]->Fill(ep_cov, tc);
					plotCats[id_idx].hists2D[2][40]->Fill(te_cov, tc);
					plotCats[id_idx].hists2D[2][41]->Fill(tp_cov, tc);
					plotCats[id_idx].hists2D[2][42]->Fill(cmsNoE_smaj, cmsNoE_t_var);
					plotCats[id_idx].hists2D[2][43]->Fill(cmsNoE_smin, cmsNoE_t_var);
					plotCats[id_idx].hists2D[2][44]->Fill(cmsLogE_smaj, cmsLogE_t_var);
					plotCats[id_idx].hists2D[2][45]->Fill(cmsLogE_smin, cmsLogE_t_var);
					plotCats[id_idx].hists2D[2][46]->Fill(ep_cov, te_cov);
					plotCats[id_idx].hists2D[2][47]->Fill(te_cov, tp_cov);
					plotCats[id_idx].hists2D[2][48]->Fill(ep_cov, tp_cov);
					plotCats[id_idx].hists2D[2][49]->Fill(ep_cov, majtime_cov);
					plotCats[id_idx].hists2D[2][50]->Fill(ep_cov, mintime_cov);
					plotCats[id_idx].hists2D[2][51]->Fill(tp_cov, majtime_cov);
					plotCats[id_idx].hists2D[2][52]->Fill(tp_cov, mintime_cov);
					plotCats[id_idx].hists2D[2][53]->Fill(te_cov, majtime_cov);
					plotCats[id_idx].hists2D[2][54]->Fill(te_cov, mintime_cov);
					plotCats[id_idx].hists2D[2][55]->Fill(ep_cov_unnorm, te_cov_unnorm);
					plotCats[id_idx].hists2D[2][56]->Fill(te_cov_unnorm, tp_cov_unnorm);
					plotCats[id_idx].hists2D[2][57]->Fill(ep_cov_unnorm, tp_cov_unnorm);
					plotCats[id_idx].hists2D[2][58]->Fill(ep_cov_unnorm, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[2][59]->Fill(ep_cov_unnorm, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[2][60]->Fill(tp_cov_unnorm, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[2][61]->Fill(tp_cov_unnorm, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[2][62]->Fill(te_cov_unnorm, majtime_cov_unnorm);
					plotCats[id_idx].hists2D[2][63]->Fill(te_cov_unnorm, mintime_cov_unnorm);
					plotCats[id_idx].hists2D[2][64]->Fill(rot2D, ep_cov);
					plotCats[id_idx].hists2D[2][65]->Fill(rot2D, ep_cov_unnorm);
				}

			}
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
		return sqrtcov(cov.at(i,j)/denom);
	}


	double Rotundity(Matrix& inmat){
		vector<Matrix> eigenvecs;
		vector<double> eigenvals;
		inmat.eigenCalc(eigenvals, eigenvecs);
		int maxd = inmat.GetDims()[0] - 1;
		double rot = 0;
		for(int i = 0; i < (int)eigenvals.size(); i++) rot += eigenvals[i];
		rot = eigenvals[maxd]/rot;
		if(rot < 0.5 || rot > 1) cout << "rot: " << rot << endl;
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

	//is a 3D rotation matrix but only does 2D rotation
	void Get2DRotationMatrix(vector<Matrix> eigenvecs, Matrix& rotmat){
		if(rotmat.GetDims()[0] != 3) return;
		if(eigenvecs.size() != 2) return;
		rotmat.reset();
		rotmat.SetEntry(eigenvecs[1].at(0,0),0,0);
		rotmat.SetEntry(eigenvecs[1].at(1,0),0,1);
		rotmat.SetEntry(eigenvecs[0].at(0,0),1,0);
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
			Point pt = pc->at(i);
			oldpt.SetEntry(pt.at(0),0,0);
			oldpt.SetEntry(pt.at(1),1,0);
			oldpt.SetEntry(pt.at(2),2,0);	

			//rotate by above transformation
			newpt.mult(rotmat,oldpt);
			//may need to preserve weight
			PointCollection tmp_newpts = newpt.MatToPoints();
			Point newpoint = tmp_newpts.at(0);
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

	Point GetCMSSWMean(PointCollection* pc, bool logw = false){
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
		Point pt;
		for(int i = 0; i < npts; i ++){
			pt = pc->at(i);
			if(logw) pt.SetWeight( log( w0 + (pc->at(i).w()/_gev)/E_tot ) );
			else pt.SetWeight( 1.0 );
			pcnew.AddPoint(pt);
		}
		//is weighted mean
		Point mean = Point(maxd);
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
		Point pt;
		for(int i = 0; i < npts; i ++){
			pt = pc->at(i);
			if(ws == 0) pt.SetWeight( 1.0 );
			//already e-weighted
			if(ws == 2) pt.SetWeight( log( w0 + (pc->at(i).w()/_gev)/E_tot ) );
			pcnew.AddPoint(pt);
		}
		//is weighted mean
		Point mean = Point(maxd);
		for(int d = 0; d < maxd; d++)
			mean.SetValue(pcnew.Centroid(d),d);
		double ent;
		for(int d1 = 0; d1 < maxd; d1++){
			for(int d2 = d1; d2 < maxd; d2++){ 
				ent = 0;
				for(int i = 0; i < npts; i++){
					ent += pcnew.at(i).w() * (pcnew.at(i).at(d1) - mean.at(d1))*(pcnew.at(i).at(d2) - mean.at(d2)) / pcnew.Sumw();
				}
				outcov.SetEntry(ent,d1,d2);
				if(d1 != d2) outcov.SetEntry(ent,d2,d1);
			}
		}
		//eta time sign convention
		if(mean.at(0) < 0){
			//time sign does NOT match eta sign
			//flip sign of eta-time entry
			outcov.SetEntry(-outcov.at(0,2),0,2);	
			outcov.SetEntry(-outcov.at(2,0),2,0);	
		}



	}






};
#endif
