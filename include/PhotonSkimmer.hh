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
                        _hists1D.push_back(azimuth_ang_2D);
                        _hists1D.push_back(etaSig_pos);
                        _hists1D.push_back(etaSig_neg);
                        _hists1D.push_back(etaphi_cov);
                        _hists1D.push_back(timeeta_cov);
                        _hists1D.push_back(timephi_cov);
			_hists1D.push_back(timemaj_cov);
			_hists1D.push_back(timemin_cov);
			
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
                	_hists2D.push_back(azAngle2D_timeMajCov);
                	_hists2D.push_back(azAngle2D_timeMinCov);
                	_hists2D.push_back(rot2D_timeMajCov);
                	_hists2D.push_back(rot2D_timeMinCov);
                	_hists2D.push_back(rot3D_timeMajCov);
                	_hists2D.push_back(rot3D_timeMinCov);

			
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
		TH1D* azimuth_ang = new TH1D("azimuth_ang","azimuth_ang",25,0.,3.5);		
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
                TH1D* timeSig = new TH1D("timeSig","timeSig",25,0,10.);
		//19 - fraction of energy in cluster
		TH1D* fracE = new TH1D("fracE","fracE",50,0.,1.1);
		//20 - azimuth angle in 2D
		TH1D* azimuth_ang_2D = new TH1D("azimuth_ang_2D","azimuth_ang_2D",50,0.,3.5);		
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
		TH1D* timemaj_cov = new TH1D("timemaj_cov","timemaj_cov",25,-0.5,0.5);
		//27 - normalized covariance - time/minor axis
		TH1D* timemin_cov = new TH1D("timemin_cov","timemin_cov",25,-0.5,0.5);


		//0 - time v subcl subcluster energy
		TH2D* time_E = new TH2D("time_subclE","time_subclE;time_center;E;a.u.", 50,-30,30,10,0,1000);
		//1 - azimuthal angle v subcl energy
		TH2D* az_E = new TH2D("az_subclE","az_subclE;azimuthal_angle;E;a.u.",50,-3.5,3.5,10,0,1000);
		//2 - rotundity (2D) v subcl energy
		TH2D* rot2D_E = new TH2D("rot2D_subclE","rot2D_subclE;rotundity2D;E;a.u.",50,0.4,1.1,10,0,1000);
		//3 - eta v phi
		TH2D* eta_phi = new TH2D("eta_phi","eta_phi;eta_center;phi_center",50,-3.5,3.5,50,-3.5,3.5);
		//4 - t v eta
		TH2D* t_eta = new TH2D("t_eta","t_eta;time_center;eta_center",50,-30,30,50,-3.5,3.5);
		//5 - t v phi
		TH2D* t_phi = new TH2D("t_phi","t_phi;time_center;phi_center;a.u.",50,-30,30,50,-3.5,3.5);
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
		TH2D* timeSig_timeCenter = new TH2D("timeSig_timeCenter","timeSig_timeCenter;timeSig;time center",25,0,10,50,-5,15);
		//18 - 2d rotundity vs 2d az angle
		TH2D* rot2D_az2D = new TH2D("rot2D_az2D","rot2D_az2D;rotundity2D;azangle2D",50,0.4,1.1,50,0.,3.5);
		//19 - time variation vs frac E
		TH2D* timeSig_fracE = new TH2D("timeSig_fracE","timeSig_fracE;timeSig;fracE",25,0,25,20,0.,1.1);
		//20 - time sigma vs total E
		TH2D* timeSig_totE = new TH2D("timeSig_totE","timeSig_totE;timeSig;totE",30,0,5,30,0,1400);
		//21 - time-eta covariance vs total E
		TH2D* timeEtaCov_totE = new TH2D("timeEtaCov_totE","timeEtaCov_totE;timeEtaCov;totE",25,-1,1,20,0,2000);
		//22 - time-phi covariance vs total E
		TH2D* timePhiCov_totE = new TH2D("timePhiCov_totE","timePhiCov_totE;timePhiCob;totE",25,-1,1,20,0,2000);
                //23 - time sigma vs. TimeMajCov
                TH2D* timeSig_timeMajCov = new TH2D("timeSig_timeMajCov","timeSig_timeMajCov;timeSig;timeMajCov",25,0,5.,25,-0.5,0.5);
                //24 - time sigma vs. TimeMinCov
                TH2D* timeSig_timeMinCov = new TH2D("timeSig_timeMinCov","timeSig_timeMinCov;timeSig;timeMinCov",25,0,5.,25,-0.5,0.5);
                //25 - az angle 2D vs. TimeMajCov
                TH2D* azAngle2D_timeMajCov = new TH2D("azAngle2D_timeMajCov","azAngle2D_timeMajCov;azAngle2D;timeMajCov",25,0,3.5,25,-0.5,0.5);
                //26 - az 2D angle vs. TimeMinCov
                TH2D* azAngle2D_timeMinCov = new TH2D("azAngle2D_timeMinCov","azAngle2D_timeMinCov;azAngle2D;timeMinCov",25,0,3.5,25,-0.5,0.5);
                //27 - rot 2D vs. TimeMajCov
                TH2D* rot2D_timeMajCov = new TH2D("rot2D_timeMajCov","rot2D_timeMajCov;rot2D;timeMajCov",25,0.4,1.1,25,-0.5,0.5);
                //28 - rot 2D angle vs. TimeMinCov
                TH2D* rot2D_timeMinCov = new TH2D("rot2D_timeMinCov","rot2D_timeMinCov;rot2D;timeMinCov",25,0.4,1.1,25,-0.5,0.5);
                //29 - rot 3D vs. TimeMajCov
                TH2D* rot3D_timeMajCov = new TH2D("rot3D_timeMajCov","rot3D_timeMajCov;rot3D;timeMajCov",25,0.4,1.1,25,-0.5,0.5);
                //30 - rot 3D angle vs. TimeMinCov
                TH2D* rot3D_timeMinCov = new TH2D("rot3D_timeMinCov","rot3D_timeMinCov;rot3D;timeMinCov",25,0.4,1.1,25,-0.5,0.5);

	
		void Skim();

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

		void MakeIDHists(){
			//total
			plotCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			plotCats.push_back(tot);	
		
		
			//notSunm
			plotCat notSunm(_hists1D, _hists2D, "notSunm","notSunm");
			//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
			notSunm.ids = {29, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
			plotCats.push_back(notSunm);
			
			//signal
			plotCat sig(_hists1D, _hists2D, "chiGam","#Chi^{0} #rightarrow #gamma");
			sig.ids = {22};
			plotCats.push_back(sig);

			//data
			plotCat data(_hists1D, _hists2D, "JetHT", "JetHT");
			data.ids = {-999};

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
					if(pc.hists2D[i][j]->GetEntries() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
					pc.hists2D[i][j]->Write();
				}
			}

		}



		void WritePlotCatStack(TFile* ofile, const vector<plotCat>& pcs){
			ofile->cd();
			string name;
			double ymax, ymin;
			vector<string> id_names;
			vector<TH1D*> hists;
			//number of histogram categories (ie leading, !leading, etc)
			int nhistCats = pcs[0].hists1D.size();
			//write 1D hists
			//variables
			for(int j = 0; j < (int)_hists1D.size(); j++){
				//cout << "var: " << _hists1D[j]->GetName() << endl;
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
					hists.clear();
					id_names.clear();
					//proc
					for(int k = 0; k < pcs.size(); k++){
						if(pcs[k].hists1D[i][j] == nullptr) continue;
						if(pcs[k].hists1D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram for proc " << pcs[k].plotName << " not filled." << endl; continue; }
			//			cout << "		adding proc " << pcs[k].plotName << " to plot with hist " << pcs[k].hists1D[i][j]->GetName() << endl;
						pcs[k].hists1D[i][j]->Write();
						//hists.push_back(pcs[k].hists1D[i][j]);			
						//id_names.push_back(pcs[k].legName);
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
			vector<double> eigenvals, eigenvals_space, norms;
			vector<Matrix> eigenvecs, eigenvecs_space;
			Matrix space_mat = Matrix(2,2);
			Matrix rotmat2D = Matrix(2,2);		
			Matrix etaphi = Matrix(2,1);	
			Matrix majmin = Matrix(2,1);	

			double npts = (double)model->GetData()->GetNPoints();
		//	cout << "FillHists - starting subcluster loop" << endl;	
			double E_k, theta, phi, r, rot2D, rot3D, vel, ec, pc, tc, pi, E_lead, phi2D;
			double v_x, v_y, v_z, ep_cov, te_cov, tp_cov, e_var, p_var, t_var;
			double majtime_cov, mintime_cov;
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

				//calculate slopes from eigenvectors
				cov.eigenCalc(eigenvals, eigenvecs);
				lead_eigenvec = eigenvecs[2];			

				e_var = cov.at(0,0);
				p_var = cov.at(1,1);
				t_var = cov.at(2,2);

				ep_cov = cov.at(1,0)/sqrt(e_var*p_var);
				te_cov = cov.at(2,0)/sqrt(t_var*e_var);
				tp_cov = cov.at(2,1)/sqrt(t_var*p_var);

				v_x = lead_eigenvec.at(0,0);	
				v_y = lead_eigenvec.at(1,0);	
				v_z = lead_eigenvec.at(2,0);	
				r = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
				//polar angle with lead eigenvector
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				theta = acos( v_z / r );
				//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
				//phi = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				phi = atan2( v_y , v_x  );
				if(signbit(v_y)) phi *= -1;
				//rotundity - 3D
				rot3D = 0;
				for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				rot3D = eigenvals[2]/rot3D;
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (lead_eigenvec.at(2,0)/sqrt(v_x*v_x + v_y*v_y)) * (acos(-1)/180.) * (1./2.2);
				vel = fabs(1./vel);
				if(isnan(vel) || isinf(vel)) vel = -999;
				
				//rotundity - 2D
				//take upper 2x2 submatrix from covariance
				space_mat.SetEntry(cov.at(0,0),0,0);	
				space_mat.SetEntry(cov.at(0,1),0,1);	
				space_mat.SetEntry(cov.at(1,0),1,0);	
				space_mat.SetEntry(cov.at(1,1),1,1);
				space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
				rot2D = 0;
				for(int i = 0; i < (int)eigenvecs_space.size(); i++) rot2D += eigenvals_space[i];
				rot2D = eigenvals_space[1]/rot2D;
				lead_eigenvec_space = eigenvecs_space[1];			
				v_x = lead_eigenvec_space.at(0,0);	
				v_y = lead_eigenvec_space.at(1,0);	
				//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
				phi2D = atan2( v_y , v_x  );
				if(signbit(v_y)) phi2D *= -1;
				
				//put into function - inputs: rotation matrix, point collection; outputs - covs
				//2D rotated with 3D covariance
				rotmat2D.SetEntry(eigenvals_space[1],0,0);
				rotmat2D.SetEntry(eigenvals_space[0],1,1);
				rotmat2D.SetEntry(0,1,0);
				rotmat2D.SetEntry(0,0,1);
			
				majtime_cov = 0;
				mintime_cov = 0;
				PointCollection majminpts;
				//calculat points rotated into major/minor axes
				for(int i = 0; i < model->GetData()->GetNPoints(); i++){
					//put eta/phi coordinates into major/minor coordinates
					Point pt = model->GetData()->at(i);
					etaphi.SetEntry(pt.at(0),0,0);
					etaphi.SetEntry(pt.at(1),1,0);

					//rotate by above transformation
					majmin.mult(rotmat2D,etaphi);	
					PointCollection newpts = majmin.MatToPoints();
					majminpts.AddPoints(newpts);	
				}
				Point majminmean = majminpts.mean();
				//calculate covariances of major/minor axes with time
				//variance of major/minor axes are eigenvalues (lead/sublead) from 2D covariance
				for(int i = 0; i < majminpts.GetNPoints(); i++){
					majtime_cov += (majminpts.at(i).at(0) - majminmean.at(0))*(model->GetData()->at(i).at(2) - model->GetData()->mean().at(2));
					mintime_cov += (majminpts.at(i).at(1) - majminmean.at(1))*(model->GetData()->at(i).at(2) - model->GetData()->mean().at(2));
	
				}
				majtime_cov /= sqrt(rotmat2D.at(0,0)*cov.at(2,2));	
				mintime_cov /= sqrt(rotmat2D.at(1,1)*cov.at(2,2));	





				//fill hists
				//centers
				plotCats[id_idx].hists1D[0][1]->Fill(tc);
				plotCats[id_idx].hists1D[0][2]->Fill(ec);
				plotCats[id_idx].hists1D[0][3]->Fill(pc);
				//4 - phoE filled in .cc
				//cluster E
				plotCats[id_idx].hists1D[0][5]->Fill(E_tot);
				//slope space - phi/eta
				plotCats[id_idx].hists1D[0][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        			//slope - eta/time
				plotCats[id_idx].hists1D[0][7]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
				//slope - phi/time
				plotCats[id_idx].hists1D[0][8]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
				//polar angle in 3D space
				plotCats[id_idx].hists1D[0][9]->Fill(theta);
				//azimuthal angle in 2D space
				plotCats[id_idx].hists1D[0][10]->Fill(phi);
				//subcluster energy
				plotCats[id_idx].hists1D[0][11]->Fill(E_k);
				//rotundity measures
				plotCats[id_idx].hists1D[0][12]->Fill(rot3D);
				plotCats[id_idx].hists1D[0][13]->Fill(rot2D);
				//velocity	
				plotCats[id_idx].hists1D[0][14]->Fill(vel);
				//2D eigenval ratio
				plotCats[id_idx].hists1D[0][15]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				//get variances
				plotCats[id_idx].hists1D[0][16]->Fill(sqrt(e_var));
				plotCats[id_idx].hists1D[0][17]->Fill(sqrt(p_var));
				plotCats[id_idx].hists1D[0][18]->Fill(sqrt(t_var));
				//fractional E
				plotCats[id_idx].hists1D[0][19]->Fill(E_k/E_tot);
				//"azimuth" angle in 2D (angle from x axis)
				plotCats[id_idx].hists1D[0][20]->Fill(phi2D);
				//pos/neg eta split sigma
				if(ec > 0) plotCats[id_idx].hists1D[0][21]->Fill(sqrt(e_var));
				else plotCats[id_idx].hists1D[0][22]->Fill(sqrt(e_var));
				//covariances
				plotCats[id_idx].hists1D[0][23]->Fill(ep_cov);
				plotCats[id_idx].hists1D[0][24]->Fill(te_cov);
				plotCats[id_idx].hists1D[0][25]->Fill(tp_cov);
				//major/minor covariances with time
				plotCats[id_idx].hists1D[0][26]->Fill(majtime_cov);
				plotCats[id_idx].hists1D[0][27]->Fill(mintime_cov);


	
				//2D hists
				plotCats[id_idx].hists2D[0][0]->Fill(tc, E_k);
				plotCats[id_idx].hists2D[0][1]->Fill(phi,E_k);
				plotCats[id_idx].hists2D[0][2]->Fill(rot2D,E_k);
				plotCats[id_idx].hists2D[0][3]->Fill(ec,pc);
				plotCats[id_idx].hists2D[0][4]->Fill(tc,ec);
				plotCats[id_idx].hists2D[0][5]->Fill(tc,pc);
				plotCats[id_idx].hists2D[0][6]->Fill(tc,pi);
				plotCats[id_idx].hists2D[0][7]->Fill(E_k,pi);
				plotCats[id_idx].hists2D[0][8]->Fill(rot3D,E_k);
				plotCats[id_idx].hists2D[0][9]->Fill(norms[k], E_k);
				plotCats[id_idx].hists2D[0][11]->Fill(nclusters, pi);
				plotCats[id_idx].hists2D[0][12]->Fill(sqrt(e_var), sqrt(p_var));
				plotCats[id_idx].hists2D[0][13]->Fill(sqrt(t_var), sqrt(e_var));
				plotCats[id_idx].hists2D[0][14]->Fill(sqrt(t_var), sqrt(p_var));
				plotCats[id_idx].hists2D[0][15]->Fill(E_k/E_tot, pi);
				plotCats[id_idx].hists2D[0][16]->Fill(nclusters, E_k/E_tot);
				plotCats[id_idx].hists2D[0][17]->Fill(sqrt(t_var), tc);
				plotCats[id_idx].hists2D[0][18]->Fill(rot2D, phi2D);
				plotCats[id_idx].hists2D[0][19]->Fill(sqrt(t_var), E_k/E_tot);
				plotCats[id_idx].hists2D[0][20]->Fill(sqrt(t_var), E_tot);
				plotCats[id_idx].hists2D[0][21]->Fill(te_cov, E_tot);
				plotCats[id_idx].hists2D[0][22]->Fill(tp_cov, E_tot);
				plotCats[id_idx].hists2D[0][23]->Fill(sqrt(t_var), majtime_cov);
				plotCats[id_idx].hists2D[0][24]->Fill(sqrt(t_var), mintime_cov);
				plotCats[id_idx].hists2D[0][25]->Fill(phi2D, majtime_cov);
				plotCats[id_idx].hists2D[0][26]->Fill(phi2D, mintime_cov);
				plotCats[id_idx].hists2D[0][27]->Fill(rot2D, majtime_cov);
				plotCats[id_idx].hists2D[0][28]->Fill(rot2D, mintime_cov);
				plotCats[id_idx].hists2D[0][29]->Fill(rot3D, majtime_cov);
				plotCats[id_idx].hists2D[0][30]->Fill(rot3D, mintime_cov);

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
					plotCats[id_idx].hists1D[1][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				//slope - eta/time
					plotCats[id_idx].hists1D[1][7]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					//slope - phi/time
					plotCats[id_idx].hists1D[1][8]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//polar angle in 3D space
					plotCats[id_idx].hists1D[1][9]->Fill(theta);
					//azimuthal angle in 2D s1ace
					plotCats[id_idx].hists1D[1][10]->Fill(phi);
					//subcluster energy
					plotCats[id_idx].hists1D[1][11]->Fill(E_k);
					//rotundity measures
					plotCats[id_idx].hists1D[1][12]->Fill(rot3D);
					plotCats[id_idx].hists1D[1][13]->Fill(rot2D);
					//velocity	
					plotCats[id_idx].hists1D[1][14]->Fill(vel);
					//2D eigenval ratio
					plotCats[id_idx].hists1D[1][15]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					//get variances
					plotCats[id_idx].hists1D[1][16]->Fill(sqrt(e_var));
					plotCats[id_idx].hists1D[1][17]->Fill(sqrt(p_var));
					plotCats[id_idx].hists1D[1][18]->Fill(sqrt(t_var));
					//fractional E
					plotCats[id_idx].hists1D[1][19]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (1ngle from x axis)
					plotCats[id_idx].hists1D[1][20]->Fill(phi2D);
					//pos/neg eta split sigma
					if(ec > 0) plotCats[id_idx].hists1D[1][21]->Fill(sqrt(e_var));
					else plotCats[id_idx].hists1D[1][22]->Fill(sqrt(e_var));
					//covariances
					plotCats[id_idx].hists1D[1][23]->Fill(ep_cov);
					plotCats[id_idx].hists1D[1][24]->Fill(te_cov);
					plotCats[id_idx].hists1D[1][25]->Fill(tp_cov);
					//major/minor covariances with time
					plotCats[id_idx].hists1D[1][26]->Fill(majtime_cov);
					plotCats[id_idx].hists1D[1][27]->Fill(mintime_cov);
			
					//2D[1] hists
					plotCats[id_idx].hists2D[1][0]->Fill(tc, E_k);
					plotCats[id_idx].hists2D[1][1]->Fill(phi,E_k);
					plotCats[id_idx].hists2D[1][2]->Fill(rot2D,E_k);
					plotCats[id_idx].hists2D[1][3]->Fill(ec,pc);
					plotCats[id_idx].hists2D[1][4]->Fill(tc,ec);
					plotCats[id_idx].hists2D[1][5]->Fill(tc,pc);
					//cout << "fill hist with time: " << tc << " and coeff: " << pi << endl;
					plotCats[id_idx].hists2D[1][6]->Fill(tc,pi);
					plotCats[id_idx].hists2D[1][7]->Fill(E_k,pi);
					plotCats[id_idx].hists2D[1][8]->Fill(rot3D,E_k);
					plotCats[id_idx].hists2D[1][9]->Fill(norms[k], E_k);
					plotCats[id_idx].hists2D[1][11]->Fill(nclusters, pi);
					plotCats[id_idx].hists2D[1][12]->Fill(e_var, p_var);
					plotCats[id_idx].hists2D[1][13]->Fill(t_var, e_var);
					plotCats[id_idx].hists2D[1][14]->Fill(t_var, p_var);
					plotCats[id_idx].hists2D[1][15]->Fill(E_k/E_tot, pi);
					plotCats[id_idx].hists2D[1][16]->Fill(nclusters, E_k/E_tot);
					plotCats[id_idx].hists2D[1][17]->Fill(sqrt(t_var), tc);
					plotCats[id_idx].hists2D[1][18]->Fill(rot2D, phi2D);
					plotCats[id_idx].hists2D[1][19]->Fill(sqrt(t_var), E_k/E_tot);
					plotCats[id_idx].hists2D[1][20]->Fill(sqrt(t_var), E_tot);
					plotCats[id_idx].hists2D[1][21]->Fill(te_cov, E_tot);
					plotCats[id_idx].hists2D[1][22]->Fill(tp_cov, E_tot);
					plotCats[id_idx].hists2D[1][23]->Fill(sqrt(t_var), majtime_cov);
					plotCats[id_idx].hists2D[1][24]->Fill(sqrt(t_var), mintime_cov);
					plotCats[id_idx].hists2D[1][25]->Fill(phi2D, majtime_cov);
					plotCats[id_idx].hists2D[1][26]->Fill(phi2D, mintime_cov);
					plotCats[id_idx].hists2D[1][27]->Fill(rot2D, majtime_cov);
					plotCats[id_idx].hists2D[1][28]->Fill(rot2D, mintime_cov);
					plotCats[id_idx].hists2D[1][29]->Fill(rot3D, majtime_cov);
					plotCats[id_idx].hists2D[1][30]->Fill(rot3D, mintime_cov);
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
					plotCats[id_idx].hists1D[2][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				//slope - eta/time
					plotCats[id_idx].hists1D[2][7]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					//slope - phi/time
					plotCats[id_idx].hists1D[2][8]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//polar angle in 3D space
					plotCats[id_idx].hists1D[2][9]->Fill(theta);
					//azimuthal angle in 2D s2ace
					plotCats[id_idx].hists1D[2][10]->Fill(phi);
					//subcluster energy
					plotCats[id_idx].hists1D[2][11]->Fill(E_k);
					//rotundity measures
					plotCats[id_idx].hists1D[2][12]->Fill(rot3D);
					plotCats[id_idx].hists1D[2][13]->Fill(rot2D);
					//velocity	
					plotCats[id_idx].hists1D[2][14]->Fill(vel);
					//2D eigenval ratio
					plotCats[id_idx].hists1D[2][15]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					//get variances
					plotCats[id_idx].hists1D[2][16]->Fill(sqrt(e_var));
					plotCats[id_idx].hists1D[2][17]->Fill(sqrt(p_var));
					plotCats[id_idx].hists1D[2][18]->Fill(sqrt(t_var));
					//fractional E
					plotCats[id_idx].hists1D[2][19]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (2ngle from x axis)
					plotCats[id_idx].hists1D[2][20]->Fill(phi2D);
					//pos/neg eta split sigma
					if(ec > 0) plotCats[id_idx].hists1D[2][21]->Fill(sqrt(e_var));
					else plotCats[id_idx].hists1D[2][22]->Fill(sqrt(e_var));
					//covariances
					plotCats[id_idx].hists1D[2][23]->Fill(ep_cov);
					plotCats[id_idx].hists1D[2][24]->Fill(te_cov);
					plotCats[id_idx].hists1D[2][25]->Fill(tp_cov);
					//major/minor covariances with time
					plotCats[id_idx].hists1D[2][26]->Fill(majtime_cov);
					plotCats[id_idx].hists1D[2][27]->Fill(mintime_cov);
					
					//2D[2] hists
					plotCats[id_idx].hists2D[2][0]->Fill(tc, E_k);
					plotCats[id_idx].hists2D[2][1]->Fill(phi,E_k);
					plotCats[id_idx].hists2D[2][2]->Fill(rot2D,E_k);
					plotCats[id_idx].hists2D[2][3]->Fill(ec,pc);
					plotCats[id_idx].hists2D[2][4]->Fill(tc,ec);
					plotCats[id_idx].hists2D[2][5]->Fill(tc,pc);
					plotCats[id_idx].hists2D[2][6]->Fill(tc,pi);
					plotCats[id_idx].hists2D[2][7]->Fill(E_k,pi);
					plotCats[id_idx].hists2D[2][8]->Fill(rot3D,E_k);
					plotCats[id_idx].hists2D[2][9]->Fill(norms[k], E_k);
					plotCats[id_idx].hists2D[2][11]->Fill(nclusters, pi);
					plotCats[id_idx].hists2D[2][12]->Fill(sqrt(e_var), sqrt(p_var));
					plotCats[id_idx].hists2D[2][13]->Fill(sqrt(t_var), sqrt(e_var));
					plotCats[id_idx].hists2D[2][14]->Fill(sqrt(t_var), sqrt(p_var));
					plotCats[id_idx].hists2D[2][15]->Fill(E_k/E_tot, pi);
					plotCats[id_idx].hists2D[2][16]->Fill(nclusters, E_k/E_tot);
					plotCats[id_idx].hists2D[2][17]->Fill(sqrt(t_var), tc);
					plotCats[id_idx].hists2D[2][18]->Fill(rot2D, phi2D);
					plotCats[id_idx].hists2D[2][19]->Fill(sqrt(t_var), E_k/E_tot);
					plotCats[id_idx].hists2D[2][20]->Fill(sqrt(t_var), E_tot);
					plotCats[id_idx].hists2D[2][21]->Fill(te_cov, E_tot);
					plotCats[id_idx].hists2D[2][22]->Fill(tp_cov, E_tot);
					plotCats[id_idx].hists2D[2][23]->Fill(sqrt(t_var), majtime_cov);
					plotCats[id_idx].hists2D[2][24]->Fill(sqrt(t_var), mintime_cov);
					plotCats[id_idx].hists2D[2][25]->Fill(phi2D, majtime_cov);
					plotCats[id_idx].hists2D[2][26]->Fill(phi2D, mintime_cov);
					plotCats[id_idx].hists2D[2][27]->Fill(rot2D, majtime_cov);
					plotCats[id_idx].hists2D[2][28]->Fill(rot2D, mintime_cov);
					plotCats[id_idx].hists2D[2][29]->Fill(rot3D, majtime_cov);
					plotCats[id_idx].hists2D[2][30]->Fill(rot3D, mintime_cov);


				}

			}
		}



};
#endif
