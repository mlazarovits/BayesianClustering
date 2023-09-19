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
		PhotonSkimmer();
		virtual ~PhotonSkimmer();

		//get rechits from file to cluster
		PhotonSkimmer(TFile* file);
		//ctor from rec hit collection - integrating into ntuplizer
		
		void CleaningSkim();
		void Skim();


		//stacked photon LLID hists
		//list of photon ids
		vector<TH1D*> llp_sig;
		vector<TH1D*> llp_ISR;
		vector<TH1D*> llp_notSunm;
		map<int, vector<double>> id_map;

		vector<vector<TH1D*>> histsID;
		vector<plotCat> plotCats;



		void MakeIDHists(){
			//signal
			plotCat sig;
			sig.legName = "#Chi^{0} #rightarrow #gamma";
			sig.plotName = "chiGam";
			sig.ids = {22, 32, 25, 35};
			plotCats.push_back(sig);
			//FSR
			plotCat FSR;
			//sparticle FSR
			FSR.legName = "sFSR";
			FSR.plotName = "sFSR";
			FSR.ids = {20, 30, 21, 31, 23, 33, 24, 34}; 
			plotCats.push_back(FSR);
			//notSunm
			plotCat notSunm;
			notSunm.legName = "notSunm";
			notSunm.plotName = "notSunm";
			notSunm.ids = {29, -1}; 
			plotCats.push_back(notSunm);

			//each category of histograms (histsID[j]) gets a list of histograms that is the same as the total list (hists1D)
			string name;
			for(int i = 0; i < (int)hists1D.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH1D* hist = (TH1D*)hists1D[i]->Clone();
					plotCats[j].hists1D.push_back(hist);
					name = hists1D[i]->GetName();
					name += "_"+plotCats[j].plotName;
					plotCats[j].hists1D[i]->SetName(name.c_str());
					plotCats[j].hists1D[i]->SetTitle("");
				}
		
			}
			for(int i = 0; i < (int)hists2D.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH2D* hist = (TH2D*)hists2D[i]->Clone();
					plotCats[j].hists2D.push_back(hist);
					name = hists2D[i]->GetName();
					name += "_"+plotCats[j].plotName;
					plotCats[j].hists2D[i]->SetName(name.c_str());
					plotCats[j].hists2D[i]->SetTitle("");
				}
		
			}

		}





		void WriteHists(TFile* ofile){
			vector<TH1D*> hists;
			vector<string> id_names;
			double ymax, ymin;
			string name;
		
			for(int i = 0; i < (int)plotCats.size(); i++)
				id_names.push_back(plotCats[i].legName);

			//normalize histograms	
			for(int i = 0; i < (int)plotCats.size(); i++){
				//relative fraction histograms
				//nSubClusters
				plotCats[i].hists1D[0]->Scale(1./plotCats[i].hists1D[0]->Integral());
				//ellipsoid center coordinates
				plotCats[i].hists1D[1]->Scale(1./plotCats[i].hists1D[1]->Integral());
				plotCats[i].hists1D[2]->Scale(1./plotCats[i].hists1D[2]->Integral());
				plotCats[i].hists1D[3]->Scale(1./plotCats[i].hists1D[3]->Integral());
				//theta + azimuthal angles
				plotCats[i].hists1D[7]->Scale(1./plotCats[i].hists1D[7]->Integral());
				plotCats[i].hists1D[8]->Scale(1./plotCats[i].hists1D[8]->Integral());
			}


			//write unnormalized and normalized time center hists
			time_center_norm = (TH1D*)time_center->Clone();
			time_center_norm->Scale(1./time_center_norm->Integral());	
			

			ofile->cd();
			//write 1D hists
			for(int i = 0; i < (int)hists1D.size(); i++){
				//write total hist to file
				name = hists1D[i]->GetName();
				if(hists1D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
				TCanvas* cv = new TCanvas((name).c_str(), "");
				TDRHist(hists1D[i], cv, name, name, "a.u.");	
				cv->Write();
				if(!_data){
					//make a vector for each type of histogram
					for(int j = 0; j < (int)plotCats.size(); j++){
						if(plotCats[j].hists1D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled for " << plotCats[j].plotName << endl; continue; }
						hists.push_back(plotCats[j].hists1D[i]);			
						//should be 3 hists in this vector
					}
					FindListHistBounds(hists, ymin, ymax);
					TCanvas* cv_stack = new TCanvas((name+"_stack").c_str(), "");
					TDRMultiHist(hists, cv_stack, name, name, "a.u.",ymin, ymax, id_names);
					//write cv to file			
			//		cv->SaveAs((fname+"/"+name+".pdf").c_str());
					cv_stack->Write();
					
					hists.clear();
				}
			}
			name = time_center_norm->GetName();
			TCanvas* cv = new TCanvas((name+"_norm").c_str(), "");
			TDRHist(time_center_norm, cv, name, name, "a.u.");	
			cv->Write();
			
			//write 2D hists
			string xname, yname;
			for(int i = 0; i < (int)hists2D.size(); i++){
				//write total hist to file
				name = hists2D[i]->GetName();
				name += "2D";
				xname = hists2D[i]->GetXaxis()->GetTitle();
				yname = hists2D[i]->GetYaxis()->GetTitle();
				if(hists2D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
				TCanvas* cv = new TCanvas((name).c_str(), "");
				TDR2DHist(hists2D[i], cv, xname, yname);
				cv->Write();
				if(!_data){
					//make a vector for each type of histogram
					for(int j = 0; j < (int)plotCats.size(); j++){
						if(plotCats[j].hists2D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled for " << plotCats[j].plotName << endl; continue; }
						TCanvas* cv_stack = new TCanvas((name+"_"+plotCats[j].plotName).c_str(), "");
						TDR2DHist(plotCats[j].hists2D[i], cv_stack, xname, yname, plotCats[j].plotName);
						cv_stack->Write();
					}
				}
			}
			
			ofile->Close();

		}

		//k = sum_n(E_n)/N
		void FillModelHists(BasePDFMixture* model, int id_idx, double transf = 1.){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es, eigenvals_space, npts_unwt;
			vector<Matrix> eigenvecs, eigenvecs_space;
			Matrix space_mat = Matrix(2,2);
			
			double npts = (double)model->GetData()->GetNPoints();
		//	cout << "FillHists - starting subcluster loop" << endl;	
			double E_k, theta, phi, r, rot2D, rot3D, vel, ec, pc, tc, pi, E_lead;
			double v_x, v_y, v_z;
			double E_tot = transf*npts;
			Matrix lead_eigenvec;
			
			int nclusters = model->GetNClusters();
			plotCats[id_idx].hists1D[0]->Fill(nclusters);
			model->GetAvgVarWeights(avg_Es);
			model->GetNormsUnwt(npts_unwt);
			
			//get leading cluster index
			vector<int> idxs;
			//sort by mixing coeffs in ascending order (smallest first)
			model->sortedIdxs(idxs);
			int leadidx = idxs[nclusters-1];
			int subleadidx = -999;
			if(nclusters > 1) subleadidx = idxs[nclusters-2];
			for(int k = 0; k < nclusters; k++){
				//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
				E_k = avg_Es[k]*transf*npts_unwt[k]; 
				
				params = model->GetPriorParameters(k);
				ec = params["mean"].at(0,0);
				pc = params["mean"].at(1,0);
				tc = params["mean"].at(2,0);
				pi = params["pi"].at(0,0);
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				lead_eigenvec = eigenvecs[2];			
	
				v_x = lead_eigenvec.at(0,0);	
				v_y = lead_eigenvec.at(1,0);	
				v_z = lead_eigenvec.at(2,0);	
				r = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
				//polar angle
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				theta = acos( lead_eigenvec.at(2,0) / r );
				//azimuthal angle
				//phi = arctan(y/x)
				phi = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				if(signbit(v_y)) phi *= -1;

				//rotundity - 3D
				rot3D = 0;
				for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				rot3D = eigenvals[2]/rot3D;
				
				//rotundity - 2D
				//take upper 2x2 submatrix from covariance
				space_mat.SetEntry(params["cov"].at(0,0),0,0);	
				space_mat.SetEntry(params["cov"].at(0,1),0,1);	
				space_mat.SetEntry(params["cov"].at(1,0),1,0);	
				space_mat.SetEntry(params["cov"].at(1,1),1,1);
				space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
				rot2D = 0;
				for(int i = 0; i < (int)eigenvecs_space.size(); i++) rot2D += eigenvals_space[i];
				rot2D = eigenvals_space[1]/rot2D;
				
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (lead_eigenvec.at(2,0)/r) * (acos(-1)/180.) * (1./2.2);
				vel = fabs(1./vel);

				//fill hists
				plotCats[id_idx].hists1D[1]->Fill(tc);
				plotCats[id_idx].hists1D[2]->Fill(ec);
				plotCats[id_idx].hists1D[3]->Fill(pc);
			
				//phi/eta
				plotCats[id_idx].hists1D[4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        			//eta/time
				plotCats[id_idx].hists1D[5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
				//phi/time
				plotCats[id_idx].hists1D[6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
				//angles in 3D space
				plotCats[id_idx].hists1D[7]->Fill(theta);
				plotCats[id_idx].hists1D[8]->Fill(phi);
				
				//average cluster energy
				//w_n = E_n/N for N pts in sample
				plotCats[id_idx].hists1D[9]->Fill(E_k);
				//rotundity measures
				plotCats[id_idx].hists1D[10]->Fill(rot3D);
				plotCats[id_idx].hists1D[11]->Fill(rot2D);
				plotCats[id_idx].hists1D[18]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				//velocity	
				plotCats[id_idx].hists1D[14]->Fill(vel);
				//2D hists
				plotCats[id_idx].hists2D[0]->Fill(tc, E_k);
				plotCats[id_idx].hists2D[3]->Fill(phi,E_k);
				plotCats[id_idx].hists2D[4]->Fill(rot2D,E_k);
				plotCats[id_idx].hists2D[5]->Fill(ec,pc);
				plotCats[id_idx].hists2D[6]->Fill(tc,ec);
				plotCats[id_idx].hists2D[7]->Fill(tc,pc);
				plotCats[id_idx].hists2D[10]->Fill(tc,pi);
				plotCats[id_idx].hists2D[11]->Fill(E_k,pi);
				plotCats[id_idx].hists2D[12]->Fill(rot3D,E_k);
				plotCats[id_idx].hists2D[15]->Fill(npts_unwt[k], E_tot);	

				//histograms for leading/subleading clusters
				if(k == leadidx){
					//histograms for leading clusters
					//leading cluster total energy
					plotCats[id_idx].hists1D[12]->Fill(E_k);
					//leading cluster npts
					plotCats[id_idx].hists1D[15]->Fill(npts_unwt[leadidx]);
					//fractional npts lead cluster
					plotCats[id_idx].hists1D[16]->Fill(npts_unwt[leadidx]/npts);	
					//npts_unwt_k = sum_n(r_nk)
					//w_n = sum_n(E_n)/npts
					E_lead = E_k;
					plotCats[id_idx].hists1D[17]->Fill(E_lead/E_tot);	
					plotCats[id_idx].hists2D[8]->Fill(nclusters,E_lead/E_tot);
					plotCats[id_idx].hists2D[1]->Fill(tc,E_k);
					plotCats[id_idx].hists1D[21]->Fill(rot2D);
					plotCats[id_idx].hists1D[23]->Fill(tc);
					plotCats[id_idx].hists2D[13]->Fill(E_k,pi);
					plotCats[id_idx].hists2D[16]->Fill(npts_unwt[k], E_lead/E_tot);	
	
				}
				if(nclusters == 1)
					plotCats[id_idx].hists2D[17]->Fill(npts_unwt[k], E_lead/E_tot);	


				if(nclusters > 1){
					if(k == subleadidx){
						//subleading cluster tot energy - if it exists
						plotCats[id_idx].hists1D[13]->Fill(E_k);
					}
					if(k != leadidx){
						//notleading cluster time v energy
						plotCats[id_idx].hists2D[2]->Fill(tc,E_k);
						plotCats[id_idx].hists1D[22]->Fill(rot2D);
						plotCats[id_idx].hists2D[14]->Fill(E_k,pi);
						plotCats[id_idx].hists1D[24]->Fill(tc);

					}
				}

			}
//			cout << "end clusters" << endl;
			plotCats[id_idx].hists1D[20]->Fill(E_tot);


		}

		void FillTotalHists(BasePDFMixture* model, double transf = 1.){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es, eigenvals_space, npts_unwt;
			vector<Matrix> eigenvecs, eigenvecs_space;
			Matrix space_mat = Matrix(2,2);
			double npts = (double)model->GetData()->GetNPoints();
			double E_k, theta, phi, r, rot2D, rot3D, vel, ec, pc, tc, pi, E_lead;
			double v_x, v_y, v_z;
			double E_tot = transf*npts;
			Matrix lead_eigenvec;
			
			int nclusters = model->GetNClusters();
			nSubClusters->Fill(nclusters);
			model->GetAvgVarWeights(avg_Es);
			model->GetNormsUnwt(npts_unwt);
			
			//get leading cluster index
			vector<int> idxs;
			//sort by mixing coeffs in ascending order (smallest first)
			model->sortedIdxs(idxs);
			int leadidx = idxs[nclusters-1];
			int subleadidx = -999;
			if(nclusters > 1) subleadidx = idxs[nclusters-2];
			for(int k = 0; k < nclusters; k++){
				//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
				E_k = avg_Es[k]*transf*npts_unwt[k]; 
				
				params = model->GetPriorParameters(k);
				ec = params["mean"].at(0,0);
				pc = params["mean"].at(1,0);
				tc = params["mean"].at(2,0);
				pi = params["pi"].at(0,0);
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				lead_eigenvec = eigenvecs[2];			

				v_x = lead_eigenvec.at(0,0);	
				v_y = lead_eigenvec.at(1,0);	
				v_z = lead_eigenvec.at(2,0);	
				r = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
				//polar angle
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				theta = acos( v_z / r );
				//azimuthal angle
				//phi = arctan(y/x)
				phi = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				if(signbit(v_y)) phi *= -1;

				//rotundity - 3D
				rot3D = 0;
				for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				rot3D = eigenvals[2]/rot3D;
				
				//rotundity - 2D
				//take upper 2x2 submatrix from covariance
				space_mat.SetEntry(params["cov"].at(0,0),0,0);	
				space_mat.SetEntry(params["cov"].at(0,1),0,1);	
				space_mat.SetEntry(params["cov"].at(1,0),1,0);	
				space_mat.SetEntry(params["cov"].at(1,1),1,1);
				space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
				rot2D = 0;
				for(int i = 0; i < (int)eigenvecs_space.size(); i++) rot2D += eigenvals_space[i];
				rot2D = eigenvals_space[1]/rot2D;
				
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (lead_eigenvec.at(2,0)/r) * (acos(-1)/180.) * (1./2.2);
				vel = fabs(1./vel);
				
				//fill histograms
				time_center->Fill(tc);
				eta_center->Fill(ec);
				phi_center->Fill(pc);
				
				//phi/eta
				slope_space->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        			//eta/time
				slope_etaT->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
				//phi/time
				slope_phiT->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
				//polar angle
				polar_ang->Fill(theta);
				//azimuthal angle
				azimuth_ang->Fill(phi);
				//velocity
				velocity->Fill(vel);
	
				//total cluster energy
				//w_n = N/W_n for N pts in sample
				e_tot->Fill(E_k);
				time_totE->Fill(tc, E_k);
	
				//rotundity - 3D
				rotundity_3D->Fill(rot3D);
				
				//rotundity - 2D
				rotundity_2D->Fill(rot2D);
				eigen2D_ratio->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				
			
				//2D hists
				time_totE->Fill(tc,E_k);			
				az_totE->Fill(phi,E_k);
				rot2D_totE->Fill(rot2D,E_k);
				eta_phi->Fill(ec,pc);
				t_eta->Fill(tc,ec);
				t_phi->Fill(tc,pc);
				t_mixcoeff->Fill(tc,pi);
				totE_mixcoeff->Fill(E_k,pi);
				rot3D_totE->Fill(rot3D,E_k);
				npts_totE->Fill(npts_unwt[k], E_tot);
			

				//leading subcluster hists
				if(k == leadidx){
					//leading cluster tot energy
					e_tot_lead->Fill(E_k);
					//leading cluster npts
					npts_lead->Fill(npts_unwt[k]);
					fracpts_lead->Fill(npts_unwt[k]/npts);	
					//npts_unwt_k = sum_n(r_nk)
					//w_n = sum_n(E_n)/npts
					E_lead = E_k;
					if(nclusters == 1 && E_lead/E_tot < 0.9) cout << "E_lead: " << E_lead << " E_tot: " << E_tot << " nptsunwt: " << npts_unwt[k] << endl;
					fracE_lead->Fill(E_lead/E_tot);	
					nsubcl_fracElead->Fill(nclusters,E_lead/E_tot);
					//leading cluster time v energy
					time_totE_lead->Fill(tc,E_lead);
					//2D rotundity - lead
					rotundity_2D_lead->Fill(rot2D);		
					time_center_lead->Fill(tc);
					totE_mixcoeff_lead->Fill(E_k,pi);
					npts_fracE_lead->Fill(npts_unwt[k], E_lead/E_tot);
				}
				if(nclusters == 1)
					npts_fracE_lead_1subcl->Fill(npts_unwt[k], E_lead/E_tot);	
				if(nclusters > 1){
					//sublead cluster
					if(k == subleadidx){
						//subleading cluster tot energy - if it exists
						e_tot_sublead->Fill(E_k);
					}


					//not lead cluster
					if(k != leadidx){
						//notleading cluster time v tot energy
						time_totE_notlead->Fill(tc,E_k);
						rotundity_2D_notlead->Fill(rot2D);		
						totE_mixcoeff_notlead->Fill(E_k,pi);
						time_center_notlead->Fill(tc);
					}


				}
	

			}
			clusterE->Fill(E_tot);

	
		}

	private:
		TH1D* time_center_norm = nullptr;





};
#endif
