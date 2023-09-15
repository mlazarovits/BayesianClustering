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
			//ISR
			plotCat ISR;
			ISR.legName = "FSR";
			ISR.plotName = "FSR";
			ISR.ids = {20, 30, 21, 31, 23, 33, 24, 34}; 
			plotCats.push_back(ISR);
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
			TH1D* time_center_norm = (TH1D*)time_center->Clone();
			time_center_norm->Scale(1./time_center_norm->Integral());	
			

			ofile->cd();
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
						if(plotCats[j].hists1D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
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
			TCanvas* cv = new TCanvas((name).c_str(), "");
			TDRHist(time_center_norm, cv, name, name, "a.u.");	
			cv->Write();
			
			ofile->Close();

		}


		void FillHists(BasePDFMixture* model, int id_idx, double w_n = 1.){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es, eigenvals_space, npts_unwt;
			vector<Matrix> eigenvecs, eigenvecs_space;
			Matrix space_mat = Matrix(2,2);

			int nclusters = model->GetNClusters();
			plotCats[id_idx].hists1D[0]->Fill(nclusters);
			model->GetAvgVarWeights(avg_Es);
			model->GetNormsUnwt(npts_unwt);
			
			double npts = (double)model->GetData()->GetNPoints();
			double E_tot = 0;
			double E_k;
			cout << "FillHists - starting subcluster loop" << endl;	
			double theta, phi, r, rot2D, rot3D, vel;
			for(int k = 0; k < nclusters; k++){
				params = model->GetParameters(k);
				plotCats[id_idx].hists1D[1]->Fill(params["mean"].at(2,0));
				plotCats[id_idx].hists1D[2]->Fill(params["mean"].at(0,0));
				plotCats[id_idx].hists1D[3]->Fill(params["mean"].at(1,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//largest eigenvalue is last
				//phi/eta
				plotCats[id_idx].hists1D[4]->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(0,0));
        			//eta/time
				plotCats[id_idx].hists1D[5]->Fill(eigenvecs[2].at(0,0)/eigenvecs[2].at(2,0));
				//phi/time
				plotCats[id_idx].hists1D[6]->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(2,0));
				//polar angle
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				r = sqrt(eigenvecs[2].at(0,0)*eigenvecs[2].at(0,0) + eigenvecs[2].at(1,0)*eigenvecs[2].at(1,0) + eigenvecs[2].at(2,0)*eigenvecs[2].at(2,0));
				theta = acos( eigenvecs[2].at(2,0) / r );
				plotCats[id_idx].hists1D[7]->Fill(theta);
				//azimuthal angle
				//phi = arctan(y/x)
				phi = atan2(eigenvecs[2].at(1,0) , eigenvecs[2].at(0,0));
				plotCats[id_idx].hists1D[8]->Fill(phi);
				
				//average cluster energy
				//w_n = E_n/N for N pts in sample
				plotCats[id_idx].hists1D[9]->Fill(avg_Es[k]/w_n);
				plotCats[id_idx].hists2D[0]->Fill(params["mean"].at(2,0), avg_Es[k]/w_n);
				E_tot += avg_Es[k]/w_n*npts_unwt[k]; 
			
				//rotundity - 3D
				for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				rot3D = eigenvals[2]/rot3D;
				plotCats[id_idx].hists1D[10]->Fill(rot3D);
				
				//rotundity - 2D
				//take upper 2x2 submatrix from covariance
				space_mat.SetEntry(params["cov"].at(0,0),0,0);	
				space_mat.SetEntry(params["cov"].at(0,1),0,1);	
				space_mat.SetEntry(params["cov"].at(1,0),1,0);	
				space_mat.SetEntry(params["cov"].at(1,1),1,1);
				space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
	
				for(int i = 0; i < (int)eigenvecs_space.size(); i++) rot2D += eigenvals_space[i];
				rot2D = eigenvals_space[1]/rot2D;
				plotCats[id_idx].hists1D[11]->Fill(rot2D);
				plotCats[id_idx].hists1D[18]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
		
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (eigenvecs[2].at(2,0)/r) * (acos(-1)/180.) * (1./2.2);
				vel = 1./vel;
				plotCats[id_idx].hists1D[14]->Fill(vel);
				cout << "filling 2d hists per cluster" << endl;
			
				//2D hists
				plotCats[id_idx].hists2D[0]->Fill(eigenvecs[2].at(2,0),E_k);			
				plotCats[id_idx].hists2D[3]->Fill(phi,E_k);
				plotCats[id_idx].hists2D[4]->Fill(rot2D,E_k);
				plotCats[id_idx].hists2D[5]->Fill(eigenvecs[2].at(0,0),eigenvecs[2].at(1,0));
				plotCats[id_idx].hists2D[6]->Fill(eigenvecs[2].at(2,0),eigenvecs[2].at(0,0));
				plotCats[id_idx].hists2D[7]->Fill(eigenvecs[2].at(2,0),eigenvecs[2].at(1,0));


			}
			cout << "end clusters" << endl;

			vector<int> idxs;
			model->sortedIdxs(idxs);
			//leading cluster avg energy
			plotCats[id_idx].hists1D[12]->Fill(avg_Es[idxs[0]]/w_n);
			//leading cluster npts
			plotCats[id_idx].hists1D[15]->Fill(npts_unwt[idxs[0]]);
			//fractional npts lead cluster
			plotCats[id_idx].hists1D[16]->Fill(npts_unwt[idxs[0]]/npts);	
			//avgE_k = sum_n(E_n*r_nk*w)/sum_n(r_nk)
			//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
			//npts_unwt_k = sum_n(r_nk)
			//E_tot = sum_k(E_k) = sum_k(sum_n(E_n*r_nk)) = sum_n(E_n*sum_k(r_nk)) = sum_n(E_n)
			double E_lead = avg_Es[idxs[0]]/w_n*npts_unwt[idxs[0]];
			plotCats[id_idx].hists1D[17]->Fill(E_lead/E_tot);	
			
			//leading cluster time v energy
			params = model->GetPriorParameters(idxs[0]);
			plotCats[id_idx].hists2D[1]->Fill(params["mean"].at(2,0),avg_Es[idxs[0]]/w_n);
			if(nclusters > 1){
				//subleading cluster avg energy - if it exists
				plotCats[id_idx].hists1D[13]->Fill(avg_Es[idxs[1]]/w_n);
				//leading cluster time v energy
				params = model->GetPriorParameters(idxs[1]);
				plotCats[id_idx].hists2D[2]->Fill(params["mean"].at(2,0),avg_Es[idxs[1]]/w_n*npts_unwt[idxs[1]]);
			}

		}

		void FillTotalHists(BasePDFMixture* model, double w_n = 1.){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es, eigenvals_space, npts_unwt;
			vector<Matrix> eigenvecs, eigenvecs_space;
			int nclusters = model->GetNClusters();
			nSubClusters->Fill(nclusters);
		
	
			model->GetAvgVarWeights(avg_Es);
			model->GetNormsUnwt(npts_unwt);
			double npts = (double)model->GetData()->GetNPoints();
			double E_tot = 0;
			double E_k;
			Matrix space_mat = Matrix(2,2);
	
			double theta, phi, r, rot2D, rot3D, vel;
			for(int k = 0; k < nclusters; k++){
				params = model->GetParameters(k);
				time_center->Fill(params["mean"].at(2,0));
				eta_center->Fill(params["mean"].at(0,0));
				phi_center->Fill(params["mean"].at(1,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//largest eigenvalue is last
				//phi/eta
				slope_space->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(0,0));
        			//eta/time
				slope_etaT->Fill(eigenvecs[2].at(0,0)/eigenvecs[2].at(2,0));
				//phi/time
				slope_phiT->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(2,0));
				//polar angle
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				r = sqrt(eigenvecs[2].at(0,0)*eigenvecs[2].at(0,0) + eigenvecs[2].at(1,0)*eigenvecs[2].at(1,0) + eigenvecs[2].at(2,0)*eigenvecs[2].at(2,0));
				theta = acos( eigenvecs[2].at(2,0) / r );
				polar_ang->Fill(theta);
		
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (eigenvecs[2].at(2,0)/r) * (acos(-1)/180.) * (1./2.2);
				vel = 1./vel;
				velocity->Fill(vel);

				//azimuthal angle
				//phi = arctan(y/x)
				phi = atan2(eigenvecs[2].at(1,0) , eigenvecs[2].at(0,0));
				azimuth_ang->Fill(phi);
				
				//average cluster energy
				//w_n = N/W_n for N pts in sample
				e_avg->Fill(avg_Es[k]/w_n);
				E_k = avg_Es[k]/w_n*npts_unwt[k]; 
				E_tot += E_k;
				time_totE->Fill(params["mean"].at(2,0), E_k);
	
				//rotundity - 3D
				for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				rot3D = eigenvals[2]/rot3D;
				rotundity_3D->Fill(rot3D);
				
				//rotundity - 2D
				//take upper 2x2 submatrix from covariance
				space_mat.SetEntry(params["cov"].at(0,0),0,0);	
				space_mat.SetEntry(params["cov"].at(0,1),0,1);	
				space_mat.SetEntry(params["cov"].at(1,0),1,0);	
				space_mat.SetEntry(params["cov"].at(1,1),1,1);
				space_mat.eigenCalc(eigenvals_space, eigenvecs_space);

				for(int i = 0; i < (int)eigenvecs_space.size(); i++) rot2D += eigenvals_space[i];
				rot2D = eigenvals_space[1]/rot2D;
				rotundity_2D->Fill(rot2D);
				eigen2D_ratio->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (eigenvecs[2].at(2,0)/r) * (acos(-1)/180.) * (1./2.2);
				vel = 1./vel;
				velocity->Fill(vel);
			
				//2D hists
				time_totE->Fill(eigenvecs[2].at(2,0),E_k);			
				az_totE->Fill(phi,E_k);
				rot2D_totE->Fill(rot2D,E_k);
				eta_phi->Fill(eigenvecs[2].at(0,0),eigenvecs[2].at(1,0));
				t_eta->Fill(eigenvecs[2].at(2,0),eigenvecs[2].at(0,0));
				t_phi->Fill(eigenvecs[2].at(2,0),eigenvecs[2].at(1,0));
				

			}
			vector<int> idxs;
			model->sortedIdxs(idxs);
			//leading cluster avg energy
			e_avg_lead->Fill(avg_Es[idxs[0]]/w_n);
			//leading cluster npts
			npts_lead->Fill(npts_unwt[idxs[0]]);
			fracpts_lead->Fill(npts_unwt[idxs[0]]/npts);	
			//avgE_k = sum_n(E_n*r_nk*w)/sum_n(r_nk)
			//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
			//npts_unwt_k = sum_n(r_nk)
			//E_tot = sum_k(E_k) = sum_k(sum_n(E_n*r_nk)) = sum_n(E_n*sum_k(r_nk)) = sum_n(E_n)
			double E_lead = avg_Es[idxs[0]]/w_n*npts_unwt[idxs[0]];
			fracE_lead->Fill(E_lead/E_tot);	

			//2D hists
			//leading cluster time v energy
			params = model->GetPriorParameters(idxs[0]);
			time_totE_lead->Fill(params["mean"].at(2,0),E_lead);
			if(nclusters > 1){
				//subleading cluster avg energy - if it exists
				e_avg_sublead->Fill(avg_Es[idxs[1]]/w_n);
				//leading cluster time v energy
				params = model->GetPriorParameters(idxs[1]);
				time_totE_sublead->Fill(params["mean"].at(2,0),avg_Es[idxs[1]]/w_n*npts_unwt[idxs[1]]);
			}
		}








};
#endif
