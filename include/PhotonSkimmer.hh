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

		void SetIsoCuts(){ _isocuts = true; }
		bool _isocuts;
		//set skip for outstream
		void SetSkip(int i){ _oskip = i; }
		int _oskip;
		void SetThresh(double t){ _thresh = t; }
		void SetBHCAlpha(double a){ _alpha = a; }
		void SetEMAlpha(double a){ _emAlpha = a; }
		double _thresh, _alpha, _emAlpha;

		//stacked photon LLID hists
		//list of photon ids
		vector<TH1D*> llp_sig;
		vector<TH1D*> llp_ISR;
		vector<TH1D*> llp_notSunm;
		map<int, vector<double>> id_map;

		vector<vector<TH1D*>> histsID;
		vector<plotCat> plotCats;

		void MakeIDHists(){
			//total
			plotCat tot;
			tot.legName = "";
			tot.plotName = "";
			tot.ids = {-999};
			plotCats.push_back(tot);	
		
			//signal
			plotCat sig;
			sig.legName = "#Chi^{0} #rightarrow #gamma";
			sig.plotName = "chiGam";
			sig.ids = {22};
			plotCats.push_back(sig);
			//FSR
			//plotCat FSR;
			//FSR.legName = "sFSR";
			//FSR.plotName = "sFSR";
			//FSR.ids = {20, 30, 21, 31, 23, 33, 24, 34}; 
			//plotCats.push_back(FSR);
			//notSunm
			plotCat notSunm;
			notSunm.legName = "bkg";
			notSunm.plotName = "bkg";
			//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
			notSunm.ids = {29, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
			//notSunm.ids = {29, -1}; 
			plotCats.push_back(notSunm);

			//each category of histograms (histsID[j]) gets a list of histograms that is the same as the total list (hists1D)
			string name;
			for(int i = 0; i < (int)hists1D.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH1D* hist = (TH1D*)hists1D[i]->Clone();
					plotCats[j].hists1D.push_back(hist);
					name = hists1D[i]->GetName();
					if(j != 0) name += "_"+plotCats[j].plotName;
					plotCats[j].hists1D[i]->SetName(name.c_str());
					if(j != 0) plotCats[j].hists1D[i]->SetTitle("");
				}
		
			}
			for(int i = 0; i < (int)hists1D_lead.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH1D* hist = (TH1D*)hists1D_lead[i]->Clone();
					plotCats[j].hists1D_lead.push_back(hist);
					name = hists1D_lead[i]->GetName();
					if(j != 0) name += "_"+plotCats[j].plotName;
					plotCats[j].hists1D_lead[i]->SetName(name.c_str());
					if(j != 0) plotCats[j].hists1D_lead[i]->SetTitle("");
				}
			}
			for(int i = 0; i < (int)hists1D_notlead.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH1D* hist = (TH1D*)hists1D_notlead[i]->Clone();
					plotCats[j].hists1D_notlead.push_back(hist);
					name = hists1D_notlead[i]->GetName();
					if(j != 0) name += "_"+plotCats[j].plotName;
					plotCats[j].hists1D_notlead[i]->SetName(name.c_str());
					if(j != 0) plotCats[j].hists1D_notlead[i]->SetTitle("");
				}
			}
			for(int i = 0; i < (int)hists2D.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH2D* hist = (TH2D*)hists2D[i]->Clone();
					plotCats[j].hists2D.push_back(hist);
					name = hists2D[i]->GetName();
					if(j != 0) name += "_"+plotCats[j].plotName;
					plotCats[j].hists2D[i]->SetName(name.c_str());
					if(j != 0) plotCats[j].hists2D[i]->SetTitle("");
				}
		
			}
			for(int i = 0; i < (int)hists2D_lead.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH2D* hist = (TH2D*)hists2D_lead[i]->Clone();
					plotCats[j].hists2D_lead.push_back(hist);
					name = hists2D_lead[i]->GetName();
					if(j != 0) name += "_"+plotCats[j].plotName;
					plotCats[j].hists2D_lead[i]->SetName(name.c_str());
					if(j != 0) plotCats[j].hists2D_lead[i]->SetTitle("");
				}
			}
			for(int i = 0; i < (int)hists2D_notlead.size(); i++){
				for(int j = 0; j < (int)plotCats.size(); j++){
					TH2D* hist = (TH2D*)hists2D_notlead[i]->Clone();
					plotCats[j].hists2D_notlead.push_back(hist);
					name = hists2D_notlead[i]->GetName();
					if(j != 0) name += "_"+plotCats[j].plotName;
					plotCats[j].hists2D_notlead[i]->SetName(name.c_str());
					if(j != 0) plotCats[j].hists2D_notlead[i]->SetTitle("");
				}
			}

		}

		void WriteHist1D(int i){
			vector<TH1D*> hists;
			double ymin, ymax;
			vector<string> id_names;
			//write total hist to file
			string name = plotCats[0].hists1D[i]->GetName();
			if(plotCats[0].hists1D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; return; }
			TCanvas* cv = new TCanvas((name).c_str(), "");
			TDRHist(plotCats[0].hists1D[i], cv, name, name, "a.u.");	
			cv->Write();
			if(!_data){
				//make a vector for each type of histogram
				for(int j = 1; j < (int)plotCats.size(); j++){
					if(plotCats[j].hists1D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled for " << plotCats[j].plotName << endl; continue; }
					hists.push_back(plotCats[j].hists1D[i]);			
					//should be 3 hists in this vector
				}
				FindListHistBounds(hists, ymin, ymax);
				TCanvas* cv_stack = new TCanvas((name+"_stack").c_str(), "");
				TDRMultiHist(hists, cv_stack, name, name, "a.u.",ymin, ymax, id_names);
				//write cv to file			
				cv_stack->Write();
				
				hists.clear();
			}
		}
			




		void WriteHists(TFile* ofile){
			vector<TH1D*> hists;
			vector<string> id_names;
			double ymax, ymin;
			string name;
		
			for(int i = 1; i < (int)plotCats.size(); i++)
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
				//rotundity
				plotCats[i].hists1D[10]->Scale(1./plotCats[i].hists1D[10]->Integral());
				plotCats[i].hists1D[11]->Scale(1./plotCats[i].hists1D[11]->Integral());
				//velocity
				plotCats[i].hists1D[14]->Scale(1./plotCats[i].hists1D[14]->Integral());
			}


			ofile->cd();
			//write 1D hists
			for(int i = 0; i < (int)hists1D.size(); i++){
				WriteHist1D(i);
			}
			
			//write 2D hists
			string xname, yname;
			for(int i = 0; i < (int)hists2D.size(); i++){
				//write total hist to file
				name = plotCats[0].hists2D[i]->GetName();
				name += "2D";
				xname = plotCats[0].hists2D[i]->GetXaxis()->GetTitle();
				yname = plotCats[0].hists2D[i]->GetYaxis()->GetTitle();
				if(plotCats[0].hists2D[i]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
				TCanvas* cv = new TCanvas((name).c_str(), "");
				TDR2DHist(plotCats[0].hists2D[i], cv, xname, yname);
				cv->Write();
				if(!_data){
					//make a vector for each type of histogram
					for(int j = 1; j < (int)plotCats.size(); j++){
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
		void FillModelHists(BasePDFMixture* model, int id_idx){
			map<string, Matrix> params;
			vector<double> eigenvals, eigenvals_space, norms;
			vector<Matrix> eigenvecs, eigenvecs_space;
			Matrix space_mat = Matrix(2,2);
			
			double npts = (double)model->GetData()->GetNPoints();
		//	cout << "FillHists - starting subcluster loop" << endl;	
			double E_k, theta, phi, r, rot2D, rot3D, vel, ec, pc, tc, pi, E_lead, phi2D;
			double v_x, v_y, v_z;
			double E_tot = 0.;
			for(int i = 0; i < npts; i++){
				E_tot += model->GetData()->at(i).w()/_gev;
			}
			Matrix lead_eigenvec, lead_eigenvec_space;
			
			int nclusters = model->GetNClusters();
			plotCats[id_idx].hists1D[0]->Fill(nclusters);
			plotCats[id_idx].hists2D[10]->Fill((double)nclusters,npts);
			model->GetNorms(norms);
			
			//get leading cluster index
			vector<int> idxs;
			//sort by mixing coeffs in ascending order (smallest first)
			model->SortIdxs(idxs);
			int leadidx = idxs[nclusters-1];
			int subleadidx = -999;
			if(nclusters > 1) subleadidx = idxs[nclusters-2];

			for(int k = 0; k < nclusters; k++){
				//E_k = sum_n(E_n*r_nk) -> avgE/w*sum_n(r_nk)
				E_k = norms[k]/_gev; 
				
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
				//polar angle with lead eigenvector
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				theta = acos( v_z / r );
				//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
				phi = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				if(signbit(v_y)) phi *= -1;

				//rotundity - 3D
				rot3D = 0;
				for(int i = 0; i < (int)eigenvecs.size(); i++) rot3D += eigenvals[i];
				rot3D = eigenvals[2]/rot3D;
				//velocity = z/r * rad/deg * deg/cm => ns/cm
				vel = (lead_eigenvec.at(2,0)/sqrt(v_x*v_x + v_y*v_y)) * (acos(-1)/180.) * (1./2.2);
				vel = fabs(1./vel);
				
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
				lead_eigenvec_space = eigenvecs_space[1];			
				v_x = lead_eigenvec_space.at(0,0);	
				v_y = lead_eigenvec_space.at(1,0);	
				r = sqrt(v_x*v_x + v_y*v_y);
				//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
				phi2D = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				if(signbit(v_y)) phi2D *= -1;
				

				//fill hists
				//centers
				plotCats[id_idx].hists1D[1]->Fill(tc);
				plotCats[id_idx].hists1D[2]->Fill(ec);
				plotCats[id_idx].hists1D[3]->Fill(pc);
				//phi/eta
				plotCats[id_idx].hists1D[4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        			//eta/time
				plotCats[id_idx].hists1D[5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
				//phi/time
				plotCats[id_idx].hists1D[6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
				//angle in 3D space
				plotCats[id_idx].hists1D[7]->Fill(theta);
				//angle in 2D space
				plotCats[id_idx].hists1D[8]->Fill(phi);
				//average cluster energy
				//w_n = E_n/N for N pts in sample
				plotCats[id_idx].hists1D[9]->Fill(E_k);
				//rotundity measures
				plotCats[id_idx].hists1D[10]->Fill(rot3D);
				plotCats[id_idx].hists1D[11]->Fill(rot2D);
				//velocity	
				plotCats[id_idx].hists1D[12]->Fill(vel);
				plotCats[id_idx].hists1D[13]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				//cluster E
				plotCats[id_idx].hists1D[14]->Fill(E_tot);
				//get variances
				plotCats[id_idx].hists1D[15]->Fill(params["cov"].at(0,0));
				plotCats[id_idx].hists1D[16]->Fill(params["cov"].at(1,1));
				plotCats[id_idx].hists1D[17]->Fill(params["cov"].at(2,2));
				//fractional E
				plotCats[id_idx].hists1D[18]->Fill(E_k/E_tot);
				//"azimuth" angle in 2D (angle from x axis)
				plotCats[id_idx].hists1D[20]->Fill(phi2D);
	
				//2D hists
				plotCats[id_idx].hists2D[0]->Fill(tc, E_k);
				plotCats[id_idx].hists2D[1]->Fill(phi,E_k);
				plotCats[id_idx].hists2D[2]->Fill(rot2D,E_k);
				plotCats[id_idx].hists2D[3]->Fill(ec,pc);
				plotCats[id_idx].hists2D[4]->Fill(tc,ec);
				plotCats[id_idx].hists2D[5]->Fill(tc,pc);
				plotCats[id_idx].hists2D[6]->Fill(tc,pi);
				plotCats[id_idx].hists2D[7]->Fill(E_k,pi);
				plotCats[id_idx].hists2D[8]->Fill(rot3D,E_k);
				plotCats[id_idx].hists2D[9]->Fill(norms[k], E_k);
				plotCats[id_idx].hists2D[11]->Fill(nclusters, pi);
				plotCats[id_idx].hists2D[12]->Fill(params["cov"].at(0,0), params["cov"].at(1,1));
				plotCats[id_idx].hists2D[13]->Fill(params["cov"].at(2,2), params["cov"].at(0,0));
				plotCats[id_idx].hists2D[14]->Fill(params["cov"].at(2,2), params["cov"].at(1,1));
				plotCats[id_idx].hists2D[15]->Fill(E_k, pi);
				plotCats[id_idx].hists2D[16]->Fill(nclusters, E_k);

				//histograms for leading/subleading clusters
				if(k == leadidx){
					//centers
					plotCats[id_idx].hists1D_lead[1]->Fill(tc);
					plotCats[id_idx].hists1D_lead[2]->Fill(ec);
					plotCats[id_idx].hists1D_lead[3]->Fill(pc);
					//phi/eta
					plotCats[id_idx].hists1D_lead[4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				//eta/time
					plotCats[id_idx].hists1D_lead[5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					//phi/time
					plotCats[id_idx].hists1D_lead[6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//angle in 3D_lead space
					plotCats[id_idx].hists1D_lead[7]->Fill(theta);
					//angle in 2D_lead space
					plotCats[id_idx].hists1D_lead[8]->Fill(phi);
					//average cluster energy
					//w_n = E_n/N for N pts in sample
					plotCats[id_idx].hists1D_lead[9]->Fill(E_k);
					//rotundity measures
					plotCats[id_idx].hists1D_lead[10]->Fill(rot3D);
					plotCats[id_idx].hists1D_lead[11]->Fill(rot2D);
					//velocity	
					plotCats[id_idx].hists1D_lead[12]->Fill(vel);
					plotCats[id_idx].hists1D_lead[13]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					
					//cluster E
					plotCats[id_idx].hists1D_lead[14]->Fill(E_tot);
					
					//get variances
					plotCats[id_idx].hists1D_lead[15]->Fill(params["cov"].at(0,0));
					plotCats[id_idx].hists1D_lead[16]->Fill(params["cov"].at(1,1));
					plotCats[id_idx].hists1D_lead[17]->Fill(params["cov"].at(2,2));
					//fractional E
					plotCats[id_idx].hists1D_lead[18]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (angle from x axis)
					plotCats[id_idx].hists1D[20]->Fill(phi2D);
			
					//2D_lead hists
					plotCats[id_idx].hists2D_lead[0]->Fill(tc, E_k);
					plotCats[id_idx].hists2D_lead[1]->Fill(phi,E_k);
					plotCats[id_idx].hists2D_lead[2]->Fill(rot2D,E_k);
					plotCats[id_idx].hists2D_lead[3]->Fill(ec,pc);
					plotCats[id_idx].hists2D_lead[4]->Fill(tc,ec);
					plotCats[id_idx].hists2D_lead[5]->Fill(tc,pc);
					plotCats[id_idx].hists2D_lead[6]->Fill(tc,pi);
					plotCats[id_idx].hists2D_lead[7]->Fill(E_k,pi);
					plotCats[id_idx].hists2D_lead[8]->Fill(rot3D,E_k);
					plotCats[id_idx].hists2D_lead[9]->Fill(norms[k], E_k);
					plotCats[id_idx].hists2D_lead[11]->Fill(nclusters, pi);
					plotCats[id_idx].hists2D_lead[12]->Fill(params["cov"].at(0,0), params["cov"].at(1,1));
					plotCats[id_idx].hists2D_lead[13]->Fill(params["cov"].at(2,2), params["cov"].at(0,0));
					plotCats[id_idx].hists2D_lead[14]->Fill(params["cov"].at(2,2), params["cov"].at(1,1));
					plotCats[id_idx].hists2D_lead[15]->Fill(E_k, pi);
					plotCats[id_idx].hists2D_lead[16]->Fill(nclusters, E_k);
				}


				else if(k != leadidx){
					//centers
					plotCats[id_idx].hists1D_notlead[1]->Fill(tc);
					plotCats[id_idx].hists1D_notlead[2]->Fill(ec);
					plotCats[id_idx].hists1D_notlead[3]->Fill(pc);
					//phi/eta
					plotCats[id_idx].hists1D_notlead[4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				//eta/time
					plotCats[id_idx].hists1D_notlead[5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					//phi/time
					plotCats[id_idx].hists1D_notlead[6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//angle in 3D_notlead space
					plotCats[id_idx].hists1D_notlead[7]->Fill(theta);
					//angle in 2D_notlead space
					plotCats[id_idx].hists1D_notlead[8]->Fill(phi);
					//average cluster energy
					//w_n = E_n/N for N pts in sample
					plotCats[id_idx].hists1D_notlead[9]->Fill(E_k);
					//rotundity measures
					plotCats[id_idx].hists1D_notlead[10]->Fill(rot3D);
					plotCats[id_idx].hists1D_notlead[11]->Fill(rot2D);
					//velocity	
					plotCats[id_idx].hists1D_notlead[12]->Fill(vel);
					plotCats[id_idx].hists1D_notlead[13]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					
					//cluster E
					plotCats[id_idx].hists1D_notlead[14]->Fill(E_tot);
					
					//get variances
					plotCats[id_idx].hists1D_notlead[15]->Fill(params["cov"].at(0,0));
					plotCats[id_idx].hists1D_notlead[16]->Fill(params["cov"].at(1,1));
					plotCats[id_idx].hists1D_notlead[17]->Fill(params["cov"].at(2,2));
					//fractional E
					plotCats[id_idx].hists1D_notlead[18]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (angle from x axis)
					plotCats[id_idx].hists1D[20]->Fill(phi2D);
			
					//2D_notlead hists
					plotCats[id_idx].hists2D_notlead[0]->Fill(tc, E_k);
					plotCats[id_idx].hists2D_notlead[1]->Fill(phi,E_k);
					plotCats[id_idx].hists2D_notlead[2]->Fill(rot2D,E_k);
					plotCats[id_idx].hists2D_notlead[3]->Fill(ec,pc);
					plotCats[id_idx].hists2D_notlead[4]->Fill(tc,ec);
					plotCats[id_idx].hists2D_notlead[5]->Fill(tc,pc);
					plotCats[id_idx].hists2D_notlead[6]->Fill(tc,pi);
					plotCats[id_idx].hists2D_notlead[7]->Fill(E_k,pi);
					plotCats[id_idx].hists2D_notlead[8]->Fill(rot3D,E_k);
					plotCats[id_idx].hists2D_notlead[9]->Fill(norms[k], E_k);
					plotCats[id_idx].hists2D_notlead[11]->Fill(nclusters, pi);
					plotCats[id_idx].hists2D_notlead[12]->Fill(params["cov"].at(0,0), params["cov"].at(1,1));
					plotCats[id_idx].hists2D_notlead[13]->Fill(params["cov"].at(2,2), params["cov"].at(0,0));
					plotCats[id_idx].hists2D_notlead[14]->Fill(params["cov"].at(2,2), params["cov"].at(1,1));
					plotCats[id_idx].hists2D_notlead[15]->Fill(E_k, pi);
					plotCats[id_idx].hists2D_notlead[16]->Fill(nclusters, E_k);

				}

			}
//			cout << "end clusters" << endl;


		}


};
#endif
