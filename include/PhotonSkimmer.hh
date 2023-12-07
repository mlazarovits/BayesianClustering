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
			plotCat tot(_hists1D, _hists2D);
			tot.ids = {-999};
			plotCats.push_back(tot);	
		
			//signal
			plotCat sig(_hists1D, _hists2D, "chiGam","#Chi^{0} #rightarrow #gamma");
			sig.ids = {22};
			plotCats.push_back(sig);
			//FSR
			//plotCat FSR;
			//FSR.legName = "sFSR";
			//FSR.plotName = "sFSR";
			//FSR.ids = {20, 30, 21, 31, 23, 33, 24, 34}; 
			//plotCats.push_back(FSR);
			//notSunm
			plotCat notSunm(_hists1D, _hists2D, "bkg","bkg");
			//bkg is id < 9 but anything other than -1 shouldn't happen but just to be safe
			notSunm.ids = {29, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8}; 
			//notSunm.ids = {29, -1}; 
			plotCats.push_back(notSunm);

		}
	


		void WritePlotCat(TFile* ofile, const plotCat& pc){
			ofile->cd();
			string name;
			vector<vector<TH1D*>> hists1D = pc.hists1D;
			//write 1D hists
			for(int i = 0; i < (int)hists1D.size(); i++){
				for(int j = 0; j < hists1D[i].size(); j++){
					if(pc.hists1D[i][j] == nullptr) continue;
					name = pc.hists1D[i][j]->GetName();
					if(pc.hists1D[i][j]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
					TCanvas* cv = new TCanvas((name).c_str(), "");
					TDRHist(pc.hists1D[i][j], cv, name, name, "a.u.");	
					cv->Write();
					//or instead write hist here
				}
			}
			
			//write 2D hists
			string xname, yname;
			vector<vector<TH2D*>> hists2D = pc.hists2D;
			for(int i = 0; i < (int)hists2D.size(); i++){
				for(int j = 0; j < hists2D[i].size(); j++){
					if(pc.hists2D[i][j] == nullptr) continue;
					//write total hist to file
					name = pc.hists2D[i][j]->GetName();
					name += "2D";
					xname = pc.hists2D[i][j]->GetXaxis()->GetTitle();
					yname = pc.hists2D[i][j]->GetYaxis()->GetTitle();
					if(pc.hists2D[i][j]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
					TCanvas* cv = new TCanvas((name).c_str(), "");
					TDR2DHist(pc.hists2D[i][j], cv, xname, yname, pc.plotName);
					cv->Write();
				}
			}

		}



		void WritePlotCats(TFile* ofile, const vector<plotCat>& pcs){
			ofile->cd();
			string name;
			double ymax, ymin;
			vector<string> id_names;
			cout << "writing " << pcs.size() << " to stack plots" << endl;
			//number of histogram categories (ie leading, !leading, etc)
			int nhistCats = pcs[0].hists1D.size();
			//write 1D hists
			for(int j = 0; j < (int)_hists1D.size(); j++){
				for(int i = 0; i < nhistCats; i++){
					vector<TH1D*> hists;
					if(pcs[0].hists1D[i][j] == nullptr) continue;
					name = pcs[0].hists1D[i][j]->GetName();
					for(int k = 0; k < pcs.size(); k++){
						if(pcs[k].hists1D[i][j] == nullptr) continue;
						if(pcs[k].hists1D[i][j]->Integral() == 0){ cout << "Histogram: " << name << " not filled." << endl; continue; }
						hists.push_back(pcs[k].hists1D[i][j]);			
						id_names.push_back(pcs[i].legName);
					}
					FindListHistBounds(hists, ymin, ymax);
					TCanvas* cv_stack = new TCanvas((name+"_stack").c_str(), "");
					TDRMultiHist(hists, cv_stack, name, name, "a.u.",ymin, ymax, id_names);
					//write cv to file			
					cv_stack->Write();
					//or write hist here	
					hists.clear();
					id_names.clear();
				}
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
				for(int j = 0; j < plotCats[i].hists1D.size(); j++){
					//relative fraction histograms
					//nSubClusters
					plotCats[i].hists1D[j][0]->Scale(1./plotCats[i].hists1D[j][0]->Integral());
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

			WritePlotCat(ofile, plotCats[0]);
			vector<plotCat> id_cats(plotCats.begin()+1, plotCats.end());
			WritePlotCats(ofile, id_cats);

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
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				lead_eigenvec = eigenvecs[2];			

				//eta - time sign convention
				//define relative sign for eta and time components
				//based on where the cluster is in the detector
				for(int i = 0; i < 3; i++){
					if(ec < 0){
						//time sign does NOT match eta sign
						//flip eta sign
						eigenvecs[i].SetEntry(-eigenvecs[i].at(0,0),0,0);	
					}
					//else time sign matches eta sign
				}

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
				//if(id_idx == 0){
				//	if(fabs(tc) > 15 && pi < 0.1) cout << "\n tail event here - cluster " << k << " has time " << tc << " ns and MM coeff " << pi << endl;
				//	if(phi > -1.8 && phi < -1.6 && model->GetData()->GetNPoints() > 15){
				//		cout << "\n az peak event here - cluster " << k << " has az angle: " << phi << " and theta: " << theta <<endl;
				//		for(int i = 0; i < 3; i++) cout << "\ni: " << i << " eigenval: " << eigenvals[i] << endl;
				//		cout << "eigenvec 0" << endl;
				//		eigenvecs[0].Print();
				//		cout << "eigenvec 1" << endl;
				//		eigenvecs[1].Print();
				//		cout << "eigenvec 2" << endl;
				//		lead_eigenvec.Print();

				//	}
				//}
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
				//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
				phi2D = acos( v_x / sqrt(v_x*v_x + v_y*v_y) );
				if(signbit(v_y)) phi2D *= -1;
				

				//fill hists
				//centers
				plotCats[id_idx].hists1D[0][1]->Fill(tc);
				plotCats[id_idx].hists1D[0][2]->Fill(ec);
				plotCats[id_idx].hists1D[0][3]->Fill(pc);
				//phi/eta
				plotCats[id_idx].hists1D[0][4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        			//eta/time
				plotCats[id_idx].hists1D[0][5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
				//phi/time
				plotCats[id_idx].hists1D[0][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
				//angle in 3D space
				plotCats[id_idx].hists1D[0][7]->Fill(theta);
				//angle in 2D space
				plotCats[id_idx].hists1D[0][8]->Fill(phi);
				//average cluster energy
				//w_n = E_n/N for N pts in sample
				plotCats[id_idx].hists1D[0][9]->Fill(E_k);
				//rotundity measures
				plotCats[id_idx].hists1D[0][10]->Fill(rot3D);
				plotCats[id_idx].hists1D[0][11]->Fill(rot2D);
				//velocity	
				plotCats[id_idx].hists1D[0][12]->Fill(vel);
				plotCats[id_idx].hists1D[0][13]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
				//cluster E
				plotCats[id_idx].hists1D[0][14]->Fill(E_tot);
				//get variances
				plotCats[id_idx].hists1D[0][15]->Fill(sqrt(params["cov"].at(0,0)));
				plotCats[id_idx].hists1D[0][16]->Fill(sqrt(params["cov"].at(1,1)));
				plotCats[id_idx].hists1D[0][17]->Fill(sqrt(params["cov"].at(2,2)));
				//fractional E
				plotCats[id_idx].hists1D[0][18]->Fill(E_k/E_tot);
				//"azimuth" angle in 2D (angle from x axis)
				plotCats[id_idx].hists1D[0][19]->Fill(phi2D);
	
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
				plotCats[id_idx].hists2D[0][12]->Fill(sqrt(params["cov"].at(0,0)), sqrt(params["cov"].at(1,1)));
				plotCats[id_idx].hists2D[0][13]->Fill(sqrt(params["cov"].at(2,2)), sqrt(params["cov"].at(0,0)));
				plotCats[id_idx].hists2D[0][14]->Fill(sqrt(params["cov"].at(2,2)), sqrt(params["cov"].at(1,1)));
				plotCats[id_idx].hists2D[0][15]->Fill(E_k/E_tot, pi);
				plotCats[id_idx].hists2D[0][16]->Fill(nclusters, E_k/E_tot);

				//histograms for leading/subleading clusters
				if(k == leadidx){
				cout << "LEAD eta sig: " << sqrt(params["cov"].at(0,0)) << " phi sig: " << sqrt(params["cov"].at(1,1)) << " time sig: " << sqrt(params["cov"].at(2,2)) << endl;
					//centers
					plotCats[id_idx].hists1D[1][1]->Fill(tc);
					plotCats[id_idx].hists1D[1][2]->Fill(ec);
					plotCats[id_idx].hists1D[1][3]->Fill(pc);
					//phi/eta
					plotCats[id_idx].hists1D[1][4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				//eta/time
					plotCats[id_idx].hists1D[1][5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					//phi/time
					plotCats[id_idx].hists1D[1][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//angle in 3D_lead space[1]
					plotCats[id_idx].hists1D[1][7]->Fill(theta);
					//angle in 2D_lead space[1]
					plotCats[id_idx].hists1D[1][8]->Fill(phi);
					//average cluster energy
					//w_n = E_n/N for N pts in sample
					plotCats[id_idx].hists1D[1][9]->Fill(E_k);
					//rotundity measures
					plotCats[id_idx].hists1D[1][10]->Fill(rot3D);
					plotCats[id_idx].hists1D[1][11]->Fill(rot2D);
					//velocity	
					plotCats[id_idx].hists1D[1][12]->Fill(vel);
					plotCats[id_idx].hists1D[1][13]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					//cluster E
					plotCats[id_idx].hists1D[1][14]->Fill(E_tot);
					//get variances
					plotCats[id_idx].hists1D[1][15]->Fill(sqrt(params["cov"].at(0,0)));
					plotCats[id_idx].hists1D[1][16]->Fill(sqrt(params["cov"].at(1,1)));
					plotCats[id_idx].hists1D[1][17]->Fill(sqrt(params["cov"].at(2,2)));
					//fractional E
					plotCats[id_idx].hists1D[1][18]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (angle from x axis)
					plotCats[id_idx].hists1D[1][19]->Fill(phi2D);
			
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
					plotCats[id_idx].hists2D[1][12]->Fill(params["cov"].at(0,0), params["cov"].at(1,1));
					plotCats[id_idx].hists2D[1][13]->Fill(params["cov"].at(2,2), params["cov"].at(0,0));
					plotCats[id_idx].hists2D[1][14]->Fill(params["cov"].at(2,2), params["cov"].at(1,1));
					plotCats[id_idx].hists2D[1][15]->Fill(E_k/E_tot, pi);
					plotCats[id_idx].hists2D[1][16]->Fill(nclusters, E_k/E_tot);
				}


				else if(k != leadidx){
				cout << "NOT LEAD eta sig: " << sqrt(params["cov"].at(0,0)) << " phi sig: " << sqrt(params["cov"].at(1,1)) << " time sig: " << sqrt(params["cov"].at(2,2)) << endl;
					//centers
					plotCats[id_idx].hists1D[2][1]->Fill(tc);
					plotCats[id_idx].hists1D[2][2]->Fill(ec);
					plotCats[id_idx].hists1D[2][3]->Fill(pc);
					//phi/eta
					plotCats[id_idx].hists1D[2][4]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(0,0));
        				//eta/time
					plotCats[id_idx].hists1D[2][5]->Fill(lead_eigenvec.at(0,0)/lead_eigenvec.at(2,0));
					//phi/time
					plotCats[id_idx].hists1D[2][6]->Fill(lead_eigenvec.at(1,0)/lead_eigenvec.at(2,0));
					//angle in 3D[2] space
					plotCats[id_idx].hists1D[2][7]->Fill(theta);
					//angle in 2D[2] space
					plotCats[id_idx].hists1D[2][8]->Fill(phi);
					//average cluster energy
					//w_n = E_n/N for N pts in sample
					plotCats[id_idx].hists1D[2][9]->Fill(E_k);
					//rotundity measures
					plotCats[id_idx].hists1D[2][10]->Fill(rot3D);
					plotCats[id_idx].hists1D[2][11]->Fill(rot2D);
					//velocity	
					plotCats[id_idx].hists1D[2][12]->Fill(vel);
					plotCats[id_idx].hists1D[2][13]->Fill(eigenvals_space[0]/eigenvals_space[1]);		
					
					//cluster E
					plotCats[id_idx].hists1D[2][14]->Fill(E_tot);
					
					//get variances
					plotCats[id_idx].hists1D[2][15]->Fill(params["cov"].at(0,0));
					plotCats[id_idx].hists1D[2][16]->Fill(params["cov"].at(1,1));
					plotCats[id_idx].hists1D[2][17]->Fill(params["cov"].at(2,2));
					//fractional E
					plotCats[id_idx].hists1D[2][18]->Fill(E_k/E_tot);
					//"azimuth" angle in 2D (angle from x axis)
					plotCats[id_idx].hists1D[2][19]->Fill(phi2D);
			
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
					plotCats[id_idx].hists2D[2][12]->Fill(params["cov"].at(0,0), params["cov"].at(1,1));
					plotCats[id_idx].hists2D[2][13]->Fill(params["cov"].at(2,2), params["cov"].at(0,0));
					plotCats[id_idx].hists2D[2][14]->Fill(params["cov"].at(2,2), params["cov"].at(1,1));
					plotCats[id_idx].hists2D[2][15]->Fill(E_k/E_tot, pi);
					plotCats[id_idx].hists2D[2][16]->Fill(nclusters, E_k/E_tot);

				}

			}
		}


};
#endif
