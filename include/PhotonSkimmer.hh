#ifndef PHOTONSKIMMER_HH
#define PHOTONSKIMMER_HH

#include "JetPoint.hh"
#include <TFile.h>
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include "PhotonProducer.hh"


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
		vector<string> id_names = {"#Chi \rightarrow #gamma", "ISR", "Not SUSY"};
		vector<TH1D*> llp_sig;
		vector<TH1D*> llp_ISR;
		vector<TH1D*> llp_notSunm;
		map<int, vector<double>> id_map;

		vector<vector<TH1D*>> histsID;

		void MakeIDHists(){
			//signal
			id_map[0] = {22, 32, 25, 35}; 
			//ISR
			id_map[1] = {20, 30, 21, 31, 23, 33, 24, 34}; 
			//not susy/unmatched
			id_map[2] = {29, -1};
			//total - combination of all of above
			id_map[3] = {};
			for(int i = 0; i < 3; i++)
				id_map[3].insert(id_map[3].end(), id_map[i].begin(), id_map[i].end());


			//each category of histograms (histsID[j]) gets a list of histograms that is the same as the total list (hists1D)
			for(int i = 0; i < (int)hists1D.size(); i++){
				for(int j = 0; (int)histsID.size(); j++){
					histsID[j].push_back((TH1D*)hists1D[i]->Clone());
				}
		
			}


		}


		void FindListHistBounds(vector<TH1D*>& hists, double& ymax, double& ymin){
			//shellsort to find max, min
			int N = (int)hists.size();
			int i, j, h;
			TH1D* v = nullptr;
			for(h = 1; h <= N/9; h = 3*h+1) ;
			for( ; h > 0; h /= 3)
				for(i = h+1; i <= N; i += 1){
					v = hists[i]; j = i;
					while(j > h && hists[j - h]->GetYaxis()->GetXmax() > v->GetYaxis()->GetXmax()){ hists[j] = hists[j - h]; j -= h; }
					hists[j] = v;
				}

			ymax = hists[0]->GetYaxis()->GetXmax();
			ymin = hists[N-1]->GetYaxis()->GetXmax();

		}


		void WriteHists(TFile* ofile){
			vector<TH1D*> hists;
			double ymax, ymin;
			string name;
			ofile->cd();
			for(int j = 0; j < (int)hists1D.size(); j++){
				//make a vector for each type of histogram
				for(int i = 0; i < (int)id_map.size() - 1; i++){
					hists.push_back(histsID[i][j]);			
					//should be 3 hists in this vector
				}
				FindListHistBounds(hists, ymax, ymin);
				name = hists1D[j]->GetName();
				TCanvas* cv = new TCanvas(name.c_str(), name.c_str());
				TDRMultiHist(hists, cv, name, name, "a.u.",ymax, ymin, id_names);
				//write cv to file			
				cv->Write();
				//write total hist to file
				histsID[j][3]->Write();

				hists.clear();
			}
			ofile->Close();

		}


		void FillHists(BasePDFMixture* model, int id_idx){
			map<string, Matrix> params;
			vector<double> eigenvals, avg_Es;
			vector<Matrix> eigenvecs;
			int nclusters = model->GetNClusters();
			histsID[id_idx][0]->Fill(nclusters);
			//e_nSubClusters->Fill(_base->Photon_energy->at(p), nclusters);
			model->GetAvgWeights(avg_Es);
			double npts = (double)model->GetData()->GetNPoints();

	
			double theta, phi, r, id;
			for(int k = 0; k < nclusters; k++){
				params = model->GetParameters(k);
				histsID[id_idx][1]->Fill(params["mean"].at(2,0));
				histsID[id_idx][2]->Fill(params["mean"].at(0,0));
				histsID[id_idx][3]->Fill(params["mean"].at(1,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//largest eigenvalue is last
				//phi/eta
				histsID[id_idx][4]->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(0,0));
        			//eta/time
				histsID[id_idx][5]->Fill(eigenvecs[2].at(0,0)/eigenvecs[2].at(2,0));
				//phi/time
				histsID[id_idx][6]->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(2,0));
				//polar angle
				//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
				r = sqrt(eigenvecs[2].at(0,0)*eigenvecs[2].at(0,0) + eigenvecs[2].at(1,0)*eigenvecs[2].at(1,0) + eigenvecs[2].at(2,0)*eigenvecs[2].at(2,0));
				theta = acos( eigenvecs[2].at(2,0) / r );
				histsID[id_idx][7]->Fill(theta);
				//azimuthal angle
				//phi = arctan(y/x)
				phi = atan(eigenvecs[2].at(1,0) / eigenvecs[2].at(0,0));
				histsID[id_idx][8]->Fill(phi);
				
				//average cluster energy
				//w_n = E_n/N for N pts in sample
				histsID[id_idx][9]->Fill(avg_Es[k]*npts);
			}
		}







};
#endif
