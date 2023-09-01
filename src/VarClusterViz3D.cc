#include "VarClusterViz3D.hh"
#include "VarEMCluster.hh"

#include "nlohmann/json.hpp"
#include <TSystem.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TF3.h>
#include <TColor.h>
#include <TROOT.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TColor.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TEllipse.h>
#include <string>
#include <iostream>
#include <fstream>
using std::string;
using json = nlohmann::json;
//check for 2D data
VarClusterViz3D::VarClusterViz3D(VarEMCluster* algo, string fname){ 
	if(algo->GetData()->Dim() != 3){
		cout << "VarClusterViz3D Error: dimensionality of data is not 2-> Dimensionality is " << algo->GetData()->Dim() << "." << endl;
		m_n = 0;
		return;
	}
	m_model = algo->GetModel();
	m_fname = fname;
	m_points = algo->GetData();
	m_n = m_points->GetNPoints();
	m_k = algo->GetNClusters();
	m_post = m_model->GetPosterior();
	m_deltaT = 0.1;
	
}

VarClusterViz3D::VarClusterViz3D(const VarClusterViz3D& viz){ 
	m_model = viz.m_model;
	if(m_model->GetData()->Dim() != 3){
		cout << "VarClusterViz3D Error: dimensionality of data is not 2-> Dimensionality is " << m_model->GetData()->Dim() << "." << endl;
		m_n = 0;
		return;
	}
	m_fname = viz.m_fname;
	m_points = viz.m_points;
	m_n = viz.m_n;
	m_k = viz.m_k;
	m_post = viz.m_post;
	m_deltaT = 0.1;
}





void VarClusterViz3D::WriteJson(string filename){
	//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
	json root = json::object();
	json clusters = json::object();
	json cluster = json::object();
	json data = json::object();
	
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> w;

	vector<double> eigenVec_0;
	vector<double> eigenVec_1;
	vector<double> eigenVec_2;
	
	if(m_n == 0){
		return;
	}
	vector<Matrix> mus, covs;
	vector<map<string, Matrix>> cluster_params;
	for(int i = 0; i < m_k; i++) cluster_params.push_back(m_model->GetParameters(i));

	for(int i = 0; i < m_n; i++){
		//eta
		x.push_back(m_points->at(i).Value(0));
		//phi
		y.push_back(m_points->at(i).Value(1));
		//time
		z.push_back(m_points->at(i).Value(2));
		//weight
		w.push_back(m_points->at(i).Weight());
	}
	data["x"] = json(x);
	data["y"] = json(y);
	data["z"] = json(z);
	data["w"] = json(w);


	root["data"] = data;
	//if no points - empty plot
//	if(x.size() == 0) return;

	//set coords for parameter circles
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;

	vector<double> avgs;
	m_model->GetAvgWeights(avgs);


	double x0, y0, z0;	
	for(int k = 0; k < m_k; k++){
		x0 = cluster_params[k]["mean"].at(0,0);
		y0 = cluster_params[k]["mean"].at(1,0);
		z0 = cluster_params[k]["mean"].at(2,0);


		cluster_params[k]["cov"].eigenCalc(eigenVals, eigenVecs);
		for(int i = 0; i < 3; i++){
			eigenVec_0.push_back(eigenVecs[0].at(i,0));
			eigenVec_1.push_back(eigenVecs[1].at(i,0));
			eigenVec_2.push_back(eigenVecs[2].at(i,0));
		}

	//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
		cluster["mixing_coeff_norm"] = cluster_params[k]["pi"].at(0,0);
		cluster["mu_x"] = x0;
		cluster["mu_y"] = y0;
		cluster["mu_z"] = z0;

		cluster["eigenVal_0"] = eigenVals[0];	
		cluster["eigenVal_1"] = eigenVals[1];	
		cluster["eigenVal_2"] = eigenVals[2];	
		
		cluster["eigenVec_0"] = json(eigenVec_0);	
		cluster["eigenVec_1"] = json(eigenVec_1);	
		cluster["eigenVec_2"] = json(eigenVec_2);	

		cluster["color"] = avgs[k];
	
		clusters[std::to_string(k)] = cluster;


		eigenVec_0.clear();
		eigenVec_1.clear();
		eigenVec_2.clear();
	}
	root["clusters"] = clusters;
	
	std::ofstream file;
	file.open(filename+".json");
	//4 space indent
	file << std::setw(4) << root << endl;
	if(_verb > 0) cout << "Writing to: " << filename << ".json" << endl;
}

void VarClusterViz3D::SeeData(){
	vector<double> x, y, z;
	for(int i = 0; i < m_n; i++){
		//eta
		x.push_back(m_points->at(i).Value(0));
		//phi
		y.push_back(m_points->at(i).Value(1));
		//time
		z.push_back(m_points->at(i).Value(2));
	}
	
if(m_n == 0 || x.size() == 0){
	cout << "Error: no data to plot." << endl;
	return;
}
	string cvName = "cv_data";//+plotName+"_tEq"+std::to_string(t).substr(0,3);
	TCanvas* cv = new TCanvas((cvName).c_str(),cvName.c_str());
	
	//sage green
	Int_t ci = TColor::GetFreeColorIndex();
	TColor* marker_color = new TColor(ci,0.61, 0.69, 0.53);  

	
//cout << "n pts: " << x.size() << endl;
	TGraph2D* gr_data = new TGraph2D((int)x.size(), &x[0], &y[0], &z[0]);
	gr_data->SetTitle("VarGMM EM Clustering - data only");
	gr_data->SetName("VarGMM EM Clustering - data only");
	gr_data->SetMarkerStyle(24);
	gr_data->SetMarkerSize(0.95);
	gr_data->SetMarkerColorAlpha(ci,1);
	//can extend x/y axes so points aren't on border
	gr_data->GetXaxis()->SetTitle("X");
	gr_data->GetYaxis()->SetTitle("Y");
	gr_data->GetZaxis()->SetTitle("Z");
	
	gr_data->Draw("p");

	m_cvs.push_back(cv);

}





