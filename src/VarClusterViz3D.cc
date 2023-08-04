#include "VarClusterViz3D.hh"
#include "VarEMCluster.hh"

#include "json/json.h"
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

void VarClusterViz3D::AddPlot(double t, string plotName){

}




void VarClusterViz3D::WriteJson(string filename){
	//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
	Json::Value root;
	Json::Value clusters;
	Json::Value cluster;
	Json::Value data;

	Json::Value x(Json::arrayValue);
	Json::Value y(Json::arrayValue);
	Json::Value z(Json::arrayValue);

	
	Json::Value eigenVec_0(Json::arrayValue);
	Json::Value eigenVec_1(Json::arrayValue);
	Json::Value eigenVec_2(Json::arrayValue);
	
	if(m_n == 0){
		return;
	}
	vector<Matrix> mus, covs;
	double pi_norm = 0;
	vector<map<string, Matrix>> cluster_params;
	vector<map<string, Matrix>> params = m_model->GetParameters();
	for(int i = 0; i < m_k; i++) cluster_params.push_back(params[i]);

	for(int i = 0; i < m_n; i++){
		//eta
		x.append(m_points->at(i).Value(0));
		//phi
		y.append(m_points->at(i).Value(1));
		z.append(m_points->at(i).Value(2));
	}
	data["x"] = x;
	data["y"] = y;
	data["z"] = z;

	root["data"] = data;
	//if no points - empty plot
//	if(x.size() == 0) return;

	for(int k = 0; k < m_k; k++)
		pi_norm += cluster_params[k]["pi"].at(0,0);
	
	//set coords for parameter circles
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;
	
	double x0, y0, z0;	
	for(int k = 0; k < m_k; k++){
		x0 = cluster_params[k]["mean"].at(0,0);
		y0 = cluster_params[k]["mean"].at(1,0);
		z0 = cluster_params[k]["mean"].at(2,0);

		cluster_params[k]["cov"].eigenCalc(eigenVals, eigenVecs);
		for(int i = 0; i < 3; i++){
			eigenVec_0.append(eigenVecs[0].at(i,0));
			eigenVec_1.append(eigenVecs[1].at(i,0));
			eigenVec_2.append(eigenVecs[2].at(i,0));
		}

	//export: data (x, y, z) in dataframe, mu (x, y, z), cov eigenvals and eigenvectors, mixing coeffs
		cluster["mixing_coeff_norm"] = cluster_params[k]["pi"].at(0,0)/pi_norm;
	
		cluster["mu_x"] = x0;
		cluster["mu_y"] = y0;
		cluster["mu_z"] = z0;

		cluster["eigenVal_0"] = eigenVals[0];	
		cluster["eigenVal_1"] = eigenVals[1];	
		cluster["eigenVal_2"] = eigenVals[2];	
		
		cluster["eigenVec_0"] = eigenVec_0;	
		cluster["eigenVec_1"] = eigenVec_1;	
		cluster["eigenVec_2"] = eigenVec_2;	
	
		clusters[std::to_string(k)] = cluster;


		eigenVec_0.clear();
		eigenVec_1.clear();
		eigenVec_2.clear();
	}
	root["clusters"] = clusters;
	Json::StreamWriterBuilder builder;
	const std::string json_file = Json::writeString(builder, root);
	std::ofstream file;
	file.open(filename+".json");
	file << json_file << endl;

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
	cout << "cvName: " << cvName << endl;
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

void VarClusterViz3D::AddAnimation(string dirname){
	/*
	//out folder does not exist
	//ie) dirname = it_0
	if(gSystem->AccessPathName((m_fname+"/"+dirname).c_str())){
		gSystem->Exec(("mkdir -p "+m_fname+"/"+dirname).c_str());
	}

	//get max + min time
	double tMax = m_points.max(2);
	double tMin = m_points.min(2);
	//nsteps = m_deltaT*(tMax - tMin)
	//loop through iterations of time
	int idx = std::to_string(m_deltaT).find("1")+1;
	double t = tMin;
	
	vector<string> tlabels;
	string label;
	while(t < tMax){
		string t_str = std::to_string(t);
		std::size_t idx1 = t_str.find(".");
		t_str.replace(idx1,1,"p");
		t_str = t_str.substr(0,idx1+3);
	//	if(std::to_string(t).find("-") != string::npos)
	//		t_str = "m"+t_str;
	//	else
	//		t_str = "p"+t_str;
		AddPlot(t, dirname+"_"+t_str);
		tlabels.push_back(t_str);
		t += m_deltaT;
	}


	//write series of plots in time
	cout << "Writing plots and gif to: ./" << m_fname << "/" << dirname << endl;
	for(int i = 0; i < m_cvs.size(); i++){
		m_cvs[i]->SaveAs((m_fname+"/"+dirname+"/cv_tEq"+tlabels[i]+".gif").c_str());
	}
	
	gSystem->cd((m_fname+"/"+dirname).c_str());
	gSystem->Exec(("rm "+dirname+".gif").c_str());
	gSystem->Exec(("convert -delay 50 -loop 1 *.gif "+dirname+".gif").c_str());
	gSystem->cd("../../../");
	m_cvs.clear();
	

	*/



}


void VarClusterViz3D::Write(){
	if(m_n == 0){
		return;
	}
	cout << "Writing plot(s) to: ./" << m_fname << ".root" << endl;
	TFile* f = TFile::Open((m_fname+".root").c_str(),"RECREATE");	
	f->cd();
	SetPalette(m_k);
	//write to root file
	for(int i = 0; i < m_cvs.size(); i++){
		m_cvs[i]->Write();
	}
	f->Close();


}




void VarClusterViz3D::SetPalette(int k){
	gStyle->cd();

	//create color palette - max number of clusters = 10
	//additional ones need to be added by hand here

	//number of gradients in the palette
	int nColors = 100;
	//set color palette
	vector<double> stops, red, green, blue;
	//number of end point colors
	int nMainColors = k;

	if(k < 2){
		cout << "Error: please give at least 2 clusters for palette creation." << endl;
		return;
	}

	//[0,1] values are R, G or B/255.
	//first color = light blue
	red.push_back(0.52);
	green.push_back(0.79);
	blue.push_back(0.96);
	//second color = light pink
	red.push_back(0.89);
	green.push_back(0.52);
	blue.push_back(0.96);

	//third color = light green
	red.push_back(0.52);
	green.push_back(0.95);
	blue.push_back(0.79);

	//fourth color = light purple
	red.push_back(0.67);
	green.push_back(0.52);
	blue.push_back(0.95);

	//fifth color = light orange
	red.push_back(0.95);
	green.push_back(0.796);
	blue.push_back(0.52);

	//sixth color = light yellow
	red.push_back(0.95);
	green.push_back(0.89);
	blue.push_back(0.52);

	
	//where to switch colors
	stops.push_back(0.0);
	stops.push_back(0.16666);
	stops.push_back(0.33333);
	stops.push_back(0.5);
	stops.push_back(0.66666);
	stops.push_back(0.83333);
	stops.push_back(1.0);

	Int_t fi = TColor::CreateGradientColorTable(nMainColors,&stops[0],&red[0],&green[0],&blue[0],nColors);
	for (int i=0;i<nColors;i++) m_palette[i] = fi+i;
	//cout << "pal - first color: " << m_palette[0] << " last color: " << m_palette[99] << endl;
}

