#include "VarClusterViz2D.hh"
#include "GaussianMixture.hh"

#include <TSystem.h>
#include <TGraph.h>
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

using std::string;

//check for 2D data
VarClusterViz2D::VarClusterViz2D(VarEMCluster* algo, string fname){ 
	if(algo->GetData()->Dim() != 2){
		cout << "VarClusterViz2D Error: dimensionality of data is not 2. Dimensionality is " << algo->GetData()->Dim() << "." << endl;
		m_n = 0;
		return;
	}
	m_algo = algo;
	m_model = algo->GetModel();
	m_fname = fname;
	m_points = algo->GetData();
	m_n = m_points->GetNPoints();
	m_k = algo->GetNClusters();
	m_post = m_model->GetPosterior();
}

VarClusterViz2D::VarClusterViz2D(const VarClusterViz2D& viz){ 
	m_model = viz.m_model;
	if(m_model->GetData()->Dim() != 2){
		cout << "VarClusterViz2D Error: dimensionality of data is not 2-> Dimensionality is " << m_model->GetData()->Dim() << "." << endl;
		m_n = 0;
		return;
	}
	m_fname = viz.m_fname;
	m_points = viz.m_points;
	m_n = viz.m_n; 
	m_k = viz.m_k;
	m_post = viz.m_post;
}
void VarClusterViz2D::AddPlot(string plotName){
	if(m_n == 0){
		return;
	}
	string cvName = "cv_"+plotName;
	TCanvas* cv = new TCanvas((cvName).c_str(),cvName.c_str());
	double pi_norm = 0;
	vector<map<string, Matrix>> clusters;
	for(int i = 0; i < m_k; i++) clusters.push_back(m_model->GetLikelihoodParameters(i));

	vector<double> x, y;	
	for(int i = 0; i < m_n; i++){
		x.push_back(m_points->at(i).Value(0.));
		y.push_back(m_points->at(i).Value(1.));
	}
	for(int k = 0; k < m_k; k++){
		pi_norm += clusters[k]["pi"].at(0,0);
	}
	//get palette colors for circles
	SetPalette(m_k);
	auto cols = TColor::GetPalette();
	cv->Update();

	//sage green
	Int_t ci = TColor::GetFreeColorIndex();
	TColor* marker_color = new TColor(ci,0.61, 0.69, 0.53);  
	
	TGraph* gr_data = new TGraph(m_n, &x[0], &y[0]);
	gr_data->SetTitle(("VarGMM EM Clustering "+plotName).c_str());
	gr_data->SetName(("VarGMM EM Clustering "+plotName).c_str());
	gr_data->SetMarkerStyle(24);
	gr_data->SetMarkerSize(0.95);
	gr_data->SetMarkerColorAlpha(ci,1);
	//can extend x/y axes so points aren't on border
	gr_data->GetXaxis()->SetTitle("x");
	gr_data->GetYaxis()->SetTitle("y");

/*
	TH1F* hist = gr_data->GetHistogram();
	//to turn off hist axes
	hist->GetXaxis()->SetLabelSize(0.);
	hist->GetXaxis()->SetTickLength(0.);
	hist->GetYaxis()->SetLabelSize(0.);
	hist->GetYaxis()->SetTickLength(0.);
	hist->SetFillStyle(4000);
	
*/
	//one pad
	TPad* graphPad = new TPad("graph pad", "graph pad",0.0,0.0,1.0,1.0);
	graphPad->Draw();
/*
	//another pad - draw on top
	TPad* circlePad = new TPad("circle pad", "circle pad",0.0,0.0,1.0,1.0);
	circlePad->SetFillStyle(4000);
	circlePad->SetFrameFillStyle(4000);
	circlePad->Draw();

*/
	//draw data	
	graphPad->cd();
	gr_data->Draw("ap"); //PCOLZ draws color palette
//	circlePad->cd();
	//draw hist to place ellispes
//	hist->Draw("axis");


	int color_idx;
	//set coords for parameter circles
	double c_x, c_y, r_x, r_y, theta; //centers, radii, and angle of each ellipse
	for(int k = 0; k < m_k; k++){
		//if mean mixing coefficient value is indistinguishable from zero, don't draw
		if(clusters[k]["pi"].at(0,0)/pi_norm < 0.01)
			continue; 

		c_x = clusters[k]["mean"].at(0,0);
		c_y = clusters[k]["mean"].at(1,0);
		
		//calculate eigenvalues + vectors for orientation (angle) of ellipse
		vector<Matrix> eigenVecs;
		vector<double> eigenVals;
		clusters[k]["cov"].eigenCalc(eigenVals, eigenVecs);
	
		r_x = sqrt(eigenVals[0]);
		r_y = sqrt(eigenVals[1]);
	
		//take direction of largest eigenvalue
		int maxValIdx = std::distance(std::begin(eigenVals),std::max_element(std::begin(eigenVals),std::end(eigenVals)));
		
		//theta = arctan(v(y)/v(x)) where v is the vector that corresponds to the largest eigenvalue
		//probably should use atan2
		theta = atan2(eigenVecs[maxValIdx].at(1,0),eigenVecs[maxValIdx].at(0,0));
	
		//convert to degrees
		theta = 180*theta/acos(-1);
	
			
		if(theta < 0) theta = -(90+theta);
		else theta += 90; //get in ROOT's weird angle definition

		TEllipse* circle = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		TEllipse* circle_bkg = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		auto col = cols[int(double(k) / double(m_k - 1)*(cols.GetSize() - 1))];
		circle->SetFillColorAlpha(col,clusters[k]["pi"].at(0,0)/pi_norm);
		circle->SetLineColor(col);
		circle->SetLineWidth(5);
		//sets transparency normalized to sum
	   	circle->SetFillStyle(1001);

		circle_bkg->SetLineColor(0); 
		circle_bkg->SetLineWidth(8);
	   	circle_bkg->SetFillStyle(0);
		circle_bkg->Draw();
		circle->Draw("f");
	}

	m_cvs.push_back(cv);

}







void VarClusterViz2D::Write(){
	if(m_n == 0){
		return;
	}
	cout << "Writing plot(s) to: " << m_fname << ".root" << endl;
	TFile* f = TFile::Open((m_fname+".root").c_str(),"RECREATE");	
	f->cd();
	SetPalette(m_k);
cout << "# cvs: " << m_cvs.size() << endl;
	//out folder does not exist
	if(gSystem->AccessPathName(m_fname.c_str())){
		gSystem->Exec(("mkdir "+m_fname).c_str());
	}
	
	for(int i = 0; i < m_cvs.size(); i++){
		m_cvs[i]->Write();
		if(i < 10)	
		m_cvs[i]->SaveAs((m_fname+"/cv_0"+std::to_string(i)+".pdf").c_str());
		else
		m_cvs[i]->SaveAs((m_fname+"/cv_"+std::to_string(i)+".pdf").c_str());
	}

	f->Close();

}




void VarClusterViz2D::SetPalette(int k){
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

