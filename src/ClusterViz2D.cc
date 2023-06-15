#include "ClusterViz2D.hh"
#include "BasePDFMixture.hh"

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
ClusterViz2D::ClusterViz2D(BasePDFMixture* model, string fname) : 
	ClusterVizBase{ model, fname }{
	if(m_model->GetData().Dim() != 2){
		cout << "ClusterViz2D Error: dimensionality of data is not 2. Dimensionality is " << m_model->GetData().Dim() << "." << endl;
		return;
	}
	m_post = m_model->GetPosterior();
	m_fname = fname;
}
/*
ClusterViz2D::ClusterViz2D(VarGaussianMixture* model, string fname) : 
	ClusterVizBase{ model, fname }{
	if(m_model->GetData().Dim() != 2){
		cout << "ClusterViz2D Error: dimensionality of data is not 2. Dimensionality is " << m_model->GetData().Dim() << "." << endl;
		return;
	}
	m_post = m_model->GetPosterior();
	m_fname = fname;
	cout << "Writing plot to: " << m_fname << ".root" << endl;
}
*/
void ClusterViz2D::AddPlot(string plotName){
	string cvName = "cv_"+plotName;
	TCanvas* cv = new TCanvas((cvName).c_str(),cvName.c_str());
	vector<Matrix> mus, covs;
	m_model->GetParameters(mus,covs);
	vector<double> x, y, z;
	for(int i = 0; i < m_n; i++){
		x.push_back(m_points.at(i).Value(0.));
		y.push_back(m_points.at(i).Value(1.));
		//just using one value of cluster weights - only works for 2 cluster models
		//weight for column 1 used because
		//w_i ~ 1 = class 1
		//w_i ~ 0 = class 0 to match class to idx
		z.push_back(m_post.at(i,1.));	
	}

	//get palette colors for circles
	SetPalette(m_model->GetNClusters());
	auto cols = TColor::GetPalette();
	
	cv->Update();

	TGraph2D* gr_data = new TGraph2D(m_n, &x[0], &y[0], &z[0]);
	gr_data->SetTitle(("GMM EM Clustering "+plotName).c_str());
	gr_data->SetName(("GMM EM Clustering "+plotName).c_str());
	gr_data->SetMarkerStyle(20);
	gr_data->SetMarkerSize(0.95);
	//can extend x/y axes so points aren't on border
	gr_data->GetXaxis()->SetTitle("x");
	gr_data->GetYaxis()->SetTitle("y");
	gr_data->GetZaxis()->SetRange(0.,1.);


	TH2D* hist = gr_data->GetHistogram();
	//to turn off hist axes
	hist->GetXaxis()->SetLabelSize(0.);
	hist->GetXaxis()->SetTickLength(0.);
	hist->GetYaxis()->SetLabelSize(0.);
	hist->GetYaxis()->SetTickLength(0.);
	hist->SetFillStyle(4000);
	

	//one pad
	TPad* graphPad = new TPad("graph pad", "graph pad",0.0,0.0,1.0,1.0);
	graphPad->Draw();

	//another pad - draw on top
	TPad* circlePad = new TPad("circle pad", "circle pad",0.0,0.0,1.0,1.0);
	circlePad->SetFillStyle(4000);
	circlePad->SetFrameFillStyle(4000);
	circlePad->Draw();


	//draw data	
	graphPad->cd();
	gr_data->GetZaxis()->SetRangeUser(0.,1.);
	gr_data->Draw("PCOLZ"); //PCOLZ draws color palette
	gr_data->GetZaxis()->SetRangeUser(0.,1.);
	graphPad->SetTheta(90);
	graphPad->SetPhi(-360);	
	circlePad->cd();
	//draw hist to place ellispes
	hist->Draw("axis");
	int color_idx;
	//set coords for parameter circles
	double c_x, c_y, r_x, r_y, theta; //centers, radii, and angle of each ellipse
	for(int k = 0; k < m_model->GetNClusters(); k++){
		c_x = mus[k].at(0,0);
		c_y = mus[k].at(1,0);
		
		//calculate eigenvalues + vectors for orientation (angle) of ellipse
		vector<Matrix> eigenVecs;
		vector<double> eigenVals;
		covs[k].eigenCalc(eigenVals, eigenVecs);
	
		r_x = sqrt(eigenVals[0]);
		r_y = sqrt(eigenVals[1]);
	
		//take direction of largest eigenvalue
		int maxValIdx = std::distance(std::begin(eigenVals),std::max_element(std::begin(eigenVals),std::end(eigenVals)));
		
		//theta = arctan(v(y)/v(x)) where v is the vector that corresponds to the largest eigenvalue
		theta = atan(eigenVecs[maxValIdx].at(1,0)/eigenVecs[maxValIdx].at(0,0));
	
	//	cout << "eigen" << endl;
	//	for(int i = 0; i < eigenVals.size(); i++){
	//		cout << "val: " << eigenVals[i] << endl;
	//		cout << "vec:" << endl;
	//		eigenVecs[i].Print();


	//	}
	//
	//	cout << "k: " << k << " theta: " << theta << " (rad)" << endl;
		//convert to degrees
		theta = 180*theta/acos(-1);
	
			
		if(theta < 0) theta = -(90+theta);
		else theta += 90; //get in ROOT's weird angle definition

		TEllipse* circle = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		TEllipse* circle_bkg = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		circle->SetLineColor(cols[int(double(k) / double(m_model->GetNClusters() - 1)*(cols.GetSize() - 1))]);
	
		circle->SetLineWidth(5);
	   	circle->SetFillStyle(0);
		circle_bkg->SetLineColor(0); 
		circle_bkg->SetLineWidth(8);
	   	circle_bkg->SetFillStyle(0);
		circle_bkg->Draw();
		circle->Draw();
	}


	m_cvs.push_back(cv);

}







void ClusterViz2D::Write(){
	cout << "Writing plot(s) to: " << m_fname << ".root" << endl;
	TFile* f = TFile::Open((m_fname+".root").c_str(),"RECREATE");	
	f->cd();
	SetPalette(m_model->GetNClusters());
	for(int i = 0; i < m_cvs.size(); i++){
		m_cvs[i]->Write();
	}

	f->Close();

}




void ClusterViz2D::SetPalette(int k){
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

	//2 clusters
	if(k == 2){
		//[0,1] values are R, G or B/255.
		//first color = light blue
		red.push_back(0.52);
		green.push_back(0.79);
		blue.push_back(0.96);
		//second color = light pink
		red.push_back(0.89);
		green.push_back(0.52);
		blue.push_back(0.96);
		

		//where to switch colors
		stops.push_back(0.0);
		stops.push_back(0.5);
		stops.push_back(1.0);

	}

	//3 clusters
	if(k == 3){
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
		
		//where to switch colors
		stops.push_back(0.0);
		stops.push_back(0.3333);
		stops.push_back(0.6666);
		stops.push_back(1.0);

	}

	Int_t fi = TColor::CreateGradientColorTable(nMainColors,&stops[0],&red[0],&green[0],&blue[0],nColors);
	for (int i=0;i<nColors;i++) m_palette[i] = fi+i;
	//cout << "pal - first color: " << m_palette[0] << " last color: " << m_palette[99] << endl;
}

