#include "ClusterViz2D.hh"
#include "GaussianMixture.hh"
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
ClusterViz2D::ClusterViz2D(GaussianMixture* model) : 
	ClusterVizBase{ model }{
	if(m_model->GetData().Dim() != 2){
		cout << "ClusterViz2D Error: dimensionality of data is not 2. Dimensionality is " << m_model->GetData().Dim() << "." << endl;
		return;
	}
	m_post = m_model->GetPosterior();
}


void ClusterViz2D::AddPlot(string plotName){
	string cvName = "cv_"+plotName;
	m_cv = new TCanvas((cvName).c_str(),cvName.c_str());
	vector<double> x, y, z;
	for(int i = 0; i < m_n; i++){
		x.push_back(m_points.at(i).Value(0.));
		y.push_back(m_points.at(i).Value(1.));
		//just using one value of cluster weights - only works for 2 cluster models
		z.push_back(m_post.at(i,0.));	
		if(i == 10)cout << "point #" << i << ": (" << x[i] << "," << y[i] << "," << z[i] << ")" << endl;
	}


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
	gr_data->Draw("PCOLZ"); //PCOLZ draws color palette
	graphPad->SetTheta(90);
	graphPad->SetPhi(-360);	
	circlePad->cd();
	//draw hist to place ellispes
	hist->Draw("axis");

	//set coords for parameter circles
	double c_x, c_y, r_x, r_y, theta; //centers, radii, and angle of each ellipse
	vector<Matrix> mus, covs;
	m_model->GetParameters(mus,covs);
	for(int k = 0; k < m_model->GetNClusters(); k++){
		c_x = mus[k].at(0,0);
		c_y = mus[k].at(1,0);
		r_x = covs[k].at(0,0);
		r_y = covs[k].at(1,1);
cout << "k: " << k << endl;	
		vector<Matrix> eigenVecs;
		vector<double> eigenVals;
		covs[k].Print();	
		covs[k].eigenCalc(eigenVals, eigenVecs);

		TEllipse* circle = new TEllipse(c_x, c_y, r_x, r_y);//0.,0,360, theta);
		TEllipse* circle_bkg = new TEllipse(c_x, c_y, r_x, r_y);//0.,0,360, theta);
		circle->SetLineColor(gStyle->GetColorPalette(k*19+10)); 
		circle->SetLineWidth(5);
	   	circle->SetFillStyle(0);
		circle_bkg->SetLineColor(0); 
		circle_bkg->SetLineWidth(8);
	   	circle_bkg->SetFillStyle(0);
		circle_bkg->Draw();
		circle->Draw();
	}




}







void ClusterViz2D::Write(string fname){
	TFile* f = TFile::Open((fname+".root").c_str(),"RECREATE");	
	cout << "Writing plot to: " << fname << ".root" << endl;
	f->cd();

		m_cv->Write();

	f->Close();






}


