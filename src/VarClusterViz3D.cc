#include "VarClusterViz3D.hh"
#include "VarGaussianMixture.hh"

#include <TSystem.h>
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

using std::string;

//check for 2D data
VarClusterViz3D::VarClusterViz3D(VarGaussianMixture* model, string fname){ 
	if(model->GetData().Dim() != 3){
		cout << "VarClusterViz3D Error: dimensionality of data is not 2. Dimensionality is " << model->GetData().Dim() << "." << endl;
		m_n = 0;
		return;
	}
	m_model = model;
	m_fname = fname;
	m_points = m_model->GetData();
	m_n = m_points.GetNPoints();
	m_k = m_model->GetNClusters(0.);
	m_post = m_model->GetPosterior();
}

VarClusterViz3D::VarClusterViz3D(const VarClusterViz3D& viz){ 
	m_model = viz.m_model;
	if(m_model->GetData().Dim() != 3){
		cout << "VarClusterViz3D Error: dimensionality of data is not 2. Dimensionality is " << m_model->GetData().Dim() << "." << endl;
		m_n = 0;
		return;
	}
	m_fname = viz.m_fname;
	m_points = m_model->GetData();
	m_n = m_points.GetNPoints();
	m_k = m_model->GetNClusters(0.);
	m_post = m_model->GetPosterior();
}
void VarClusterViz3D::AddPlot(string plotName){
	if(m_n == 0){
		return;
	}
	string cvName = "cv_"+plotName;
	TCanvas* cv = new TCanvas((cvName).c_str(),cvName.c_str());
	vector<Matrix> mus, covs;
	vector<double> pis;
	double pi_norm = 0;
	m_model->GetGausParameters(mus,covs);
	m_model->GetMixingCoeffs(pis);
	vector<double> x, y, z;
	for(int i = 0; i < m_n; i++){
		x.push_back(m_points.at(i).Value(0.));
		y.push_back(m_points.at(i).Value(1.));
		z.push_back(m_points.at(i).Value(2.));
	}
	for(int k = 0; k < m_k; k++){
		cout << "k: " << k << " pi: " << pis[k] << endl;
		pi_norm += pis[k];
	}
	//get palette colors for circles
	SetPalette(m_k);
	auto cols = TColor::GetPalette();
	cv->Update();

	//sage green
	Int_t ci = TColor::GetFreeColorIndex();
	TColor* marker_color = new TColor(ci,0.61, 0.69, 0.53);  
	
	TGraph2D* gr_data = new TGraph2D(m_n, &x[0], &y[0], &z[0]);
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
	double c_x, c_y, c_z, r_x, r_y, r_z, alpha, beta, gamma;
	string eqn, func_x, func_y, func_z;
	Matrix A = Matrix(3,3);
	Matrix R = Matrix(3,3);
	Matrix R_T = Matrix(3,3);
	Matrix Rx = Matrix(3,3);
	Matrix Ry = Matrix(3,3);
	Matrix Rz = Matrix(3,3);
	for(int k = 0; k < m_k; k++){
		//if mean mixing coefficient value is indistinguishable from zero, don't draw
		if(pis[k]/pi_norm < 0.01)
			continue; 

		c_x = mus[k].at(0,0);
		c_y = mus[k].at(1,0);
		c_z = mus[k].at(2,0);

		func_x = "(x - "+std::to_string(c_x)+")";
		func_y = "(y - "+std::to_string(c_y)+")";
		func_z = "(z - "+std::to_string(c_z)+")";
		
		//calculate eigenvalues + vectors for orientation (angle) of ellipse
		vector<Matrix> eigenVecs;
		vector<double> eigenVals;
		covs[k].eigenCalc(eigenVals, eigenVecs);
	
		r_x = sqrt(eigenVals[0]);
		r_y = sqrt(eigenVals[1]);
		r_z = sqrt(eigenVals[2]);
		
		//alpha = arccos(-z2/sqrt(1 - z3^2))
		alpha = acos(eigenVecs[2].at(1,0)/sqrt(1 - eigenVecs[2].at(2,0)));
		//beta = arccos(z3)
		beta = acos(eigenVecs[2].at(2,0));
		//gamma = arccos(y3/sqrt(1 - z3^2))
		gamma = acos(eigenVecs[1].at(2,0)/sqrt(1 - eigenVecs[2].at(2,0)));


		//rotate radius matrix A => R_T*A*R
		A.reset();
		R.reset();
		Rx.reset();
		Ry.reset();
		Rz.reset();

		//construct A matrix - diagonal with entries being inverse squares of axes
		A.SetEntry(1./(r_x*r_x), 0, 0);
		A.SetEntry(1./(r_y*r_y), 1, 1);
		A.SetEntry(1./(r_z*r_z), 2, 2);
		
		//construct rotation matrix from euler angles - R = Rz(gam)*Ry(beta)*Rx(alpha)
		Rz.SetEntry(cos(gamma), 0, 0);
		Rz.SetEntry(-sin(gamma), 0, 1);
		Rz.SetEntry(sin(gamma), 1, 0);
		Rz.SetEntry(cos(gamma), 1, 1);
		Rz.SetEntry(1, 2, 2);

		Ry.SetEntry(cos(beta), 0, 0);
		Ry.SetEntry(-sin(beta), 2, 0);
		Ry.SetEntry(sin(beta), 0, 2);
		Ry.SetEntry(cos(beta), 2, 2);
		Ry.SetEntry(1, 1, 1);

		Rx.SetEntry(cos(alpha), 1, 1);
		Rx.SetEntry(-sin(alpha), 1, 2);
		Rx.SetEntry(sin(alpha), 2, 1);
		Rx.SetEntry(cos(alpha), 2, 2);
		Rx.SetEntry(1, 0, 0);
	
	
		R.mult(Rz,Ry);
		R.mult(R,Rx);

		A.mult(R_T,A);
		A.mult(A,R);	

		eqn  = func_x+"*( "+std::to_string(A.at(0,0))+"*"+func_x+" + "+std::to_string(A.at(0,1))+"*"+func_y+" + "+std::to_string(A.at(0,2))+"*"+func_z+" ) + "; 
		eqn += func_y+"*( "+std::to_string(A.at(1,0))+"*"+func_x+" + "+std::to_string(A.at(1,1))+"*"+func_y+" + "+std::to_string(A.at(1,2))+"*"+func_z+" ) + ";
		eqn += func_z+"*( "+std::to_string(A.at(2,0))+"*"+func_x+" + "+std::to_string(A.at(2,1))+"*"+func_y+" + "+std::to_string(A.at(2,2))+"*"+func_z+" ) = 1";

		TF3* circle = new TF3(("cluster_"+std::to_string(k)).c_str(),eqn.c_str());
		//TEllipse* circle = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		//TEllipse* circle_bkg = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		auto col = cols[int(double(k) / double(m_k - 1)*(cols.GetSize() - 1))];
		circle->SetFillColorAlpha(col,pis[k]/pi_norm);
		circle->SetLineColor(col);
		circle->SetLineWidth(5);
		//sets transparency normalized to sum
		cout << "k: " << k << " transparency: " << pis[k]/pi_norm << endl;
	   	circle->SetFillStyle(1001);

		//circle_bkg->SetLineColor(0); 
		//circle_bkg->SetLineWidth(8);
	   	//circle_bkg->SetFillStyle(0);
		//circle_bkg->Draw();
		circle->Draw("f");
	}

	m_cvs.push_back(cv);

}







void VarClusterViz3D::Write(){
	if(m_n == 0){
		return;
	}
	cout << "Writing plot(s) to: " << m_fname << ".root" << endl;
	TFile* f = TFile::Open((m_fname+".root").c_str(),"RECREATE");	
	f->cd();
	SetPalette(m_k);

	//out folder does not exist
	if(gSystem->AccessPathName(m_fname.c_str())){
		gSystem->Exec(("mkdir "+m_fname).c_str());
	}
	
	for(int i = 0; i < m_cvs.size(); i++){
		m_cvs[i]->Write();
		m_cvs[i]->SaveAs((m_fname+"/cv_"+std::to_string(i)+".pdf").c_str());
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

