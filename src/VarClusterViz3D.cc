#include "VarClusterViz3D.hh"
#include "VarGaussianMixture.hh"

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
	m_deltaT = 0.1;

	//normalize data
	m_shift = m_points.Center();
	//want to apply same transformation to points: (x - shift)/scale
	m_shift.Scale(-1);
	m_scale = m_points.Normalize();
	cout << "max - min (1/scale)" << endl;
	m_scale.Print();
	m_scale.Invert();
	

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
	m_deltaT = 0.1;
	
	//normalize data
	m_shift = m_points.Center();
	//want to apply same transformation to points: (x - shift)/scale
	m_shift.Scale(-1);
	m_scale = m_points.Normalize();
	m_scale.Invert();
}

void VarClusterViz3D::AddPlot(double t, string plotName){
gErrorIgnoreLevel = kWarning;
	if(m_n == 0){
		return;
	}
	vector<Matrix> mus, covs;
	vector<double> pis;
	double pi_norm = 0;
	m_model->GetGausParameters(mus,covs);
	//normalize params

	m_model->GetMixingCoeffs(pis);
	vector<double> x, y;
	for(int i = 0; i < m_n; i++){
		//round time to nearest deltaT decimal place and plot if 0 < round(time) < t 
		if(std::ceil(m_points.at(i).Value(2) / m_deltaT) * m_deltaT > t) continue;
		//eta
		x.push_back(m_points.at(i).Value(0));
		//phi
		y.push_back(m_points.at(i).Value(1));
	}
	//if no points - empty plot
	if(x.size() == 0) return;

	for(int k = 0; k < m_k; k++){
	//	cout << "k: " << k << " pi: " << pis[k] << endl;
		pi_norm += pis[k];	
	}
	string cvName = "cv_"+plotName;
	TCanvas* cv = new TCanvas((cvName).c_str(),cvName.c_str());
	
	//get palette colors for circles
	SetPalette(m_k);
	auto cols = TColor::GetPalette();
	cv->Update();

	//sage green
	Int_t ci = TColor::GetFreeColorIndex();
	TColor* marker_color = new TColor(ci,0.61, 0.69, 0.53);  

	TGraph* gr_data = new TGraph((int)x.size(), &x[0], &y[0]);
	gr_data->SetTitle(("VarGMM EM Clustering "+plotName).c_str());
	gr_data->SetName(("VarGMM EM Clustering "+plotName).c_str());
	gr_data->SetMarkerStyle(24);
	gr_data->SetMarkerSize(0.95);
	gr_data->SetMarkerColorAlpha(ci,1);
	//can extend x/y axes so points aren't on border
	gr_data->GetXaxis()->SetTitle("eta-norm");
	gr_data->GetYaxis()->SetTitle("phi-norm");


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
	gr_data->Draw("ap");
	gr_data->GetYaxis()->SetRangeUser(-0.1,1.1);
	gr_data->GetXaxis()->SetLimits(-0.1,1.1);


//	circlePad->cd();
	//draw hist to place ellispes
//	hist->Draw("axis");
	int color_idx;
	//set coords for parameter circles
	double x0, y0, z0, c_x, c_y, r_x, r_y, theta;
	double rx0;
        double ry0;
        double rz0;
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;
	//scale first
	Matrix mat_scale = Matrix(3, 3); //3D points
	Matrix mat_shift = Matrix(3, 1); //only for mean
	
	
	mat_scale.PointToScale(m_scale);
	mat_shift.PointToShift(m_shift);
	for(int k = 0; k < m_k; k++){
		//if mean mixing coefficient value is indistinguishable from zero, don't draw
		if(pis[k]/pi_norm < 0.01)
			continue; 

		//do same normalization transformation on parameters
		//then shift -> mu'' = mu' - shift = scale*mu - shift
		//shift mus[k]
		mus[k].add(mat_shift);
		//scale mus[k] => mu' = scale*mu
		mus[k].mult(mat_scale,mus[k]);
		x0 = mus[k].at(0,0);
		y0 = mus[k].at(1,0);
		z0 = mus[k].at(2,0);


		covs[k].eigenCalc(eigenVals, eigenVecs);
	
		//scale = 1/(max - min)
		rx0 = sqrt(eigenVals[0]*m_scale.at(0));
		ry0 = sqrt(eigenVals[1]*m_scale.at(1));
		rz0 = sqrt(eigenVals[2]*m_scale.at(2));
	
		//project onto x-y plane
		r_x = rx0*sqrt(1 - ((t - z0)*(t - z0))/(rz0*rz0) );
		r_y = ry0*sqrt(1 - ((t - z0)*(t - z0))/(rz0*rz0) );
	
		
		c_x = x0;
		c_y = y0;
		//if rx or ry is nan -> ellipse doesn't exist at this value of t
		if(isnan(r_x) || isnan(r_y)) continue;
		//take direction of largest eigenvalue
		int maxValIdx = std::distance(std::begin(eigenVals),std::max_element(std::begin(eigenVals),std::end(eigenVals)));
		
		//theta = arctan(v(y)/v(x)) where v is the vector that corresponds to the largest eigenvalue
		theta = atan2(eigenVecs[maxValIdx].at(1,0),eigenVecs[maxValIdx].at(0,0));
if(k == 1 && t == 0.6){
cout << "k: " << k <<  " t: " << t << " c_x: " << c_x << " c_y: " << c_y << " r_x: " << r_x << " r_y: " << r_y << endl;
cout << "x0: " << x0 << " y0: " << y0 << " z0: " << z0 << " pi: " << pis[k] << endl;	
cout << "rx0: " << rx0 << " ry0: " << ry0 << " rz0: " << rz0 << endl;
cout << "eigenx: " << eigenVals[0] << " eigeny: " << eigenVals[1] << " eigenz: " << eigenVals[2] << endl; 
cout << "scale x: " << m_scale.at(0) << endl;
cout << "t-z0: " << (t - z0) << endl;
cout << "(t-z0)^2: " << (t - z0)*(t - z0) << endl;
cout << "(t-z0)^2/rz0^2: " << (t - z0)*(t - z0)/(rz0*rz0)  << endl;
cout << "1 - (t-z0)^2/rz0^2: " << 1 - (t - z0)*(t - z0)/(rz0*rz0)  << endl;
cout << "sqrt: " << sqrt(1 - (t - z0)*(t - z0)/(rz0*rz0) ) << endl;
covs[k].Print();
cout << "theta: " << theta << endl;
cout << "\n" << endl;
}
		//convert to degrees
		theta = 180*theta/acos(-1);	
	
		TEllipse* circle = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		TEllipse* circle_bkg = new TEllipse(c_x, c_y, r_x, r_y,0,360, theta);
		auto col = cols[int(double(k) / double(m_k - 1)*(cols.GetSize() - 1))];
		circle->SetFillColorAlpha(col,pis[k]/pi_norm);
		circle->SetLineColor(col);
		circle->SetLineWidth(5);
		//sets transparency normalized to sum
	//	cout << "k: " << k << " transparency: " << pis[k]/pi_norm << endl;
	   	circle->SetFillStyle(1001);

	///cout << "\n" << endl;

		circle_bkg->SetLineColor(0); 
		circle_bkg->SetLineWidth(8);
	   	circle_bkg->SetFillStyle(0);
		//circle_bkg->Draw();
		circle->Draw("f");
	}
	m_cvs.push_back(cv);

}




void VarClusterViz3D::SeeData(){
	vector<double> x, y, z;
	for(int i = 0; i < m_n; i++){
		//eta
		x.push_back(m_points.at(i).Value(0));
		//phi
		y.push_back(m_points.at(i).Value(1));
		//time
		z.push_back(m_points.at(i).Value(2));
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
		t_str = t_str.substr(0,idx1+2);
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
	gSystem->Exec(("convert -delay 50 -loop 1 *.gif "+dirname+".gif").c_str());
	gSystem->cd("../../../");
	m_cvs.clear();
	

	



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

