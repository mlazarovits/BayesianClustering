#include "GaussianMixture.hh"
#include "RandomSample.hh"
#include "ClusterViz2D.hh"
#include <iostream>

#include <TStyle.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>

using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	

/*
//cmd line i/o
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i],"--output", 8) == 0){
     			i++;
    	 		fname = string(argv[i]);
   		}
	}

*/


cout << "Free sha-va-ca-doo!" << endl;


//dimensionality
int N = 2;
//n data points
int Nsample = 500;



//create symmetric matrix
Matrix sigma = Matrix(N,N);
sigma.InitRandomSymPosDef();
Matrix mu = Matrix(N,1);
mu.InitRandom();
////sample points from an n-dim gaussian for one cluster
Matrix mat;
mat.SampleNDimGaussian(mu,sigma,Nsample);
PointCollection pc = mat.MatToPoints();
//pc.Print();


/*
int m_n = pc.GetNPoints();
vector<double> x, y;
for(int i = 0; i < m_n; i++){
	x.push_back(pc.at(i).Value(0.));
	y.push_back(pc.at(i).Value(1.));
	//just using one value of cluster weights
//	cout << "point #" << i << ": (" << x[i] << "," << y[i] << "," << z[i] << ")" << endl;
}	

TGraph* gr_data = new TGraph(m_n, &x[0], &y[0]);
gr_data->SetMarkerStyle(20);
gr_data->SetMarkerSize(0.95);
gr_data->SetTitle("GMM EM Clustering");
gr_data->GetXaxis()->SetTitle("x");
gr_data->GetYaxis()->SetTitle("y");
TCanvas* cv = new TCanvas("cv");
cv->cd();
gr_data->Draw("AP");

std::string fname = "test.root";
TFile* f = TFile::Open(fname.c_str(),"RECREATE");	
cout << "Writing plot to: " << fname << endl;
f->cd();
cv->Write();
f->Close();



cout << "x max: " << *std::max_element(std::begin(x),std::end(x)) << " x min: " << *std::min_element(std::begin(x),std::end(x)) << endl;
cout << "y max: " << *std::max_element(std::begin(y),std::end(y)) << " y min: " << *std::min_element(std::begin(y),std::end(y)) << endl;
*/

//sample points for another cluster
Matrix sigma2 = Matrix(N,N);
sigma2.InitRandomSymPosDef(0.,1.,111);
Matrix mu2 = Matrix(N,1);
mu2.InitRandom(0.,1.,112);

Matrix mat2;
mat2.SampleNDimGaussian(mu2,sigma2,Nsample);
pc += mat2.MatToPoints();

cout << "mean of cluster 1" << endl;
mu.Print();
cout << "mean of cluster 2" << endl;
mu2.Print();



//create GMM model
int k = 2; //number of clusters for GMM (may or may not be true irl)
GaussianMixture gmm = GaussianMixture(k);
gmm.AddData(pc);
//run EM algo

//Initialize - randomize parameters 
gmm.Initialize();

//////////ITERATION 1//////////
cout << "iteration #1" << endl;
//E-step: calculate posterior 
gmm.CalculatePosterior();

//M-step: get new GMM parameters
gmm.UpdateParameters();

//calculate LogL for convergence test
cout << "log-likelihood: " << gmm.EvalLogL() << endl;
//viz stuff
//1D(?), 2D and 3D
ClusterViz2D cv2D = ClusterViz2D(&gmm);
cv2D.AddPlot("it1");

/*
//////////ITERATION 2//////////
cout << "\n" << endl;
cout << "iteration #2" << endl;
//repeat
gmm.CalculatePosterior();
gmm.UpdateParameters();
cout << "log-likelihood: " << gmm.EvalLogL() << endl;

gStyle->SetPalette(kCandy);
cv2D.UpdatePosterior();
cv2D.AddPlot("it2");
*/

cv2D.Write();










}
