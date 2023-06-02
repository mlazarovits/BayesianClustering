#include "GaussianMixture.hh"
#include "RandomSample.hh"
#include <iostream>

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



//sample gaussians to simulate data
int N = 3;
////n data points
int Nsample = 5;

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

//sample points for another cluster
Matrix sigma2 = Matrix(N,N);
sigma2.InitRandomSymPosDef(0.,1.,111);
Matrix mu2 = Matrix(N,1);
mu2.InitRandom(0.,1.,111);

Matrix mat2;
mat2.SampleNDimGaussian(mu2,sigma2,Nsample);
pc += mat2.MatToPoints();
//pc.Print();

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

cout << "iteration #1" << endl;
//E-step: calculate posterior 
gmm.CalculatePosterior();

//M-step: get new GMM parameters
gmm.UpdateParameters();

//calculate LogL for convergence test
cout << "log-likelihood: " << gmm.EvalLogL() << endl;

cout << "\n" << endl;
cout << "iteration #2" << endl;
//repeat
gmm.CalculatePosterior();
gmm.UpdateParameters();
cout << "log-likelihood: " << gmm.EvalLogL() << endl;


//need to write viz stuff
//2D and 3D



}
