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

//cout << "mean of cluster 1" << endl;
//mu.Print();
//cout << "mean of cluster 2" << endl;
//mu2.Print();

//create GMM model
GaussianMixture gmm = GaussianMixture(2);
gmm.AddData(pc);
//run EM algo
Initialize - randomize parameters 
gmm.Initialize();
//E-step 
gmm.CalculatePosterior();

//get new GMM parameters




















}
