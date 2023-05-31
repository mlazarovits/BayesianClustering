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
RandomSample rs;
rs.SetRange(0.,1.);
//n data points
int Nsample = 10;

//create symmetric matrix
Matrix sigma = Matrix(N,N);
sigma.InitRandomSymPosDef();
Matrix mu = Matrix(N,1);
mu.InitRandom();


Matrix mat;
mat.SampleNDimGaussian(mu,sigma,Nsample);
////eventually make a PointCollection


cout << "end" << endl;



//create GMM model

//run EM algo


//get new GMM parameters




















}
