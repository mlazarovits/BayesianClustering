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
int Nsample = 10;
//RandomSample rs;
//rs.SetRange(0.,1.);

//create symmetric matrix
cout << "init sigma" << endl;
Matrix sigma = Matrix(N,N);
cout << "fill sigma" << endl;
sigma.InitRandomSymPosDef();
cout << "init mu" << endl;
Matrix mu = Matrix(N,1);
cout << "fill mu" << endl;
//mu.InitRandom();


Matrix mat;
mat.SampleNDimGaussian(mu,sigma,Nsample);

//crases with SampleNDimGaussian + Matrix sigma = Matrix(N,N); + sigma.InitRandomSymPosDef();



////eventually make a PointCollection



cout << "end" << endl;



//create GMM model

//run EM algo


//get new GMM parameters




















}
