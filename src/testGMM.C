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
RandomSample rs(123);
rs.SetRange(0.,1.);
//2D Gauss
vector<double> mu = {0.1, 0.4};
Matrix cov = Matrix(2,2);
//with diagonal covariance matrix
cov.SetEntry(0.5,0,0);
cov.SetEntry(0.6,1,1);
//n data points
int n = 100;

//test Cholesky decomp
//create symmetric matrix
int N = 3;
Matrix mat = Matrix(N,N);
mat.InitRandomSym();
Matrix L = Matrix(N,N);
L.cholesky(mat);

//eventually make a PointCollection
//vector<vector<double>> pts = rs.SampleNDimGaussian(





//create GMM model

//run EM algo


//get new GMM parameters




















}
