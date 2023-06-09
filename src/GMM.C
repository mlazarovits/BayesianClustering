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
	
	string fname;// = "test";
	int nSamples = 100;
	bool hprint = false;
	//dimensionality
	int N = 2;
	//n data points
	int Nsample = 500;
	int k = 2; //number of clusters for GMM (may or may not be true irl)
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i],"--help", 6) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"-h", 2) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"--output", 8) == 0){
     			i++;
    	 		fname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-o", 2) == 0){
     			i++;
    	 		fname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-n", 2) == 0){
     			i++;
    	 		Nsample = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nSamples", 10) == 0){
     			i++;
    	 		Nsample = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-n", 2) == 0){
     			i++;
    	 		nSamples = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nSamples", 10) == 0){
     			i++;
    	 		nSamples = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-d", 2) == 0){
     			i++;
    	 		N = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nDims", 6) == 0){
     			i++;
    	 		N = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-k", 2) == 0){
     			i++;
    	 		k = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nClusters", 11) == 0){
     			i++;
    	 		k = std::atoi(argv[i]);
   		}
	
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)            print options" << endl;
   		cout << "   --ouput(-o) [file]    output root file (in test/)" << endl;
   		cout << "   --nSamples(-n) [n]    sets number of data points to simulate per cluster (default = 500)" << endl;
   		cout << "   --nDims(-d) [d]       sets dimensionality of data (default = 2)" << endl;
   		cout << "   --nClusters(-k) [k]   sets number of clusters in GMM (default = 2)" << endl;
   		cout << "Example: ./runGMM_EM.x -n 100 -o testViz.root" << endl;

   		return 0;
  	}

	fname = "plots/"+fname
	cout << "Free sha-va-ca-doo!" << endl;
	
	
	
	
	/////SIMULATE DATA//////
	//create symmetric matrix
	Matrix sigma = Matrix(N,N);
	sigma.InitRandomSymPosDef();
	Matrix mu = Matrix(N,1);
	mu.InitRandom();
	////sample points from an n-dim gaussian for one cluster
	Matrix mat;
	mat.SampleNDimGaussian(mu,sigma,Nsample);
	PointCollection pc = mat.MatToPoints();
	
	//sample points for another cluster
	Matrix sigma2 = Matrix(N,N);
	sigma2.InitRandomSymPosDef(0.,1.,111);
	Matrix mu2 = Matrix(N,1);
	mu2.InitRandom(0.,1.,112);
	
	Matrix mat2;
	mat2.SampleNDimGaussian(mu2,sigma2,Nsample);
	pc += mat2.MatToPoints();
	
	cout << "Original parameters" << endl;
	cout << "mean 1" << endl;
	mu.Print();
	cout << "cov 1" << endl;
	sigma.Print();
	cout << "mean 2" << endl;
	mu2.Print();
	cout << "cov 1" << endl;
	sigma2.Print();
	
	
	
	//create GMM model
	GaussianMixture gmm = GaussianMixture(k);
	gmm.AddData(pc);
	ClusterViz2D cv2D = ClusterViz2D(&gmm);
	////////run EM algo////////
	
	//Initialize - randomize parameters 
	gmm.Initialize();
	
	//loop
	int nIts = 50;
	double dLogL, newLogL, oldLogL;
	double LogLThresh = 0.0005;
	double it = 0;
	for(int it = 0; it < nIts; it++){
		oldLogL = gmm.EvalLogL();
		
		//E step
		gmm.CalculatePosterior();
		//M step
		gmm.UpdateParameters();
		
		//Plot
		cv2D.UpdatePosterior();
		cv2D.AddPlot("it"+std::to_string(it));
		cv2D.Write();
		
		//Check for convergence
		cout << "iteration #" << it << " log-likelihood: " << gmm.EvalLogL() << endl;
		newLogL = gmm.EvalLogL();
		dLogL = fabs(oldLogL - newLogL);
		cout << "dLogL: " << dLogL << endl;
		if(dLogL < LogLThresh){
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
		
	}
	vector<Matrix> mus, covs;
	gmm.GetParameters(mus,covs);
	
	cout << "Estimated parameters" << endl;
	cout << "mean 1" << endl;
	mus[0].Print();
	cout << "covs 1" << endl;
	covs[0].Print();
	
	
	cout << "mean 2" << endl;
	mus[1].Print();
	cout << "covs" << endl;
	covs[1].Print();


}
