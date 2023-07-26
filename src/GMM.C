#include "EMCluster.hh"
#include "RandomSample.hh"
#include "GaussianMixture.hh"
#include "ClusterViz2D.hh"
#include <iostream>


#include <TStyle.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <string>
using std::cout;
using std::endl;
using std::string;

int main(int argc, char *argv[]){
	
	string fname;// = "test";
	bool hprint = false;
	//dimensionality
	int N = 2;
	//n data points
	int Nsample = 500;
	int k = 2; //number of clusters for GMM (may or may not be true irl)
	int nIts = 1;
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
		if(strncmp(argv[i],"--nIterations", 13) == 0){
     			i++;
    	 		nIts = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-it", 5) == 0){
			i++;
    	 		nIts = std::atoi(argv[i]);
   		}
	
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --ouput(-o) [file]            output root file (in test/)" << endl;
   		cout << "   --nSamples(-n) [n]            sets number of data points to simulate per cluster (default = 500)" << endl;
   		cout << "   --nDims(-d) [d]               sets dimensionality of data (default = 2)" << endl;
   		cout << "   --nClusters(-k) [k]           sets number of clusters in GMM (default = 2)" << endl;
   		cout << "   --nIterations(-it) [nIts]     sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "Example: ./runGMM_EM.x -n 100 -o testViz.root" << endl;

   		return 0;
  	}

	fname = "plots/"+fname;
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
	GaussianMixture* gmm = new GaussianMixture(k);
	gmm->SetData(&pc);
	gmm->InitParameters();

	//create EM algo
	EMCluster* algo = new EMCluster(gmm,k);

	//viz object
	ClusterViz2D cv2D = ClusterViz2D(algo);
	cv2D.SeeData();
	
	map<string, vector<Matrix>> params;
	
	//loop
	double dLogL, newLogL, oldLogL;
	double LogLThresh = 0.0005;
	////////run EM algo////////
	for(int it = 0; it < nIts; it++){
		oldLogL = algo->EvalLogL();
		////run EM in for loop
		algo->Estimate();
		algo->Update();
		
		//Plot
		cv2D.UpdatePosterior();
		cv2D.AddPlot("it"+std::to_string(it));
		
		//Check for convergence
		newLogL = algo->EvalLogL();
		dLogL = fabs(oldLogL - newLogL);
		cout << "iteration #" << it << " dLogL: " << dLogL << endl;
		if(dLogL < LogLThresh){
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
	}
	cv2D.Write();
	
	params = gmm->GetParameters();	
	vector<Matrix> mus = params["means"];
	vector<Matrix> covs = params["covs"];

	cout << "Estimated parameters" << endl;
	cout << "mean 1" << endl;
	mus[0].Print();
	cout << "cov 1" << endl;
	covs[0].Print();
	cout << "mean 2" << endl;
	mus[1].Print();
	cout << "cov 2" << endl;
	covs[1].Print();

	cout << "Original parameters" << endl;
	cout << "mean 1" << endl;
	mu.Print();
	cout << "cov 1" << endl;
	sigma.Print();
	cout << "mean 2" << endl;
	mu2.Print();
	cout << "cov 2" << endl;
	sigma2.Print();
}
