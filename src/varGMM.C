#include "GaussianMixture.hh"
#include "RandomSample.hh"
#include "VarClusterViz2D.hh"
#include "VarClusterViz3D.hh"
#include <iostream>
#include <cmath>
#include <TSystem.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>

using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	
	string fname = "test";
	bool hprint = false;
	//dimensionality
	int N = 3;
	//n data points
	int Nsample = 50;
	int k = 2; //number of clusters for GMM (may or may not be true irl)
	int nIts = 50; //number of iterations to run EM algorithm
	bool viz = false;
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
		if(strncmp(argv[i],"-n", 3) == 0){
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
		if(strncmp(argv[i],"-v", 2) == 0){
    	 		viz = true;
   		}
		if(strncmp(argv[i],"--viz", 5) == 0){
    	 		viz = true;
   		}
	
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --ouput(-o) [file]            output root file (in plots/)" << endl;
   		cout << "   --nSamples(-n) [n]            sets number of data points to simulate per cluster (default = 500)" << endl;
   		cout << "   --nDims(-d) [d]               sets dimensionality of data (default = 2)" << endl;
   		cout << "   --nClusters(-k) [k]           sets number of clusters in GMM (default = 2)" << endl;
   		cout << "   --nIterations(-it) [nIts]     sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "   --viz(-v)                     makes plots (and gifs if N == 3)" << endl;
   		cout << "Example: ./runGMM_EM.x -n 100 -o testViz.root" << endl;

   		return 0;
  	}

	fname = "plots/"+fname;
	if(gSystem->AccessPathName((fname).c_str())){
		gSystem->Exec(("mkdir -p "+fname).c_str());
	}
	cout << "Free sha-va-ca-doo!" << endl;
	
	
	
	
	/////SIMULATE DATA//////
	PointCollection pc;
	cout << "Original parameters" << endl;
	for(int i = 0; i < 2; i++){
		//create symmetric matrix
		Matrix sigma = Matrix(N,N);
		sigma.InitRandomSymPosDef(0., 1., 111+i);
		Matrix mu = Matrix(N,1);
		mu.InitRandom(0., 10., 1121+i);
		cout << "mean " << i << endl;
		mu.Print();
		cout << "cov " << i << endl;
		sigma.Print();	
		cout << "\n" << endl;	
		////sample points from an n-dim gaussian for one cluster
		Matrix mat;
		mat.SampleNDimGaussian(mu,sigma,Nsample);
		pc += mat.MatToPoints();
	}	
	
	//create GMM model
	GaussianMixture* gmm = new GaussianMixture(k);
	gmm->SetData(&pc);
	gmm->InitParameters();
	gmm->InitPriorParameters();
	
	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,k);
	algo->SetThresh(1.);

	//viz object
	VarClusterViz3D cv3D = VarClusterViz3D(algo);
	cv3D.SeeData();

	
	//loop
	double dLogL, newLogL;
	double LogLThresh = 0.001;
	double oldLogL = algo->EvalLogL();
	////////run EM algo////////
	for(int it = 0; it < nIts; it++){
		//cout << "------------- it #" << it << " BEGIN -------------" << endl;
		//E step
		algo->Estimate();
		//M step
		algo->Update();
		
		//Plot
		cv3D.UpdatePosterior();
		if(viz) cv3D.WriteJson(fname+"/it"+std::to_string(it));

	
		//Check for convergence
		newLogL = algo->EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return -1;
		}
		//ELBO should not decrease with iterations, dLogL should always be negative
		dLogL = oldLogL - newLogL;
		cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << " old ELBO: " << oldLogL << " new ELBO: " << newLogL << endl;
		oldLogL = newLogL;
		if(fabs(dLogL) < LogLThresh){
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
	}
	cv3D.Write();

	cout << "Estimated parameters" << endl;
	map<string, Matrix> params;
	for(int i = 0; i < gmm->GetNClusters(); i++){
		params = gmm->GetPriorParameters(i);	
		cout << "weight " << i << ": " << params["pi"].at(0,0) << endl;
		cout << "mean " << i << endl;
		params["mean"].Print();
		cout << "cov " << i << endl;
		params["cov"].Print();
		params.clear();
	}







}
