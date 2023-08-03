#include "GaussianMixture.hh"
#include "RandomSample.hh"
#include "VarClusterViz2D.hh"
#include "VarClusterViz3D.hh"
#include "BayesHierCluster.hh"
#include "Gaussian.hh"
#include "NormalInvWishart.hh"


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
	//for dirichlet prior in BHC
	double alpha = 1;
	unsigned long long seed = 112;
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
		if(strncmp(argv[i],"-a", 2) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--alpha", 7) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
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
   		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 1)" << endl;
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
	PointCollection pc = PointCollection();
	for(int i = 0; i < k; i++){
		//create symmetric matrix
		Matrix sigma = Matrix(N,N);
		sigma.InitRandomSymPosDef(0.,1.,seed+i);
		Matrix mu = Matrix(N,1);
		mu.InitRandom(0.,1.,seed+i);
		////sample points from an n-dim gaussian for one cluster
		Matrix mat;
		mat.SampleNDimGaussian(mu,sigma,Nsample);
		pc += mat.MatToPoints();
	}
	
	//sample points for another cluster
	//Matrix sigma2 = Matrix(N,N);
	//sigma2.InitRandomSymPosDef(0.,1.,111);
	//Matrix mu2 = Matrix(N,1);
	//mu2.InitRandom(0.,1.,112);
	//Matrix mat2;
	//mat2.SampleNDimGaussian(mu2,sigma2,Nsample);
	//pc += mat2.MatToPoints();


	
	//create Gaussian model with conjugate prior
	Gaussian* gaus = new Gaussian(N);
	//set prior parameters - default for now
	Matrix W = Matrix(N,N);
	W.InitIdentity();
	Matrix m = Matrix(N,1);
	m.InitEmpty();
	double dof = N + 0.1; //dof > D - 1
	double scale = 0.5;
	NormalInvWishart* prior = new NormalInvWishart(W, m, dof, scale);
	gaus->SetPrior(prior);
	
	//Bayesian Hierarchical Clustering algo
	BayesHierCluster* bhc = new BayesHierCluster(gaus);

	bhc->SetAlpha(alpha);
	bhc->AddData(&pc);
	bhc->Cluster();


cout << "end" << endl;	


/*	
	//create GMM model
	GaussianMixture* gmm = new GaussianMixture(k);
	gmm->SetData(&pc);
	gmm->InitParameters();
	gmm->InitPriorParameters();

	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,k);

	//viz object
	VarClusterViz3D cv3D = VarClusterViz3D(algo);
	cv3D.SeeData();

	map<string, vector<Matrix>> params;
	
	//loop
	double dLogL, newLogL, oldLogL;
	double LogLThresh = 0.0001;
	////////run EM algo////////
	for(int it = 0; it < nIts; it++){
		oldLogL = algo->EvalLogL();
		
		//E step
		algo->Estimate();
		//M step
		algo->Update();
		
		//Plot
		//cv3D.UpdatePosterior();
		//if(viz) cv3D.WriteJson(fname+"/it"+std::to_string(it));

	
		//Check for convergence
		newLogL = algo->EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return -1;
		}
		//ELBO should not decrease with iterations, dLogL should always be negative
		dLogL = oldLogL - newLogL;
		cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(fabs(dLogL) < LogLThresh){
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
	}
	cv3D.Write();

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
	cout << "\n" << endl;	

	cout << "Original parameters" << endl;
	cout << "mean 1" << endl;
	mu.Print();
	cout << "cov 1" << endl;
	sigma.Print();
	
	cout << "mean 2" << endl;
	mu2.Print();
	cout << "cov 2" << endl;
	sigma2.Print();
*/







}
