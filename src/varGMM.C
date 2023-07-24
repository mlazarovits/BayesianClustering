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
	//create symmetric matrix
	Matrix sigma = Matrix(N,N);
	sigma.InitRandomSymPosDef();
	Matrix mu = Matrix(N,1);
	mu.InitRandom(1121);
	////sample points from an n-dim gaussian for one cluster
	Matrix mat;
	mat.SampleNDimGaussian(mu,sigma,Nsample);
	PointCollection pc = mat.MatToPoints();
	
	//sample points for another cluster
	Matrix sigma2 = Matrix(N,N);
	sigma2.InitRandomSymPosDef(0.,1.,111);
//	sigma2.SetEntry(0.5,0,0);
//	sigma2.SetEntry(0.4,1,1);
//	sigma2.SetEntry(0.6,2,2);
	Matrix mu2 = Matrix(N,1);
	mu2.InitRandom(0.,1.,112);

	cout << "original sigma2" << endl;
	sigma2.Print();
	
	Matrix mat2;
	mat2.SampleNDimGaussian(mu2,sigma2,Nsample);
	pc += mat2.MatToPoints();
	cout << "min z: " << pc.min(2) << " max z: " << pc.max(2) << endl;

	cout << "original: " << pc.at(0).Value(0) << endl;
	//get normalization + centering factors
//	Point shift = pc.Center();
//	Point scale = pc.Normalize();

//	cout << "transformed: " << pc.at(0).Value(0) << endl;
//	cout << "point shift" << endl;
//	shift.Print();
//	cout << "point scale" << endl;
//	scale.Print();

	vector<Matrix> mus, covs;
	vector<double> pis;
	vector<Matrix> eigenVecs;
	vector<double> eigenVals;
	double theta;
	
	//create GMM model
	GaussianMixture gmm = GaussianMixture(k);
	gmm.SetData(&pc);
/*
	//Initialize - randomize parameters 
	vgmm.Initialize();
	VarClusterViz3D cv3D = VarClusterViz3D(&vgmm, fname);
	


	//loop
	double dLogL, newLogL, oldLogL;
	double LogLThresh = 0.0001;
	vgmm.GetGausParameters(mus,covs);
	cout << "Original parameters" << endl;
	cout << "mu 1" << endl;
	mus[0].Print();
	cout << "cov 1" << endl;
	covs[0].Print();
	////////run EM algo////////
	for(int it = 0; it < nIts; it++){
		oldLogL = vgmm.EvalLogL();
		
		//E step
		vgmm.Estimate();
		//M step
		vgmm.Update();
		
		//Plot
		cv3D.UpdatePosterior();
		if(viz) cv3D.WriteJson(fname+"/it"+std::to_string(it));
		//if(viz) cv3D.AddAnimation("it"+std::to_string(it));
		
		mus.clear(); covs.clear();
		pis.clear();

		vgmm.GetGausParameters(mus,covs);
		vgmm.GetMixingCoeffs(pis);
		cout << "iteration #" << it << ": Estimated parameters" << endl;
		for(int i = 0; i < k; i++){
			cout << "mean " << i+1 << endl;
			mus[i].Print();
			cout << "covs " << i+1 << endl;
			covs[i].Print();
			
		}

	
		//Check for convergence
		newLogL = vgmm.EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return -1;
		}
		//ELBO should not decrease with iterations, dLogL should always be negative
		dLogL = oldLogL - newLogL;
		cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(dLogL < LogLThresh){
			cout << "Reached convergence at iteration " << it << endl;
	//		break;
		}
		cout << "\n" << endl;
		
	}
//	cv3D.SeeData();
//	cv3D.Write();

	
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
cout << "\n" << endl;	


	cout << "min z: " << pc.min(2) << " max z: " << pc.max(2) << endl;


*/




}
