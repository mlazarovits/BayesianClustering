#include "PhotonProducer.hh"
#include "VarClusterViz3D.hh"
#include "Jet.hh"
#include "GaussianMixture.hh"
#include "VarEMCluster.hh"

#include "TSystem.h"
#include <TFile.h>
#include <iostream>
#include <cmath>
#include<string>

using std::string;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	
	string fname = "testphoton";
	bool hprint = false;
	//dimensionality
	int N = 2;
	//n data points
	int Nsample = 50;
	int k = 2; //number of clusters for GMM (may or may not be true irl)
	int nIts = 50; //number of iterations to run EM algorithm
	double thresh = 1.;
	double alpha = 0.1;
	bool viz = false;
	int verb = 0;
	int npho = 0;
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
		if(strncmp(argv[i],"--viz", 5) == 0){
    	 		viz = true;
   		}
		if(strncmp(argv[i],"--verbosity", 11) == 0){
    	 		i++;
			verb = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-v", 2) == 0){
    	 		i++;
			verb = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-t", 2) == 0){
			i++;
    	 		thresh = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--thresh", 8) == 0){
			i++;
    	 		thresh = std::stod(argv[i]);
   		}
	
		if(strncmp(argv[i],"-a", 2) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--alpha", 7) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"-p", 2) == 0){
    	 		i++;
			npho = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--photon", 8) == 0){
    	 		i++;
			npho = std::atoi(argv[i]);
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
   		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 1)" << endl;
   		cout << "   --thresh(-t) [t]              sets threshold for cluster cutoff" << endl;
		cout << "   --nIterations(-it) [nIts]     sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "   --viz                         makes plots (and gifs if N == 3)" << endl;
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --photon(-p) [npho   ]        set photon number to analyze (default = 0)" << endl;
   		cout << "Example: ./runGMM_EM.x -n 100 -o testViz.root" << endl;

   		return 0;
  	}

	fname = "plots/"+fname;
	cout << "Free sha-va-ca-doo!" << endl;
	if(gSystem->AccessPathName((fname).c_str())){
		gSystem->Exec(("mkdir -p "+fname).c_str());
	}
	
	
	
	/////GET DATA FROM NTUPLE//////
	string in_file = "gmsb_AODSIM_KUCMSNtuplizer_v4.root";
	TFile* file = TFile::Open(in_file.c_str());
	PhotonProducer prod(file);
	vector<JetPoint> rhs;
	int evt = 0;
	//get corresponding PV information - TODO: assuming jet is coming from interation point or PV or somewhere else?
	Point vtx;
	prod.GetPrimaryVertex(vtx, evt);
	prod.GetRecHits(rhs,evt,npho);
	cout << rhs.size() << " rechits in first photon in first event" << endl;
	if(rhs.size() < 1) return -1;

	fname += "evt"+std::to_string(evt)+"_pho"+std::to_string(npho)+"_kmax"+std::to_string(k);

	//combine rechits in eta-phi area to simulate merged jet to find subjets
	Jet testpho;
	//set PV for momentum direction calculations
	testpho.SetVertex(vtx);
	for(int i = 0; i < rhs.size(); i++){
		testpho.add(rhs[i]);
	}

	cout << testpho.GetNConstituents() << " constituents in test photon" << endl;


	PointCollection* pc = new PointCollection();
	testpho.GetEtaPhiConstituents(*pc);
	
	//create GMM model
	GaussianMixture* gmm = new GaussianMixture(k);
	gmm->SetData(pc);
	gmm->InitParameters();
	gmm->InitPriorParameters();
	
	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,k);
	algo->SetThresh(thresh);

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

	cout << "Estimated parameters" << endl;
	map<string, Matrix> params;
	for(int i = 0; i < gmm->GetNClusters(); i++){
		params = gmm->GetPriorParameters(i);	
		cout << "weight " << i << ": " << params["pi"].at(0,0) << " alpha: " << params["alpha"].at(0,0) << endl;
		cout << "mean " << i << endl;
		params["mean"].Print();
		cout << "cov " << i << endl;
		params["cov"].Print();
		cout << "scale " << i << ": " << params["scale"].at(0,0) << " dof " << i << ": " << params["dof"].at(0,0) << endl;
		cout << "m " << i << endl;
		params["m"].Print();
		cout << "scalemat " << i << endl;
		params["scalemat"].Print();
		params.clear();
	}

}
