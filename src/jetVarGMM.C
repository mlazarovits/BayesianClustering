#include "JetProducer.hh"
#include "JetClusterizer.hh"


#include <TFile.h>
#include <iostream>
#include <cmath>
#include<string>

using std::string;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	
	string fname = "test";
	bool hprint = false;
	//dimensionality
	int N = 2;
	//n data points
	int Nsample = 50;
	int k = 2; //number of clusters for GMM (may or may not be true irl)
	int nIts = 50; //number of iterations to run EM algorithm
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
	
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --ouput(-o) [file]            output root file (in test/)" << endl;
   		cout << "   --nSamples(-n) [n]            sets number of data points to simulate per cluster (default = 500)" << endl;
   		cout << "   --nDims(-d) [d]               sets dimensionality of data (default = 2)" << endl;
   		cout << "   --nClusters(-k) [k]           sets number of clusters in GMM (default = 2)" << endl;
   		cout << "   --nIterations(-it) [nIts]   sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "Example: ./runGMM_EM.x -n 100 -o testViz.root" << endl;

   		return 0;
  	}

	fname = "plots/"+fname;
	cout << "Free sha-va-ca-doo!" << endl;
	
	
	/////GET DATA FROM NTUPLE//////
	string in_file = "gmsb_AODSIM_KUCMSNtuplizer_v4.root";
	TFile* file = TFile::Open(in_file.c_str());
	JetProducer prod(file);
	vector<JetPoint> rhs;
	int evt = 0;
	//get corresponding PV information - TODO: assuming jet is coming from interation point or PV or somewhere else?
	Point vtx;
	prod.GetPrimaryVertex(vtx, evt);
	prod.GetRecHits(rhs,0);
	cout << rhs.size() << " rechits in first event" << endl;

	//combine rechits in eta-phi area to simulate merged jet to find subjets
	Jet testjet;
	//set PV for momentum direction calculations
	testjet.SetVertex(vtx);
	double etaMax = 0.5;
	double etaMin = -etaMax;  
	double phiMax = 2.;
	double phiMin = -2.8;
	int nRhs = 0;
	for(int i = 0; i < rhs.size(); i++){
//		if(nRhs > 10) break;
		if(rhs[i].eta() > etaMax || rhs[i].eta() < etaMin)
			continue;
		if(rhs[i].phi() > phiMax || rhs[i].phi() < phiMin)
			continue;
			testjet.add(rhs[i]);
		nRhs++;
	}

	cout << testjet.GetNConstituents() << " constituents in testjet" << endl;

	//cluster jets for 1 event
	JetClusterizer jc;
	//calculate subjets for all rechits in a eta-phi area - pretend they have been merged into a jet
	double logLthresh = 0.001;
	int maxIt = 1;
	int maxK = 3;
	bool viz = true;
	jc.FindSubjets_etaPhi(testjet, logLthresh, maxIt, maxK, viz);


		
	//vector<Jet> finalJets = clusterTree.GetJets(depth = d)
	//write finalJets to same file (space, time, momentum, energy)




/*	
	//put jets into algorithm friendly format - point collection
	PointCollection points;
	VarGaussianMixture vgmm(m_maxK);
	//match points to jets by idx
	int n_pts = jet.GetConstituents();
	for(int i = 0; i < n_pts; i++){
		points += jets[i].four_space();
		jets[i].SetUserIdx(i);
	}

	

	
	//create GMM model
	VarGaussianMixture vgmm = VarGaussianMixture(k);
	vgmm.AddData(pc);
	VarClusterViz2D cv2D = VarClusterViz2D(&vgmm, fname);

	//Initialize - randomize parameters 
	vgmm.Initialize();
	//loop
	double dLogL, newLogL, oldLogL;
	double LogLThresh = 0.0001;
	double it = 0;
	////////run EM algo////////
	for(int it = 0; it < nIts; it++){
		oldLogL = vgmm.EvalLogL();
		
		//E step
		vgmm.CalculatePosterior();
		//M step
		vgmm.UpdateParameters();
		
		//Plot
		cv2D.UpdatePosterior();
		cv2D.AddPlot("it"+std::to_string(it));
		
		//Check for convergence
		newLogL = vgmm.EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it << " log-likelihood: " << newLogL << endl;
			return -1;
		}
		dLogL = fabs(oldLogL - newLogL);
		cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(dLogL < LogLThresh){
			cout << "Reached convergence at iteration " << it << endl;
			break;
		}
		
	}
	cv2D.Write();
	vector<Matrix> mus, covs;
	vector<double> pis;
	vgmm.GetGausParameters(mus,covs);
	vgmm.GetMixingCoeffs(pis);
	cout << "Estimated parameters" << endl;
	for(int i = 0; i < k; i++){
		cout << "mean " << i+1 << endl;
		mus[i].Print();
		cout << "covs " << i+1 << endl;
		covs[i].Print();
		cout << "mixing coeff " << i+1 << " " << pis[i] << endl;
	}

	cout << "Original parameters" << endl;
	cout << "mean 1" << endl;
	mu.Print();
	cout << "cov 1" << endl;
	sigma.Print();
	cout << "mean 2" << endl;
	mu2.Print();
	cout << "cov 1" << endl;
	sigma2.Print();
cout << "\n" << endl;	
*/
}
