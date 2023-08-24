#include "JetProducer.hh"
#include "JetClusterizer.hh"
#include "VarClusterViz3D.hh"

#include <TFile.h>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	
	string fname = "testjet";
	bool hprint = false;
	int k = 2; //number of clusters for GMM (may or may not be true irl)
	int nIts = 50; //number of iterations to run EM algorithm
	double thresh = 1.;
	double alpha = 0.1;
	bool viz = false;
	int verb = 0;
	int evt = 0;
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
		if(strncmp(argv[i],"-e", 2) == 0){
    	 		i++;
			evt = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--event", 7) == 0){
    	 		i++;
			evt = std::atoi(argv[i]);
   		}
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --ouput(-o) [file]            output root file (in test/)" << endl;
   		cout << "   --nClusters(-k) [k]           sets number of clusters in GMM (default = 2)" << endl;
   		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 1)" << endl;
   		cout << "   --thresh(-t) [t]              sets threshold for cluster cutoff" << endl;
		cout << "   --nIterations(-it) [nIts]     sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "   --viz                         makes plots (and gifs if N == 3)" << endl;
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --event(-e) [evt]             set event number to analyze (default = 0)" << endl;
   		cout << "Example: ./jetAlgo.x -a 0.5 -t 1.6 --viz" << endl;

   		return 0;
  	}

	fname = "plots/"+fname;
	string a_string;
	std::stringstream stream;
	stream << std::fixed << std::setprecision(3) << alpha;
	a_string = stream.str();
	int idx = a_string.find(".");
	a_string.replace(idx,1,"p");	

	string t_string;
	stream.str("");
	stream << std::fixed << std::setprecision(3) << thresh;
	t_string = stream.str();
	idx = t_string.find(".");
	t_string.replace(idx,1,"p");	


	fname += "_evt"+std::to_string(evt)+"_kmax"+std::to_string(k)+"_alpha"+a_string+"_thresh"+t_string;
	cout << "Free sha-va-ca-doo!" << endl;
	
	
	/////GET DATA FROM NTUPLE//////
	string in_file = "gmsb_AODSIM_KUCMSNtuplizer_v4.root";
	TFile* file = TFile::Open(in_file.c_str());
	JetProducer prod(file);
	vector<JetPoint> rhs;
	//get corresponding PV information - TODO: assuming jet is coming from interation point or PV or somewhere else?
	Point vtx;
	prod.GetPrimaryVertex(vtx, evt);
	prod.GetRecHits(rhs,evt);
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

	cout << "Clustering with alpha = " << alpha << " and cutoff threshold = " << thresh << endl;
	//cluster jets for 1 event
	JetClusterizer jc;
	//calculate subjets for all rechits in a eta-phi area - pretend they have been merged into a jet
//	jc.FindSubjets_etaPhi(testjet, thresh, nIts, k, viz, alpha);
	jc.Cluster(testjet, alpha, thresh, viz, verb, fname);

		
	//vector<Jet> finalJets = clusterTree.GetJets(depth = d)
	//write finalJets to same file (space, time, momentum, energy)


}
