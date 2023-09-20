#include "JetProducer.hh"
#include "JetSkimmer.hh"
#include "Clusterizer.hh"
#include "VarClusterViz3D.hh"
#include <TSystem.h>
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
	
	string fname = "jet";
	bool hprint = false;
	int nIts = 50; //number of iterations to run EM algorithm
	double thresh = 1.;
	double alpha = 0.5;
	bool viz = false;
	int verb = 0;
	int evt = 0;
	bool weighted = false;
	bool smeared = false;
	bool skim = false;
	bool distconst = false;
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
		if(strncmp(argv[i],"--weight", 8) == 0){
    	 		weighted = true;
   		}
		if(strncmp(argv[i],"--smear", 7) == 0){
    	 		smeared = true;
   		}
		if(strncmp(argv[i],"--dist", 6) == 0){
    	 		distconst = true;
   		}
		if(strncmp(argv[i],"--skim", 6) == 0){
    	 		skim = true;
   		}
	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --ouput(-o) [file]            output root file (in plots/)" << endl;
   		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 1)" << endl;
   		cout << "   --thresh(-t) [t]              sets threshold for cluster cutoff" << endl;
		cout << "   --nIterations(-it) [nIts]     sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --event(-e) [evt]             set event number to analyze (default = 0)" << endl;
   		cout << "   --viz                         makes plots (and gifs if N == 3)" << endl;
   		cout << "   --smear                       smears data according to preset covariance (default = false)" << endl;
   		cout << "   --weight                      weights data points (default = false)" << endl;
   		cout << "   --dist                        clusters must be within pi/2 in phi (default = false)" << endl;
   		cout << "   --skim                        skim over all photons to make distributions (default = false)" << endl;
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


	string extra = "";

	fname += "_evt"+std::to_string(evt)+"_alpha"+a_string+"_thresh"+t_string+extra;
	cout << "Free sha-va-ca-doo!" << endl;
	
	if(weighted) fname += "_Eweighted";
	if(smeared) fname += "_EtaPhiSmear";
	if(distconst) fname += "_distanceConstrained";
	
	/////GET DATA FROM NTUPLE//////
	string in_file = "GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";
	string cmslab = "GMSB_L-350TeV_Ctau-200cm_2017_v9";
	fname += "_v9";
	if(viz){
	//directory is made in python script
	//	if(gSystem->AccessPathName((fname).c_str())){
	//		gSystem->Exec(("mkdir -p "+fname).c_str());
	//	}
	//	else{
	//		gSystem->Exec(("rm -rf "+fname).c_str());
	//		gSystem->Exec(("mkdir -p "+fname).c_str());

	//	}
		cout << "Writing to directory: " << fname << endl;
	}
	TFile* file = TFile::Open(in_file.c_str());
	if(skim){
		cout << "Skimming jets" << endl;
		JetSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		skimmer.Skim();
		return 0;
	}
	JetProducer prod(file);


	vector<JetPoint> rhs;
	//get corresponding PV information - TODO: assuming jet is coming from interation point or PV or somewhere else?
	Point vtx;
	prod.GetPrimaryVertex(vtx, evt);
	prod.GetRecHits(rhs,evt);

	//combine rechits in eta-phi area to simulate merged jet to find subjets
	Jet testjet;
	//set PV for momentum direction calculations
	testjet.SetVertex(vtx);
	for(int i = 0; i < rhs.size(); i++){
		testjet.add(rhs[i]);
	}

	cout << testjet.GetNConstituents() << " constituents in event " << evt << endl;

	/*
	map<string, Matrix> params;
	params["scale"] = Matrix(0.001);
	params["dof"] = Matrix(3 - 1 + 10.);
	Matrix mat = Matrix(3, 3);
	mat.InitIdentity();
	Matrix mean = Matrix(3, 1);
	mat.InitEmpty();
	params["scalemat"] = mat;
	params["m"] = mean;
	*/

	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi; //-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	

	
	cout << "Clustering with alpha = " << alpha << " and cutoff threshold = " << thresh << endl;
	//cluster jets for 1 event
	Clusterizer* algo = new Clusterizer();
	algo->SetAlpha(alpha);
	algo->SetThresh(thresh);
	algo->SetVerbosity(verb);
	if(weighted) algo->SetWeighted(weighted);
	if(smeared) algo->SetDataSmear(smear);
	algo->SetDistanceConstraint(distconst);
	
	if(viz)	algo->Cluster(testjet, fname); 
	else algo->Cluster(testjet);
	
		
	//vector<Jet> finalJets = clusterTree.GetJets(depth = d)
	//write finalJets to same file (space, time, momentum, energy)


}
