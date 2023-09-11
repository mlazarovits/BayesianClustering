#include "PhotonProducer.hh"
#include "Jet.hh"
#include "Clusterizer.hh"
#include "PhotonSkimmer.hh"
#include "TSystem.h"
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
	
	string fname = "photon";
	bool hprint = false;
	int k = 5; //number of clusters for GMM (may or may not be true irl)
	int nIts = 50; //number of iterations to run EM algorithm
	double thresh = 1.;
	double alpha = 0.1;
	bool viz = false;
	int verb = 0;
	int npho = 0;
	int evt = 0;
	bool weighted = false;
	bool smeared = false;
	bool skim = false;
	bool data = false;
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
		if(strncmp(argv[i],"-p", 2) == 0){
    	 		i++;
			npho = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--photon", 8) == 0){
    	 		i++;
			npho = std::atoi(argv[i]);
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
		if(strncmp(argv[i],"--skim", 6) == 0){
    	 		skim = true;
   		}
		if(strncmp(argv[i],"--data", 6) == 0){
    	 		data = true;
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
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --photon(-p) [npho]           set photon number to analyze (default = 0)" << endl;
   		cout << "   --event(-e) [evt]             set event number to analyze (default = 0)" << endl;
   		cout << "   --viz                         makes plots (and gifs if N == 3)" << endl;
   		cout << "   --smear                       smears data according to preset covariance (default = false)" << endl;
   		cout << "   --weight                      weights data points (default = false)" << endl;
   		cout << "   --skim                        skim over all photons to make distributions (default = false)" << endl;
   		cout << "   --data                        run over data (JetHT sample) (default = false)" << endl;
   		cout << "Example: ./photonAlgo.x -a 0.5 -t 1.6 --viz -o photonViz" << endl;

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


	fname += "_evt"+std::to_string(evt)+"_pho"+std::to_string(npho)+"_kmax"+std::to_string(k)+"_alpha"+a_string+"_thresh"+t_string;
	cout << "Free sha-va-ca-doo!" << endl;

	
	if(weighted) fname += "_Eweighted";
	if(smeared) fname += "_EtaPhiSmear";

	/////GET DATA FROM NTUPLE//////
	string in_file, cmslab;

	//local file
	//in-file = "GMSB_AOD_v6_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";//"gmsb_AODSIM_KUCMSNtuplizer_v4.root";
	

	if(data){
		in_file = "/uscms/home/mlazarov/nobackup/CMSSW_13_0_0/src/KUCMSNtupleizer/KUCMSNtupleizer/JetHT_Met50_AOD_v2_JetHT_AOD_Run2018A-15Feb2022_UL2018-v2.root";
		cmslab = "JetHT_Met50_2018_v2";
		fname += "_JetHT_v2";
	}
	
	else{
		in_file = "/uscms/home/mlazarov/nobackup/CMSSW_13_0_0/src/KUCMSNtupleizer/KUCMSNtupleizer/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X.root";
		cmslab = "GMSB_L-350TeV_Ctau-200cm_2017_v9";	
		fname += "_GMSB_v9";
	}

	if(viz){
		if(gSystem->AccessPathName((fname).c_str())){
			gSystem->Exec(("mkdir -p "+fname).c_str());
		}
		else{
			gSystem->Exec(("rm -rf "+fname).c_str());
			gSystem->Exec(("mkdir -p "+fname).c_str());

		}
		cout << "Writing to directory: " << fname << endl;
	}
	

	TFile* file = TFile::Open(in_file.c_str());
	if(skim){
		cout << "Skimming photons + subclusters" << endl;
		PhotonSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		skimmer.Skim();
		return 0;
	}


	vector<JetPoint> rhs;
	//get corresponding PV information - TODO: assuming jet is coming from interation point or PV or somewhere else?
	PhotonProducer prod(file);
	prod.GetRecHits(rhs,evt,npho);
	cout << rhs.size() << " rechits in photon " << npho << " in event " << evt << endl;
	if(rhs.size() < 1) return -1;





	//combine rechits in eta-phi area to simulate merged jet to find subjets
	Jet testpho;
	//set PV for momentum direction calculations
	Point vtx;
	prod.GetPrimaryVertex(vtx, evt);
	testpho.SetVertex(vtx);
	for(int i = 0; i < rhs.size(); i++){
		testpho.add(rhs[i]);
	}

	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi; //=sigma
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	


	Clusterizer* algo = new Clusterizer();
	algo->SetAlpha(alpha);
	algo->SetThresh(thresh);
	algo->SetVerbosity(verb);
	algo->SetMaxNClusters(k);
	algo->SetWeighted(weighted);
	if(smeared) algo->SetDataSmear(smear);
	
	if(viz)	algo->FindSubjets(testpho, fname); 
	else algo->FindSubjets(testpho);
}
