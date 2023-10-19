#include "JetProducer.hh"
#include "JetSkimmer.hh"
#include "PhotonSkimmer.hh"
#include "BayesCluster.hh"
#include "FullViz3D.hh"
#include "Clusterizer.hh"
//#include "VarClusterViz3D.hh"
#include <TSystem.h>
#include <TFile.h>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <regex>

using std::string;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	//to track computation time from beginning of program
	clock_t t = clock();

	
	string fname = "fullcluster";
	string in_file = "GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";
	bool hprint = false;
	int nIts = 50; //number of iterations to run EM algorithm
	double thresh = 1.;
	double alpha = 0.1;
	double emAlpha = 0.5;
	bool viz = false;
	int verb = 0;
	bool weighted = false;
	bool smeared = false;
	bool skim = false;
	//by default in BayesCluster
	bool distconst = true;
	//clustering strategy for skimmer
	int strat = 0; //0 is NlnN
	int evti = 0; //for skimming from evti to evtj
	int evtj = 0;
	int obj = 0; //object to cluster (0 : jets, 1 : photons)
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i],"--help", 6) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"-h", 2) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"--input", 7) == 0){
     			i++;
    	 		in_file = string(argv[i]);
   		}
		if(strncmp(argv[i],"-i", 2) == 0){
     			i++;
    	 		in_file = string(argv[i]);
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
	
		if(strncmp(argv[i],"-EMa", 3) == 0){
			i++;
    	 		emAlpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--EMalpha", 7) == 0){
			i++;
    	 		emAlpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"-a", 2) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--alpha", 7) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
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
		if(strncmp(argv[i],"--strategy", 10) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-s", 2) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--evtFirst", 6) == 0){
    	 		i++;
			evti = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--evtLast", 6) == 0){
    	 		i++;
			evtj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--object", 8) == 0){
    	 		i++;
			obj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-obj", 4) == 0){
    	 		i++;
			obj = std::atoi(argv[i]);
   		}

	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --input(-i) [file]            input root file" << endl;
   		cout << "   --output(-o) [file]           output root file (in plots/)" << endl;
   		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 0.1)" << endl;
   		cout << "   --EMalpha(-EMa) [a]           sets concentration parameter alpha for variational EM GMM (default = 0.5)" << endl;
   		cout << "   --thresh(-t) [t]              sets threshold for cluster cutoff" << endl;
		cout << "   --nIterations(-it) [nIts]     sets number of iterations for EM algorithm (default = 50)" << endl;
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --strategy(-s) [strat]        set clustering strategy for skimmer (default = NlnN, does not apply to photons)" << endl;
   		cout << "   --object(-obj) [obj]          set object to cluster (0 : jets, default; 1 : photons)" << endl;
   		cout << "   --evtFirst [i] --evtLast [j]  skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
   		cout << "   --viz                         makes plots (and gifs if N == 3)" << endl;
   		cout << "   --smear                       smears data according to preset covariance (default = false)" << endl;
   		cout << "   --weight                      weights data points (default = false)" << endl;
   		cout << "   --dist                        clusters must be within pi/2 in phi (default = false)" << endl;
   		cout << "   --skim                        skim over all photons to make distributions (default = false)" << endl;
   		cout << "Example: ./jetAlgo.x -a 0.5 -t 1.6 --viz" << endl;

   		return 0;
  	}

	if(obj == 0)
		fname = "jets";
	else if(obj == 1)
		fname = "pho";
	else{
		cout << "Object number " << obj << " not supported. Only 0 : jets, 1 : photons." << endl;
		return -1;
	}

	fname = "plots/"+fname;
	string a_string;
	std::stringstream stream;
	stream << std::fixed << std::setprecision(3) << alpha;
	a_string = stream.str();
	int idx = a_string.find(".");
	a_string.replace(idx,1,"p");	
	
	string ema_string;
	stream.str("");
	stream << std::fixed << std::setprecision(3) << emAlpha;
	ema_string = stream.str();
	idx = ema_string.find(".");
	ema_string.replace(idx,1,"p");	

	string t_string;
	stream.str("");
	stream << std::fixed << std::setprecision(3) << thresh;
	t_string = stream.str();
	idx = t_string.find(".");
	t_string.replace	(idx,1,"p");	


	string extra = "";

	//make sure if analyzing one event, not skimming + no range set
	if(evti != evtj && skim == false){
		cout << "Only going to analyze evti = " << evti << endl;
	}
	//make sure evti < evtj
	if(evti > evtj){
		int evt = evtj;
		evtj = evti;
		evti = evt;
	}

	fname += "_evt"+std::to_string(evti)+"_bhcAlpha"+a_string+"_emAlpha"+ema_string+"_thresh"+t_string+extra;
	cout << "Free sha-va-ca-doo!" << endl;
	
	if(weighted) fname += "_Eweighted";
	if(smeared) fname += "_EtaPhiSmear";
	if(distconst) fname += "_distanceConstrained";
	
	/////GET DATA FROM NTUPLE//////
	//get version
	std::smatch m;
	std::regex re("_v[0-9]_");
	string version = "";
	std::regex_search(in_file,m,re);
	for(auto x : m) version += x;
	string cmslab = in_file.substr(in_file.find(version)+4,in_file.find("_AODSIM") - in_file.find(version)-4);//"GMSB_L-350TeV_Ctau-200cm_2017_v9";
	cmslab += version.substr(0,3);
	fname += version.substr(0,3); //"_v9"
	if(viz){
		cout << "Writing to directory: " << fname << endl;
	}
	TFile* file = TFile::Open(in_file.c_str());
	if(skim){
		cout << "Skimming ";
		if(obj == 0){
			cout << "jets" << endl;
			JetSkimmer skimmer(file);
			skimmer.SetCMSLabel(cmslab);
			skimmer.SetStrategy(strat);
			//set alpha, EMalpha
			skimmer.SetEventRange(evti,evtj);
			skimmer.Skim();
		}
		else if(obj == 1){
			cout << "photons" << endl;
			PhotonSkimmer skimmer(file);
			bool data;
			if(in_file.find("SIM") != string::npos)
				data = false;
                	else
				data = true;
			skimmer.SetData(data);
                	//skimmer.SetDebug(debug);
			//set EMalpha
                	skimmer.SetCMSLabel(cmslab);
			skimmer.SetEventRange(evti,evtj);
                	skimmer.Skim();
		}
                return 0;

	}
	vector<node*> trees;
	vector<Jet> jets; 
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi; //-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	


	if(obj == 0){
		JetProducer prod(file);
		vector<Jet> rhs;
		Point vtx;
		prod.GetPrimaryVertex(vtx, evti);
		prod.GetRecHits(rhs,evti);	

		cout << rhs.size() << " rechits in event " << evti << endl;
		double gev;
		if(weighted){
			//need to transfer from GeV (energy) -> unitless (number of points)
			//transfer factor is over all points in event
			gev = 0;
			for(int i = 0; i < (int)rhs.size(); i++) gev += rhs[i].E();
			gev = gev/(double)rhs.size(); //gev = k = sum_n E_n/n pts
//			cout << "gev: " << gev << endl;
			for(int i = 0; i < (int)rhs.size(); i++){ rhs[i].SetWeight(rhs[i].E()/gev); }//weights[i] /= gev; } //sums to n pts, w_n = E_n/k 
		}

		
		//cluster jets for 1 event

//		int nrhs = 75;
//		cout << "clustering first " << nrhs << " rechits" << endl;
//		rhs.resize(nrhs);
//		rhs.shrink_to_fit();

		BayesCluster *algo = new BayesCluster(rhs);
		if(smeared) algo->SetDataSmear(smear);
		algo->SetThresh(thresh);
		algo->SetAlpha(alpha);
		algo->SetSubclusterAlpha(emAlpha);
		algo->SetVerbosity(verb);
		trees = algo->Cluster();
		
		if(viz){
			//plotting stuff here
			FullViz3D plots = FullViz3D(trees);
			plots.SetVerbosity(verb);
			plots.SetTransfFactor(gev);
			//add info of true jets
			prod.GetTrueJets(jets, evti);
			plots.AddTrueJets(jets);
			plots.Write(fname);
		}
	}


	else if(obj == 1){
		vector<JetPoint> rhs;
        	//get corresponding PV information - TODO: assuming jet is coming from interation point or PV or somewhere else?
        	PhotonProducer prod(file);
		int npho = 0; //which photon to analyze
        	prod.GetRecHits(rhs,evti,npho);
        	cout << rhs.size() << " rechits in photon " << npho << " in event " << evti << endl;
        	if(rhs.size() < 1) return -1;
        
		//combine rechits in eta-phi area to simulate merged jet to find subjets
		Jet testpho;
        	//set PV for momentum direction calculations
		Point vtx;
        	prod.GetPrimaryVertex(vtx, evti);
        	testpho.SetVertex(vtx);
        	for(int i = 0; i < rhs.size(); i++){
        	        testpho.add(rhs[i]);
        	}

		Clusterizer* algo = new Clusterizer();
	        algo->SetSubclusterAlpha(alpha);
	        algo->SetThresh(thresh);
	        algo->SetVerbosity(verb);
	        int k = 5; //max number of clusters
		algo->SetMaxNClusters(k);
	        algo->SetWeighted(weighted);
	        if(smeared) algo->SetDataSmear(smear);
	
	        if(viz) algo->FindSubjets(testpho, fname);
	        else algo->FindSubjets(testpho);

	}

	

// old implementation
//	Clusterizer* algo = new Clusterizer();
//	algo->SetClusterAlpha(alpha);
//	algo->SetSubclusterAlpha(emAlpha);
//	algo->SetThresh(thresh);
//	algo->SetVerbosity(verb);
//	if(weighted) algo->SetWeighted(weighted);
//	if(smeared) algo->SetDataSmear(smear);
//	algo->SetDistanceConstraint(distconst);
//	
//	if(viz)	algo->Cluster(testjet, fname); 
//	else algo->Cluster(testjet);
	
		
	//vector<Jet> finalJets = clusterTree.GetJets(depth = d)
	//write finalJets to same file (space, time, momentum, energy)


	t = clock() - t;
	cout << "Total program time elapsed: " << (float)t/CLOCKS_PER_SEC << " seconds." << endl;
}
