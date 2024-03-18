#include "JetProducer.hh"
#include "JetSkimmer.hh"
#include "PhotonSkimmer.hh"
#include "BayesCluster.hh"
#include "FullViz3D.hh"
//#include "Clusterizer.hh"
#include "VarClusterViz3D.hh"
#include "BasicDetectorSim.hh"

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
	string fname = "";
	string oname = "";
	string in_file = "rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";
	bool hprint = false;
	double thresh = 1.;
	double alpha = 0.1;
	double emAlpha = 0.5;
	bool viz = false;
	int verb = 0;
	bool weighted = true;
	bool smeared = true;
	bool timesmeared = false;
	//by default in BayesCluster
	bool distconst = true;
	//clustering strategy for skimmer
	int strat = 0; //0 is NlnN
	int evt = 0; //for event evt
	int obj = 0; //object to cluster (0 : jets, 1 : photons)
	int nobj = 0; //number of object in event to cluster (only for photons)
	//this should be in N/GeV
	//at least at 1 GeV but 1 GeV rh shouldn’t be able to seed a cluster so 1 GeV should be a fraction of entries
	double gev = 1./30.;
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
    	 		oname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-o", 2) == 0){
     			i++;
    	 		oname = string(argv[i]);
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
		if(strncmp(argv[i],"--noWeight", 10) == 0){
    	 		weighted = false;
   		}
		if(strncmp(argv[i],"--noSmear", 9) == 0){
    	 		smeared = false;
   		}
		if(strncmp(argv[i],"--timeSmear", 11) == 0){
    	 		timesmeared = true;
   		}
		if(strncmp(argv[i],"--noDist", 8) == 0){
    	 		distconst = false;
   		}
		if(strncmp(argv[i],"--strategy", 10) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-s", 2) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--evt", 5) == 0){
    	 		i++;
			evt = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-e", 2) == 0){
    	 		i++;
			evt = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--object", 8) == 0){
    	 		i++;
			obj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nobj", 6) == 0){
    	 		i++;
			nobj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--gev", 5) == 0){
			i++;
    	 		gev = std::stod(argv[i]);
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
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --strategy(-s) [strat]        set clustering strategy for skimmer (0 : NlnN (default), 1 : N^2,  does not apply to photons)" << endl;
   		cout << "   --object [obj]                set object to cluster (0 : jets, default; 1 : photons, 2 : detector sim)" << endl;
   		cout << "   --nobj [n]                    set number of object in event to cluster (photons only, default = 0)" << endl;
   		cout << "   --evt(-e) [evt]               get plots for event e (default = 0)" << endl;
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/30 GeV)" << endl;
   		cout << "   --viz                         makes plots (and gifs if N == 3)" << endl;
   		cout << "   --noSmear                     turns off smearing data according to preset covariance (default = true)" << endl;
   		cout << "   --timeSmear                   turns on time smearing data according to preset covariance (default = false)" << endl;
   		cout << "   --noWeight                    turns off weighting data points (default = true)" << endl;
   		cout << "   --noDist                      turns off - clusters must be within pi/2 in phi (default = true)" << endl;
   		cout << "Example: ./FullClusterSingle.x -a 0.5 -t 1.6 --viz" << endl;

   		return 0;
  	}
	cout << "Free sha-va-ca-doo!" << endl;
	string cmslab, version;	
	/////GET DATA FROM NTUPLE//////
	if(!in_file.empty()){
		//get version
		std::smatch m;
		std::regex re("_v[0-9]+_");
		version = "";
		std::regex_search(in_file,m,re);
		for(auto x : m) version += x;
		cmslab = in_file.substr(in_file.find(version)+version.size(),in_file.find(".root") - in_file.find(version)-version.size());//"GMSB_L-350TeV_Ctau-200cm_2017_v9";
		version.pop_back();
		cmslab += version;
	}
	/////MAKE DATA WITH PYTHIA + BASIC DETECTOR SIM//////
	else{
		in_file = "BDSIM";
	} 

	if(obj == 0 || obj == 2){
		if(obj == 0) fname += "jets";
		else fname += "jetsSim"; 
		if(strat == 0)
			fname += "NlnN";
		else if(strat == 1)
			fname += "N2";
		else if(strat == 2)
			fname += "GMMonly";
		else{
			cout << "Strategy number " << strat << " not supported. Only 0 : NlnN, 1 : N^2, 2 : GMM only." << endl;
			return -1;
		}
	}
	else if(obj == 1)
		fname += "photon";
	else{
		cout << "Object number " << obj << " not supported. Only 0 : jets, 1 : photons." << endl;
		return -1;
	}




	if(oname != "")
		fname = "skims/"+fname+"Skim_"+oname;
	else fname = "skims/"+fname+"Skim";
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
	t_string.replace(idx,1,"p");	
	
	string gev_string;
	stream.str("");
	stream << std::fixed << std::setprecision(3) << gev;
	gev_string = stream.str();
	idx = gev_string.find(".");
	gev_string.replace(idx,1,"p");	

	
	if(obj != 1) fname += "_bhcAlpha"+a_string+"_emAlpha"+ema_string+"_thresh"+t_string+"_";
	else fname += "_emAlpha"+ema_string+"_thresh"+t_string+"_";
	fname += "NperGeV"+gev_string+"_";
	fname += cmslab; //long sample name


	fname = "plots/"+fname;




	if(viz){
		cout << "Writing to directory: " << fname << endl;
		if(gSystem->AccessPathName(fname.c_str())){
			gSystem->Exec(("mkdir -p "+fname).c_str());
		}
	}
	TFile* file = TFile::Open(in_file.c_str());
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = 2*acos(-1)/360.; //1 degree in radians
	double deta = dphi; //-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(0.,2,2); //no smear in time	
	double tres_c = 0.2; //ns
	//double tres_n = 30*sqrt(tres_c - tres_c*tres_c); //v1//ns*E (set s.t. 30 GeV gives sig_t = 1 ns)
	double tres_n = 30*sqrt(1 - tres_c*tres_c); //v1//ns*E (set s.t. 30 GeV gives sig_t = 1 ns)

	vector<Jet> rhs, jets, phos;
	vector<node*> trees;
	//get rhs (as Jets) for event
	if(obj == 0){
		JetProducer prod(file);
		prod.SetTransferFactor(gev);
		prod.SetMinPt(30);
		prod.SetMinNrhs(15);
		prod.SetMinEmE(20);
		if(strat != 2){
			cout << "Getting rec hits for jets at event " << evt << endl;
			prod.GetRecHits(rhs,evt);	
			cout << rhs.size() << " rechits in event " << evt << endl;
		}
		else{
			cout << "Getting rec hits for jet " << nobj << " at event " << evt << endl;
			prod.GetTrueJets(jets, evt);
			if(jets.size() < 1){ cout << "No jets passing selection found for event " << evt << endl; return -1; }
			if(nobj > jets.size() - 1){ cout << "Only " << jets.size() << " jets passing selection found for event " << evt << endl; return -1; }
			jets[nobj].GetJets(rhs);
			cout << rhs.size() << " rechits in jet " << nobj << " in event " << evt << endl;
		}
	}

	else if(obj == 1){
        	PhotonProducer prod(file);
		prod.SetTransferFactor(gev);
		int npho = nobj; //which photon to analyze
		if(viz) cout << "Making plots for photon " << npho << " at event " << evt << endl;
        	//prod.GetRecHits(rhs,evt,npho);
        	prod.SetIsoCut();
		prod.GetTruePhotons(phos, evt);
		if(phos.size() < 1){ cout << "No photons passing selection found for event " << evt << endl; return -1; }
		if(nobj > phos.size() - 1){ cout << "Only " << phos.size() << " photons passing selection found for event " << evt << endl; return -1; }
		phos[npho].GetJets(rhs);
		cout << rhs.size() << " rechits in photon " << npho << " in event " << evt << endl;
	}
	else if(obj == 2){
		BasicDetectorSim det;
		det.SetNEvents(1);
		det.SetVerbosity(verb);
		det.SetEnergyThreshold(1.); //set to 1 GeV
		//could set time resolution cte to MTD specs here
		
		det.SimTTbar();
		det.TurnOnPileup();
		det.TurnOnSpikes(0.05);
		//default arg is all events
		det.SimulateEvents(evt);
		det.GetRecHits(rhs);
		det.GetTrueJets(jets);
		cout << rhs.size() << " rechits in event " << evt << endl;
	}
	
        if(rhs.size() < 1) return -1;
	vector<double> ws;
	rhs[0].GetWeights(ws);

	//to debug - use less rechits
	//int nrhs = 100;
	//rhs.resize(nrhs);
	//rhs.shrink_to_fit();


	
	BayesCluster *algo = new BayesCluster(rhs);
	if(smeared) algo->SetDataSmear(smear);
	//set time resolution smear: c^2 + n^2/e^2
	//remember time is already in ns
	//e = w/gev
	if(timesmeared) algo->SetTimeResSmear(tres_c, tres_n*gev);
	algo->SetThresh(thresh);
	algo->SetAlpha(alpha);
	algo->SetSubclusterAlpha(emAlpha);
	algo->SetVerbosity(verb);
	clock_t t = clock();
	//jets
	if(obj == 0 || obj == 2){
		cout << "Using clustering strategy " << strat << ": ";
		if(strat == 0 || strat == 1){
			if(strat == 0){
				cout << "Delauney (NlnN)" << endl;
				//to track computation time from beginning of program
				trees = algo->NlnNCluster();
			}
			else if(strat == 1){
				cout << "N^2 (naive)" << endl;
				//to track computation time from beginning of program
				if(obj == 0) trees = algo->N2Cluster();
			}
			if(viz){
				//plotting stuff here
				FullViz3D plots = FullViz3D(trees);
				plots.SetVerbosity(verb);
				plots.SetTransfFactor(gev);
				//add info of true jets
				plots.AddTrueJets(jets);
				plots.Write(fname);
			}
			cout << jets.size() << " true jets" << endl;
		}
		else if(strat == 2){
			cout << "mixture model with pre-clustered AK4 jets" << endl;
			string oname = "";
			if(viz) oname = fname;
			GaussianMixture* gmm = algo->SubCluster(oname);
			int nk = gmm->GetNClusters();
			cout << nk << " clusters in model." << endl;
			map<string, Matrix> params;
			for(int k = 0; k < nk; k++){
				params = gmm->GetPriorParameters(k);
				cout << "cluster " << k << " has time center " << params["mean"].at(2,0) << " with mixing coeff " << params["pi"].at(0,0) << endl;
			}
			

		}
		else{ cout << " undefined. Please use --strategy(-s) [strat] with options 0 (NlnN) or 1 (N^2)." << endl; return -1; }
	}
	//photons
	else if(obj == 1){
		string oname = "";
		if(viz) oname = fname;
		GaussianMixture* gmm = algo->SubCluster(oname);
		int nk = gmm->GetNClusters();
		cout << nk << " clusters in model." << endl;
		map<string, Matrix> params;
		for(int k = 0; k < nk; k++){
			params = gmm->GetPriorParameters(k);
			cout << "cluster " << k << " has time center " << params["mean"].at(2,0) << " with mixing coeff " << params["pi"].at(0,0) << endl;
		}
	}	



	else{ cout << "Object undefined. Please use --object [obj] with options 0 for jets or 1 for photons." << endl; return -1; }


	//causing crashes - //file->Close();
	//delete file;
	

	t = clock() - t;
	cout << "Total program time elapsed: " << (float)t/CLOCKS_PER_SEC << " seconds." << endl;
}
