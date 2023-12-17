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
	string oname = "";
	string in_file = "";//rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";
	bool hprint = false;
	double thresh = 1.;
	double alpha = 0.1;
	double emAlpha = 0.5;
	int verb = 0;
	bool weighted = true;
	bool nosmear = false;
	//by default in BayesCluster
	bool distconst = true;
	//clustering strategy for skimmer
	int strat = 0; //0 is NlnN
	int evti = 0; //for skimming from evti to evtj
	int evtj = 0;
	int obj = 0; //object to cluster (0 : jets, 1 : photons)
	//this should be in N/GeV
	double gev = 1./30.;
	//put cuts on jets (ie min pt) here
	double minpt = 30;
	double minnrhs = 2;
	double minEmE = 10;
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
    	 		weighted = false;
			cout << "Turning off energy weighting." << endl;
   		}
		if(strncmp(argv[i],"--noTimeSmear", 11) == 0){
    	 		nosmear = true;
			cout << "Turning off time smearing." << endl;
   		}
		if(strncmp(argv[i],"--dist", 6) == 0){
    	 		distconst = false;
			cout << "Turning off distance constraint." << endl;
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
		if(strncmp(argv[i],"--gev", 5) == 0){
			i++;
    	 		gev = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minpt", 7) == 0){
			i++;
    	 		minpt = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minNrhs", 9) == 0){
			i++;
    	 		minnrhs = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minemE", 8) == 0){
			i++;
    	 		minEmE = std::stod(argv[i]);
   		}

	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --input(-i) [file]            input root file" << endl;
   		cout << "   --output(-o) [file]           output root file" << endl;
   		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 0.1)" << endl;
   		cout << "   --EMalpha(-EMa) [a]           sets concentration parameter alpha for variational EM GMM (default = 0.5)" << endl;
   		cout << "   --thresh(-t) [t]              sets threshold for cluster cutoff" << endl;
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --strategy(-s) [strat]        set clustering strategy for skimmer (0 : NlnN (default), 1 : N^2, 2 : GMM only, does not apply to photons)" << endl;
   		cout << "   --object [obj]                set object to cluster (0 : jets, default; 1 : photons)" << endl;
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/30 GeV)" << endl;
   		cout << "   --minpt [minpt]               set minimum pt (default = 30 GeV)" << endl;
   		cout << "   --minNrhs [minnrhs]           set minimum # of rhs (default = 2)" << endl;
   		cout << "   --minemE [mineme]             set minimum ECAL energy (default = 10 GeV)" << endl;
   		cout << "   --evtFirst [i] --evtLast [j]  skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
   		cout << "   --NoTimeSmear                 turns off smearing data in time dimension (default = false)" << endl;
   		cout << "   --NoWeight                    turns off weighting data points (default = false)" << endl;
   		cout << "   --NoDist                      turns off distance constraint: clusters must be within pi/2 in phi (default = false)" << endl;
   		cout << "Example: ./jetAlgo.x -a 0.5 -t 1.6 --viz" << endl;

   		return 0;
  	}


	cout << "Free sha-va-ca-doo!" << endl;
	
	string cmslab, version;	
	TFile* file = nullptr;
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
		file = TFile::Open(in_file.c_str());
		if(evti != evtj) cout << "Skimming events " << evti << " to " << evtj << " for ";
		else cout << "Skimming all events for ";
	}
	/////MAKE DATA WITH PYTHIA + BASIC DETECTOR SIM//////
	else{
		in_file = "BDSIM";
	} 

	if(gSystem->AccessPathName(in_file.c_str())){
		cout << "Error: file " << in_file << " not found." << endl;
		return -1;
	}	

	string fname;
	if(obj == 0)
		fname = "jets";
	else if(obj == 1)
		fname = "photon";
	else{
		cout << "Object number " << obj << " not supported. Only 0 : jets, 1 : photons." << endl;
		return -1;
	}

	//condor is in outname
	if(oname.find("condor") != string::npos){
		fname = oname;	
	}
	else{
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

		
		if(evti != evtj) fname += "_evt"+std::to_string(evti)+"to"+std::to_string(evtj);
		if(obj != 1) fname += "_bhcAlpha"+a_string+"_emAlpha"+ema_string+"_thresh"+t_string+"_";
		else fname += "_emAlpha"+ema_string+"_thresh"+t_string+"_";
		fname += "NperGeV"+gev_string+"_";
		fname += cmslab; //long sample name
		//fname += cmslab.substr(0,cmslab.find("_")); //short sample name
		//fname += version; //"_v9"

	}
	fname = fname+".root";
	
	//make sure evti < evtj
	if(evti > evtj){
		int evt = evtj;
		evtj = evti;
		evti = evt;
	}

	if(obj == 0){
		cout << "jets" << endl;
		JetSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		skimmer.SetMinPt(minpt);
		skimmer.SetMinNrhs(minnrhs);
		skimmer.SetMinEmE(minEmE);
		if(in_file.find("JetHT") != string::npos)
			skimmer.SetData(true);
		if(strat == 2) 
			skimmer.GMMOnly();
		else
			skimmer.SetStrategy(strat);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
		//set alpha, EMalpha
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetTimeSmear(!nosmear);
		//do only mm/true jet pv times
		skimmer.Skim();
	}
	else if(obj == 1){
		cout << "photons" << endl;
		PhotonSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		bool data;
		if(in_file.find("SIM") != string::npos)
			data = false;
        	else
			data = true;
		skimmer.SetIsoCuts();
		skimmer.SetData(data);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
        	//skimmer.SetDebug(debug);
		//set EMalpha
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetTimeSmear(false);
        	skimmer.Skim();
	}
        return 0;

}
