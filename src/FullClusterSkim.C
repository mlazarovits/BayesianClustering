#include "JetProducer.hh"
#include "JetSkimmer.hh"
#include "PhotonSkimmer.hh"
#include "SuperClusterSkimmer.hh"
#include "BayesCluster.hh"
#include "FullViz3D.hh"
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
	bool smear = true;
	bool timesmear = false;
	//by default in BayesCluster
	bool distconst = true;
	int evti = 0; //for skimming from evti to evtj
	int evtj = 0;
	int obj = 0; //object to cluster (0 : jets, 1 : photons)
	//this should be in N/GeV
	double gev = -999;
	//put cuts on jets (ie min pt) here
	double minpt = 30;
	double minnrhs = 15;
	double minEmE = 20;
	double minRhE = 0.5;
	double maxRhE = -999;
	bool frac = false;
	bool calib = true;
	bool spikes = false;
	int skip = 1;
	int bh = 1;
	bool iso = true;
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
		if(strncmp(argv[i],"--noWeight", 10) == 0){
    	 		weighted = false;
			cout << "Turning off energy weighting." << endl;
   		}
		if(strncmp(argv[i],"--noSmear", 9) == 0){
    	 		smear = false;
			cout << "Turning off smearing." << endl;
   		}
		if(strncmp(argv[i],"--timeSmear", 11) == 0){
    	 		timesmear = true;
			cout << "Turning on time smearing." << endl;
   		}
		if(strncmp(argv[i],"--noDist", 8) == 0){
    	 		distconst = false;
			cout << "Turning off distance constraint." << endl;
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
		if(strncmp(argv[i],"--minRhE", 8) == 0){
			i++;
    	 		minRhE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--maxRhE", 8) == 0){
			i++;
    	 		maxRhE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--applyFrac", 11) == 0){
			frac = true;
			cout << "Apply fractions for rec hits." << endl;
   		}
		if(strncmp(argv[i],"--noCalibrate", 13) == 0){
			calib = false;
			cout << "Not applying time calibration to rhs" << endl;
   		}
		if(strncmp(argv[i],"--rejectSpikes", 14) == 0){
			spikes = true;
			cout << "Rejecting spikes with swiss cross cut" << endl;
   		}
		if(strncmp(argv[i],"--BHFilter", 10) == 0){
			i++;
    	 		bh = std::stoi(argv[i]);
   		}
		if(strncmp(argv[i],"--skip", 6) == 0){
			i++;
    	 		skip = std::stoi(argv[i]);
   		}
		if(strncmp(argv[i],"--noIso", 7) == 0){
			iso = false;
			cout << "Not applying isolation to photons" << endl;
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
   		cout << "   --object [obj]                set object to cluster (0 : jets, default; 1 : superclusters, 2 : photons)" << endl;
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/30 GeV)" << endl;
   		cout << "   --minpt [minpt]               set minimum pt (default = 30 GeV)" << endl;
   		cout << "   --minNrhs [minnrhs]           set minimum # of rhs (default = 2)" << endl;
   		cout << "   --minemE [mineme]             set minimum ECAL energy (default = 10 GeV)" << endl;
   		cout << "   --minRhE [minRhe]             set minimum ECAL rechit energy (default = 0.5 GeV)" << endl;
   		cout << "   --maxRhE [maxRhe]             set maximum ECAL rechit energy (default = -999, off)" << endl;
   		cout << "   --BHFilter [bh]               set how beam halo filter is applied (0 : not applied, 1 : applied (default), 2 : inversely applied)" << endl;
   		cout << "   --skip [skip]                 set skip for event loop (default = 1)" << endl;
   		cout << "   --evtFirst [i] --evtLast [j]  skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
   		cout << "   --noSmear                     turns off smearing data (default = true, on)" << endl;
   		cout << "   --timeSmear                   turns on time smearing data (default = false, off)" << endl;
   		cout << "   --noWeight                    turns off weighting data points (default = false, on)" << endl;
   		cout << "   --noDist                      turns off distance constraint: clusters must be within pi/2 in phi (default = false, on)" << endl;
   		cout << "   --applyFrac                   applying fractions for rec hits PHOTONS ONLY (default = false, off)" << endl;
   		cout << "   --noCalibrate                 turn off channel-by-channel calibration for rh time (default = false, on)" << endl;
   		cout << "   --rejectSpikes                reject spikes based on swiss cross cut (default = false, off)" << endl;
   		cout << "   --noIso                       turn off isolation in preselection (photons only, default = true, on)" << endl;
   		cout << "Example: ./jetAlgo.x -a 0.5 -t 1.6 --viz" << endl;

   		return 0;
  	}


	cout << "Free sha-va-ca-doo!" << endl;
	
	if(gev == -999){
		if(obj == 0) gev = 1./10.;
		else gev = 1./30.;
	}

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

	if(gSystem->AccessPathName(in_file.c_str())){
		cout << "Error: file " << in_file << " not found." << endl;
		return -1;
	}	

	string fname;
	if(obj == 0)
		fname = "jets";
	else if(obj == 1)
		fname = "supercluster";
	else if(obj == 2)
		fname = "photon";
	else{
		cout << "Object number " << obj << " not supported. Only 0 : jets, 1 : superclusters, 2 : photons." << endl;
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
cout << "fname " << fname << endl;
	//choose time calibration file
	string calibfile = "";
        if(fname.find("GJets") != string::npos)
                calibfile = "info/KUCMS_GJets_R17_v16_rhE5_mo_Cali.root";
        else if(fname.find("JetHT") != string::npos)
                calibfile = "info/KUCMS_JetHT_R17_v18_rhE5_Cali.root";
        else if(fname.find("DEG") != string::npos || fname.find("DoubleEG") != string::npos || fname.find("EGamma") != string::npos)
                calibfile = "info/KUCMS_DoubleEG_R17_v18_rhE5_Cali.root";
        else if(fname.find("QCD") != string::npos)
                calibfile = "info/KUCMS_QCD_R17_v16_rhE5_mo_Cali.root";
        //else default to GJets
        else
                calibfile = "info/KUCMS_GJets_R17_v16_rhE5_mo_Cali.root";
	if(obj == 0){
		cout << "jets" << endl;
		JetSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		skimmer.SetMinPt(minpt);
		skimmer.SetMinNrhs(minnrhs);
		skimmer.SetMinRhE(0.5);
		skimmer.SetMinRhE_PV(minRhE);
		if(maxRhE != -999) skimmer.SetMaxRhE(maxRhE);
		skimmer.SetMinEmE(minEmE);
		skimmer.SetSpikeRejection(spikes); //if true, reject spikes
		skimmer.SetSkip(skip);
		if(calib) skimmer.SetTimeCalibrationMap(calibfile);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
		//set alpha, EMalpha
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetSmear(smear);
		skimmer.SetTimeSmear(timesmear); 
		skimmer.SetBeamHaloFilter(bh);
		//do only mm/true jet pv times
		skimmer.Skim();
	}
	else if(obj == 1){
		cout << "superclusters" << endl;
		SuperClusterSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		bool data;
		if(in_file.find("SIM") != string::npos)
			data = false;
        	else
			data = true;
		if(calib) skimmer.SetTimeCalibrationMap(calibfile);
		if(iso){
			cout << "Applying isolation preselection for training labels." << endl;
			skimmer.SetIsoCuts();
		}
		else cout << "Not applying isolation preselection for training labels." << endl;
		skimmer.SetMinPt(minpt);
		if(iso) skimmer.SetIsoCuts();
		skimmer.SetMinRhE(minRhE);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
        	skimmer.ApplyFractions(frac);
		//skimmer.SetDebug(debug);
		//set EMalpha
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetSmear(smear);
		skimmer.SetTimeSmear(timesmear); 
		skimmer.SetBeamHaloFilter(bh);
        	skimmer.Skim();
	}
	else if(obj == 2){
		cout << "photons" << endl;
		PhotonSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		bool data;
		if(in_file.find("SIM") != string::npos)
			data = false;
        	else
			data = true;
		if(calib) skimmer.SetTimeCalibrationMap(calibfile);
		if(iso) skimmer.SetIsoCuts();
		skimmer.SetMinRhE(minRhE);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
        	skimmer.ApplyFractions(frac);
		//skimmer.SetDebug(debug);
		//set EMalpha
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetSmear(smear);
		skimmer.SetTimeSmear(timesmear); 
        	skimmer.Skim();
	}
	if(calib) cout << "Using calibration file " << calibfile << endl;
	else cout << "No timing calibration applied" << endl;
        return 0;

}
