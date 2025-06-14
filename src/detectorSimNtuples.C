#include "BasicDetectorSim.hh"
#include "Jet.hh"
#include <string>
#include <iostream>
#include <TFile.h>
#include <Math/Vector4D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TBranch.h>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using PtEtaPhiEVector = ROOT::Math::PtEtaPhiEVector;
using XYZTVector = ROOT::Math::XYZTVector;

int main(int argc, char *argv[]){
	vector<Jet> rhs, jets;
	vector<PtEtaPhiEVector> genmoms;
	vector<PtEtaPhiEVector> recomoms;
	vector<XYZTVector> genpos;
	vector<XYZTVector> recopos;
	vector<vector<JetPoint>> ems;
	int nevts = 1;
	int evt = 0;
	int evti = 0;
	int evtj = 0;
	int verb = 0;
	bool hprint = false;
	bool pu = false;
	bool skim = false;
	bool ntuple = false;	
	string oname = "";
	bool ttbar = false;
	bool qcd = false;
	bool sig_delayed = false;
	bool sig_boosted = false;
	double spikeProb = 0.;
	double energy_c = 0.26;
	double tres_cte = 0.1727 * 1e-9;
	double tres_stoch = 0.5109 * 1e-9;
	double tres_noise = 2.106 * 1e-9;
	double ethresh = 0.1; //zero suppression threshold for rechit reconstruction (and inclusion in AK4 jet)
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i],"--help", 6) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"-h", 2) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"--verbosity", 11) == 0){
    	 		i++;
			verb = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-v", 2) == 0){
    	 		i++;
			verb = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--pileup", 8) == 0){
    	 		pu = true;
   		}
		if(strncmp(argv[i],"-pu", 3) == 0){
    	 		pu = true;
   		}
		if(strncmp(argv[i],"--skim", 6) == 0){
    	 		skim = true;
   		}
		if(strncmp(argv[i],"--ntuple", 8) == 0){
    	 		ntuple = true;
   		}
		if(strncmp(argv[i],"--nevts", 7) == 0){
    	 		i++;
			nevts = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--output", 8) == 0){
     			i++;
    	 		oname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-o", 2) == 0){
     			i++;
    	 		oname = string(argv[i]);
   		}
		if(strncmp(argv[i],"--ttbar", 7) == 0){
    	 		ttbar = true;
   		}
		if(strncmp(argv[i],"--QCD", 5) == 0){
    	 		qcd = true;
   		}
		if(strncmp(argv[i],"--sigDelayed", 12) == 0){
    	 		sig_delayed = true;
   		}
		if(strncmp(argv[i],"--sigBoosted", 12) == 0){
    	 		sig_boosted = true;
   		}
		if(strncmp(argv[i],"--spikeProb", 11) == 0){
			i++;
    	 		spikeProb = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--eThresh", 9) == 0){
			i++;
    	 		ethresh = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--energyCte", 11) == 0){
			i++;
    	 		energy_c = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--evtFirst", 6) == 0){
                        i++;
                        evti = std::atoi(argv[i]);
                }
                if(strncmp(argv[i],"--evtLast", 6) == 0){
                        i++;
                        evtj = std::atoi(argv[i]);
                }
		if(strncmp(argv[i],"--tResCte", 9) == 0){
			i++;
    	 		tres_cte = std::stod(argv[i])*1e-9; //need to convert to s
   		}
		if(strncmp(argv[i],"--tResNoise", 11) == 0){
			i++;
    	 		tres_noise = std::stod(argv[i])*1e-9; //need to convert to s
   		}
		if(strncmp(argv[i],"--tResStoch", 11) == 0){
			i++;
    	 		tres_stoch = std::stod(argv[i])*1e-9; //need to convert to s
   		}

	}
	if(hprint){
		cout << "Making simulated ntuples" << endl;
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --pileup(-pu)                 simulate pileup" << endl;
		cout << "   --ttbar                       simulate ttbar" << endl;
		cout << "   --QCD                         simulate QCD" << endl;
		cout << "   --sigDelayed                  simulate delayed signal" << endl;
		cout << "   --sigBoosted                  simulate boosted signal" << endl;
		cout << "   --output(-o) [ofile]          set output file name" << endl; 
   		cout << "   --nevts [nevts]               set number of events to simulate (default = 1)" << endl;
   		cout << "   --eThresh [ethresh]           set energy threshold for rechit reco (default = 0.1)" << endl;
   		cout << "   --spikeProb [p]               set probability of spike occuring (default = 0, off)" << endl;
   		cout << "   --energyCte [c]               set energy smearing constant (default = 0.26)" << endl;
   		cout << "   --tResCte [t]                 set time smearing constant parameter in ns (default = 0.133913 ns)" << endl;
   		cout << "   --tResNoise [t]               set time smearing noise (n*n/(e*e)) parameter in ns (default = 0.00691415 ns)" << endl;
   		cout << "   --tResStoch [t]               set time smearing stochastic (s*s/e) parameter in ns (default = 1.60666 ns)" << endl;
		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
		cout << "   --evtFirst [i] --evtLast [j]  skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;

		return -1;	
	}


	if(!ttbar && !qcd && !sig_delayed && !sig_boosted){
		cout << "No process specified to simulate. Exiting..." << endl;
		return -1;
	}
	else{
		if(!oname.empty()){
			if(oname.find("condor") == string::npos)
				oname = "simNtuples_"+oname;
		}
		else
			oname = "simNtuples";
	}
	if(spikeProb < 0 || spikeProb > 1){
		cout << "Invalid spike probability " << spikeProb << ". Must be [0,1]" << endl;
		return -1;
	}	

	cout << "Simulating events from ";	
	if(ttbar){
		cout << "ttbar ";
		if(oname.find("ttbar") == string::npos) oname += "_ttbar";	
	}
	if(qcd){
		cout << "QCD ";
		if(oname.find("QCD") == string::npos) oname += "_QCD";	
	}
	cout << endl;

	//TODO: change onames when processes are decided
	if(sig_delayed){
		cout << "delayed signal " << endl;
		oname += "sig_delayed";
	}

	if(sig_boosted){
		cout << "boosted signal " << endl;
		oname += "sig_boosted";
	}

	if(pu){
		cout << "and pileup " << endl;
		oname += "_PU";
	}
       //make sure evti < evtj
       if(evti > evtj){
       	int evt = evtj;
       	evtj = evti;
       	evti = evt;
       }

	//consider doing det from pythia cmnd card
	BasicDetectorSim det;
	det.SetNEvents(nevts);
	//for reconstructing rechits
	det.SetEnergyThreshold(ethresh);
	det.SetEventRange(evti,evtj);
	det.SetVerbosity(verb);
	det.SetEnergySmear(energy_c);
	det.SetTimeResCts(tres_cte, tres_stoch, tres_noise);
	if(ttbar) det.SimTTbar();
	if(qcd) det.SimQCD();
	//if(sig_delayed)
	//if(sig_boosted)
	if(pu) det.TurnOnPileup();
	if(spikeProb > 0) det.TurnOnSpikes(0.01);
	
	///////make ntuple///////
	det.InitTree(oname+".root");
	det.SimulateEvents();
	det.WriteTree();
	

}
