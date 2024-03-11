#include "BasicDetectorSim.hh"
#include "Jet.hh"
#include <string>
#include <iostream>
#include <TFile.h>
#include <Math/Vector4D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TBranch.h>

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
	int verb = 0;
	bool hprint = false;
	bool pu = false;
	bool skim = false;
	bool ntuple = false;	
	double gev = 1./10.;
	string oname = "ntuples";
	bool ttbar = false;
	bool qcd = false;
	bool sig_delayed = false;
	bool sig_boosted = false;
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
		if(strncmp(argv[i],"--gev", 5) == 0){
			i++;
    	 		gev = std::stod(argv[i]);
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
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/10 GeV)" << endl;
		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
		return -1;	
	}


	if(!ttbar && !qcd && !sig_delayed && !sig_boosted){
		cout << "No process specified to simulate. Exiting..." << endl;
		return -1;
	}	
	//consider doing det from pythia cmnd card

	BasicDetectorSim det;
	det.SetNEvents(nevts);
	det.SetEnergyThreshold(1.); //set to 1 GeV
	//set energy transfer factor in N/GeV
	det.SetTransferFactor(gev);
	det.SetVerbosity(verb);
	if(ttbar) det.SimTTbar();
	//if(qcd)
	//if(sig_delayed)
	//if(sig_boosted)
	if(pu) det.TurnOnPileup();
	//det.TurnOnSpikes(0.01);
	///////make ntuple///////
	det.InitTree("rootfiles/"+oname+".root");
	det.SimulateEvents();
	det.WriteTree();
	

}
