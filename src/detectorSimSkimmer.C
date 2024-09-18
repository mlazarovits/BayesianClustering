#include "BasicDetectorSim.hh"
#include "Jet.hh"
#include "BHCJetSkimmer.hh"
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
	int evti = 0; //for skimming from evti to evtj
	int evtj = 0;
	int verb = 0;
	bool hprint = false;
	bool pu = false;
	bool skim = false;
	bool ntuple = false;	
	double gev = 1./10.;
	//set clustering strategy
	//0 = NlnN
	//1 = N2
	int strat = 0;
	string infile = "";
	string oname = "";
	double thresh = 1.;
	double emAlpha = 0.5;
	double alpha = 0.1;
	double minpt = 30.;
	double minnrhs = 15;
	double minRhE = 0.5;
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
		if(strncmp(argv[i],"--gev", 5) == 0){
			i++;
    	 		gev = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--strategy", 10) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-s", 2) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--input", 7) == 0){
     			i++;
    	 		infile = string(argv[i]);
   		}
		if(strncmp(argv[i],"-i", 2) == 0){
     			i++;
    	 		infile = string(argv[i]);
   		}
		if(strncmp(argv[i],"--output", 8) == 0){
     			i++;
    	 		oname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-o", 2) == 0){
     			i++;
    	 		oname = string(argv[i]);
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
		if(strncmp(argv[i],"--evtFirst", 6) == 0){
    	 		i++;
			evti = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--evtLast", 6) == 0){
    	 		i++;
			evtj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--minpt", 7) == 0){
			i++;
    	 		minpt = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minNrhs", 9) == 0){
			i++;
    	 		minnrhs = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minRhE", 8) == 0){
			i++;
    	 		minRhE = std::stod(argv[i]);
   		}



	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --input(-i) [file]            input root file" << endl;
   		cout << "   --output(-o) [file]           output root file" << endl;
   		cout << "   --strategy(-s) [strat]        sets clustering strategy (0 = NlnN, default; 1 = N2)" << endl;
		cout << "   --alpha(-a) [a]               sets concentration parameter alpha for DPM in BHC (default = 0.1)" << endl;
   		cout << "   --EMalpha(-EMa) [a]           sets concentration parameter alpha for variational EM GMM (default = 0.5)" << endl;
   		cout << "   --thresh(-t) [t]              sets threshold for cluster cutoff" << endl;
   		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/10 GeV)" << endl;
   		cout << "   --minpt [minpt]               set minimum pt (default = 30 GeV)" << endl;
   		cout << "   --minNrhs [minnrhs]           set minimum # of rhs (default = 2)" << endl;
   		cout << "   --minRhE [minRhe]             set minimum ECAL rechit energy (default = 0.5 GeV)" << endl;
   		cout << "   --evtFirst [i] --evtLast [j]  skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
   		cout << "Example: ./detectorSimSkimmer.x -i rootfiles/simNtuples_ttbar.root -a 0.5 -t 1.6" << endl;
		return 0;	
	}

	if(gSystem->AccessPathName(infile.c_str())){
		cout << "Error: file " << infile << " not found." << endl;
		return -1;
	}

	if(evti != evtj) cout << "Skimming events " << evti << " to " << evtj << " for ";
	else cout << "Skimming all events for ";

	//make sure evti < evtj
	if(evti > evtj){
		int evt = evtj;
		evtj = evti;
		evti = evt;
	}
	TFile* file = TFile::Open(infile.c_str());
	if(oname.empty()){
		oname = file->GetName();
		string match = "simNtuples_";
		oname = oname.substr(oname.find(match)+match.size(),oname.find(".root")-(oname.find(match)+match.size()));
		oname = "simSkim_"+oname;
		
		std::stringstream stream;
		string gev_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3) << gev;
		gev_string = stream.str();
		int idx = gev_string.find(".");
		gev_string.replace(idx,1,"p");	
		oname += "_NperGeV"+gev_string;

		if(strat == 0) oname += "_NlnN";
		else oname += "_N2";

		oname += ".root";
	}
	else{
		if(oname.find("condor") != string::npos){ 
			oname = oname+".root";
		}
		else{
			string oname_extra = file->GetName();
			string match = "simNtuples_";
			oname_extra = oname_extra.substr(oname_extra.find(match)+match.size(),oname_extra.find(".root")-(oname_extra.find(match)+match.size()))+".root";
			oname = "simSkim_"+oname+"_"+oname_extra;
		}
	}	

	cout << "Energy transfer factor: " << gev << endl;
	BHCJetSkimmer skimmer(file);
	skimmer.SetOutfile(oname);
	skimmer.SetMinRhE(minRhE);
	skimmer.SetStrategy(strat);
	skimmer.SetVerbosity(verb);
	skimmer.SetTransferFactor(gev);
	skimmer.SetAlpha(alpha);
	skimmer.SetSubclusterAlpha(emAlpha);
	skimmer.SetThreshold(thresh);
	skimmer.SetEventRange(evti,evtj);
	skimmer.Skim();

}
