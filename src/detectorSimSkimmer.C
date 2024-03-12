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
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/30 GeV)" << endl;
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
		cout << "1 ctor oname " << oname << endl;
		string match = "simNtuples_";
		oname = oname.substr(oname.find(match)+match.size(),oname.find(".root")-(oname.find(match)+match.size()))+".root";
		oname = "simSkim_"+oname;
		cout << "ctor oname " << oname << endl;
	}
	else{
		string oname_extra = file->GetName();
		cout << "1 ctor oname_extra " << oname_extra << endl;
		string match = "simNtuples_";
		oname_extra = oname_extra.substr(oname_extra.find(match)+match.size(),oname_extra.find(".root")-(oname_extra.find(match)+match.size()))+".root";
		oname = "simSkim_"+oname+"_"+oname_extra;
		cout << "ctor oname " << oname << endl;

	}	

	BHCJetSkimmer skimmer(file);
	//skimmer.SetMinPt(minpt);
	//skimmer.SetMinNrhs(minnrhs);
	skimmer.SetOutfile(oname);
	skimmer.SetTransferFactor(gev);
	skimmer.SetEventRange(evti,evtj);
	skimmer.Skim();
	/*
	TFile f(fname.c_str(),"RECREATE");
	int nhists = 7;
	TH1D* hists[nhists];
	hists[0] = new TH1D("rhtimes","rhtimes",100,-25,25);
	hists[1] = new TH1D("rhE","rhE",100,0,100);
	hists[2] = new TH1D("rheta","rheta",50,-1.5,1.5);
	hists[3] = new TH1D("rhphi","rhphi",50,-3.2,3.2);
	//gen hists - need for "calibration" (validation of detector effects)
	//ratio between reco and gen energy
	hists[4] = new TH1D("recoE_genE","recoE_genE",10,-0.5,10.5);
	hists[5] = new TH1D("recophi","recophi",50,-3.2,3.2);
	hists[6] = new TH1D("genphi","genphi",50,-3.2,3.2);


	int n2dhists = 6;
	TH2D* hists2d[n2dhists];
	hists2d[0] = new TH2D("rheta_rhphi","rheta_rhphi", 25, -1.5, 1.5, 25, -3.2, 3.2);
	hists2d[1] = new TH2D("recoeta_recophi","recoeta_recophi", 25, -1.5, 1.5, 25, -3.2, 3.2);
	hists2d[2] = new TH2D("geneta_genphi","geneta_genphi", 25, -1.5, 1.5, 25, -3.2, 3.2);
	hists2d[3] = new TH2D("recoMingenPhi_recopt","recoMingenPhi_recopt",50,0.,3.2,50,0,90);	
	hists2d[4] = new TH2D("recoTogenE_genE","recoTogenE_genE",20,-0.5,10.5,100,0,250.);
	hists2d[5] = new TH2D("rht_rhE","rht_rhE",100,-25,25,100,0,100);

	//fill rh histograms
	for(int i = 0; i < rhs.size(); i++){
		hists[0]->Fill(rhs[i].t());
		hists[1]->Fill(rhs[i].E());
		hists[2]->Fill(rhs[i].eta());
		hists[3]->Fill(rhs[i].phi_std());

		hists2d[0]->Fill(rhs[i].eta(), rhs[i].phi_std());
		hists2d[5]->Fill(rhs[i].t(),rhs[i].E());
	} 

	//fill reco/gen hists
	cout << genmoms.size() << " gen particles and " << recomoms.size() << " reco particles in event " << evt << endl;
	for(int i = 0; i < genmoms.size(); i++){
		hists[4]->Fill(recomoms[i].e()/genmoms[i].e());
		hists[5]->Fill(recopos[i].phi());
		hists[6]->Fill(genpos[i].phi());

		hists2d[1]->Fill(recopos[i].eta(), recopos[i].phi());
		hists2d[2]->Fill(genpos[i].eta(), genpos[i].phi());
		hists2d[3]->Fill(fabs(recopos[i].phi() - genpos[i].phi()), genmoms[i].pt());
		hists2d[4]->Fill(recomoms[i].e()/genmoms[i].e(),genmoms[i].e());

	}


	f.cd();
	for(int i = 0; i < nhists; i++){
		hists[i]->Write();
	}
	for(int i = 0; i < n2dhists; i++){
		hists2d[i]->Write();
	}
	f.Close();
	*/

}
