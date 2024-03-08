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


	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                    print options" << endl;
   		cout << "   --pileup(-pu)                 simulate pileup" << endl;
   		cout << "   --nevts [nevts]               set number of events to simulate (default = 1)" << endl;
   		cout << "   --gev [gev]                   set energy weight transfer factor in N/GeV (default = 1/10 GeV)" << endl;
		cout << "   --verbosity(-v) [verb]        set verbosity (default = 0)" << endl;
   		cout << "   --skim                        make histograms (default = false)" << endl;
   		cout << "   --ntuple                      save to file (default = false)" << endl;
		return -1;	
	}


	
	//consider doing det from pythia cmnd card

	BasicDetectorSim det;
	det.SetNEvents(nevts);
	det.SetEnergyThreshold(1.); //set to 1 GeV
	//set energy transfer factor in N/GeV
	det.SetTransferFactor(gev);
	det.SetVerbosity(verb);
	det.SimTTbar();
	if(pu) det.TurnOnPileup();
	//det.TurnOnSpikes(0.01);

	



	///////make ntuple///////
	if(ntuple){
		det.InitTree();
		det.SimulateEvents();
		det.WriteTree();
	}
	///////make histograms///////
	if(skim){
		cout << "skimming for evt " << evt << endl;
		//default arg is all events
		det.SimulateEvents(evt);
		det.GetRecHits(rhs);
		det.GetTrueJets(jets);
		det.GetParticlesMom(genmoms, recomoms);
		det.GetParticlesPos(genpos, recopos);
		det.GetEmissions(ems);
		cout << rhs.size() << " rechits and " << jets.size() << " true jets in event " << evt << endl;
		return -1;
		string fname = "plots/detectorSimSkim.root";
		cout << "Writing to " << fname << endl;
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
	}

}
