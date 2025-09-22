#include <iostream>
#include <fstream>
#include "../include/BaseProducer.hh"
#include "../include/Jet.hh"
#include "../include/JetProducer.hh"
#include "../include/PhotonProducer.hh"
#include "../include/SampleWeight.hh"
#include "TSystem.h"
//#include "../include/BaseSkimmer.hh"
using std::cout;
using std::endl;

//make tchain from filelist
TChain* MakeTChain(string flist){
        if(gSystem->AccessPathName(flist.c_str())){
                cout << "Error: file " << flist << " doesn't exist." << endl;
                return nullptr;
        }
        std::ifstream infile(flist);
        TChain* ch = new TChain("tree/llpgtree");
        ch->SetTitle(flist.c_str());
        string file;
        cout << "TChaining files in " << flist << "..." << endl;
        while(std::getline(infile,file)){
                if( file[0] == '#' ) continue;
                if(gSystem->AccessPathName(file.c_str())){
                        cout << "Skipping file " << file << " - not found." << endl;
                        continue;
                }
                //std::cout << "--  adding file: " << file << std::endl;
                ch->Add(file.c_str());
        }
        cout << "Done TChaining" << endl;
        return ch;
}

//returns vector sum of given objects
Jet VectorSum(vector<Jet>& jets){
        Jet ret;
        for(auto j : jets){
                ret.add(j);
        }
        return ret;
}

int EventWeightCalc_FileList(string selection = "", string year = "18", bool isoBkgSel = true){
	gSystem->Load("lib/libBayesianClustering.so");
	double gev_jet = 0.1;
	double gev_pho = 1./30.;

	double minht = 50;
	double maxmet = 50;
	double pi = 4*atan(1);

	map<string, string> recodates;
	recodates["17"] = "RunIIFall17DRPremix";
	recodates["18"] = "RunIISummer20UL18";

	//just for MCs - data weight = 1
	vector<string> filelists;
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT50to100_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT100to200_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT200to300_AODSIM_"+recodates[year]+"_list");
        filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT300to500_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT500to700_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT700to1000_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT1000to1500_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT1500to2000_AODSIM_"+recodates[year]+"_list");
	filelists.push_back("filelists/kucmsntuple_QCD_R"+year+"_"+selection+"_v31_QCD_HT2000toInf_AODSIM_"+recodates[year]+"_list");
	if(year == "18")
		recodates[year] += "RECO";
//	filelists.push_back("filelists/kucmsntuple_GJets_R"+year+"_"+selection+"_v31_GJets_HT-40To100_TuneCP5_AODSIM_"+recodates[year]+"_list");
//	filelists.push_back("filelists/kucmsntuple_GJets_R"+year+"_"+selection+"_v31_GJets_HT-100To200_TuneCP5_AODSIM_"+recodates[year]+"_list");
//	filelists.push_back("filelists/kucmsntuple_GJets_R"+year+"_"+selection+"_v31_GJets_HT-200To400_TuneCP5_AODSIM_"+recodates[year]+"_list");
//	filelists.push_back("filelists/kucmsntuple_GJets_R"+year+"_"+selection+"_v31_GJets_HT-400To600_TuneCP5_AODSIM_"+recodates[year]+"_list");
//	filelists.push_back("filelists/kucmsntuple_GJets_R"+year+"_"+selection+"_v31_GJets_HT-600ToInf_TuneCP5_AODSIM_"+recodates[year]+"_list");



	ofstream ofile;
	string ofilename = "info/EventWeights_"+selection+"_R"+year;
	if(isoBkgSel)
		ofilename += "_isoBkgSel";
	ofilename += ".txt";
	ofile.open(ofilename);
	//dont write just for remembering what's being written
	//ofile << "file	jet_weight	pho_weight" << endl;
	for(int f = 0; f < filelists.size(); f++){
		cout << "Sample " << filelists[f] << endl;
		//loop over all files in a list
		TChain* chj = MakeTChain(filelists[f]+".txt");
		JetProducer* jet_prod = new JetProducer(chj);
		
		TChain* chp = MakeTChain(filelists[f]+".txt");
		PhotonProducer* pho_prod = new PhotonProducer(chp);
	
		ReducedBase* base = jet_prod->GetBase();
		int nEvts = base->fChain->GetEntries();

		pho_prod->SetMinPt(30.); //50 for iso bkg (same for jets + photons), 30 for nominal
		if(isoBkgSel) pho_prod->SetMinPt(50.); //50 for iso bkg (same for jets + photons), 30 for nominal
		pho_prod->SetMinNrhs(15);
                pho_prod->SetMinEmE(30); 
                pho_prod->SetMinRhE(0.5);
	       	pho_prod->SetTransferFactor(1/30.);
		pho_prod->SetIsoCut();	

		jet_prod->SetMinPt(50.); //50 for iso bkg (same for jets + photons), 30 for nominal
		jet_prod->SetMinNrhs(15);
                jet_prod->SetMinEmE(10); 
                jet_prod->SetMinRhE(0.5); 
	       	jet_prod->SetTransferFactor(1/10.);	
		

		SampleWeight swts;
		swts.Init();
		double scale, xsec;
		swts.GetWeights(filelists[f],scale,xsec);
cout << "scale " << scale << " xsec " << xsec << endl;
		int nSelEvts_jet = 0;
		int nSelEvts_pho = 0;
		vector<Jet> jets, phos;
		cout << " nEvts " << nEvts << endl;
		//get total number of selected events for weighting
		for(int i = 0; i < nEvts; i++){
		        base->GetEntry(i);
		        cout << "\33[2K\r"<< "evt: " << i << " of " << nEvts << flush;
			jet_prod->GetTrueJets(jets, i);
		
			pho_prod->GetTruePhotons(phos, i);
			//add in event level selection (ie min ht, etc.)
			//do iso bkg evt selection in photons to compare data/MC
                	if(isoBkgSel){
		        	if(jets.size() >= 1){ nSelEvts_jet++; }
                        	//L1 seed
                        	if(!base->Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS) continue;
                        	//cout << "passed L1 seed" << endl;     
                        	//L1 to HLT Regional EGM matching leg
                        	if(!base->Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter) continue;
                        	//cout << "passed L1 to HLT" << endl;   
                        	//photon pt > 60
                        	if(!base->Trigger_hltEG60EtFilter) continue;
                        	//cout << "passed photon pt > 60" << endl;      
                        	//jet ht > 175 && jet pt > 10 && |eta jet| < 3
                        	if(!base->Trigger_hltHT175Jet10) continue;
                        	//cout << "passed HT > 175" << endl;    
                        	//jet ht > 350 && jet pt > 15 && |eta jet| < 3
                        	if(!base->Trigger_hltPFHT350Jet15) continue;
                        	//cout << "passed HT > 135" << endl;    

                        	//min photon multiplicity
                        	if(phos.size() < 1) continue;
                        	//cout << "passed pho mult" << endl;    
                        	//min jet multiplicity
                        	if(jets.size() < 1) continue;
                        	//cout << "passed jet mult" << endl;
				//ht - scalar sum
                        	double ht = 0;
                        	for(auto j : jets) ht += j.pt();
                        	//dphi bw photon and jet systems (vector sum of objects)
                        	Jet jet_sys = VectorSum(jets);
                        	Jet pho_sys = VectorSum(phos);
                        	double pho_phi = pho_sys.phi_02pi();
                        	//cout << "pho system E " << pho_sys.E() << " phi " << pho_sys.phi_02pi() << " jet system E " << jet_sys.E() << " phi " << jet_sys.phi_02pi() << endl;
                        	double jet_phi = jet_sys.phi_02pi();
                        	double dphi_phojet = pho_phi - jet_phi;
                        	dphi_phojet = acos(cos(dphi_phojet)); //wraparound - will always be < pi
                        	//MET upper limit - orthogonal to signal MET selection
                        	if(base->Met_pt > maxmet) continue;
                        	//cout << "passed max met" << endl;
                        	//min jet ht
                        	if(ht < minht) continue;
                        	//cout << "passed min ht" << endl;
                        	//az angle bw hardest presel photon + jet system
                        	//cout << "dphi " << dphi_phojet << " met " << base->Met_pt << endl;
                        	if(dphi_phojet < pi-0.3) continue; //want dphi ~ phi - implies less MET in event
                        	//cout << "passed dphi " << endl;
                        //trigger req - take baseline, photon pt leg + jet ht legs from HLT Photon60 R9Id90 CaloIdL IsoL DisplacedIdL PFHT350MinPFJet15
				nSelEvts_pho++;
			}
			else{
		        	if(jets.size() >= 1){ nSelEvts_jet++; }
		        	if(phos.size() >= 1){ nSelEvts_pho++; }
			}
			jets.clear();
			phos.clear();
		}
		cout << endl;
		//divide by number of selected events
		//initial set to total events in ctor
		double jet_weight;
		double pho_weight; 
       		if(nSelEvts_jet == 0) jet_weight = 0;
		else jet_weight = (scale * (xsec * 1000) * base->Evt_genWgt) / (double)nSelEvts_jet;
		if(nSelEvts_pho == 0) pho_weight = 0;
		else pho_weight = (scale * (xsec * 1000) * base->Evt_genWgt) / (double)nSelEvts_pho;
		
cout << "scale " << scale << " xsec " << xsec << " genwt " << base->Evt_genWgt << " n sel evts pho " << nSelEvts_pho << " n sel evts jet " << nSelEvts_jet << " out of " << nEvts << " total events" << endl;

		ofile << filelists[f] << "	" << jet_weight << "	" << pho_weight << endl;

	}
	cout << "Wrote event weights to " << ofilename << endl;
	ofile.close();
	return 1;
};
