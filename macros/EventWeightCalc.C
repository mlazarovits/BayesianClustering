#include <iostream>
#include <fstream>
#include "../include/BaseProducer.hh"
#include "../include/Jet.hh"
#include "../include/JetProducer.hh"
#include "../include/SampleWeight.hh"

using std::cout;
using std::endl;

//configurable parameters
//gev

int EventWeightCalc(string selection){
	double gev_jet = 0.1;
	double gev_pho = 1./30.;

	//just for MCs - data weight = 1
	vector<string> files;
	files.push_back("GJets_R17_"+selection+"_v24_GJets_HT-40To100_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GJets_R17_"+selection+"_v24_GJets_HT-100To200_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GJets_R17_"+selection+"_v24_GJets_HT-200To400_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GJets_R17_"+selection+"_v24_GJets_HT-400To600_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GJets_R17_"+selection+"_v24_GJets_HT-600ToInf_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GMSB_R17_"+selection+"_v24_GMSB_L-150TeV_Ctau-0_1cm_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GMSB_R17_"+selection+"_v24_GMSB_L-150TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GMSB_R17_"+selection+"_v24_GMSB_L-350TeV_Ctau-0_1cm_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GMSB_R17_"+selection+"_v24_GMSB_L-350TeV_Ctau-1000cm_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("GMSB_R17_"+selection+"_v24_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT1000to1500_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT100to200_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT1500to2000_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT2000toInf_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT200to300_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT300to500_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT500to700_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT50to100_AODSIM_RunIIFall17DRPremix.root");
	files.push_back("QCD_R17_"+selection+"_v24_QCD_HT700to1000_AODSIM_RunIIFall17DRPremix.root");
	ofstream ofile;
	string ofilename = "info/EventWeights"+selection".txt";
	ofile.open(ofilename);
	//dont write just for remembering what's being written
	//ofile << "file	jet_weight	pho_weight" << endl;
	for(int f = 0; f < files.size(); f++){
		cout << "File " << files[f];
		//loop over all files in a list
		TFile* file = TFile::Open(("root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSNtuples/"+files[f]).c_str());
		BaseProducer* prod = new JetProducer(file);
		ReducedBase* base = prod->GetBase();
		int nEvts = base->fChain->GetEntries();

		SampleWeight swts;
		swts.Init();
		double scale, xsec;
		swts.GetWeights(file,scale,xsec);

		int nSelEvts_jet = 0;
		int nSelEvts_pho = 0;
		vector<Jet> jets, phos;
		cout << " nEvts " << nEvts << endl;
		//get total number of selected events for weighting
		for(int i = 0; i < nEvts; i++){
		        base->GetEntry(i);
		        cout << "\33[2K\r"<< "evt: " << i << " of " << nEvts << flush;
		        prod->GetTrueJets(jets, i, gev_jet);
			prod->GetTruePhotons(phos, i, gev_pho);
		        if(jets.size() >= 1){ nSelEvts_jet++; }
		        if(phos.size() >= 1){ nSelEvts_pho++; }
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
		


		ofile << files[f] << "	" << jet_weight << "	" << pho_weight << endl;

	}
	ofile.close();
	return 1;
};
