#include <iostream>
#include <filesystem>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RCsvDS.hxx>
#include "TTreeInterface.h"

enum JetType{
	bhc = 0,
	recoAK4 = 1,
	genAK4 = 2, 
	recoAK8 = 3,
	genAK8 = 4,
	recoAK15 = 3,
	genAK15 = 5,
	bhc_noPU = 6
};

enum SelType{
	singleW = 0,
	boostTop = 1,
	QCDdijets = 2
};


map<string, pair<double,double>> bins;
map<JetType,string> jetnames;

unordered_map<string, char> coltypes;

void BuildMaps(SelType sel){
	bins["Mass"] = make_pair(0,250.);
	bins["Energy"] = make_pair(0,1000.);
	bins["Pt"] = make_pair(0,1000.);
	bins["EtaCenter"] = make_pair(-3.1,3.1);
	bins["PhiCenter"] = make_pair(0,8*atan(1));
	bins["TimeCenter"] = make_pair(-1,1);
	bins["JetSize"] = make_pair(0,2);

	jetnames[bhc] = "BHC";
	jetnames[bhc_noPU] = "BHCnoPU";
	jetnames[recoAK4] = "recoAK4";
	jetnames[recoAK8] = "recoAK8";
	jetnames[recoAK15] = "recoAK15";
	jetnames[genAK4] = "genAK4";
	jetnames[genAK8] = "genAK8";
	jetnames[genAK15] = "genAK15";

	if(sel == singleW){
		coltypes["W_Energy"] = 'D';
		coltypes["W_EtaCenter"] = 'D';
		coltypes["W_PhiCenter"] = 'D';
		coltypes["Wq_Energy"] = 'D';
		coltypes["Wq_EtaCenter"] = 'D';
		coltypes["Wq_PhiCenter"] = 'D';
	}
	else if(sel == boostTop){
		coltypes["Top_Energy"] = 'D';
		coltypes["Top_EtaCenter"] = 'D';
		coltypes["Top_PhiCenter"] = 'D';
		coltypes["parton_Energy"] = 'D';
		coltypes["parton_EtaCenter"] = 'D';
		coltypes["parton_PhiCenter"] = 'D';
		//coltypes["TopW_Energy"] = 'D';
		//coltypes["TopW_EtaCenter"] = 'D';
		//coltypes["TopW_PhiCenter"] = 'D';
		//coltypes["Topb_Energy"] = 'D';
		//coltypes["Topb_EtaCenter"] = 'D';
		//coltypes["Topb_PhiCenter"] = 'D';
		//coltypes["Wq_Energy"] = 'D';
		//coltypes["Wq_EtaCenter"] = 'D';
		//coltypes["Wq_PhiCenter"] = 'D';
	}
	else if(sel == QCDdijets){
		coltypes["q_Energy"] = 'D';
		coltypes["q_EtaCenter"] = 'D';
		coltypes["q_PhiCenter"] = 'D';
	}
	else { }
}
double dR(double eta1, double eta2, double phi1, double phi2){
	double dphi = (phi1-phi2);
	dphi = acos(cos(dphi));
	return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
}

void removefile(string filename){
    	if (std::filesystem::exists(filename)) {
    	    try {
    	        std::filesystem::remove(filename);
    	        std::cout << "File " << filename << " successfully removed." << std::endl;
    	    } catch (const std::filesystem::filesystem_error& e) {
    	        std::cerr << "Error removing file " << filename << ": " << e.what() << std::endl;
    	    }
    	} else {
    	    std::cout << "File " << filename << " does not exist." << std::endl;
    	}
}	

void MakeHistsSim(string file = "", string proc = ""){ 
// TODO -make also want to have this script do hist formatting so it spits out pretty hists (not just ugly ones)
// see macros/HistFormatSim.C for formatting code
	file = "simSkims/condorSim_singleW_defaultv9p13_ptHatMin200_PFCand_defaultv9_bhcAlpha-1e-300_emAlpha-1e-5_NperGeV0p250_beta0-1e-5_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NlnN_singleW.root";	
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	//TODO - turn on output info
	//TFile* f = TFile::Open(file.c_str(),"READ");
	TTreeInterface TI(file,"tree");
	//do for each jet type
	vector<int> jettypes = {1, 2, 3, 4, 5 ,6};

	string pt_thresh = "175";
	vector<string> obs = {"Energy","Mass","Pt","EtaCenter","PhiCenter","TimeCenter","EtaVar","PhiVar","TimeVar"};

	vector<string> genparts;
	SelType sel = singleW; 
	if(file.find("singleW") != string::npos){
		sel = singleW;
	}
	else if(file.find("boostTop") != string::npos){
		sel = boostTop;
		pt_thresh = "450"; //higher pthatmin
	}
	else if(file.find("QCDdijets") != string::npos){
		sel = QCDdijets;
	}
	BuildMaps(sel);
	

	TFile* ofile = new TFile("testHist.root","RECREATE");
	ofile->cd();

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;

	for(int jt = 0; jt < 1; jt++){//jettypes.size(); jt++){
		//create unrolled dataframe
		string jetname = jetnames.at(JetType(jt))+"Jet"; 

		vector<string> branchlist;
		for(int o = 0; o < obs.size(); o++){
			branchlist.push_back(jetname+"_"+obs[o]);
		}
		branchlist.push_back(jetname+"_JetSize");
		if(sel == boostTop){
			TI.MapIdx(jetname+"_TopMatchedIdx","genpart_Energy","Top_Energy");
			TI.MapIdx(jetname+"_TopMatchedIdx","genpart_EtaCenter","Top_EtaCenter");
			TI.MapIdx(jetname+"_TopMatchedIdx","genpart_PhiCenter","Top_PhiCenter");
		}
		else if(sel == QCDdijets){
			TI.MapIdx(jetname+"_qMatchedIdx","genpart_Energy","q_Energy");
			TI.MapIdx(jetname+"_qMatchedIdx","genpart_EtaCenter","q_EtaCenter");
			TI.MapIdx(jetname+"_qMatchedIdx","genpart_PhiCenter","q_PhiCenter");
		}
		else{	
			TI.MapIdx(jetname+"_WMatchedIdx","genpart_Energy","W_Energy");
			TI.MapIdx(jetname+"_WMatchedIdx","genpart_EtaCenter","W_EtaCenter");
			TI.MapIdx(jetname+"_WMatchedIdx","genpart_PhiCenter","W_PhiCenter");
		}
		
		//checked that these are matched properly
		//TI.MapIdx(jetname+"_WMatchedIdx","genpart_ID",jetname+"W_genpartID");
		vector<string> subbranchlist;
		//only do subcluster unrolling for bhc jets
		if(jt == bhc){
			TI.SetNSubBranch(jetname+"_nSubclustersJet");
			for(int o = 0; o < obs.size(); o++){
				subbranchlist.push_back(jetname+"_subcluster"+obs[o]);
			}
			if(sel == boostTop){
				//see which 'partons' (ie W's, b's, q's) are which
				TI.MapIdxToSubIdx(jetname+"Top_subclusterPartonMatchedIdx","genpart_ID","parton_GenPartID");
				TI.MapIdxToSubIdx(jetname+"Top_subclusterPartonMatchedIdx","genpart_Energy","parton_Energy");
				TI.MapIdxToSubIdx(jetname+"Top_subclusterPartonMatchedIdx","genpart_EtaCenter","parton_EtaCenter");
				TI.MapIdxToSubIdx(jetname+"Top_subclusterPartonMatchedIdx","genpart_PhiCenter","parton_PhiCenter");
			}
			else if(sel == singleW){	
				TI.MapIdxToSubIdx(jetname+"W_subclusterPartonMatchedIdx","genpart_ID","Wq_GenPartID");
				TI.MapIdxToSubIdx(jetname+"W_subclusterPartonMatchedIdx","genpart_Energy","Wq_Energy");
				TI.MapIdxToSubIdx(jetname+"W_subclusterPartonMatchedIdx","genpart_EtaCenter","Wq_EtaCenter");
				TI.MapIdxToSubIdx(jetname+"W_subclusterPartonMatchedIdx","genpart_PhiCenter","Wq_PhiCenter");
			}
			else{ }
		}
		string csvname = jetname+"Unrolled.csv";
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname);
		ROOT::RDataFrame df = ROOT::RDF::FromCSV(csvname.c_str(),true,' ',-1LL,std::move(coltypes));

		//auto cols = df.GetColumnNames();
		//for(auto s : cols) cout << "col " << s << " col type " << df.GetColumnType(s) << endl;
		//df.Display({"evtidx","jetidx","subclidx","BHCJet_Energy", "W_Energy", "W_EtaCenter"},15,15)->Print();

		ROOT::RDataFrame df_njets("tree",file,{jetname+"_Pt"});
		//remove intermediate csv
		removefile(csvname);
		

		//df.Display({"evtidx","jetidx","subclidx","BHCJet_subclusterTimeVar"},10,10)->Print();
		string leadcut = jetname+"_Pt > "+pt_thresh;
		//take only 1 jet row (first subcluster - all subcl rows have same jet info)
		string jetcut =  "subclidx == 0";
		string leadjetcut = leadcut + " && " + jetcut;

		//njets hists
		auto njets = df_njets.Define(jetname+"_nJets_lead",jetname+"_Pt["+jetname+"_Pt > "+pt_thresh+"].size()")
				.Histo1D({(jetname+"_nJets_lead").c_str(), (jetname+"_nJets_lead").c_str(), 10,0,10},jetname+"_nJets_lead");
		hists1d.push_back(njets);

		//n subclusters hists
		string highmassjetcut = jetname+"_Mass > 100 && "+jetcut+" && "+leadcut;
		string lowmassjetcut = jetname+"_Mass < 50 && "+jetcut+" && "+leadcut;
		string Wmassjetcut = jetname+"_Mass > 70 && "+jetname+"_Mass < 90 && "+jetcut+" && "+leadcut;
		if(sel == singleW){
			auto nsubcls_highmass = df.Filter(highmassjetcut)
			   		.Histo1D({(jetname+"_nSubclusters_lead_highMass").c_str(), (jetname+"_nSubclusters_lead_highMass").c_str(), 10,0,10},jetname+"_nSubclustersJet");
			hists1d.push_back(nsubcls_highmass);
			auto nsubcls_lowmass = df.Filter(lowmassjetcut)
			   		.Histo1D({(jetname+"_nSubclusters_lead_lowMass").c_str(), (jetname+"_nSubclusters_lead_lowMass").c_str(), 10,0,10},jetname+"_nSubclustersJet");
			hists1d.push_back(nsubcls_lowmass);
			auto nsubcls_Wmass = df.Filter(Wmassjetcut)
			   		.Histo1D({(jetname+"_nSubclusters_lead_WMass").c_str(), (jetname+"_nSubclusters_lead_WMass").c_str(), 10,0,10},jetname+"_nSubclustersJet");
			hists1d.push_back(nsubcls_Wmass);
		}

		//jet-gen particle matching hists
		if(sel == singleW){
			//string genmatch_ratio = jetname+"_Energy / W_Energy";
			//auto hist_jetGen_energyRatio = df.Filter(leadjetcut+" && W_Energy != -1")
			//			.Define(jetname+"MatchedToW_EnergyRatio",genmatch_ratio)
			//			.Histo1D({(jetname+"_genW_Eratio_lead").c_str(),(jetname+"_genW_Eratio_lead;genq_Eratio;a.u.").c_str(), 50,0,2},jetname+"MatchedToW_EnergyRatio");
			//hists1d.push_back(hist_jetGen_energyRatio);
			
			string genmatch_ratio = jetname+"_Energy / W_Energy";
			auto hist_jetGen_match = df.Filter(leadjetcut+" && W_Energy != -1")
						.Define(jetname+"MatchedToW_EnergyRatio",genmatch_ratio)
						.Define(jetname+"MatchedToW_deltaR", dR,{jetname+"_EtaCenter","W_EtaCenter",jetname+"_PhiCenter","W_PhiCenter"})
						.Histo2D({(jetname+"_genW_Eratio_dR_lead").c_str(),(jetname+"_genW_Eratio_dR_lead;genW_Eratio;genW_dR;a.u.").c_str(), 50,0,2, 50, 0, 1.4},jetname+"MatchedToW_EnergyRatio",jetname+"MatchedToW_deltaR");
			hists2d.push_back(hist_jetGen_match);

		}
		
		//subcluster-gen particle matching hists
		if(sel == singleW){
			//string genmatch_ratio = jetname+"_subclusterEnergy / Wq_Energy";
			//auto hist_jetGen_energyRatio = df.Filter(leadcut+" && Wq_Energy != -1")
			//			.Define(jetname+"MatchedToWq_EnergyRatio",genmatch_ratio)
			//			.Histo1D({(jetname+"_genq_subclusterEratio_lead").c_str(),(jetname+"_genq_subclusterEratio_lead;genq_subclusterEratio;a.u.").c_str(), 50,0,2},jetname+"MatchedToWq_EnergyRatio");
			//hists1d.push_back(hist_jetGen_energyRatio);
			
			string genmatch_ratio = jetname+"_subclusterEnergy / Wq_Energy";
			auto hist_jetGen_match = df.Filter(leadcut+" && Wq_Energy != -1")
						.Define(jetname+"MatchedToWq_subclusterEnergyRatio",genmatch_ratio)
						.Define(jetname+"MatchedToWq_subclusterdeltaR", dR,{jetname+"_subclusterEtaCenter","Wq_EtaCenter",jetname+"_subclusterPhiCenter","Wq_PhiCenter"})
						.Histo2D({(jetname+"_genWq_subclEratio_dR_lead").c_str(),(jetname+"_genWq_subclEratio_dR_lead;genq_subclEratio;genq_subcldR;a.u.").c_str(), 50,0,2, 50, 0, 1.4},jetname+"MatchedToWq_subclusterEnergyRatio",jetname+"MatchedToWq_subclusterdeltaR");
			hists2d.push_back(hist_jetGen_match);

		}

		//jet size vs jet mass for high pt jets
		auto hist2d_jetsize_jetmass = df.Filter(leadjetcut)
				.Histo2D({(jetname+"_Mass_JetSize_lead").c_str(), (jetname+"_Mass_JetSize_lead;mass;jet size;a.u.").c_str(), 50,0,250, 50,0,2}, jetname+"_Mass",jetname+"_JetSize"); 
		hists2d.push_back(hist2d_jetsize_jetmass);

		//PU cleaning histograms
		string wmatched_cut = " W_Energy != -1 && "+jetname+"_nSubclustersJet > 1";
		string relEtaVar = jetname+"_subclusterEtaVar / "+jetname+"_EtaVar";
		string relPhiVar = jetname+"_subclusterPhiVar / "+jetname+"_PhiVar";
		string relTimeVar = jetname+"_subclusterTimeVar / "+jetname+"_TimeVar";
		string relGeoAvgFuncStr = "pow( "+relEtaVar+" * "+relPhiVar+" * "+relTimeVar+", 1./3.)";
		auto hist2d_pu_cleaning = df.Filter(leadjetcut+" && "+wmatched_cut)
				.Define("subclRelEnergy",jetname+"_subclusterEnergy / "+jetname+"_Energy") //relative energy
				.Define("subclRelGeoAvgVar",relGeoAvgFuncStr) //relative geo avg var
				.Histo2D({(jetname+"_subclRelGeoAvgVar_subclRelEnergy_lead_WMatched_ge2Subcls").c_str(),(jetname+"_subclRelGeoAvgVar_subclRelEnergy_lead_WMatched_ge2Subcls;subclRelGeoAvgVar;subclRelEnergy:a.u.").c_str(),50,0,1.2,50,0,1.2},"subclRelGeoAvgVar","subclRelEnergy");
			//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
		hists2d.push_back(hist2d_pu_cleaning);
	
		//write jet 1D histograms
		for(int b = 0; b < branchlist.size(); b++){
			string obs = branchlist[b].substr(branchlist[b].find("_")+1);
			//skip covariance branches
			if(obs.find("Var") != string::npos || obs.find("Cov") != string::npos) continue;
			auto hist1d = df.Filter(leadjetcut)
				.Histo1D({(branchlist[b]+"_lead").c_str(),(branchlist[b]+"_lead").c_str(),50,bins.at(obs).first, bins.at(obs).second},branchlist[b]);
			hists1d.push_back(hist1d);
		}
	}

	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
cout << "Wrote all histograms to " << ofile->GetName() << endl;
	ofile->Close();	


}
