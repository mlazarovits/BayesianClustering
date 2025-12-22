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
	recoAK15 = 5,
	genAK15 = 6,
	bhc_noPU = 7
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

void MakeHistsSim2(string file = "", string proc = ""){ 
// TODO -make also want to have this script do hist formatting so it spits out pretty hists (not just ugly ones)
// see macros/HistFormatSim.C for formatting code
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}
	string pt_thresh = "175";
	vector<string> obs = {"Energy","Mass","Pt","EtaCenter","PhiCenter","TimeCenter","EtaVar","PhiVar","TimeVar"};

	vector<string> genparts;
	SelType sel = singleW; 
	if(file.find("singleW") != string::npos){
		sel = singleW;
	}
	else if(file.find("boostTop") != string::npos){
		sel = boostTop;
		pt_thresh = "250"; //higher pthatmin
	}
	else if(file.find("QCDdijets") != string::npos){
		sel = QCDdijets;
	}
	BuildMaps(sel);
	
	string oname = file.substr(0,file.find(".root"));
	oname = oname+"_paperPlots.root";
	TFile* ofile = new TFile(oname.c_str(),"RECREATE");
	ofile->cd();

	vector<ROOT::RDF::RResultPtr<TH1D>> hists1d;
	vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;

	for(int jt = 0; jt < 8; jt++){//jt < 1; jt++){//jettypes.size(); jt++){
		string jetname = jetnames.at(JetType(jt))+"Jet";
		//skip gen jets
		if(jetname.find("gen") != string::npos)
			continue;
		//create unrolled dataframe
		TTreeInterface TI(file,"tree");

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
		string csvname = jetname+"Unrolled.csv";
		TI.CreateFlattenedCSV(branchlist, subbranchlist, csvname);
		ROOT::RDataFrame df0 = ROOT::RDF::FromCSV(csvname.c_str(),true,' ',-1LL,std::move(coltypes));

		//define sigma columns
		auto df = df0.Define(jetname+"_EtaSig","sqrt("+jetname+"_EtaVar)").Define(jetname+"_PhiSig","sqrt("+jetname+"_PhiVar)").Define(jetname+"_TimeSig","sqrt("+jetname+"_TimeVar)");

		//auto cols = df.GetColumnNames();
		//for(auto s : cols) cout << "col " << s << " col type " << df.GetColumnType(s) << endl;
		//df.Display({"evtidx","jetidx","subclidx","BHCJet_Energy", "W_Energy", "W_EtaCenter"},15,15)->Print();

		ROOT::RDataFrame df_njets("tree",file,{jetname+"_Pt"});
		//remove intermediate csv
		removefile(csvname);
		

		//df.Display({"evtidx","jetidx","subclidx","BHCJet_subclusterTimeVar"},10,10)->Print();
		string leadcut = jetname+"_Pt > "+pt_thresh;
		string notleadcut = jetname+"_Pt <= "+pt_thresh;
		string nocut = jetname+"_Pt >= 0";
		map<string, string> ptcuts;
		ptcuts["lead"] = leadcut;
		ptcuts["notlead"] = notleadcut;
		ptcuts[""] = nocut;
		//take only 1 jet row (first subcluster - all subcl rows have same jet info)
		string jetcut =  "subclidx == 0";
		for(auto pit = ptcuts.begin(); pit != ptcuts.end(); pit++){
			string ptcut = pit->second;
			string ptjetcut = ptcut + " && " + jetcut;
			string ptname = pit->first;
			if(ptname != ""){
				ptname = "_"+ptname;	
				//njets hists
				auto njets = df_njets.Define(jetname+"_nJets"+ptname,jetname+"_Pt["+ptcut+"].size()")					.Histo1D({(jetname+"_nJets"+ptname).c_str(), (jetname+"_nJets"+ptname).c_str(), 10,0,10},jetname+"_nJets"+ptname);
				hists1d.push_back(njets);
			}
			else{
				auto njets = df_njets.Histo1D({(jetname+"_nJets"+ptname).c_str(), (jetname+"_nJets"+ptname).c_str(), 10,0,10},jetname+"_nJets"+ptname);
				hists1d.push_back(njets);
	
			}
			//write jet 1D histograms
			for(int b = 0; b < branchlist.size(); b++){
				string obs = branchlist[b].substr(branchlist[b].find("_")+1);
				//skip covariance branches
				if(obs.find("Var") != string::npos || obs.find("Cov") != string::npos) continue;
				auto hist1d = df.Filter(ptcut)
					.Histo1D({(branchlist[b]+ptname).c_str(),(branchlist[b]+ptname).c_str(),25,bins.at(obs).first, bins.at(obs).second},branchlist[b]);
				hists1d.push_back(hist1d);
			}
			auto h_jet_etasig = df.Filter(ptcut)
				.Histo1D({(jetname+"_EtaSig"+ptname).c_str(),(jetname+"_EtaSig"+ptname).c_str(),25,0,2.5},jetname+"_EtaSig");
			hists1d.push_back(h_jet_etasig);
			auto h_jet_phisig = df.Filter(ptcut)
				.Histo1D({(jetname+"_PhiSig"+ptname).c_str(),(jetname+"_PhiSig"+ptname).c_str(),25,0,1.2},jetname+"_PhiSig");
			hists1d.push_back(h_jet_phisig);
			auto h_jet_timesig = df.Filter(ptcut)
				.Histo1D({(jetname+"_TimeSig"+ptname).c_str(),(jetname+"_TimeSig"+ptname).c_str(),25,0, 10},jetname+"_TimeSig");
			hists1d.push_back(h_jet_timesig);

cout << "making plot " << jetname+"_nSubclusters"+ptname << endl;
			auto h_jet_nsubclusters = df.Filter(ptcut)
				.Histo1D({(jetname+"_nSubclusters"+ptname).c_str(),(jetname+"_nSubclusters"+ptname).c_str(),10,0,10},jetname+"_nSubclustersJet");
			hists1d.push_back(h_jet_nsubclusters);

			//n subclusters hists
			string highmassjetcut = jetname+"_Mass > 100 && "+jetcut+" && "+ptcut;
			string lowmassjetcut = jetname+"_Mass < 50 && "+jetcut+" && "+ptcut;
			string Wmassjetcut = jetname+"_Mass > 70 && "+jetname+"_Mass < 90 && "+jetcut+" && "+ptcut;
			string topmassjetcut = jetname+"_Mass > 165 && "+jetname+"_Mass < 185 && "+jetcut+" && "+ptcut;
			if(sel == singleW){
				auto nsubcls_highmass = df.Filter(highmassjetcut)
				   		.Histo1D({(jetname+"_nSubclusters"+ptname+"_highMass").c_str(), (jetname+"_nSubclusters"+ptname+"_highMass").c_str(), 10,0,10},jetname+"_nSubclustersJet");
				hists1d.push_back(nsubcls_highmass);
				auto nsubcls_lowmass = df.Filter(lowmassjetcut)
				   		.Histo1D({(jetname+"_nSubclusters"+ptname+"_lowMass").c_str(), (jetname+"_nSubclusters"+ptname+"_lowMass").c_str(), 10,0,10},jetname+"_nSubclustersJet");
				hists1d.push_back(nsubcls_lowmass);
				auto nsubcls_Wmass = df.Filter(Wmassjetcut)
				   		.Histo1D({(jetname+"_nSubclusters"+ptname+"_WMass").c_str(), (jetname+"_nSubclusters"+ptname+"_WMass").c_str(), 10,0,10},jetname+"_nSubclustersJet");
				hists1d.push_back(nsubcls_Wmass);
			}

			//jet size vs jet mass
			auto hist2d_jetsize_jetmass = df.Filter(ptcut)
					.Histo2D({(jetname+"_Mass_JetSize"+ptname).c_str(), (jetname+"_Mass_JetSize"+ptname+";mass;jet size;a.u.").c_str(), 50,0,250, 50,0,2}, jetname+"_Mass",jetname+"_JetSize"); 
			hists2d.push_back(hist2d_jetsize_jetmass);
			//jet pt vs jet mass
			auto hist2d_jetpt_jetmass = df.Filter(ptcut)
					.Histo2D({(jetname+"_Mass_JetPt"+ptname).c_str(), (jetname+"_Mass_JetPt"+ptname+";mass;pt;a.u.").c_str(), 50,0,250,50,0,std::stod(pt_thresh)+50}, jetname+"_Mass",jetname+"_Pt"); 
			hists2d.push_back(hist2d_jetpt_jetmass);
			
			if(jetname.find("BHC") == string::npos)
				continue;

			//PU cleaning histograms
			string matched_cut;
			string matched_name;
			if(sel == singleW){
				matched_cut = " W_Energy != -1 && "+jetname+"_nSubclustersJet > 1";
				matched_name = "_WMatched";
			}
			if(sel == boostTop){
				matched_cut = " Top_Energy != -1 && "+jetname+"_nSubclustersJet > 1";
				//matched_cut = " "+topmassjetcut+" && "+jetname+"_nSubclustersJet > 1";
				matched_name = "_TopMatched";
			}
			string relEtaVar = jetname+"_subclusterEtaVar / "+jetname+"_EtaVar";
			string relPhiVar = jetname+"_subclusterPhiVar / "+jetname+"_PhiVar";
			string relTimeVar = jetname+"_subclusterTimeVar / "+jetname+"_TimeVar";
			string relGeoAvgFuncStr = "pow( "+relEtaVar+" * "+relPhiVar+" * "+relTimeVar+", 1./3.)";
			auto df_relvar = df.Define("subclRelEnergy",jetname+"_subclusterEnergy / "+jetname+"_Energy").Define("subclRelGeoAvgVar",relGeoAvgFuncStr)
					.Define("subclRelEtaVar",relEtaVar).Define("subclRelPhiVar",relPhiVar).Define("subclRelTimeVar",relTimeVar);
			auto hist2d_pu_cleaning = df_relvar.Filter(ptcut+" && "+matched_cut)
					.Histo2D({(jetname+"_subclRelGeoAvgVar_subclRelEnergy"+ptname+matched_name+"_ge2Subcls").c_str(),(jetname+"_subclRelGeoAvgVar_subclRelEnergy"+ptname+matched_name+"_ge2Subcls;subclRelGeoAvgVar;subclRelEnergy;a.u.").c_str(),50,0,5.,50,0,1.1},"subclRelGeoAvgVar","subclRelEnergy");
				//subcluster relative energy vs geo avg of relative variances for jets matched to relevant gen particles (ie W for single W, top for boostTop, etc) 
			hists2d.push_back(hist2d_pu_cleaning);

			//1D PU plots
			string relEcut_ge = "subclRelEnergy >= 0.5";
			auto h_reletavar_relEge0p5 = df_relvar.Filter(ptcut+" && "+matched_cut+" && "+relEcut_ge)
						.Histo1D({(jetname+"_subclRelEtaVar"+ptname+matched_name+"_ge2Subcls_subclRelEge0p5").c_str(),(jetname+"_subclRelEtaVar"+ptname+matched_name+"_ge2Subcls_subclRelEge0p5;subclRelEtaVar;a.u.").c_str(),50,0,5},"subclRelEtaVar");
			hists1d.push_back(h_reletavar_relEge0p5);
			auto h_relphivar_relEge0p5 = df_relvar.Filter(ptcut+" && "+matched_cut+" && "+relEcut_ge)
						.Histo1D({(jetname+"_subclRelPhiVar"+ptname+matched_name+"_ge2Subcls_subclRelEge0p5").c_str(),(jetname+"_subclRelPhiVar"+ptname+matched_name+"_ge2Subcls_subclRelEge0p5;subclRelPhiVar;a.u.").c_str(),50,0,5},"subclRelPhiVar");
			hists1d.push_back(h_relphivar_relEge0p5);
			auto h_reltimevar_relEge0p5 = df_relvar.Filter(ptcut+" && "+matched_cut+" && "+relEcut_ge)
						.Histo1D({(jetname+"_subclRelTimeVar"+ptname+matched_name+"_ge2Subcls_subclRelEge0p5").c_str(),(jetname+"_subclRelTimeVar"+ptname+matched_name+"_ge2Subcls_subclRelEge0p5;subclRelTimeVar;a.u.").c_str(),50,0,5},"subclRelTimeVar");
			hists1d.push_back(h_reltimevar_relEge0p5);

			string relEcut_lt = "subclRelEnergy < 0.5";
			auto h_reletavar_relElt0p5 = df_relvar.Filter(ptcut+" && "+matched_cut+" && "+relEcut_lt)
						.Histo1D({(jetname+"_subclRelEtaVar"+ptname+matched_name+"_ge2Subcls_subclRelElt0p5").c_str(),(jetname+"_subclRelEtaVar"+ptname+matched_name+"_ge2Subcls_subclRelElt0p5;subclRelEtaVar;a.u.").c_str(),50,0,5},"subclRelEtaVar");
			hists1d.push_back(h_reletavar_relElt0p5);
			auto h_relphivar_relElt0p5 = df_relvar.Filter(ptcut+" && "+matched_cut+" && "+relEcut_lt)
						.Histo1D({(jetname+"_subclRelPhiVar"+ptname+matched_name+"_ge2Subcls_subclRelElt0p5").c_str(),(jetname+"_subclRelPhiVar"+ptname+matched_name+"_ge2Subcls_subclRelElt0p5;subclRelPhiVar;a.u.").c_str(),50,0,5},"subclRelPhiVar");
			hists1d.push_back(h_relphivar_relElt0p5);
			auto h_reltimevar_relElt0p5 = df_relvar.Filter(ptcut+" && "+matched_cut+" && "+relEcut_lt)
						.Histo1D({(jetname+"_subclRelTimeVar"+ptname+matched_name+"_ge2Subcls_subclRelElt0p5").c_str(),(jetname+"_subclRelTimeVar"+ptname+matched_name+"_ge2Subcls_subclRelElt0p5;subclRelTimeVar;a.u.").c_str(),50,0,5},"subclRelTimeVar");
			hists1d.push_back(h_reltimevar_relElt0p5);
	
		}
	}

	ofile->cd();
	for(auto hist : hists1d) hist->Write();
	for(auto hist : hists2d) hist->Write();
cout << "Wrote all histograms to " << ofile->GetName() << endl;
	ofile->Close();	


}
