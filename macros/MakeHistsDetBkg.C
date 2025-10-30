#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RCsvDS.hxx>

using std::cout;
using std::endl;
using std::vector;
using std::string;

void MakeHistsDetBkg(string file = ""){
	file = "skims/condor_superclusters_defaultv7_beta0-1e-5_m0-0p0-0p0-0p0_W0diag-0p013-0p013-33p333_nu0-3_NperGeV-0p0333333_emAlpha-1e-5_MET_R17_AL1NpSC_nolumimask_v31_MET.root";
	if(gSystem->AccessPathName(file.c_str())){
		cout << "File " << file << " does not exist." << endl;
		return;
	}

	ROOT::RDataFrame df("tree",file);

	string CRdef = "MET < 100";

	auto sc_time_eta = df.Filter(CRdef)
		.Histo2D({"SC_TimeCenter_EtaCenter_METlt100","SC_TimeCenter_EtaCenter_METlt100;time [ns];eta;a.u.",50,-20,20,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");
	
	//auto sc_time_eta_bhFlag = df.Filter(CRdef+" && Flag_globalSuperTightHalo2016Filter == 1")
		//.Histo2D({"SC_TimeCenter_EtaCenter_METlt100_BHFilter","SC_TimeCenter_EtaCenter_METlt100_BHFilter;time [ns];eta;a.u.",50,-20,20,50,-1.6,1.6},"SC_TimeCenter","SC_EtaCenter");


	string ofilename = file;
	ofilename = ofilename.substr(0,ofilename.find(".root"));
	ofilename += "_hists.root";
	cout << "Writing hists to " << ofilename << endl;	
	TFile* ofile = new TFile(ofilename.c_str(),"RECREATE");
	ofile->cd();
	sc_time_eta->Write();


}
