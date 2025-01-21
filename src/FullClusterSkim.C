#include "JetSkimmer.hh"
#include "PhotonSkimmer.hh"
#include "SuperClusterSkimmer.hh"
#include "BayesCluster.hh"
#include "FullViz3D.hh"
//#include "VarClusterViz3D.hh"
#include <TSystem.h>
#include <TFile.h>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <regex>

using std::string;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
	string oname = "";
	string in_file = "";//rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";
	bool hprint = false;
	double thresh = 1.;
	
	//prior parameters
	double emAlpha = 0.5;
	Matrix scale = Matrix(1e-3);
	Matrix dof = Matrix(3);
	Matrix W(3,3);
	W.InitIdentity();
	W.mult(W,1./3.);
	Matrix m(3,1);

	int verb = 0;
	bool weighted = true;
	bool smear = true;
	bool timesmear = false;
	//by default in BayesCluster
	bool distconst = true;
	int evti = 0; //for skimming from evti to evtj
	int evtj = 0;
	int obj = 0; //object to cluster (0 : jets, 1 : photons)
	//this should be in N/GeV
	double gev = -999;
	//put cuts on jets (ie min pt) here
	double minpt = 30;
	double minnrhs = 15;
	double minEmE = 20;
	double minRhE = 0.5;
	double maxRhE = -999;
	bool frac = false;
	bool calib = true;
	bool spikes = false;
	int skip = 1;
	int bh = 1;
	bool iso = true;

	double minpt_isobkg = 50;
	double minht_isobkg = 50;
	double minjetpt_isobkg = 50;
	double maxmet_isobkg = 150;
	bool isobkg = false; 
	for(int i = 0; i < argc; i++){
		if(strncmp(argv[i],"--help", 6) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"-h", 2) == 0){
    	 		hprint = true;
   		}
		if(strncmp(argv[i],"--input", 7) == 0){
     			i++;
    	 		in_file = string(argv[i]);
   		}
		if(strncmp(argv[i],"-i", 2) == 0){
     			i++;
    	 		in_file = string(argv[i]);
   		}
		if(strncmp(argv[i],"--output", 8) == 0){
     			i++;
    	 		oname = string(argv[i]);
   		}
		if(strncmp(argv[i],"-o", 2) == 0){
     			i++;
    	 		oname = string(argv[i]);
   		}
		if(strncmp(argv[i],"--verbosity", 11) == 0){
    	 		i++;
			verb = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-v", 2) == 0){
    	 		i++;
			verb = std::atoi(argv[i]);
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
		if(strncmp(argv[i],"--beta0", 7) == 0){
			i++;
    	 		double beta0 = std::stod(argv[i]);
			scale = Matrix(beta0);
   		}
		if(strncmp(argv[i],"--nu0", 5) == 0){
			i++;
    	 		double nu0 = std::stod(argv[i]);
			dof = Matrix(nu0);
   		}
		if(strncmp(argv[i],"--W0diag", 8) == 0){
			i++;
    	 		double W_ee = std::stod(argv[i]);
			i++;
    	 		double W_pp = std::stod(argv[i]);
			i++;
    	 		double W_tt = std::stod(argv[i]);
			//no covariance bw dimensions
			W = Matrix(3,3);
			W.SetEntry(W_ee,0,0);
			W.SetEntry(W_pp,1,1);
			W.SetEntry(W_tt,2,2);
   		}
		if(strncmp(argv[i],"--m0", 4) == 0){
			i++;
    	 		double m_e = std::stod(argv[i]);
			i++;
    	 		double m_p = std::stod(argv[i]);
			i++;
    	 		double m_t = std::stod(argv[i]);
			m = Matrix(3,1);
			m.SetEntry(m_e,0,0);
			m.SetEntry(m_p,1,0);
			m.SetEntry(m_t,2,0);
   		}
		if(strncmp(argv[i],"--noWeight", 10) == 0){
    	 		weighted = false;
			cout << "Turning off energy weighting." << endl;
   		}
		if(strncmp(argv[i],"--noSmear", 9) == 0){
    	 		smear = false;
			cout << "Turning off smearing." << endl;
   		}
		if(strncmp(argv[i],"--timeSmear", 11) == 0){
    	 		timesmear = true;
			cout << "Turning on time smearing." << endl;
   		}
		if(strncmp(argv[i],"--noDist", 8) == 0){
    	 		distconst = false;
			cout << "Turning off distance constraint." << endl;
   		}
		if(strncmp(argv[i],"--evtFirst", 6) == 0){
    	 		i++;
			evti = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--evtLast", 6) == 0){
    	 		i++;
			evtj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--object", 8) == 0){
    	 		i++;
			obj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--gev", 5) == 0){
			i++;
    	 		gev = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minpt", 7) == 0){
			i++;
    	 		minpt = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minNrhs", 9) == 0){
			i++;
    	 		minnrhs = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minemE", 8) == 0){
			i++;
    	 		minEmE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minRhE", 8) == 0){
			i++;
    	 		minRhE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--maxRhE", 8) == 0){
			i++;
    	 		maxRhE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--applyFrac", 11) == 0){
			frac = true;
			cout << "Apply fractions for rec hits." << endl;
   		}
		if(strncmp(argv[i],"--noCalibrate", 13) == 0){
			calib = false;
			cout << "Not applying time calibration to rhs" << endl;
   		}
		if(strncmp(argv[i],"--rejectSpikes", 14) == 0){
			spikes = true;
			cout << "Rejecting spikes with swiss cross cut" << endl;
   		}
		if(strncmp(argv[i],"--BHFilter", 10) == 0){
			i++;
    	 		bh = std::stoi(argv[i]);
   		}
		if(strncmp(argv[i],"--skip", 6) == 0){
			i++;
    	 		skip = std::stoi(argv[i]);
   		}
		if(strncmp(argv[i],"--noIso", 7) == 0){
			iso = false;
			cout << "Not applying isolation to photons" << endl;
   		}
		if(strncmp(argv[i],"--isoBkg", 8) == 0){
			isobkg = true;
			cout << "Selecting photons with isolated background cuts" << endl;
   		}
		if(strncmp(argv[i],"--minphopt_isobkg", 17) == 0){
			i++;
    	 		minpt_isobkg = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minht_isobkg", 14) == 0){
			i++;
    	 		minht_isobkg = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minjetpt_isobkg", 17) == 0){
			i++;
    	 		minjetpt_isobkg = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--maxmet_isobkg", 15) == 0){
			i++;
    	 		maxmet_isobkg = std::stod(argv[i]);
   		}

	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                           print options" << endl;
   		cout << "   --input(-i) [file]                   input root file" << endl;
   		cout << "   --output(-o) [file]                  output root file" << endl;
   		cout << "   --thresh(-t) [t]                     set threshold for cluster cutoff" << endl;
   		cout << "   --EMalpha(-EMa) [a]                  set concentration parameter alpha for variational EM GMM (default = 0.5)" << endl;
   		cout << "   --beta0 [beta0]                      set scale parameter on covariance for prior on mu (N(mu | m0, (beta0*Lambda)^-1) (default = 0.001)" << endl;
   		cout << "   --m0 [m0_eta] [m0_phi] [m0_time]     set mean parameter for prior on mu (N(mu | m0, (beta0*Lambda)^-1) (default = [0,0,0])" << endl;
   		cout << "   --W0diag [W0_ee] [W0_pp] [W0_tt]     set *diagonal elements* in covariance parameter for prior on lambda (InverseWishart(Lambda | W0, nu0) (default = [1/3,1/3,1/3])" << endl;
   		cout << "   --nu0 [nu0]                          set dof parameter for prior on lambda (InverseWishart(Lambda | W0, nu0) (default = 3 = dim)" << endl;
   		cout << "   --verbosity(-v) [verb]               set verbosity (default = 0)" << endl;
   		cout << "   --object [obj]                       set object to cluster (0 : jets, default; 1 : superclusters, 2 : photons)" << endl;
   		cout << "   --gev [gev]                          set energy weight transfer factor in N/GeV (default = 1/30 GeV)" << endl;
   		cout << "   --minpt [minpt]                      set minimum pt (default = 30 GeV)" << endl;
   		cout << "   --minNrhs [minnrhs]                  set minimum # of rhs (default = 2)" << endl;
   		cout << "   --minemE [mineme]                    set minimum ECAL energy (default = 10 GeV)" << endl;
   		cout << "   --minRhE [minRhe]                    set minimum ECAL rechit energy (default = 0.5 GeV)" << endl;
   		cout << "   --maxRhE [maxRhe]                    set maximum ECAL rechit energy (default = -999, off)" << endl;
   		cout << "   --BHFilter [bh]                      set how beam halo filter is applied (0 : not applied, 1 : applied (default), 2 : inversely applied)" << endl;
   		cout << "   --skip [skip]                        set skip for event loop (default = 1)" << endl;
   		cout << "   --evtFirst [i] --evtLast [j]         skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
   		cout << "   --noSmear                            turns off smearing data (default = true, on)" << endl;
   		cout << "   --timeSmear                          turns on time smearing data (default = false, off)" << endl;
   		cout << "   --noWeight                           turns off weighting data points (default = false, on)" << endl;
   		cout << "   --noDist                             turns off distance constraint: clusters must be within pi/2 in phi (default = false, on)" << endl;
   		cout << "   --applyFrac                          applying fractions for rec hits PHOTONS ONLY (default = false, off)" << endl;
   		cout << "   --noCalibrate                        turn off channel-by-channel calibration for rh time (default = false, on)" << endl;
   		cout << "   --rejectSpikes                       reject spikes based on swiss cross cut (default = false, off)" << endl;
   		cout << "   --noIso                              turn off isolation in preselection (photons only, default = true, on)" << endl;
   		cout << "   --isoBkg                             apply isolated background selection (photons only, default = false, off)" << endl;
   		cout << "   --minphopt_isobkg [minpt]               set mininum photon pt for iso bkg selection (photons only, default = 70)" << endl;
   		cout << "   --minht_isobkg [minht]               set minimum jet ht for iso bkg selection (photons only, default = 50)" << endl;
   		cout << "   --minjetpt_isobkg [minjetpt]         set minimum jet pt for iso bkg selection (photons only, default = 50)" << endl;
   		cout << "   --maxmet_isobkg [maxmet]             set maximum met for iso bkg selection (photons only, default = 150)" << endl;
   		cout << "Example: ./FullClusterSkim.x -a 0.5 -t 1.6 -o test" << endl;

   		return 0;
  	}


	cout << "Free sha-va-ca-doo!" << endl;
	if(gev == -999){
		if(obj == 0) gev = 1./10.;
		else gev = 1./30.;
	}

	map<string, Matrix> prior_params;
	prior_params["scale"] = scale;
	prior_params["dof"] = dof;
	prior_params["scalemat"] = W;
	prior_params["mean"] = m;


	string cmslab, version;	
	TFile* file = nullptr;
	/////GET DATA FROM NTUPLE//////
	if(!in_file.empty()){
		//get version
		std::smatch m;
		std::regex re("_v[0-9]+_");
		version = "";
		std::regex_search(in_file,m,re);
		for(auto x : m) version += x;
		cmslab = in_file.substr(in_file.find(version)+version.size(),in_file.find(".root") - in_file.find(version)-version.size());//"GMSB_L-350TeV_Ctau-200cm_2017_v9";
		version.pop_back();
		cmslab += version;
		file = TFile::Open(in_file.c_str());
	}

	if(gSystem->AccessPathName(in_file.c_str())){
		cout << "Error: file " << in_file << " not found." << endl;
		return -1;
	}	

	string fname;
	if(obj == 0)
		fname = "jets";
	else if(obj == 1)
		fname = "supercluster";
	else if(obj == 2)
		fname = "photon";
	else{
		cout << "Object number " << obj << " not supported. Only 0 : jets, 1 : superclusters, 2 : photons." << endl;
		return -1;
	}

	//condor is in outname
	if(oname.find("condor") != string::npos){
		fname = oname;	
	}
	else{
		if(oname != "")
			fname = "skims/"+fname+"Skim_"+oname;
		else fname = "skims/"+fname+"Skim";
		if(evti != evtj) fname += "_evt"+std::to_string(evti)+"to"+std::to_string(evtj);
		
		std::stringstream stream;
		int idx;

		string ema_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3) << emAlpha;
		ema_string = stream.str();
		idx = ema_string.find(".");
		ema_string.replace(idx,1,"p");	
		fname += "_emAlpha"+ema_string+"_";
		
		string scale_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3) << scale.at(0,0);
		scale_string = stream.str();
		idx = scale_string.find(".");
		scale_string.replace(idx,1,"p");	
		fname += "beta0-"+scale_string+"_";
		
		string mean_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3);
		for(int i = 0; i < m.GetDims()[0]; i++){
			stream << m.at(i,0);
			if(i < m.GetDims()[0]-1) stream << "-";
		}
		mean_string = stream.str();
		idx = 0;
		while(idx != string::npos){
			idx = mean_string.find(".");
			mean_string.replace(idx,1,"p");	
			idx = mean_string.find(".");
		}
		fname += "m0-"+mean_string+"_";
		
		string W_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3);
		for(int i = 0; i < W.GetDims()[0]; i++){
			stream << W.at(i,i);
			if(i < W.GetDims()[0]-1) stream << "-";
		}
		W_string = stream.str();
		idx = 0;
		while(idx != string::npos){
			idx = W_string.find(".");
			W_string.replace(idx,1,"p");	
			idx = W_string.find(".");
		}
		//do replacement of . to p
		fname += "W0diag-"+W_string+"_";
		
		string dof_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3) << dof.at(0,0);
		dof_string = stream.str();
		idx = dof_string.find(".");
		dof_string.replace(idx,1,"p");	
		fname += "nu0-"+dof_string+"_";

		string t_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3) << thresh;
		t_string = stream.str();
		idx = t_string.find(".");
		t_string.replace(idx,1,"p");	
		fname += "thresh"+t_string+"_";
		
		string gev_string;
		stream.str("");
		stream << std::fixed << std::setprecision(3) << gev;
		gev_string = stream.str();
		idx = gev_string.find(".");
		gev_string.replace(idx,1,"p");	
		fname += "NperGeV"+gev_string+"_";

		
		fname += cmslab; //long sample name
		//fname += cmslab.substr(0,cmslab.find("_")); //short sample name
		//fname += version; //"_v9"

	}
	fname = fname+".root";

	cout << "transfer factor (gev) N/Energy " << gev << endl;	
	cout << "Prior Parameters" << endl;
	cout << "beta0" << endl;
	scale.Print();
	cout << "mean0" << endl;
	m.Print();
	cout << "nu0" << endl;
	dof.Print();
	cout << "W0" << endl;
	W.Print(); 
	
cout << "fname " << fname << endl;
	//make sure evti < evtj
	if(evti > evtj){
		int evt = evtj;
		evtj = evti;
		evti = evt;
	}
	if(evti != evtj) cout << "Skimming events " << evti << " to " << evtj << " for ";
	else cout << "Skimming all events for ";
	//choose time calibration file
	string calibfile = "";
        if(fname.find("GJets") != string::npos)
                calibfile = "info/KUCMS_GJets_R17_v16_rhE5_mo_Cali.root";
        else if(fname.find("JetHT") != string::npos)
                calibfile = "info/KUCMS_JetHT_R17_v18_rhE5_Cali.root";
        else if(fname.find("DEG") != string::npos || fname.find("DoubleEG") != string::npos || fname.find("EGamma") != string::npos)
                calibfile = "info/KUCMS_DoubleEG_R17_v18_rhE5_Cali.root";
        else if(fname.find("QCD") != string::npos)
                calibfile = "info/KUCMS_QCD_R17_v16_rhE5_mo_Cali.root";
        //else default to GJets
        else
                calibfile = "info/KUCMS_GJets_R17_v16_rhE5_mo_Cali.root";
	if(obj == 0){
		cout << "jets" << endl;
		JetSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		skimmer.SetMinPt(minpt);
		skimmer.SetMinNrhs(minnrhs);
		skimmer.SetMinRhE(0.5);
		skimmer.SetMinRhE_PV(minRhE);
		if(maxRhE != -999) skimmer.SetMaxRhE(maxRhE);
		skimmer.SetMinEmE(minEmE);
		skimmer.SetSpikeRejection(spikes); //if true, reject spikes
		skimmer.SetSkip(skip);
		if(calib) skimmer.SetTimeCalibrationMap(calibfile);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
		//set EMalpha, thresh
		skimmer.SetThresh(thresh);
		skimmer.SetPriorParameters(prior_params);
		skimmer.SetEMAlpha(emAlpha);
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetSmear(smear);
		skimmer.SetTimeSmear(timesmear); 
		skimmer.SetBeamHaloFilter(bh);
		//do only mm/true jet pv times
		skimmer.Skim();
	}
	else if(obj == 1){
		cout << "superclusters" << endl;
		SuperClusterSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		bool data;
		if(in_file.find("SIM") != string::npos)
			data = false;
        	else
			data = true;
		if(calib) skimmer.SetTimeCalibrationMap(calibfile);
		if(iso){
			cout << "Applying isolation preselection for training labels and object selection." << endl;
			skimmer.SetIsoCuts();
		}
		else cout << "Not applying isolation preselection for training labels and object selection." << endl;
		skimmer.SetMinPt(minpt);
		skimmer.SetMinRhE(minRhE);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
        	skimmer.ApplyFractions(frac);
		//skimmer.SetDebug(debug);
		//set EMalpha, thresh
		skimmer.SetThresh(thresh);
		skimmer.SetEMAlpha(emAlpha);
		skimmer.SetPriorParameters(prior_params);
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetSmear(smear);
		skimmer.SetTimeSmear(timesmear); 
		skimmer.SetBeamHaloFilter(bh);
        	skimmer.Skim();
	}
	else if(obj == 2){
		cout << "photons" << endl;
		PhotonSkimmer skimmer(file);
		skimmer.SetCMSLabel(cmslab);
		bool data;
		if(in_file.find("SIM") != string::npos)
			data = false;
        	else
			data = true;
		if(calib) skimmer.SetTimeCalibrationMap(calibfile);
		if(iso) skimmer.SetIsoCuts();
		skimmer.SetMinRhE(minRhE);
		skimmer.SetOutfile(fname);
		skimmer.SetTransferFactor(gev);
        	skimmer.ApplyFractions(frac);
		//skimmer.SetDebug(debug);
		//set EMalpha, thresh
		skimmer.SetThresh(thresh);
		skimmer.SetEMAlpha(emAlpha);
		skimmer.SetPriorParameters(prior_params);
		skimmer.SetEventRange(evti,evtj);
		skimmer.SetSmear(smear);
		skimmer.SetTimeSmear(timesmear);
		skimmer.SetBeamHaloFilter(bh);
		
		skimmer.SetMinPt_IsoBkg(minpt_isobkg);
		skimmer.SetMinHt_IsoBkg(minht_isobkg);
		skimmer.SetMinJetPt_IsoBkg(minjetpt_isobkg);
		skimmer.SetMaxMet_IsoBkg(maxmet_isobkg);
		skimmer.SetIsoBkgSel(isobkg);

        	skimmer.Skim();
	}
	if(calib) cout << "Using calibration file " << calibfile << endl;
	else cout << "No timing calibration applied" << endl;
        return 0;

}
