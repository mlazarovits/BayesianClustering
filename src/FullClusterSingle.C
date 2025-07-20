#include "JetProducer.hh"
#include "PhotonProducer.hh"
#include "JetSimProducer.hh"
#include "BayesCluster.hh"
#include "FullViz3D.hh"
#include "VarClusterViz3D.hh"
#include "BasicDetectorSim.hh"

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
	string fname = "";
	string oname = "";
	string in_file = "rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root";
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
	double alpha = 0.1;

	bool viz = false;
	int verb = 0;
	bool weighted = true;
	bool smeared = true;
	bool timesmeared = false;
	//by default in BayesCluster
	bool distconst = true;
	//clustering strategy for skimmer
	int strat = 0; //0 is NlnN
	int evt = 0; //for event evt
	int obj = 0; //object to cluster (0 : jets, 1 : photons)
	int nobj = 0; //number of object in event to cluster (only for photons)
	//this should be in N/GeV
	//at least at 1 GeV but 1 GeV rh shouldn’t be able to seed a cluster so 1 GeV should be a fraction of entries
	double gev = 1./30.;
	double minRhE = 0.5;
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
		if(strncmp(argv[i],"--viz", 5) == 0){
    	 		viz = true;
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
		if(strncmp(argv[i],"-a", 2) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--alpha", 7) == 0){
			i++;
    	 		alpha = std::stod(argv[i]);
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
   		}
		if(strncmp(argv[i],"--noSmear", 9) == 0){
    	 		smeared = false;
   		}
		if(strncmp(argv[i],"--timeSmear", 11) == 0){
    	 		timesmeared = true;
   		}
		if(strncmp(argv[i],"--noDist", 8) == 0){
    	 		distconst = false;
   		}
		if(strncmp(argv[i],"--strategy", 10) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-s", 2) == 0){
    	 		i++;
			strat = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--evt", 5) == 0){
    	 		i++;
			evt = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"-e", 2) == 0){
    	 		i++;
			evt = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--object", 8) == 0){
    	 		i++;
			obj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--nobj", 6) == 0){
    	 		i++;
			nobj = std::atoi(argv[i]);
   		}
		if(strncmp(argv[i],"--gev", 5) == 0){
			i++;
    	 		gev = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minRhE", 8) == 0){
			i++;
    	 		minRhE = std::stod(argv[i]);
   		}
		

	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                           print options" << endl;
   		cout << "   --input(-i) [file]                   input root file" << endl;
   		cout << "   --output(-o) [file]                  output root file (in plots/)" << endl;
   		cout << "   --alpha(-a) [a]                      sets concentration parameter alpha for DPM in BHC (default = 0.1)" << endl;
   		cout << "   --EMalpha(-EMa) [a]                  sets concentration parameter alpha for variational EM GMM (default = 0.5)" << endl;
   		cout << "   --thresh(-t) [t]                     sets threshold for cluster cutoff" << endl;
   		cout << "   --beta0 [beta0]                      set scale parameter on covariance for prior on mu (N(mu | m0, (beta0*Lambda)^-1) (default = 0.001)" << endl;
   		cout << "   --m0 [m0_eta] [m0_phi] [m0_time]     set mean parameter for prior on mu (N(mu | m0, (beta0*Lambda)^-1) (default = [0,0,0])" << endl;
   		cout << "   --W0diag [W0_ee] [W0_pp] [W0_tt]     set *diagonal elements* in covariance parameter for prior on lambda (InverseWishart(Lambda | W0, nu0) (default = [1/3,1/3,1/3])" << endl;
   		cout << "   --nu0 [nu0]                          set dof parameter for prior on lambda (InverseWishart(Lambda | W0, nu0) (default = 3 = dim)" << endl;
   		cout << "   --verbosity(-v) [verb]               set verbosity (default = 0)" << endl;
   		cout << "   --strategy(-s) [strat]               set clustering strategy for skimmer (0 : NlnN (default), 1 : N^2, 2 : GMM only, 0+1 do not apply to photons)" << endl;
   		cout << "   --object [obj]                       set object to cluster (0 : jets, default; 1 : photons, 2 : detector sim)" << endl;
   		cout << "   --nobj [n]                           set number of object in event to cluster (photons only, default = 0)" << endl;
   		cout << "   --evt(-e) [evt]                      get plots for event e (default = 0)" << endl;
   		cout << "   --gev [gev]                          set energy weight transfer factor in N/GeV (default = 1/30 GeV)" << endl;
   		cout << "   --minRhE [minRhe]                    set minimum ECAL rechit energy (default = 0.5 GeV)" << endl;
   		cout << "   --viz                                makes plots (and gifs if N == 3)" << endl;
   		cout << "   --noSmear                            turns off smearing data according to preset covariance (default = true)" << endl;
   		cout << "   --timeSmear                          turns on time smearing data according to preset covariance (default = false)" << endl;
   		cout << "   --noWeight                           turns off weighting data points (default = true)" << endl;
   		cout << "   --noDist                             turns off - clusters must be within pi/2 in phi (default = true)" << endl;
   		cout << "Example: ./FullClusterSingle.x -a 0.5 -t 1.6 --viz" << endl;

   		return 0;
  	}
	cout << "Free sha-va-ca-doo!" << endl;
	/////GET DATA FROM NTUPLE//////
	string cmslab = in_file.substr(0,in_file.find(".root"));
	cmslab = cmslab.substr(cmslab.find("/")+1);
	if(cmslab.find("cmseos") != string::npos)
		cmslab = "";
	
	map<string, Matrix> prior_params;
	prior_params["scale"] = scale;
	prior_params["dof"] = dof;
	prior_params["scalemat"] = W;
	prior_params["mean"] = m;

	if(obj == 0 || obj == 2){
		if(obj == 0) fname += "jets";
		else fname += "jetsSim"; 
		if(strat == 0)
			fname += "NlnN";
		else if(strat == 1)
			fname += "N2";
		else if(strat == 2)
			fname += "GMMonly";
		else{
			cout << "Strategy number " << strat << " not supported. Only 0 : NlnN, 1 : N^2, 2 : GMM only." << endl;
			return -1;
		}
	}
	else if(obj == 1)
		fname += "photon";
	else{
		cout << "Object number " << obj << " not supported. Only 0 : jets, 1 : photons." << endl;
		return -1;
	}




	if(oname != "")
		fname = "skims/"+fname+"Skim_"+oname;
	else fname = "skims/"+fname+"Skim";
	string a_string;
	std::stringstream stream;
	stream << alpha;
	a_string = stream.str();
	int idx = a_string.find(".");
	if(idx != -1)
		a_string.replace(idx,1,"p");	
	if(strat != 2) fname += "_bhcAlpha-"+a_string;

	string ema_string;
	stream.str("");
	stream << emAlpha;
	ema_string = stream.str();
	idx = ema_string.find(".");
	if(idx != -1)
		ema_string.replace(idx,1,"p");	
	fname += "_emAlpha-"+ema_string+"_";

		
	string scale_string;
	stream.str("");
	stream << scale.at(0,0);
	scale_string = stream.str();
	idx = scale_string.find(".");
	if(idx != -1)
		scale_string.replace(idx,1,"p");	
	fname += "beta0-"+scale_string+"_";
	
	string mean_string;
	stream.str("");
	for(int i = 0; i < m.GetDims()[0]; i++){
		stream << m.at(i,0);
		if(i < m.GetDims()[0]-1) stream << "-";
	}
	mean_string = stream.str();
	idx = 0;
	while(idx != string::npos){
		idx = mean_string.find(".");
		if(idx == -1) break;
		mean_string.replace(idx,1,"p");	
		idx = mean_string.find(".");
	}
	fname += "m0-"+mean_string+"_";
	
	string W_string;
	stream.str("");
	for(int i = 0; i < W.GetDims()[0]; i++){
		stream << W.at(i,i);
		if(i < W.GetDims()[0]-1) stream << "-";
	}
	W_string = stream.str();
	idx = 0;
	while(idx != string::npos){
		idx = W_string.find(".");
		if(idx == -1) break;
		W_string.replace(idx,1,"p");	
		idx = W_string.find(".");
	}
	//do replacement of . to p
	fname += "W0diag-"+W_string+"_";
	
	string dof_string;
	stream.str("");
	stream << dof.at(0,0);
	dof_string = stream.str();
	idx = dof_string.find(".");
	if(idx != -1)
		dof_string.replace(idx,1,"p");	
	fname += "nu0-"+dof_string+"_";


	string t_string;
	stream.str("");
	stream << thresh;
	t_string = stream.str();
	idx = t_string.find(".");
	if(idx != -1)
		t_string.replace(idx,1,"p");
	fname += "thresh"+t_string+"_";

	
	string gev_string;
	stream.str("");
	stream << gev;
	gev_string = stream.str();
	idx = gev_string.find(".");
	if(idx != -1)
		gev_string.replace(idx,1,"p");	
	fname += "NperGeV"+gev_string+"_";
	
	fname += cmslab; //long sample name
	fname += "evt"+to_string(evt);

	string root_fname = fname;
	fname = "plots/"+fname;

	cout << "Prior Parameters" << endl;
	cout << "EMalpha0 " << emAlpha << endl;
	cout << "Energy transfer factor: " << gev << endl;
	cout << "beta0" << endl;
	scale.Print();
	cout << "mean0" << endl;
	m.Print();
	cout << "nu0" << endl;
	dof.Print();
	cout << "W0" << endl;
	W.Print(); 
	cout << "fname " << fname << endl;
	if(viz){
		cout << "Writing viz output to directory: " << fname << endl;
		if(gSystem->AccessPathName(fname.c_str())){
			gSystem->Exec(("mkdir -p "+fname).c_str());
		}
	}
	
	TFile* file = TFile::Open(in_file.c_str());
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = 2*acos(-1)/360.; //1 degree in radians
	double deta = dphi; //-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(0.,2,2); //no smear in time	
	double tres_c = 0.2; //ns
	//double tres_n = 30*sqrt(tres_c - tres_c*tres_c); //v1//ns*E (set s.t. 30 GeV gives sig_t = 1 ns)
	double tres_n = 30*sqrt(1 - tres_c*tres_c); //v1//ns*E (set s.t. 30 GeV gives sig_t = 1 ns)

	vector<Jet> rhs, jets, phos, jetsAK8, jetsAK15;
	vector<node*> trees;
	//get rhs (as Jets) for event
	if(obj == 0){
		JetProducer prod(file);
		prod.SetTransferFactor(gev);
		prod.SetMinPt(30);
		prod.SetMinNrhs(15);
		prod.SetMinEmE(30);
		//calibrate
		prod.SetTimeCalibration(true);
		if(strat != 2){
			cout << "Getting rec hits for jets at event " << evt << endl;
			prod.GetRecHits(rhs,evt);	
			cout << rhs.size() << " rechits in event " << evt << endl;
		}
		else{
			cout << "Getting rec hits for jet " << nobj << " at event " << evt << endl;
			prod.GetTrueJets(jets, evt, gev);
			if(jets.size() < 1){ cout << "No jets passing selection found for event " << evt << endl; return -1; }
			if(nobj > jets.size() - 1){ cout << "Only " << jets.size() << " jets passing selection found for event " << evt << endl; return -1; }
			jets[nobj].GetJets(rhs);
			cout << rhs.size() << " rechits in jet " << nobj << " in event " << evt << endl;
		}
	}

	else if(obj == 1){
        	PhotonProducer prod(file);
		prod.SetTransferFactor(gev);
		int npho = nobj; //which photon to analyze
		if(viz) cout << "Making plots for photon " << npho << " at event " << evt << endl;
        	//prod.GetRecHits(rhs,evt,npho);
        	prod.SetIsoCut();
		//calibrate
		prod.SetTimeCalibration(true);
		prod.GetTruePhotons(phos, evt);
		if(phos.size() < 1){ cout << "No photons passing selection found for event " << evt << endl; return -1; }
		if(nobj > phos.size() - 1){ cout << "Only " << phos.size() << " photons passing selection found for event " << evt << endl; return -1; }
		phos[npho].GetJets(rhs);
		cout << rhs.size() << " rechits in photon " << npho << " in event " << evt << endl;
	}
	//sim jets
	else if(obj == 2){
		JetSimProducer prod(file);
		//JetSimProducer prod(in_file);
		prod.SetTransferFactor(gev);
		prod.SetRecoMinPt(50);
		prod.SetRecoMinE(100);
		prod.SetMinRhE(minRhE); 
		//no need to calibrate
		prod.PrintPreselection();
		if(strat != 2){
		      cout << "Getting all rec hits for event " << evt << endl;
		      prod.GetRecHits(rhs,evt);	
		      cout << rhs.size() << " rechits in event " << evt << endl;
		}
		else{
			cout << "Getting rec hits for jet " << nobj << " at event " << evt << endl;
			prod.GetRecoJets(jets, jetsAK8, jetsAK15, evt);
			if(jets.size() < 1){ cout << "No jets passing selection found for event " << evt << endl; return -1; }
			if(nobj > jets.size() - 1){ cout << "Only " << jets.size() << " jets passing selection found for event " << evt << endl; return -1; }
			jets[nobj].GetJets(rhs);
			root_fname += ".root";
			cout << "Writing eta phi map to " << root_fname << endl;
			prod.EtaPhiMap(root_fname, rhs);
			cout << rhs.size() << " rechits in jet " << nobj << " in event " << evt << endl;
		}
	}
	
        if(rhs.size() < 1) return -1;
	vector<double> ws;
	rhs[0].GetWeights(ws);

	//to debug - use less rechits
	//int nrhs = 100;
	//rhs.resize(nrhs);
	//rhs.shrink_to_fit();


	
	BayesCluster *algo = new BayesCluster(rhs);
	if(smeared) algo->SetDataSmear(smear);
	//set time resolution smear: c^2 + n^2/e^2
	//remember time is already in ns
	//e = w/gev
	algo->SetThresh(thresh);
	algo->SetAlpha(alpha);
	algo->SetSubclusterAlpha(emAlpha);
	algo->SetVerbosity(verb);
	algo->SetPriorParameters(prior_params);
	clock_t t = clock();
	//jets
	if(obj == 0 || obj == 2){
		cout << "Using clustering strategy " << strat << ": ";
		if(strat == 0 || strat == 1){
			if(strat == 0){
				cout << "Delauney (NlnN)" << endl;
				//to track computation time from beginning of program
				trees = algo->NlnNCluster();
			}
			else if(strat == 1){
				cout << "N^2 (naive)" << endl;
				//to track computation time from beginning of program
				if(obj == 0) trees = algo->N2Cluster();
			}
			//cout << trees.size() << " true trees" << endl;
			if(viz){
				//plotting stuff here
				FullViz3D plots = FullViz3D(trees);
				plots.SetVerbosity(verb);
				plots.SetTransfFactor(gev);
				//add info of true jets
				plots.AddTrueJets(jets);
				plots.Write(fname);
			}
		}
		else if(strat == 2){
			cout << "mixture model with pre-clustered AK4 jets" << endl;
			string oname = "";
			if(viz) oname = fname;
			GaussianMixture* gmm = algo->SubCluster(oname);
			int nk = gmm->GetNClusters();
			cout << nk << " clusters in model." << endl;
			//check pis sum to 1
			double pisum = 0;
			map<string, Matrix> params;
			vector<double> norms;
			gmm->GetNorms(norms);
			for(int k = 0; k < nk; k++){
				params = gmm->GetLHPosteriorParameters(k);
				cout << "cluster k " << k << " pi " << params["pi"].at(0,0) << " energy " << norms[k]/gev << " center" << endl;
				params["mean"].Print();
				pisum += params["pi"].at(0,0);
			}
			cout << "pisum " << pisum << endl;
			

		}
		else{ cout << " undefined. Please use --strategy(-s) [strat] with options 0 (NlnN) or 1 (N^2)." << endl; return -1; }
	}
	//photons
	else if(obj == 1){
		string oname = "";
		if(viz) oname = fname;
		GaussianMixture* gmm = algo->SubCluster(oname);
		int nk = gmm->GetNClusters();
		cout << nk << " clusters in model." << endl;
		map<string, Matrix> params;
		for(int k = 0; k < nk; k++){
			params = gmm->GetLHPosteriorParameters(k);
			cout << "cluster " << k << " has time center " << params["mean"].at(2,0) << " with mixing coeff " << params["pi"].at(0,0) << endl;
		}
	}	



	else{ cout << "Object undefined. Please use --object [obj] with options 0 for jets or 1 for photons." << endl; return -1; }


	//causing crashes - //file->Close();
	//delete file;
	

	t = clock() - t;
	cout << "Total program time elapsed: " << (float)t/CLOCKS_PER_SEC << " seconds." << endl;
}
