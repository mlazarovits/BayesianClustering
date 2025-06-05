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
	double gev = 1./10.;
	bool smear = false;
	double tres_cte = 0.1727;//0.133913;
	double tres_stoch = 0.5109;//1.60666;
	double tres_noise = 2.106;//0.00691415;

	//set clustering strategy
	//0 = NlnN
	//1 = N2
	//2 = gmm only
	int strat = 0;
	string infile = "";
	string oname = "";
	double thresh = 1.;
	//prior parameters
	double emAlpha = 0.5;
	double alpha = 0.1;
	Matrix scale = Matrix(1e-3);
	Matrix dof = Matrix(3);
	Matrix W(3,3);
	W.InitIdentity();
	W.mult(W,1./3.);
	Matrix m(3,1);
	
	double minpt = 30.;
	double minE = 30.;
	int minnrhs = 1;
	double minRhE = 0.5;
	double minNconsts = 5;
	double mintoppt = 0;
	double minwpt = 0;
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
		if(strncmp(argv[i],"--minE", 6) == 0){
			i++;
    	 		minE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minNrhs", 9) == 0){
			i++;
    	 		minnrhs = std::stoi(argv[i]);
   		}
		if(strncmp(argv[i],"--minRhE", 8) == 0){
			i++;
    	 		minRhE = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minNconsts", 12) == 0){
			i++;
    	 		minNconsts = std::stoi(argv[i]);
   		}
		if(strncmp(argv[i],"--smear", 7) == 0){
    	 		smear = true;
			cout << "Turning on smearing (no measurement error)." << endl;
   		}
		if(strncmp(argv[i],"--tResCte", 9) == 0){
			i++;
    	 		tres_cte = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--tResStoch", 11) == 0){
			i++;
    	 		tres_stoch = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--tResNoise", 11) == 0){
			i++;
    	 		tres_noise = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minTopPt", 10) == 0){
			i++;
    	 		mintoppt = std::stod(argv[i]);
   		}
		if(strncmp(argv[i],"--minWPt", 8) == 0){
			i++;
    	 		minwpt = std::stod(argv[i]);
   		}



	}
	if(hprint){
		cout << "Usage: " << argv[0] << " [options]" << endl;
   		cout << "  options:" << endl;
   		cout << "   --help(-h)                           print options" << endl;
   		cout << "   --input(-i) [file]                   input root file" << endl;
   		cout << "   --output(-o) [file]                  output root file" << endl;
   		cout << "   --strategy(-s) [strat]               sets clustering strategy (0 = NlnN, default; 1 = N2, 2 = gmm only, 3 = NlnN with reco AK4 rhs)" << endl;
		cout << "   --alpha(-a) [a]                      sets concentration parameter alpha for DPM in BHC (default = 0.1)" << endl;
   		cout << "   --EMalpha(-EMa) [a]                  sets concentration parameter alpha for variational EM GMM (default = 0.5)" << endl;
   		cout << "   --beta0 [beta0]                      set scale parameter on covariance for prior on mu (N(mu | m0, (beta0*Lambda)^-1) (default = 0.001)" << endl;
   		cout << "   --m0 [m0_eta] [m0_phi] [m0_time]     set mean parameter for prior on mu (N(mu | m0, (beta0*Lambda)^-1) (default = [0,0,0])" << endl;
   		cout << "   --W0diag [W0_ee] [W0_pp] [W0_tt]     set *diagonal elements* in covariance parameter for prior on lambda (InverseWishart(Lambda | W0, nu0) (default = [1/3,1/3,1/3])" << endl;
   		cout << "   --nu0 [nu0]                          set dof parameter for prior on lambda (InverseWishart(Lambda | W0, nu0) (default = 3 = dim)" << endl;
   		cout << "   --thresh(-t) [t]                     sets threshold for cluster cutoff" << endl;
   		cout << "   --verbosity(-v) [verb]               set verbosity (default = 0)" << endl;
   		cout << "   --gev [gev]                          set energy weight transfer factor in N/GeV (default = 1/10 GeV)" << endl;
   		cout << "   --minpt [minpt]                      set gen minimum pt (default = 30 GeV)" << endl;
   		cout << "   --minE [minE]                        set gen minimum E (default = 30 GeV)" << endl;
   		cout << "   --minTopPt [mintoppt]                set gen top minimum pt (default = 0 GeV)" << endl;
   		cout << "   --minWPt [minwpt]                    set gen W minimum pt (default = 0 GeV)" << endl;
   		cout << "   --minNrhs [minnrhs]                  set minimum # of rhs (default = 2)" << endl;
   		cout << "   --minRhE [minRhE]                    set minimum rechit energy (default = 0.5 GeV)" << endl;
   		cout << "   --minNconsts [minNconsts]            set minimum number of constituents for gen jets (default = 5)" << endl;
   		cout << "   --smear                              smear cov (spatial only, turns off meas error)" << endl;
   		cout << "   --tResCte [t]                        set time smearing constant parameter in ns (default = 0.1727 ns)" << endl;
   		cout << "   --tResNoise [t]                      set time smearing noise (n*n/(e*e)) parameter in ns (default = 2.106 ns)" << endl;
   		cout << "   --tResStoch [t]                      set time smearing stochastic (s*s/e) parameter in ns (default = 0.5109 ns)" << endl;
   		cout << "   --evtFirst [i] --evtLast [j]         skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
   		cout << "Example: ./detectorSimSkimmer.x -i rootfiles/simNtuples_ttbar.root -a 0.5 -t 1.6" << endl;
		return 0;	
	}

	if(gSystem->AccessPathName(infile.c_str())){
		cout << "Error: file " << infile << " not found." << endl;
		return -1;
	}

	if(evti != evtj) cout << "Skimming events " << evti << " to " << evtj << endl;
	else cout << "Skimming all events" << endl;

	map<string, Matrix> prior_params;
	prior_params["scale"] = scale;
	prior_params["dof"] = dof;
	prior_params["scalemat"] = W;
	prior_params["mean"] = m;
	
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
		if(oname.find("/") != string::npos)
			oname = oname.substr(oname.find("/")+1);	
		oname = "simSkims/simSkim_"+oname;

	}
	else{
		if(oname.find("condor") == string::npos){
			string oname_extra = file->GetName();
			string match = "simNtuples_";
			oname_extra = oname_extra.substr(oname_extra.find(match)+match.size(),oname_extra.find(".root")-(oname_extra.find(match)+match.size()))+".root";
			oname = "simSkim_"+oname+"_"+oname_extra;
			if(strat == 0) oname += "_NlnN";
			else if(strat == 1) oname += "_N2";
			else if(strat == 2) oname += "_gmmOnly";
			else oname += "_NlnNonAK4";
			string a_string;
			std::stringstream stream;
			stream << alpha;
			a_string = stream.str();
			int idx = a_string.find(".");
			if(idx != -1)
				a_string.replace(idx,1,"p");	
			if(strat != 2) oname += "_bhcAlpha-"+a_string;

			string ema_string;
			stream.str("");
			stream << emAlpha;
			ema_string = stream.str();
			idx = ema_string.find(".");
			if(idx != -1)
				ema_string.replace(idx,1,"p");	
			oname += "_emAlpha-"+ema_string+"_";

				
			string scale_string;
			stream.str("");
			stream << scale.at(0,0);
			scale_string = stream.str();
			idx = scale_string.find(".");
			if(idx != -1)
				scale_string.replace(idx,1,"p");	
			oname += "beta0-"+scale_string+"_";
			
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
			oname += "m0-"+mean_string+"_";
			
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
			oname += "W0diag-"+W_string+"_";
			
			string dof_string;
			stream.str("");
			stream << dof.at(0,0);
			dof_string = stream.str();
			idx = dof_string.find(".");
			if(idx != -1)
				dof_string.replace(idx,1,"p");	
			oname += "nu0-"+dof_string+"_";


			string t_string;
			stream.str("");
			stream << thresh;
			t_string = stream.str();
			idx = t_string.find(".");
			if(idx != -1)
				t_string.replace(idx,1,"p");
			oname += "thresh"+t_string+"_";

			
			string gev_string;
			stream.str("");
			stream << gev;
			gev_string = stream.str();
			idx = gev_string.find(".");
			if(idx != -1)
				gev_string.replace(idx,1,"p");	
			oname += "NperGeV"+gev_string;
		}
		if(oname.find("/") != string::npos)
			oname = oname.substr(0,oname.find("/")) + oname.substr(oname.find("/")+1);	
	}	
	
	oname += ".root";
	if(oname.find("condor") == string::npos) oname = "plots/"+oname;

	cout << "Prior Parameters" << endl;
	cout << "alpha0 " << alpha << " EMalpha0 " << emAlpha << endl;
	cout << "Energy transfer factor: " << gev << endl;
	cout << "beta0" << endl;
	scale.Print();
	cout << "mean0" << endl;
	m.Print();
	cout << "nu0" << endl;
	dof.Print();
	cout << "W0" << endl;
	W.Print();
	cout << "Using tres_cte = " << tres_cte << " ns, tres_stoch = " << tres_stoch << " ns and tres_noise = " << tres_noise << " ns" << endl;
 
	BHCJetSkimmer skimmer(file);
	skimmer.SetOutfile(oname);
	skimmer.SetMinRhE(minRhE);
	skimmer.SetMinNrhs(minnrhs);
	skimmer.SetGenMinE(minE);
	skimmer.SetGenMinPt(minpt);
	skimmer.SetRecoMinPt(0);
	skimmer.SetRecoMinE(0);
	skimmer.SetGenTopMinPt(mintoppt);
	skimmer.SetGenWMinPt(minwpt);
	skimmer.SetMinNGenConsts(minNconsts);
	skimmer.SetStrategy(strat);
	skimmer.SetVerbosity(verb);
	skimmer.SetTransferFactor(gev);
	skimmer.SetAlpha(alpha);
	skimmer.SetSubclusterAlpha(emAlpha);
	skimmer.SetPriorParameters(prior_params);
	skimmer.SetThreshold(thresh);
	skimmer.SetEventRange(evti,evtj);
	skimmer.SetSmear(smear);
	skimmer.SetMeasErrParams(acos(-1)/180, tres_cte, tres_stoch, tres_noise);

	skimmer.Skim();

}
