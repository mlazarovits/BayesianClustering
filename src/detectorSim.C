#include "BasicDetectorSim.hh"
#include "Jet.hh"
#include <string>
#include <iostream>
#include <TFile.h>
#include <Math/Vector4D.h>

using std::string;
using std::cout;
using std::endl;

using PtEtaPhiEVector = ROOT::Math::PtEtaPhiEVector;
using XYZTVector = ROOT::Math::XYZTVector;

int main(int argc, char *argv[]){

	vector<Jet> rhs;
	vector<PtEtaPhiEVector> genmoms;
	vector<PtEtaPhiEVector> recomoms;
	int nevts = 10;
	int evt = 0;
	int verb = 0;

	BasicDetectorSim det;
	det.SetNEvents(nevts);
	det.SetVerbosity(verb);
	det.SimTTbar();
	//default arg is all events
	det.SimulateEvents(evt);
	det.GetRecHits(rhs);
	cout << rhs.size() << " rechits in event " << evt << endl;
	
	

	string fname = "plots/detectorSimSkim.root";
	cout << "Writing to " << fname << endl;
	TFile f(fname.c_str(),"RECREATE");
	int nhists = 5;
	TH1D* hists[nhists];
	hists[0] = new TH1D("rhtimes","rhtimes",100,-10,10);
	hists[1] = new TH1D("rhE","rhE",100,0,500);
	hists[2] = new TH1D("rheta","rheta",50,-1.5,1.5);
	hists[3] = new TH1D("rhphi","rhphi",50,-0.1,6.4);
	//gen hists - need for "calibration" (validation of detector effects)
	//ratio between reco and gen energy
	hists[4] = new TH1D("recoE_genE","recoE_genE",100,0,5);

	//fill histograms
	for(int i = 0; i < rhs.size(); i++){
		hists[0]->Fill(rhs[i].t());
		hists[1]->Fill(rhs[i].E());
		hists[2]->Fill(rhs[i].eta());
		hists[3]->Fill(rhs[i].phi());

	} 

	f.cd();
	for(int i = 0; i < nhists; i++){
		hists[i]->Write();
	}
	f.Close();


}
