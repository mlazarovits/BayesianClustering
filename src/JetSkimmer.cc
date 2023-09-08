#include "JetSkimmer.hh"

#include "Clusterizer.hh"
#include "Matrix.hh"
#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>

JetSkimmer::JetSkimmer(){ };



JetSkimmer::~JetSkimmer(){ 
	_file->Close();
	delete _base;
	delete _file;
}


JetSkimmer::JetSkimmer(TFile* file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	_prod = new JetProducer(_file);
	_base = _prod->GetBase();
		//	_base = new ReducedBase(tree);
	_nEvts = _base->fChain->GetEntries();

	hists1D.push_back(nClusters);
	hists1D.push_back(nTrueJets);
	

}



//make cluster param histograms
void JetSkimmer::Skim(){
	TFile* ofile = new TFile("plots/jet_skims_v6.root","RECREATE");
	
	

	//create data smear matrix - smear in eta/phi
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	
	
	Clusterizer* algo = new Clusterizer();
	algo->SetAlpha(0.1);
	algo->SetThresh(1.);
	algo->SetMaxNClusters(5);
	algo->SetWeighted(true);
	algo->SetVerbosity(0);
	algo->SetDataSmear(smear);


	GaussianMixture* gmm = new GaussianMixture();
	
	map<string, Matrix> params;

	vector<JetPoint> rhs;
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		//find subclusters for each photon
		_prod->GetRecHits(rhs, i);
		cout << "evt: " << i << " of " << _nEvts << " nrhs: " << rhs.size() << "\r" << flush;

		//need to cluster - full algo	
		gmm = algo->FindSubjets(Jet(rhs));
		FillTotalHists(gmm);
		
		//jet specific hists
		nTrueJets->Fill((double)_base->Jet_energy->size());	


	}

}





