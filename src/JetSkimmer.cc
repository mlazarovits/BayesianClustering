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
	_file = file;
	TTree* tree = (TTree*)file->Get("tree/llpgtree");
	_base = new ReducedBase(tree);
	_nEvts = _base->fChain->GetEntries();
}

//make cluster param histograms
void JetSkimmer::Skim(){
	TFile* ofile = new TFile("plots/jet_skims_v6.root","RECREATE");
	

	
	//# of subclusters
	TH1I* nSubClusters = new TH1I("nSubClusters","nSubClusters",7,0,7.);
	//# of subclusters vs. photon reco energy
	TH2D* e_nSubClusters = new TH2D("e_nSubClusters","e_nSubClusters",50,0.,1000.,7,0.,7.);	


	

	int nPho;
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = -log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta,0,0);
	smear.SetEntry(dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	
	
	Clusterizer* algo = new Clusterizer();
	algo->SetAlpha(0.1);
	algo->SetThresh(1.);
	algo->SetMaxNClusters(5);
	algo->SetWeighted(true);
	algo->SetVerbosity(0);
//	algo->SetDataSmear(smear);


	GaussianMixture* gmm = new GaussianMixture();
	
	map<string, Matrix> params;

	vector<JetPoint> rhs;
	int nclusters;
	double theta, phi, r;
	vector<double> eigenvals, avg_Es;
	vector<Matrix> eigenvecs;
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		//find subclusters for each photon
		_prod->GetRecHits(rhs, i);
		cout << "evt: " << i << " of " << _nEvts << " nrhs: " << rhs.size() << "\r" << flush;

		//need to cluster - full algo	
		gmm = algo->FindSubjets(Jet(rhs));
		
		nclusters = gmm->GetNClusters();
		nSubClusters->Fill(nclusters);
		//e_nSubClusters->Fill(_base->Photon_energy->at(p), nclusters);
		
		gmm->GetAvgWeights(avg_Es);		


		for(int k = 0; k < nclusters; k++){
			params = gmm->GetParameters(k);
			eta_center->Fill(params["mean"].at(0,0));
			phi_center->Fill(params["mean"].at(1,0));
			time_center->Fill(params["mean"].at(2,0));
		
			//calculate slopes from eigenvectors
			params["cov"].eigenCalc(eigenvals, eigenvecs);
			
			//largest eigenvalue is last
			//phi/eta
			slope_space->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(0,0));
        		//eta/time
			slope_etaT->Fill(eigenvecs[2].at(0,0)/eigenvecs[2].at(2,0));
			//phi/time
			slope_phiT->Fill(eigenvecs[2].at(1,0)/eigenvecs[2].at(2,0));
			//polar angle
			//theta = arccos(z/r), r = sqrt(x2 + y2 + z2)
			r = sqrt(eigenvecs[2].at(0,0)*eigenvecs[2].at(0,0) + eigenvecs[2].at(1,0)*eigenvecs[2].at(1,0) + eigenvecs[2].at(2,0)*eigenvecs[2].at(2,0));
			theta = acos( eigenvecs[2].at(2,0) / r );
			polar_ang->Fill(theta);
			//azimuthal angle
			//phi = arctan(y/x)
			phi = atan(eigenvecs[2].at(1,0) / eigenvecs[2].at(0,0));
			azimuth_ang->Fill(phi);

			//average cluster energy
			e_avg->Fill(avg_Es[k]);
	
		}
		rhs.clear();
		eigenvals.clear();
		eigenvecs.clear();

	}
	ofile->cd();
	for(int i = 0; i < (int)hists1D.size(); i++) hists1D[i]->Write();
	nSubClusters->Write();
	//e_nSubClusters->Write();

}





