#include "JetSkimmer.hh"
#include "JetProducer.hh"
#include "Clusterizer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"
#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>

JetSkimmer::JetSkimmer(){ };



JetSkimmer::~JetSkimmer(){ 
}


JetSkimmer::JetSkimmer(TFile* file) : BaseSkimmer(file){
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
	string fname = "plots/jet_skims_"+_cms_label+".root";
	cout << "Writing skim to: " << fname << endl;
	cout << "Using clustering strategy";
	if(_strategy == NlnN)
		cout << " NlnN (Delauney)" << endl;
	else if(_strategy == N2)
		cout << " N2 (naive)" << endl;
	else
		cout << " undefined. Please use SetStrategy(i) with i == 0 (NlnN), 1 (N2)" << endl;
	TFile* ofile = new TFile(fname.c_str(),"RECREATE");

	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	

	double alpha = 0.1;
	double emAlpha = 0.5;
	double thresh = 1.;
	
	
	//nominal N^2 version
	Clusterizer* algo = new Clusterizer();
	algo->SetClusterAlpha(alpha);
	algo->SetSubclusterAlpha(emAlpha);
	algo->SetThresh(thresh);
	algo->SetMaxNClusters(5);
	algo->SetWeighted(true);
	algo->SetVerbosity(0);
	algo->SetDataSmear(smear);


	BayesCluster* algoDelauney = nullptr;

	map<string, Matrix> params;
	vector<node*> trees;
	vector<node*> delauneytrees;
	vector<JetPoint> rhs;
	vector<Jet> rhs_jet;
	double k;
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		_prod->GetRecHits(rhs_jet, i);
		_prod->GetRecHits(rhs, i);
		
		cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;

		for(int r = 0; r < rhs.size(); r++){
			t_rhs->Fill(rhs[r].t());
		}

		int njets = _base->Jet_energy->size();	
		nTrueJets->Fill((double)_base->Jet_energy->size());	
		if(njets < 1) continue;


		//run clustering
		//delauney NlnN version
		if(_strategy == NlnN){
			//need to transfer from GeV (energy) -> unitless (number of points)
			//transfer factor is over all points in event
			double gev = 0;
			for(int i = 0; i < (int)rhs_jet.size(); i++) gev += rhs_jet[i].E();
			gev = gev/(double)rhs_jet.size(); //gev = k = sum_n E_n/n pts
//			cout << "gev: " << gev << endl;
			for(int i = 0; i < (int)rhs_jet.size(); i++){ rhs_jet[i].SetWeight(rhs_jet[i].E()/gev); }//weights[i] /= gev; } //sums to n pts, w_n = E_n/k  
			algoDelauney = new BayesCluster(rhs_jet);
			algoDelauney->SetDataSmear(smear);
			algoDelauney->SetThresh(thresh);
			algoDelauney->SetAlpha(alpha);
			algoDelauney->SetSubclusterAlpha(emAlpha);
			algoDelauney->SetVerbosity(0);
			delauneytrees = algoDelauney->Cluster();
		}
		//N^2 version
		else if(_strategy == N2)
			trees = algo->Cluster(Jet(rhs));
		


		//calculate transfer factor
		k = rhs[0].E()/algo->GetData()->at(0).w();


		FillPVHists(trees);
		
		for(int i = 0; i < (int)trees.size(); i++){	
			BasePDFMixture* model = trees[i]->model;
			FillModelHists(model, k);
		}
		
		//jet specific hists
		nClusters->Fill((double)trees.size());
	
		for(int j = 0; j < njets; j++){	
			e_nRhs->Fill(_base->Jet_energy->at(j),(double)_base->Jet_drRhIds->at(j).size());
		}

	}
	WriteHists(ofile);
}





