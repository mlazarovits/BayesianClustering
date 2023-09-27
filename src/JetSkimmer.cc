#include "JetSkimmer.hh"
#include "JetProducer.hh"
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
	TFile* ofile = new TFile(fname.c_str(),"RECREATE");

	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	
	
	Clusterizer* algo = new Clusterizer();
	algo->SetClusterAlpha(0.1);
	algo->SetSubclusterAlpha(0.1);
	algo->SetThresh(1.);
	algo->SetMaxNClusters(5);
	algo->SetWeighted(true);
	algo->SetVerbosity(0);
	algo->SetDataSmear(smear);

	map<string, Matrix> params;
	vector<node*> tree;
	vector<JetPoint> rhs;
	double k;
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		_prod->GetRecHits(rhs, i);
		cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << flush;

		//calculate transfer factor
		k = 0;
		for(int r = 0; r < (int)rhs.size(); i++) k += rhs[i].e();
		k /= (double)rhs.size();


		tree = algo->Cluster(Jet(rhs));

		FillPVHists(tree);
		
		for(int i = 0; i < (int)tree.size(); i++){	
			BasePDFMixture* model = tree[i]->model;
			FillModelHists(model, k);
		}
		
		//jet specific hists
		int njets = _base->Jet_energy->size();	
		nTrueJets->Fill((double)_base->Jet_energy->size());	
		nClusters->Fill((double)tree.size());
	
		for(int j = 0; j < njets; j++){	
			e_nRhs->Fill(_base->Jet_energy->at(j),(double)_base->Jet_drRhIds->at(j).size());
		}

	}

}





