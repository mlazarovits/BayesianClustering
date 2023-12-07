#include "PhotonSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
PhotonSkimmer::PhotonSkimmer(){ 
	_evti = 0;
	_evtj = 0;
	_isocuts = false;
	_oskip = 10;
	_thresh = 1.;
	_alpha = 0.1;
	_emAlpha = 0.5;
};



PhotonSkimmer::~PhotonSkimmer(){ 
}


PhotonSkimmer::PhotonSkimmer(TFile* file) : BaseSkimmer(file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time
	_prod = new PhotonProducer(file);
	_base = _prod->GetBase();
		//	_base = new ReducedBase(tree);
	_nEvts = _base->fChain->GetEntries();
	_evti = 0;
	_evtj = _nEvts;
	_oname = "plots/photon_skims_"+_cms_label+".root";

	_isocuts = false;
	_oskip = 10;
	_thresh = 1.;
	_alpha = 0.1;
	_emAlpha = 0.5;
	objE->SetTitle("phoE");
	objE->SetName("phoE");

	objE_clusterE->SetTitle("phoE_clusterE");
	objE_clusterE->SetName("phoE_clusterE");
}

//make cluster param histograms
void PhotonSkimmer::Skim(){

	cout << "Writing skim to: " << _oname << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");

	//create histograms to be filled
	MakeIDHists();
	
	//set energy weight transfer factor
	_prod->SetTransferFactor(_gev);
	
	int nPho;
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = 2*acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(0.,2,2); //no smear in time	
	double tres_c = 0.2;
	double tres_n = 30*sqrt(1 - tres_c*tres_c)*_gev;	

	
	vector<Jet> rhs;
	vector<Jet> phos;
	double phoid, k;
	if(_debug){ _oskip = 1000; }
	double sumE;


	//set iso cuts
	_prod->SetIsoCut();

	//loop over events
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	for(int e = _evti; e < _evtj; e++){
		_base->GetEntry(e);
		_prod->GetTruePhotons(phos, e);
		int nPho = phos.size();
		//loop over selected photons
		for(int p = 0; p < nPho; p++){
			sumE = 0;
			//if(e % _oskip == 0) cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " of " << nPho << " nrhs: " << rhs.size()  << endl;
			phos[p].GetJets(rhs);
		cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " pho: " << p << " nrhs: " << rhs.size()  << flush;
			//cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " of " << nPho << " nrhs: " << rhs.size()  << endl;

			BayesCluster *algo = new BayesCluster(rhs);
			algo->SetDataSmear(smear);
			//set time resolution smearing
			algo->SetTimeResSmear(tres_c, tres_n);
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(0);
			
			if(rhs.size() < 1){ continue; }
			GaussianMixture* gmm = algo->SubCluster();
			for(int r = 0; r < rhs.size(); r++) sumE += rhs[r].E();
			
			if(!_data){
				//find corresponding histogram category (signal, ISR, notSunm)	
				//split by LLP ID
				//0 = signal: chi_any -> gamma (22, 32, 25, 35)
				//1 = ISR: chi -> W/Z -> gamma, gino/sq -> q -> gamma (23, 33, 24, 34, 21, 31, 20, 30)
				//2 = not susy: p -> gamma, not matched (29, -1)
				phoid = _base->Photon_genLlpId->at(p);
				for(int i = 0; i < (int)plotCats.size(); i++){ //exclude total category - overlaps with above categories
					vector<double> ids = plotCats[i].ids;
					if(std::any_of(ids.begin(), ids.end(), [&](double iid){return iid == phoid;}) || i == 0){
						FillModelHists(gmm, i);
						if(i != 0) break;
					}
				}
			}
			objE->Fill(_base->Photon_energy->at(p));
			objE_clusterE->Fill(_base->Photon_energy->at(p), sumE);
		}
	}
	cout << "\n" << endl;
	ofile->WriteTObject(objE);
	ofile->WriteTObject(objE_clusterE);
	WriteHists(ofile);

}

