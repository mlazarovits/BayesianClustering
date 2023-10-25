#include "PhotonSkimmer.hh"
#include "Clusterizer.hh"
#include "Matrix.hh"

#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>
PhotonSkimmer::PhotonSkimmer(){ 
	_evti = 0;
	_evtj = 0;
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

	objE->SetTitle("phoE");
	objE->SetName("phoE");
}

//make cluster param histograms
void PhotonSkimmer::Skim(){

	cout << "Writing skim to: " << _oname << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");

	//create histograms to be filled
	MakeIDHists();
	
	int nPho;
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
	algo->SetWeighted(true);
	algo->SetVerbosity(0);
	algo->SetDataSmear(smear);


	GaussianMixture* gmm = new GaussianMixture();

	
	vector<JetPoint> rhs;
	double phoid, k;
	int eSkip = 10;
	if(_debug){ eSkip = 1000; }
	double sumE;



	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	for(int e = _evti; e < _evtj; e++){
		_base->GetEntry(e);
		nPho = (int)_base->Photon_energy->size();
		//if(i % 10 == 0) cout << "Event " << i << " events of " << _nEvts << endl;
		for(int p = 0; p < nPho; p++){
			sumE = 0;
			//find subclusters for each photon
			_prod->GetRecHits(rhs, e, p);
			//if(_debug) cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " nrhs: " << rhs.size()  << endl;
			if(e % eSkip == 0) cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " nrhs: " << rhs.size()  << endl;
			//else cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << "  pho: " << p << " nrhs: " << rhs.size()  << flush;
			
			if(rhs.size() < 1){ continue; }
			algo->SetMaxNClusters(int(rhs.size()));
			gmm = algo->FindSubjets(Jet(rhs));
			//get weight transfer factor - k = sum_n(E_n)/N = E_n/w_n = E_n*k/E_n for N rhs in a photon supercluster
			//w_n = E_n/k
			k = rhs[0].E()/gmm->GetData()->at(0).w();
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
					if(std::any_of(ids.begin(), ids.end(), [&](double iid){return iid == phoid;})){
						FillModelHists(gmm, i, k);
						plotCats[i].hists1D[19]->Fill(_base->Photon_energy->at(p));
						plotCats[i].hists2D[9]->Fill(_base->Photon_energy->at(p), sumE);
						break;
					}
				}
			}
			objE->Fill(_base->Photon_energy->at(p));
			objE_clusterE->Fill(_base->Photon_energy->at(p), sumE);
			FillTotalHists(gmm, k);
			rhs.clear();
		}
	}
	cout << "\n" << endl;
	WriteHists(ofile);

}

void PhotonSkimmer::CleaningSkim(){
	string _oname = "plots/photon_clean_skims_"+_cms_label+".root";
	cout << "Writing skim to: " << _oname << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	//rh time
	TH1D* t_rh = new TH1D("t_rh","t_rh",50,-150, 150);
	//rh time vs. rh e
	TH2D* TvErh = new TH2D("TvErh","TvErh",50, -150, 150., 1000, 0., 1000);
	TH2D* TvErh_lowE = new TH2D("TvErh_lowE","TvErh_lowE",50, -150, 150., 1000, 0., 10.);
	TH2D* TvErh_bx = new TH2D("TvErh_bx","TvErh_bx",50, -150, 150., 100, 0., 1000);
	TH2D* TvErh_cut = new TH2D("TvErh_cut","TvErh_cut",50, -100, 100., 1000, 0., 1000);

	TH1D* phoE = new TH1D("phoE","phoE",1000, 0, 1000);


	int nPho;
	vector<JetPoint> rhs;
	int eSkip = 10;
	if(_debug){ eSkip = 1000; }

	for(int i = 0; i < _nEvts; i+=eSkip){
		_base->GetEntry(i);
		nPho = (int)_base->Photon_energy->size();
		for(int p = 0; p < nPho; p++){
			phoE->Fill(_base->Photon_energy->at(p));
			//find subclusters for each photon
			_prod->GetRecHits(rhs, i, p);
			cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << "  pho: " << p << " nrhs: " << rhs.size()  << flush;
			for(int r = 0; r < (int)rhs.size(); r++){
				t_rh->Fill(rhs[r].t());
				TvErh->Fill(rhs[r].t(), rhs[r].E());
				if(rhs[r].E() <= 10.0) TvErh_lowE->Fill(rhs[r].t(), rhs[r].E());
				if(fabs(rhs[r].t()) < 25) TvErh_bx->Fill(rhs[r].t(), rhs[r].E());
				if(rhs[r].E() < 3.0 && fabs(rhs[r].t()) > 50) continue;
				TvErh_cut->Fill(rhs[r].t(), rhs[r].E());

			}
		}
	}
	
	ofile->cd();
	ofile->WriteTObject(t_rh);
	ofile->WriteTObject(TvErh);
	ofile->WriteTObject(TvErh_lowE);
	ofile->WriteTObject(TvErh_bx);
	ofile->WriteTObject(TvErh_cut);
	ofile->WriteTObject(phoE);
	ofile->Close();

}
