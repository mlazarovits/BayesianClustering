#include "PhotonSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
//make cluster param histograms
void PhotonSkimmer::Skim(){

	cout << "Writing skim to: " << _oname << endl;
	cout << "Using clustering strategy mixture model with pre-clustered AK4 jets (time calculated using MM components + naive methods)" << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");

	MakeProcCats(_oname);
	
	//make output csv file
	//write to csv dir if not condor job else write to current dir
	if(_oname.find("condor") == string::npos){
		_csvname = _oname.substr(_oname.find("/"));
		_csvname = _csvname.substr(0,_csvname.find(".root"));
		_csvname = "csv"+_csvname+".csv";
	}
	else{
		_csvname = _oname.substr(0,_oname.find(".root"));
		_csvname = _csvname+".csv";
	}
	_csvfile.open(_csvname);
	//write header
	SetObs();
	
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
	//if(_timesmear) cout << "Smearing covariance in time with energy dependence." << endl;	

	
	vector<Jet> rhs;
	vector<Jet> phos;
	int phoid, genidx;
	if(_debug){ _oskip = 1000; }
	double sumE;


	//to check histogram indices
	//for(int i = 0; i < _hists2D.size(); i++)
	//	cout << "i: " << i << " hist: " << _hists2D[i]->GetName() << endl;
	//return;	
	//set iso cuts
	_prod->SetIsoCut();
	//set energy weight transfer factor
	_prod->SetTransferFactor(_gev);
	_prod->ApplyFractions(_applyFrac);

	_prod->PrintPreselection();
	//loop over events
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	double pvx, pvy, pvz;
	_timeoffset = 0;
	_swcross = 0;
	int phoidx, scidx;
	for(int e = _evti; e < _evtj; e++){
		_base->GetEntry(e);
		_prod->GetTruePhotons(phos, e, _gev);
		//PV info
		pvx = _base->PV_x;
		pvy = _base->PV_y;
		pvz = _base->PV_z;
		
		int nPho = phos.size();
		//loop over selected photons
		for(int p = 0; p < nPho; p++){
			sumE = 0;
			//if(e % _oskip == 0) cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " of " << nPho << " nrhs: " << rhs.size()  << endl;
			phos[p].GetJets(rhs);
			phoidx = phos[p].GetUserIdx();
			if(rhs.size() < 1){ continue; }
			cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " of " << nPho << " nrhs: " << rhs.size()  << endl;
		//cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " pho: " << p << " nrhs: " << rhs.size()  << flush;

			BayesCluster *algo = new BayesCluster(rhs);
			if(_smear) algo->SetDataSmear(smear);
			//set time resolution smearing
			if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n);
			algo->SetThresh(_thresh);
			algo->SetAlpha(_alpha);
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetVerbosity(0);
			
			GaussianMixture* gmm = algo->SubCluster();
			for(int r = 0; r < rhs.size(); r++) sumE += rhs[r].E();
	
			_swcross = swissCross(rhs);
			vector<double> obs;				
			if(!_data){
				//find corresponding histogram category (signal, ISR, notSunm)	
				//split by LLP ID
				//0 = all
				//1 = signal: chi_any -> gamma (22, 32, 25, 35)
				//2 = not susy or not matched: p -> gamma, not matched (29, -1)
				genidx = _base->Photon_genIdx->at(p);
				if(genidx == -1) phoid = -1;
				else phoid = _base->Gen_susId->at(genidx);
				for(int i = 0; i < (int)_procCats.size(); i++){ //exclude total category - overlaps with above categories
					vector<double> ids = _procCats[i].ids;
					if(std::any_of(ids.begin(), ids.end(), [&](double iid){return (iid == double(phoid)) || (iid == -999);})){
						FillModelHists(gmm, i, obs);
						FillCMSHists(rhs,i);
						_procCats[i].hists1D[0][4]->Fill(_base->Photon_energy->at(phoidx));
						_procCats[i].hists1D[0][226]->Fill(_base->Photon_sieie->at(phoidx));
						_procCats[i].hists1D[0][227]->Fill(_base->Photon_sipip->at(phoidx));
		
						scidx = _base->Photon_scIndex->at(phoidx);
						_procCats[i].hists1D[0][224]->Fill(_base->SuperCluster_smaj->at(scidx));
						_procCats[i].hists1D[0][225]->Fill(_base->SuperCluster_smin->at(scidx));
						_procCats[i].hists1D[0][228]->Fill(double(rhs.size()));
						if(_base->Photon_energy->at(phoidx) >= 0 && _base->Photon_energy->at(phoidx) < 200)
							_procCats[i].hists1D[0][229]->Fill(rhs.size());
						if(_base->Photon_energy->at(phoidx) >= 200 && _base->Photon_energy->at(phoidx) < 400)
							_procCats[i].hists1D[0][230]->Fill(rhs.size());
						if(_base->Photon_energy->at(phoidx) >= 400 && _base->Photon_energy->at(phoidx) < 600)
							_procCats[i].hists1D[0][231]->Fill(rhs.size());
						if(_base->Photon_energy->at(phoidx) >= 600 && _base->Photon_energy->at(phoidx) < 1000)
							_procCats[i].hists1D[0][232]->Fill(rhs.size());

					}
				}
			}
			else{
				for(int i = 0; i < (int)_procCats.size(); i++){ //exclude total category - overlaps with above categories
					FillModelHists(gmm, i, obs);
					FillCMSHists(rhs,i);
					_procCats[i].hists1D[0][4]->Fill(_base->Photon_energy->at(phoidx));
					_procCats[i].hists1D[0][226]->Fill(_base->Photon_sieie->at(phoidx));
					_procCats[i].hists1D[0][227]->Fill(_base->Photon_sipip->at(phoidx));
		
					scidx = _base->Photon_scIndex->at(phoidx);
					_procCats[i].hists1D[0][224]->Fill(_base->SuperCluster_smaj->at(scidx));
					_procCats[i].hists1D[0][225]->Fill(_base->SuperCluster_smin->at(scidx));
					_procCats[i].hists1D[0][228]->Fill(rhs.size());
				}
				

			}
			int label = GetTrainingLabel(p,gmm);
			WriteObs(e,p,obs,label);
			objE_clusterE->Fill(_base->Photon_energy->at(p), sumE);
		}
	}
	cout << "\n" << endl;
	ofile->WriteTObject(objE_clusterE);
	WriteHists(ofile);
	cout << "Wrote skim to: " << _oname << endl;
	cout << "Wrote MVA inputs to " << _csvname << endl;
	_csvfile.close();

}




void PhotonSkimmer::AddSample(TFile* file){


	PhotonProducer* prod = new PhotonProducer(file);
	ReducedBase* base = prod->GetBase();
	

}
