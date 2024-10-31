#include "SuperClusterSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
//make cluster param histograms
void SuperClusterSkimmer::Skim(){

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
	vector<Jet> scs;
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
	
	_prod->SetTimeSmear(_timesmear);
	_prod->PrintPreselection();
	//loop over events
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	double pvx, pvy, pvz;
	_timeoffset = 0;
	_swcross = 0;
	int scidx;
	vector<JetPoint> bhRhs;
	double BHclusterPass = 0;
	double BHclusterFail = 0;
	double BHPass = 0;
	double BHFail = 0;
	int nscran = 0;	
	for(int e = _evti; e < _evtj; e++){
		_base->GetEntry(e);


		_prod->GetTrueSuperClusters(scs, e, _gev);
		//PV info
		pvx = _base->PV_x;
		pvy = _base->PV_y;
		pvz = _base->PV_z;
	
			
		int nSC = scs.size();
		int npho = _base->Photon_energy->size();
		//loop over selected scs
		for(int s = 0; s < nSC; s++){
			sumE = 0;
			//if(e % _oskip == 0) cout << "evt: " << e << " of " << _nEvts << "  sc: " << s << " of " << nPho << " nrhs: " << rhs.size()  << endl;
			scs[s].GetJets(rhs);
			//index in ntuples (before preselection)
			scidx = scs[s].GetUserIdx();
			if(rhs.size() < 1){ continue; }
			cout << "evt: " << e << " of " << _nEvts << "  sc: " << s << " of " << nSC << " nrhs: " << rhs.size() << endl;
		//cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " sc: " << p << " nrhs: " << rhs.size()  << flush;



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
	
			//make swiss cross
			_swcross = swissCross(rhs);
			//make BH filter
			_clusterSize = MakeBHFilterCluster(scs[s],bhRhs);
			if(_clusterSize <= 3) _BHcluster = false;
			else if(_clusterSize >= 6) _BHcluster = true;
			else{
				double td = BHTimeDiscriminant(bhRhs);
				if(_clusterSize == 4){
					//make sure seed > 10 GeV
					double maxE = 0;
					for(auto rh : bhRhs){
						if(rh.E() > maxE) maxE = rh.E();
					}
					if(maxE <= 10) _BHcluster = false;
					else{
						if(!BHclusterIso(bhRhs)) _BHcluster = false;
						else{
							if(td < 0) _BHcluster = true;
							else _BHcluster = false;
						}
					}
					//cout << "clustersize 4 - maxE " << maxE << " iso " << BHclusterIso(bhRhs) << " td " << td << " ~BH " << _BHcluster << endl;
				}
				if(_clusterSize == 5){
					if(td < 0) _BHcluster = true;
					else _BHcluster = false;
				}
			}
			if(!_base->Flag_globalSuperTightHalo2016Filter){
				BHFail++;
				if(_BHcluster) BHclusterPass++;
				else BHclusterFail++;
			}
			else BHPass++;
		//if(_clusterSize > 3) cout << "~BH " << _BHcluster << " cluster size " << _clusterSize << " true BH filter " << _base->Flag_globalSuperTightHalo2016Filter << endl;
		
			//one map per subcluster
			vector<map<string,double>> mapobs;			
			//get id_idx of procCat that matches sample - still 1/sample but with correct labels now
			int id_idx = -999;
			//skip "total" procCat for always separated hists
			FillModelHists(gmm, 1, mapobs);
			FillCMSHists(rhs,1);
			nscran++;
			//no separate categories for SCs - only fill "total" category
			_procCats[1].hists1D[0][4]->Fill(_base->SuperCluster_energy->at(scidx));
			_procCats[1].hists1D[0][226]->Fill(_base->SuperCluster_covEtaEta->at(scidx));
			_procCats[1].hists1D[0][227]->Fill(_base->SuperCluster_covPhiPhi->at(scidx));
		
			_procCats[1].hists1D[0][224]->Fill(_base->SuperCluster_smaj->at(scidx));
			_procCats[1].hists1D[0][225]->Fill(_base->SuperCluster_smin->at(scidx));
			_procCats[1].hists1D[0][228]->Fill(double(rhs.size()));
			if(_base->SuperCluster_energy->at(scidx) >= 0 && _base->SuperCluster_energy->at(scidx) < 200)
				_procCats[1].hists1D[0][229]->Fill(rhs.size());
			if(_base->SuperCluster_energy->at(scidx) >= 200 && _base->SuperCluster_energy->at(scidx) < 400)
				_procCats[1].hists1D[0][230]->Fill(rhs.size());
			if(_base->SuperCluster_energy->at(scidx) >= 400 && _base->SuperCluster_energy->at(scidx) < 600)
				_procCats[1].hists1D[0][231]->Fill(rhs.size());
			if(_base->SuperCluster_energy->at(scidx) >= 600 && _base->SuperCluster_energy->at(scidx) < 1000)
				_procCats[1].hists1D[0][232]->Fill(rhs.size());



			//add CMS benchmark variables - R9, Sietaieta, Siphiiphi, Smajor, Sminor
			//add CMS benchmark variable - isolation information
			//need to find associated photon
			//SuperCluster_photonIndx
			int np = _base->SuperCluster_PhotonIndx->at(scidx);
			
			
			for(int k = 0; k < mapobs.size(); k++){
				if(np != -1){
					//do 2017 preselection
					double r9 = _base->Photon_r9->at(np);
					double HoE = _base->Photon_hadOverEM->at(np);
					double Sieie = _base->Photon_SigmaIEtaIEta->at(np);
					double eIso = _base->Photon_ecalRHSumEtConeDR04->at(np);//_base->Photon_ecalPFClusterIso->at(np);
					double hIso = _base->Photon_hcalTowerSumEtConeDR04->at(np);//_base->Photon_hcalPFClusterIso->at(np);
					double tIso = _base->Photon_trkSumPtHollowConeDR03->at(np);
					double pt = _base->Photon_pt->at(np);
					//cout << "sc " << scidx << " pho " << np << " eIso/pt " << eIso/pt << " hIso/pt " << hIso/pt << " tIso/pt " << tIso/pt << " eIso " << eIso << " hIso " << hIso << " tIso " << tIso << " pt " << pt << endl;	
					if(r9 >= 0.9 && HoE <= 0.15 && Sieie <= 0.014 && eIso <= 5.0 + 0.01*pt && hIso <= 12.5 + 0.03*pt + 3.0e-5*pt*pt && tIso <= 6.0 + 0.002*pt && pt > 40){
						mapobs[k]["R9"] = r9;
						mapobs[k]["Sietaieta"] = Sieie;
						mapobs[k]["Siphiiphi"] = _base->SuperCluster_covPhiPhi->at(scidx);
						mapobs[k]["Smajor"] = _base->SuperCluster_smaj->at(scidx);
						mapobs[k]["Sminor"] = _base->SuperCluster_smin->at(scidx);
						//iso/pT
						mapobs[k]["ecalPFClusterIsoOvPt"] = eIso/pt;
						mapobs[k]["hcalPFClusterIsoOvPt"] = hIso/pt;
						mapobs[k]["trkSumPtHollowConeDR03OvPt"] = tIso/pt;
				
					}
					else{ //failed preselection
						mapobs[k]["R9"] = -999;
						mapobs[k]["Sietaieta"] = -999;
						mapobs[k]["Siphiiphi"] = -999;
						mapobs[k]["Smajor"] = -999;
						mapobs[k]["Sminor"] = -999;
						//iso/pT
						mapobs[k]["ecalPFClusterIsoOvPt"] = -999;
						mapobs[k]["hcalPFClusterIsoOvPt"] = -999;
						mapobs[k]["trkSumPtHollowConeDR03OvPt"] = -999;
					}
				}
				else{ //not matched to photon
					mapobs[k]["R9"] = -999;
					mapobs[k]["Sietaieta"] = -999;
					mapobs[k]["Siphiiphi"] = -999;
					mapobs[k]["Smajor"] = -999;
					mapobs[k]["Sminor"] = -999;
					//iso/pT
					mapobs[k]["ecalPFClusterIsoOvPt"] = -999;
					mapobs[k]["hcalPFClusterIsoOvPt"] = -999;
					mapobs[k]["trkSumPtHollowConeDR03OvPt"] = -999;
				}
				mapobs[k]["event"] = e;
				mapobs[k]["object"] = scidx;
				mapobs[k]["subcl"] = k;
				int label = GetTrainingLabel(scidx, k, gmm);
				mapobs[k]["label"] = label;

				BaseSkimmer::WriteObs(mapobs[k],"superclusters");
			}
		//cout << "n obs " << mapobs.size() << " " << _inputs.size() << endl;			
			
			
			objE_clusterE->Fill(_base->SuperCluster_energy->at(scidx), sumE);
		}
	}
	cout << "\n" << endl;
	ofile->WriteTObject(objE_clusterE);
	WriteHists(ofile);
	cout << "Ran over " << nscran << " superclusters" << endl;
	cout << "Wrote skim to: " << _oname << endl;
	cout << "Wrote MVA inputs to " << _csvname << endl;
	_csvfile.close();
	cout << "Total events that passed " << BHPass << " failed BH Filter " << BHFail << " \% events that passed cluster selection " << BHclusterPass/BHFail << "\% and \% of events that failed cluster selection given (for both) that the event failed the BH filter " << BHclusterFail/BHFail << "\%" << endl;
}




void SuperClusterSkimmer::AddSample(TFile* file){

	

}
