#include "PhotonSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
//make cluster param histograms
void PhotonSkimmer::Skim(){
	//set histogram weights for HT slices, etc
	_weight = 1;
	if(_data){ _weight = 1.; }
	else if(_fname.find("QCD") != string::npos && !_isoBkgSel){
		cout << "Getting weights from info/EventWeights_QCD_SVIPM100_R18.txt for QCD" << endl;
	        ifstream weights("info/EventWeights_QCD_SVIPM100_R18.txt", std::ios::in);
	        string filein;
	        double jet_weight, pho_weight;
	        while( weights >> filein >> jet_weight >> pho_weight){
			if(_fname.find(filein) == string::npos) continue;
	                _weight = pho_weight;
	                break;
	        }
	}		
	else if(_fname.find("GJets") != string::npos && _isoBkgSel){
		cout << "Getting weights from info/EventWeights_GJets_SVIPM100_R18_isoBkgSel.txt for GJets with isolated bkg selection" << endl;
	        ifstream weights("info/EventWeights_GJets_SVIPM100_R18_isoBkgSel.txt", std::ios::in);
	        string filein;
	        double jet_weight, pho_weight;
	        while( weights >> filein >> jet_weight >> pho_weight){
			if(_fname.find(filein) == string::npos) continue;
	                _weight = pho_weight;
	                break;
	        }
	}
	else _weight = 1;	



	cout << "Writing skim to: " << _oname << endl;
	cout << "Using clustering strategy mixture model with pre-clustered photons" << endl;
	cout << setprecision(10) << "weight " << _weight << endl; 
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
	vector<JetPoint> rh_pts;
	int phoid, genidx;
	if(_debug){ _oskip = 1000; }
	double sumE;

	//set kin reqs for jets
	if(_isoBkgSel){
		_jetprod->SetTransferFactor(0.1);
		_jetprod->SetMinPt(_minJetPt_isoBkg);
		_jetprod->SetMinNrhs(15);
		_jetprod->SetMinEmE(10);
		_jetprod->SetMinRhE(0.5);
	}
	double ht, pho_phi, jet_phi, dphi_phojet;
	Jet jet_sys, pho_sys;
	double pi = 4*atan(1);
	bool l1seed, l1tohlt, hlt_phopt, hlt_jetht;

	double nIsoBkgPass = 0;
	double totEvt = 0;

	cout << "transfer factor (gev) N/Energy " << _gev << " EM alpha " << _emAlpha << " BHC alpha " << _alpha << endl;	
	cout << "Prior Parameters" << endl;
	cout << "beta0" << endl;
	_prior_params["scale"].Print();
	cout << "mean0" << endl;
	_prior_params["mean"].Print();
	cout << "nu0" << endl;
	_prior_params["dof"].Print();
	cout << "W0" << endl;
	_prior_params["W"].Print(); 
	
	//set iso cuts
	if(_isocuts){
		_prod->SetIsoCut();
		cout << "Applying isolation preselection for photons." << endl;
	}
	//set energy weight transfer factor
	_prod->SetTransferFactor(_gev);
	_prod->ApplyFractions(_applyFrac);

	cout << "Photon preselection" << endl;
	_prod->PrintPreselection();
	if(_isoBkgSel){
		cout << "\nSelecting photons with isolated background preselection." << endl;
		cout << "Jet preselection" << endl;
		_jetprod->PrintPreselection();
		cout << "Minimum ht: " << _minHt_isoBkg << endl;
		cout << "Maximum met: " << _maxMet_isoBkg << endl;
	}
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
		cout << "evt: " << e << " of " << _nEvts;
		if(_BHFilter != notApplied){
		cout << "BH filter applied - flag " << _base->Flag_globalSuperTightHalo2016Filter << " BHFilter " << _BHFilter << endl;
                        if(_BHFilter == applied){
                                //apply beam halo filter - other noise filters needed for full Run2 recommendations
                                if(!_base->Flag_globalSuperTightHalo2016Filter) continue;
                        }
                        else{
                                //inversely apply beam halo filter - other noise filters needed for full Run2 recommendations
                                if(_base->Flag_globalSuperTightHalo2016Filter) continue;
                        }
                }  
		totEvt++;

		_prod->GetTruePhotons(phos, e, _gev);
		if(phos.size() < 1){ cout << endl; continue; }
		//PV info
		pvx = _base->PV_x;
		pvy = _base->PV_y;
		pvz = _base->PV_z;
		BayesPoint PV({pvx, pvy, pvz});	
		//do iso bkg evt selection to compare data/MC
		if(_isoBkgSel){
			//L1 seed
			if(!_base->Trigger_hltL1sSingleEGNonIsoOrWithJetAndTauNoPS) continue;
			//cout << "passed L1 seed" << endl;	
			//L1 to HLT Regional EGM matching leg
			if(!_base->Trigger_hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter) continue;
			//cout << "passed L1 to HLT" << endl;	
			//photon pt > 60
			if(!_base->Trigger_hltEG60EtFilter) continue;
			//cout << "passed photon pt > 60" << endl;	
			//jet ht > 175 && jet pt > 10 && |eta jet| < 3
			if(!_base->Trigger_hltHT175Jet10) continue;
			//cout << "passed HT > 175" << endl;	
			//jet ht > 350 && jet pt > 15 && |eta jet| < 3
			if(!_base->Trigger_hltPFHT350Jet15) continue;
			//cout << "passed HT > 135" << endl;	
			
			//gev = 0.1 for jets
			_jetprod->GetTrueJets(_jets, e);
			//min photon multiplicity
			if(phos.size() < 1) continue;
			//cout << "passed pho mult" << endl;	
			//min jet multiplicity
			if(_jets.size() < 1) continue;
			//cout << "passed jet mult" << endl;	
		
			//ht - scalar sum
			ht = 0;
			for(auto j : _jets) ht += j.pt();
			//dphi bw photon and jet systems (vector sum of objects)
			jet_sys = VectorSum(_jets);
			
			//MET upper limit - orthogonal to signal MET selection
			if(_base->Met_pt > _maxMet_isoBkg) continue;
			//cout << "passed max met" << endl;	
			//min jet ht
			if(ht < _minHt_isoBkg) continue; 
			//cout << "passed min ht" << endl;	
			//az angle bw hardest presel photon + jet system
			//cout << "dphi " << dphi_phojet << " met " << _base->Met_pt << endl;	
	
			nIsoBkgPass++;
		}
		
		int nPho = phos.size();
		//loop over selected photons
		for(int p = 0; p < nPho; p++){
			sumE = 0;
			//if(e % _oskip == 0) cout << "evt: " << e << " of " << _nEvts << "  pho: " << p << " of " << nPho << " nrhs: " << rhs.size()  << endl;
			phos[p].GetJets(rhs);
			double pho_rhE = 0;
			for(auto rh : rhs)
				pho_rhE += rh.E();
			phoidx = phos[p].GetUserIdx();
			scidx = _base->Photon_scIndex->at(phoidx);
			if(rhs.size() < 1){ cout << endl; continue; }
			cout << "  pho: " << p << " of " << nPho << " nrhs: " << rhs.size()  << " pt " << phos[p].pt() << " E " << phos[p].E() << " rh E " << pho_rhE << endl;
			if(_isoBkgSel){
				pho_phi = phos[p].phi_02pi(); 
				cout << "\npho system E " << phos[p].E() << " phi " << phos[p].phi_02pi() << " jet system E " << jet_sys.E() << " phi " << jet_sys.phi_02pi() << endl;
				cout << "pho system pt " << phos[p].pt() << " phi " << phos[p].phi_02pi() << " jet system pt " << jet_sys.pt() << " phi " << jet_sys.phi_02pi() << endl;
				jet_phi = jet_sys.phi_02pi();
				dphi_phojet = pho_phi - jet_phi;
				dphi_phojet = acos(cos(dphi_phojet)); //wraparound - will always be < pi
				cout << "# jets " << _jets.size() << " ht " << ht << " met " << _base->Met_pt << " dphi " << dphi_phojet << endl;
				if(dphi_phojet < pi-0.3) continue; //want dphi ~ phi - implies less MET in event
				//cout << "passed dphi " << endl;	
				//trigger req - take baseline, photon pt leg + jet ht legs from HLT Photon60 R9Id90 CaloIdL IsoL DisplacedIdL PFHT350MinPFJet15
			}
		//cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " pho: " << p << " nrhs: " << rhs.size()  << flush;
			Jet bhc_pho;
			int ret = RunClustering(rhs, PV, bhc_pho);
			if(ret == -1){
				continue;
			}
			rhs.clear();
			bhc_pho.GetJets(rhs); 
			//get parameters for model
	
			_procCats[1].hists1D[0][266]->Fill(bhc_pho.E());

			_procCats[1].hists1D[0][259]->Fill(phos[p].pt());
			
			vector<pair<int,int>> icoords;
			vector<double> neighborEs;
			rh_pts = phos[p].GetJetPoints();
			GetNeighborE(rh_pts, -1, icoords, neighborEs,false,9);
			for(int e = 0; e < (int)neighborEs.size(); e++){
				_procCats[1].hists2D[0][237]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
			}
			double maxE = 0;
			Jet maxE_rh;
			PointCollection* points = new PointCollection();
			for(int r = 0; r < rhs.size(); r++){
				sumE += rhs[r].E();
				if(rhs[r].E() > maxE){
					maxE = rhs[r].E();
					maxE_rh = rhs[r];
				}
				_procCats[0].hists1D[0][258]->Fill(rhs[r].t());
				_procCats[1].hists1D[0][258]->Fill(rhs[r].t());
				BayesPoint pt({rhs[r].eta(), rhs[r].phi(), rhs[r].t()});
				pt.SetWeight(rhs[r].E()*_gev);
				points->AddPoint(pt);
			}
			for(int r = 0; r < rhs.size(); r++){
				_procCats[0].hists2D[0][238]->Fill(rhs[r].eta() - maxE_rh.eta(), acos(cos(rhs[r].phi() - maxE_rh.phi())), rhs[r].E()*_weight);
				_procCats[1].hists2D[0][238]->Fill(rhs[r].eta() - maxE_rh.eta(), acos(cos(rhs[r].phi() - maxE_rh.phi())), rhs[r].E()*_weight);

			}


			_swcross = swissCross(rhs);
			//vector<double> obs;				
			map<string,double> obs; //init obs map
			for(int d = 0; d < _inputs.size(); d++){
				obs[_inputs[d]] = -999;
			}

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
						FillModelHists(bhc_pho, i, obs);
						//FillModelHists(gmm, i, obs);
						//FillCMSHists(rhs,i);
						_procCats[i].hists1D[0][4]->Fill(_base->Photon_energy->at(phoidx));
						_procCats[i].hists1D[0][226]->Fill(_base->Photon_sieie->at(phoidx));
						_procCats[i].hists1D[0][227]->Fill(_base->Photon_sipip->at(phoidx));
		
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
					FillModelHists(bhc_pho, i, obs);
					//FillModelHists(gmm, i, obs);
					//FillCMSHists(rhs,i);
					_procCats[i].hists1D[0][4]->Fill(_base->Photon_energy->at(phoidx));
					_procCats[i].hists1D[0][226]->Fill(_base->Photon_sieie->at(phoidx));
					_procCats[i].hists1D[0][227]->Fill(_base->Photon_sipip->at(phoidx));
		
					scidx = _base->Photon_scIndex->at(phoidx);
					_procCats[i].hists1D[0][224]->Fill(_base->SuperCluster_smaj->at(scidx));
					_procCats[i].hists1D[0][225]->Fill(_base->SuperCluster_smin->at(scidx));
					_procCats[i].hists1D[0][228]->Fill(rhs.size());
				}
				

			}
			//add CMS benchmark variables - R9, Sietaieta, Siphiiphi, Smajor, Sminor
			//add CMS benchmark variable - isolation information
			//need to find associated photon
			//do 2017 preselection
			double r9 = _base->Photon_r9->at(phoidx);
			double HoE = _base->Photon_hadOverEM->at(phoidx);
			double Sieie = _base->Photon_sieie->at(phoidx);
			double Sipip = _base->Photon_sipip->at(phoidx);
			double eIso = _base->Photon_ecalRHSumEtConeDR04->at(phoidx);//_base->Photon_ecalPFClusterIso->at(phoidx);
			double hIso = _base->Photon_hcalTowerSumEtConeDR04->at(phoidx);//_base->Photon_hcalPFClusterIso->at(phoidx);
			double tIso = _base->Photon_trkSumPtHollowConeDR03->at(phoidx);
			double pt = _base->Photon_pt->at(phoidx);
			double ptmin;
			if(p == 0)
				ptmin = 70;
			else
				ptmin = 40;
			//cout << "sc " << scidx << " pho " << phoidx << " eIso/pt " << eIso/pt << " hIso/pt " << hIso/pt << " tIso/pt " << tIso/pt << " eIso " << eIso << " hIso " << hIso << " tIso " << tIso << " pt " << pt << endl;	
			if(r9 >= 0.9 && HoE <= 0.15 && Sieie <= 0.014 && eIso <= 5.0 + 0.01*pt && hIso <= 12.5 + 0.03*pt + 3.0e-5*pt*pt && tIso <= 6.0 + 0.002*pt && pt > ptmin){
                                        obs["R9"] = r9;
                                        obs["Sietaieta"] = Sieie;
					obs["Siphiiphi"] = Sipip;
                                        obs["Smajor"] = _base->SuperCluster_smaj->at(scidx);
                                        obs["Sminor"] = _base->SuperCluster_smin->at(scidx);
                                        //iso/pT
                                        obs["ecalPFClusterIsoOvPt"] = eIso/pt;
                                        obs["hcalPFClusterIsoOvPt"] = hIso/pt;
                                        obs["trkSumPtHollowConeDR03OvPt"] = tIso/pt;	
			
			}
			else{ //failed preselection
                                        obs["R9"] = -999;
                                        obs["Sietaieta"] = -999;
					obs["Siphiiphi"] = -999;
                                        obs["Smajor"] = -999;
                                        obs["Sminor"] = -999;
                                        //iso/pT
                                        obs["ecalPFClusterIsoOvPt"] = -999;
                                        obs["hcalPFClusterIsoOvPt"] = -999;
                                        obs["trkSumPtHollowConeDR03OvPt"] = -999;	
			}
			
			int label = GetTrainingLabel(phoidx,bhc_pho);
		cout << " label for photon " << p << ", " << phoidx << " : " << label << endl; 
			//cout << "label: " << label << endl;
				obs["event"] = e;
                                obs["object"] = phoidx;
				obs["event_weight"] = _weight;
			//only get lead subcluster -> ncl = 0
                                obs["subcl"] = 0;
				obs["lead"] = 0;
				if(p == 0)
					obs["lead"] = 1;
                                obs["label"] = label;
				BaseSkimmer::WriteObs(obs,"photons");
				objE_clusterE->Fill(_base->SuperCluster_energy->at(scidx), sumE);
		//cout << "n obs " << obs.size() << " " << _inputs.size() << endl;
		//int no = 0;
		//for(auto o : obs)
		//	cout << o.first << ": " << o.second << endl;
		//cout << endl;
		//for(int o = 0; o < _inputs.size(); o++)
		//	cout << o << ": " << _inputs[o] << endl;	
		}
	}
	cout << "\n" << endl;
	ofile->WriteTObject(objE_clusterE);
	WriteHists(ofile);
	cout << "Wrote skim to: " << _oname << endl;
	cout << "Wrote MVA inputs to " << _csvname << endl;
	_csvfile.close();

	cout << "Total number of events ran over: " << totEvt << " events that passed isolated bkg selection: " << nIsoBkgPass << " fraction: " << nIsoBkgPass/totEvt << endl;
}




void PhotonSkimmer::AddSample(TFile* file){


	PhotonProducer* prod = new PhotonProducer(file);
	ReducedBase* base = prod->GetBase();
	

}
