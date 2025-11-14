#include "PhotonSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
//make cluster param histograms
void PhotonSkimmer::Skim(){
	if(_jsonfile != "" && _applyLumiMask && _data){
		_jsonfile = "config/json/"+_jsonfile;
		cout << "Applying lumi mask " << _jsonfile << endl;
		_jsonTool.BuildMap(_jsonfile);

	}
	//set histogram weights for HT slices, etc
	_weight = 1;
	if(_data){ _weight = 1.; }
	else if(_fname.find("QCD") != string::npos){
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
	else if(_fname.find("GJets") != string::npos){
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
	WriteHeader();

	int nPho;
	
	vector<Jet> rhs;
	vector<Jet> phos;
	vector<JetPoint> rh_pts;
	int phoid, genidx;
	if(_debug){ _oskip = 1000; }
	double sumE;

	_jetprod->SetTransferFactor(0.0333333);
	_jetprod->SetMinPt(_minJetPt_CRsel);
	_jetprod->SetMinNrhs(15);
	_jetprod->SetMinEmE(10);
	_jetprod->SetMinRhE(0.5);

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
	cout << "Jet preselection for isolation enforcement" << endl;
	_jetprod->PrintPreselection();
	cout << "Minimum ht: " << _minHt_CRsel << endl;
	cout << "Maximum met: " << _maxMet_CRsel << endl;
	//loop over events
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	double pvx, pvy, pvz;
	int phoidx, scidx;
	int nEvts_tot = 0;
	int nEvts_EvtFilterPass = 0;
	int nEvts_GJetsPass = 0;
	for(int e = _evti; e < _evtj; e++){
		_base->GetEntry(e);
		//apply lumi mask
		if(_applyLumiMask){
			if(!_jsonTool.IsGood(_base->Evt_run, _base->Evt_luminosityBlock) && _jsonfile != "" && _data){
				cout << "Skipping event " << e << " in run " << _base->Evt_run << " and lumi section " << _base->Evt_luminosityBlock << " due to lumi mask." << endl;
				continue;
			}
		}
		cout << "evt: " << e << " of " << _nEvts << " for run " << _base->Evt_run;
		nEvts_tot++;
		_prod->GetTruePhotons(phos, e, _gev);
		if(phos.size() < 1){ cout << endl; continue; }
		//PV info
		pvx = _base->PV_x;
		pvy = _base->PV_y;
		pvz = _base->PV_z;
		BayesPoint PV({pvx, pvy, pvz});	
		
		int nPho = phos.size();
		FillBranch(_base->Flag_BadChargedCandidateFilter ,"Flag_BadChargedCandidateFilter");
		FillBranch(_base->Flag_BadPFMuonDzFilter ,"Flag_BadPFMuonDzFilter");
		FillBranch(_base->Flag_BadPFMuonFilter ,"Flag_BadPFMuonFilter");
		FillBranch(_base->Flag_EcalDeadCellTriggerPrimitiveFilter ,"Flag_EcalDeadCellTriggerPrimitiveFilter");
		FillBranch(_base->Flag_HBHENoiseFilter ,"Flag_HBHENoiseFilter");
		FillBranch(_base->Flag_HBHENoiseIsoFilter ,"Flag_HBHENoiseIsoFilter");
		FillBranch(_base->Flag_ecalBadCalibFilter ,"Flag_ecalBadCalibFilter");
		FillBranch(_base->Flag_goodVertices ,"Flag_goodVertices");
		FillBranch(_base->Flag_hfNoisyHitsFilter ,"Flag_hfNoisyHitsFilter");
		FillBranch(_base->Flag_globalSuperTightHalo2016Filter ,"Flag_globalSuperTightHalo2016Filter");

		FillBranch((double)e,"evt");
		FillBranch(_weight,"evt_wt");
		FillBranch(_base->Met_pt,"MET");
		_jetprod->GetTrueJets(_jets, e);
		double ht = 0;
		for(auto j : _jets) ht += j.pt();
		FillBranch(ht,"ht");
		FillBranch((double)_jets.size(),"nSelJets");
		SetGJetsCR_EvtSel(e);
		FillBranch(_passGJetsEvtSel,"PassGJetsCR");	
		SetDijetsCR_EvtSel(e);
		FillBranch(_passDijetsEvtSel,"PassDijetsCR");	
		bool evtfilters = _base->Flag_BadChargedCandidateFilter && _base->Flag_BadPFMuonDzFilter && _base->Flag_BadPFMuonFilter && _base->Flag_EcalDeadCellTriggerPrimitiveFilter && _base->Flag_HBHENoiseFilter && _base->Flag_HBHENoiseIsoFilter && _base->Flag_ecalBadCalibFilter && _base->Flag_goodVertices && _base->Flag_hfNoisyHitsFilter;
		if((_oname.find("EGamma") != string::npos || _oname.find("DoubleEG") != string::npos || _oname.find("GJets") != string::npos)){
			if(!evtfilters){
				cout << "skipping event - failed event filters" << endl; 
				_tree->Fill();
				_reset();
				continue;
			} else nEvts_EvtFilterPass++;

			if(!_passGJetsEvtSel){
				cout << "skipping event - failed GJets CR selection" << endl; 
				_tree->Fill();
				_reset();
				continue;
			} else nEvts_GJetsPass++;
		}	
		int bhc_pho_idx = 0;
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
			vFillBranch(GJetsCR_ObjSel(phos[p]),"PassGJetsCR_Obj");

		//cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " pho: " << p << " nrhs: " << rhs.size()  << flush;
			Jet bhc_pho;
			int ret = RunClustering(phos[p], bhc_pho, false, bhc_pho_idx); //downweighting rhs
			if(ret < 0){
				continue;
			}
			bhc_pho_idx++;
			rhs.clear();
			bhc_pho.GetJets(rhs);
			//get parameters for model
			
			//vector<double> obs;				
			map<string,double> obs; //init obs map
			InitObs(obs);
			obs.at("event") = e;
			obs.at("event_weight") = _weight;
			obs.at("object") = phoidx;

			FillJetObs(bhc_pho, obs);
			//add CMS benchmark variables - R9, Sietaieta, Siphiiphi, Smajor, Sminor
			//add CMS benchmark variable - isolation information
			//need to find associated photon
			double r9 = _base->Photon_r9->at(phoidx);
			double HoE = _base->Photon_hadOverEM->at(phoidx);
			double Sieie = _base->Photon_sieie->at(phoidx);
			double Sipip = _base->Photon_sipip->at(phoidx);

			//isolation variables
			double eIso = _base->Photon_ecalRHSumEtConeDR04->at(phoidx);//_base->Photon_ecalPFClusterIso->at(phoidx);
			double hIso = _base->Photon_hcalTowerSumEtConeDR04->at(phoidx);//_base->Photon_hcalPFClusterIso->at(phoidx);
			double tIso = _base->Photon_trkSumPtHollowConeDR03->at(phoidx);
			double pt = _base->Photon_pt->at(phoidx);
			obs.at("R9") = r9;
                        obs.at("Sietaieta") = Sieie;
			obs.at("Siphiiphi") = Sipip;
                        obs.at("Smajor") = _base->SuperCluster_smaj->at(scidx);
                        obs.at("Sminor") = _base->SuperCluster_smin->at(scidx);
                        //iso/pT
                        obs.at("Pt") = pt;
                        obs.at("hcalTowerSumEtConeDR04") = hIso;
                        obs.at("trkSumPtHollowConeDR04") = _base->Photon_trkSumPtHollowConeDR04->at(phoidx);
                        obs.at("trkSumPtSolidConeDR04") = _base->Photon_trkSumPtSolidConeDR04->at(phoidx);
                        obs.at("hadTowOverEM") = _base->Photon_hadTowOverEM->at(phoidx);
                        obs.at("ecalRHSumEtConeDR04") = eIso;
			double ptmin;
			if(p == 0)
				ptmin = 70;
			else
				ptmin = 40;
			//cout << "sc " << scidx << " pho " << phoidx << " eIso " << eIso << " hIso " << hIso << " tIso " << tIso << " pt " << pt << " trkSumPtHollowConeDR04 "  << _base->Photon_trkSumPtHollowConeDR04->at(phoidx) << " trkSumPtSolidConeDR04 " << _base->Photon_trkSumPtSolidConeDR04->at(phoidx) << " hadTowOverEM " << _base->Photon_hadTowOverEM->at(phoidx) << endl;
			//do 2017 preselection
			if(r9 >= 0.9 && HoE <= 0.15 && Sieie <= 0.014 && eIso <= 5.0 + 0.01*pt && hIso <= 12.5 + 0.03*pt + 3.0e-5*pt*pt && tIso <= 6.0 + 0.002*pt && pt > ptmin){
                               		obs.at("2017_presel") = 1;
			
			}
			else{ //failed preselection
                               		obs.at("2017_presel") = 0;
			}
			
			int label = GetTrainingLabel(phoidx,bhc_pho,phos[p]);
		cout << " label for photon " << p << ", " << phoidx << " : " << label << endl; 
			//cout << "label: " << label << endl;
				obs.at("event") = e;
                                obs.at("object") = phoidx;
				obs.at("event_weight") = _weight;
			//only get lead subcluster -> ncl = 0
                                obs.at("label") = label;
				BaseSkimmer::WriteObs(obs,"photons");
				vFillBranch(label,"trueLabel");

			//do DNN prediction
			vector<double> dnn_scores;
			map<string, double> dnn_obs;
			MakeDNNInputs(bhc_pho, phoidx, dnn_obs);

			DNNPredict(dnn_obs, dnn_scores);
			//TODO - update with discriminator score cut
			auto max_el = max_element(dnn_scores.begin(), dnn_scores.end());
			double predval = *max_el;
			//labeling starts from 1
			int nclass = std::distance(dnn_scores.begin(), max_el) + 1;
			nclass = (nclass == 0) ? 4 : 6;
			cout << "class " << nclass << " predval " << predval << " for photon " << phoidx << " with label " << label << endl;	
			vFillBranch(nclass,"predLabel");
			vFillBranch(dnn_scores[0],"predScore_isoBkg");
			vFillBranch(dnn_scores[1],"predScore_nonIsoBkg");	

			FillBranches(obs);	
		}
		_tree->Fill();
		_reset();
	}
	cout << "\n" << endl;
	//ofile->WriteTObject(objE_clusterE);
	//WriteHists(ofile);
	_csvfile.close();
	cout << "\n" << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	_tree->SetDirectory(ofile);
	ofile->cd();
	_tree->Write();
	ofile->Close();
	cout << "Wrote skim to: " << _oname << endl;
	cout << "Wrote MVA inputs to " << _csvname << endl;
	cout << "Total events ran over " << nEvts_tot << " events that passed event filters " << nEvts_EvtFilterPass << " events that pass event filters and GJets selection " << nEvts_GJetsPass << endl;
}




void PhotonSkimmer::AddSample(TFile* file){


	PhotonProducer* prod = new PhotonProducer(file);
	ReducedBase* base = prod->GetBase();
	

}
