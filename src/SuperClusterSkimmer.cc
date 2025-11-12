#include "SuperClusterSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
//make cluster param histograms
void SuperClusterSkimmer::Skim(){
	if(_jsonfile != "" && _applyLumiMask && _data){
		_jsonfile = "config/json/"+_jsonfile;
		cout << "Applying lumi mask " << _jsonfile << endl;
		_jsonTool.BuildMap(_jsonfile);

	}

	cout << "Writing skim to: " << _oname << endl;
	cout << "Using clustering strategy NlnN full BHC with supercluster objects" << endl;

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
	//set NN model + features - can move to .C for more flexibility
	SetNNModel("config/json/small3CNN_EMultr_2017and2018.json");
	
	int nPho;
	
	vector<Jet> rhs;
	vector<JetPoint> rh_pts;
	vector<Jet> scs;
	if(_debug){ _oskip = 1000; }
	double sumE;

	if(_isocuts) cout << "Applying isolation for SCs matched to photons." << endl;	
	//set iso cuts
	if(_isocuts) _prod->SetIsoCut();
	//set energy weight transfer factor
	_prod->SetTransferFactor(_gev);
	_prod->ApplyFractions(_applyFrac);
	
	_prod->PrintPreselection();
	cout << "Jet selection for GJets CR" << endl;
	_jetprod->PrintPreselection();
	cout << "Minimum ht: " << _minHt_CRsel << endl;
	cout << "Maximum met: " << _maxMet_CRsel << endl;
	//loop over events
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	double pvx, pvy, pvz;
	_timeoffset = 0;
	vector<JetPoint> bhRhs;
	int nscran = 0;
	vector<pair<int,int>> icoords, icoords_nocenter;
	vector<double> Es, Es_nocenter;	
	//genpt of photons whose mom is ~40 (ISR)
	//genpt of photons whose mom is ~50 (meson decay)
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
		cout << "evt " << e << " ntuple event " << _base->Evt_event << endl;//" base is nullptr? " << (_base == nullptr) << endl;
		nEvts_tot++;

		//event filters - true is good event
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
		SetGJetsCR_EvtSel(e);
		FillBranch(_passGJetsEvtSel,"PassGJetsCR");	

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
		

		_prod->GetTrueSuperClusters(scs, e, _gev);
		//PV info
		pvx = _base->PV_x;
		pvy = _base->PV_y;
		pvz = _base->PV_z;
		BayesPoint PV({pvx, pvy, pvz});	

			
		int nSC = scs.size();
		int npho = _base->Photon_energy->size();
cout << "event " << e << " has " << nSC << " scs" << endl;
		//loop over selected scs
		int jet_scIdx = 0;
		int scidx;
		if(nSC < 1) continue;	
		for(int s = 0; s < nSC; s++){
			sumE = 0;
			//if(e % _oskip == 0) cout << "evt: " << e << " of " << _nEvts << "  sc: " << s << " of " << nPho << " nrhs: " << rhs.size()  << endl;
			scs[s].GetJets(rhs);

			//index in ntuples (before preselection)
			scidx = scs[s].GetUserIdx();
			if(rhs.size() < 1){ cout << "sc #" << s << " has " << rhs.size() << " rhs " << endl; continue; }
			cout << "evt: " << e << " of " << _nEvts << "  sc: " << s << " of " << nSC << " nrhs: " << rhs.size() << " E: " << scs[s].E() << endl;
		//cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " sc: " << p << " nrhs: " << rhs.size()  << flush;

			//fill eta phi map
			vector<double> neighborEs; 
			std::vector<std::pair<int, int>> icoords;
			vector<JetPoint> rh_pts; scs[s].GetJetPoints(rh_pts);
			GetNeighborE(rh_pts, -1, icoords, neighborEs,false,61);
			for(int e = 0; e < (int)neighborEs.size(); e++){
				ENeighbors->Fill(icoords[e].first, icoords[e].second, neighborEs[e]);
			}

			Jet bhc_sc;
			int ret = RunClustering(scs[s], bhc_sc, true, jet_scIdx); //fully remove PU clusters
			if(ret < 0){
				continue;
			}
		
			map<string,double> mapobs;
			//get id_idx of procCat that matches sample - still 1/sample but with correct labels now
			int id_idx = -999;
			//skip "total" procCat for always separated hists (id_idx == 1)
			

			vFillBranch((double)rhs.size(), "nRHs_ogSC");
			_prod->GetRhIdResMap(_rhIdToRes);

			nscran++;
			//add CMS benchmark variables - R9, Sietaieta, Siphiiphi, Smajor, Sminor
			//add CMS benchmark variable - isolation information
			//need to find associated photon
			//SuperCluster_photonIndx
			int np = _base->SuperCluster_PhotonIndx->at(scidx);
			//get jet points (rhs with IDs) for CNN grid
			bhc_sc.GetJetPoints(rh_pts);
			
			FillBranches(bhc_sc);
			int label = GetTrainingLabel(scidx, bhc_sc, scs[s]);
			//make CNN training grid
			addVector("rh_iEta",false);
			addVector("rh_iPhi",false);
			addVector("rh_energy",false);
			MakeCNNInputGrid(rh_pts, mapobs, jet_scIdx);
			PointCollection* points = GetPointsFromJet(scs[s]);
			double timesig_seed = CalcTimeSignificance(points,_base->SuperCluster_XtalSeedID->at(scs[s].GetUserIdx()));
			vFillBranch(timesig_seed, "timeSignificanceSeed");			
			
			mapobs["event"] = e;
			mapobs["event_weight"] = 1.;
			mapobs["object"] = scidx;
			mapobs["label"] = label;
			
			vector<float> ovalues; //discriminator output value, pass-by-ref
			CNNPredict(mapobs,ovalues);
			//TODO - update with discriminator score cut
			auto max_el = max_element(ovalues.begin(), ovalues.end());
			double predval = *max_el;
			//labeling starts from 1
			int nclass = std::distance(ovalues.begin(), max_el) + 1;
			cout << "class " << nclass << " predval " << predval << " for SC " << scidx << " with label " << label << endl;	
			vFillBranch((double)label, "trueLabel");
			vFillBranch((double)nclass, "predLabel");
			vFillBranch(ovalues[0], "predScore_physBkg");
			vFillBranch(ovalues[1], "predScore_BH");
			vFillBranch(ovalues[2], "predScore_spike");
			//write good SCs to CSV for training
			if(label != -1){
				BaseSkimmer::WriteObs(mapobs,"superclusters");
			}

			//make rh maps
			vector<JetPoint> rrhs; JetPoint rrh; 
			for(int r = 0; r < rhs.size(); r++){
				rhs[r].GetJetPointAt(0,rrh);
				rrhs.push_back(rrh);
			}
			cout << "label " << label << endl;
			for(int r = 0; r < rrhs.size(); r++){
				GetNeighborE(rrhs,r,icoords,Es);
				GetNeighborE(rrhs,r,icoords_nocenter,Es_nocenter,true);
				for(int ee = 0; ee < Es.size(); ee++){
					if(label == 1) //phys bkg
						ENeighbors_physBkg->Fill(icoords[ee].first,icoords[ee].second,Es[ee]);
					else if(label == 2) //BH
						ENeighbors_BH->Fill(icoords[ee].first,icoords[ee].second,Es[ee]);
					else if(label == 3) //spike
						ENeighbors_spikes->Fill(icoords[ee].first,icoords[ee].second,Es[ee]);
					else{ }
				}
				for(int ee = 0; ee < Es_nocenter.size(); ee++){
					if(label == 1) //phys bkg
						ENeighborsSkipCenter_physBkg->Fill(icoords_nocenter[ee].first,icoords_nocenter[ee].second,Es_nocenter[ee]);
					else if(label == 2) //BH
						ENeighborsSkipCenter_BH->Fill(icoords_nocenter[ee].first,icoords_nocenter[ee].second,Es_nocenter[ee]);
					else if(label == 3) //spike
						ENeighborsSkipCenter_spikes->Fill(icoords_nocenter[ee].first,icoords_nocenter[ee].second,Es_nocenter[ee]);
					else{ }
				}
			}
	
			jet_scIdx++;	
		}
		_tree->Fill();
		_reset();
	}
	_csvfile.close();
	cout << "\n" << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	_tree->SetDirectory(ofile);
	ofile->cd();
	_tree->Write();
	ofile->WriteTObject(ENeighbors_physBkg);	
	ofile->WriteTObject(ENeighbors_BH);	
	ofile->WriteTObject(ENeighbors_spikes);	
	ofile->WriteTObject(ENeighborsSkipCenter_physBkg);	
	ofile->WriteTObject(ENeighborsSkipCenter_BH);	
	ofile->WriteTObject(ENeighborsSkipCenter_spikes);	
	ofile->WriteTObject(ENeighbors);
	ofile->Close();
	cout << "Ran over " << nscran << " superclusters" << endl;
	cout << "Wrote skim to: " << _oname << endl;
	cout << "Wrote MVA inputs to " << _csvname << endl;
	cout << "Total events ran over " << nEvts_tot << " events that passed event filters " << nEvts_EvtFilterPass << " events that pass event filters and GJets selection " << nEvts_GJetsPass << endl;
}




void SuperClusterSkimmer::AddSample(TFile* file){

	

}
