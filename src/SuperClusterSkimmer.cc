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
	//_csvfile.open(_csvname);
	//write header
	//WriteHeader();
	
	int nPho;
	
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
		if(_applyLumiMask && _data){
			if(!_jsonTool.IsGood(_base->Evt_run, _base->Evt_luminosityBlock) && _jsonfile != "" && _data){
				cout << "Skipping event " << e << " in run " << _base->Evt_run << " and lumi section " << _base->Evt_luminosityBlock << " due to lumi mask." << endl;
				continue;
			}
		}
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
		_jetprod->GetTrueJets(_jets, e);
		SetGJetsCR_EvtSel(e);
		FillBranch(_passGJetsEvtSel,"PassGJetsCR");	

		SetEventFilterPass();
		if((_oname.find("EGamma") != string::npos || _oname.find("DoubleEG") != string::npos || _oname.find("GJets") != string::npos)){
			if(!_evtfilters){
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
		//loop over selected scs
		int jet_scIdx = 0;
		int scidx;
		if(nSC < 1) continue;
		map<string, int> nSCs_type;
		nSCs_type[SCtypes[0]] = 0;
		nSCs_type[SCtypes[1]] = 0;
		nSCs_type[SCtypes[2]] = 0;
			
		for(int s = 0; s < nSC; s++){
			sumE = 0;
			//if(e % _oskip == 0) cout << "evt: " << e << " of " << _nEvts << "  sc: " << s << " of " << nPho << " nrhs: " << rhs.size()  << endl;
			//index in ntuples (before preselection)
			scidx = scs[s].GetUserIdx();
			if(scs[s].GetNPoints() < 1){ cout << "sc #" << s << " has " << scs[s].GetNPoints() << " rhs " << endl; continue; }
			if(_verb > -1) cout << "evt: " << e << " ntuple event " << _base->Evt_event << " of " << _nEvts << "  sc: " << s << " of " << nSC << " nrhs: " << scs[s].GetNPoints() << " E: " << scs[s].E() << endl;
		//cout << "\33[2K\r"<< "evt: " << e << " of " << _nEvts << " sc: " << p << " nrhs: " << rhs.size()  << flush;
			
			//per event
			_prod->GetRhIdResMap(_rhIdToRes);
			
			nscran++;
			
			map<string,Jet> SCtypes_map; //keys need to match SCtypes member var
			SCtypes_map[SCtypes[2]] = scs[s];
			
			unique_ptr<PointCollection> points = GetPointsFromJet(scs[s]);
			double timesig_seed = CalcTimeSignificance(points,_base->SuperCluster_XtalSeedID->at(scs[s].GetUserIdx()),false);
			vFillBranch(timesig_seed, "seedTimeSignificance_CMS");
			vFillBranch(0,"nSubclusters_CMS");


			Jet bhc_sc, bhc_sc_pucleaned;
			addVectors();
			int ret = RunClustering(scs[s], bhc_sc, bhc_sc_pucleaned, false, jet_scIdx,"BHCPUCleaned"); //downweight rhs according to good PU clusters
			SCtypes_map[SCtypes[0]] = bhc_sc;
			SCtypes_map[SCtypes[1]] = bhc_sc_pucleaned;
	
			int cmssc_label = -1;
			for(auto jt = SCtypes_map.begin(); jt != SCtypes_map.end(); jt++){
				Jet sc = jt->second;
				string tag = jt->first;
				//don't save invalid BHC clusterings
				if(_verb > 0) cout << "SC type " << tag << endl;
				//get id_idx of procCat that matches sample - still 1/sample but with correct labels now
				int id_idx = -999;
				//skip "total" procCat for always separated hists (id_idx == 1)
				//get SC points (rhs with IDs) for CNN grid
				sc.GetJetPoints(rh_pts);
				nSCs_type[tag]++;
			
				vFillBranch((double)s,"object_"+tag);	
				vFillBranch((double)rh_pts.size(), "nRHs_"+tag);
				
				FillBranches(sc,tag);
				int label = GetTrainingLabel(scidx, sc, scs[s], tag);
				if(tag == "CMS"){
					cmssc_label = label;
				}
				//make CNN training grid
				vector<vector<double>> rh_grid(_ngrid, vector<double>(_ngrid, 0.0));
				MakeCNNInputGrid(rh_pts, rh_grid, jet_scIdx, tag);
			
	
				vector<float> ovalues; //discriminator output value, pass-by-ref
				CNNPredict(rh_grid, ovalues);
				//TODO - update with discriminator score cut
				auto max_el = max_element(ovalues.begin(), ovalues.end());
				double predval = *max_el;
				//labeling starts from 1
				int nclass = std::distance(ovalues.begin(), max_el) + 1;
				if(_verb > 0) cout << "class " << nclass << " predval " << predval << " for SC " << scidx << " with label " << label << endl;	
				//only do for CMS superclusters since that's what the network was trained on
				
				if(tag != "CMS"){
					vFillBranch(-999, "predScore_physBkg_"+tag);
					vFillBranch(-999, "predScore_BH_"+tag);
					vFillBranch(-999, "predScore_spike_"+tag);
					if(ret < 0){
						vFillBranch(-999, "trueLabel_"+tag);
					}
					else{
						vFillBranch((double)label, "trueLabel_"+tag);
					}
				}
				else{
					vFillBranch((double)label, "trueLabel_"+tag);
					vFillBranch(ovalues[0], "predScore_physBkg_"+tag);
					vFillBranch(ovalues[1], "predScore_BH_"+tag);
					vFillBranch(ovalues[2], "predScore_spike_"+tag);
				}
			}
			jet_scIdx++;	
		}
		//cout << "event " << e << " has amount of following types of SCs " << endl;
		//for(auto it = nSCs_type.begin(); it != nSCs_type.end(); it++)
		//	cout << it->first << " " << it->second << endl;
		//if(e == 3){
		//	PrintTreeEntry();
		//}

		//cout << "============================================================================" << endl;
		//cout << endl;
		_tree->Fill();
		_reset();
	}
	//_csvfile.close();
	cout << "\n" << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	_tree->SetDirectory(ofile);
	ofile->cd();
	_tree->Write();
	ofile->Close();
	delete ofile;
	cout << "Ran over " << nscran << " superclusters" << endl;
	cout << "Wrote skim to: " << _oname << endl;
	//cout << "Wrote MVA inputs to " << _csvname << endl;
	cout << "Total events ran over " << nEvts_tot << " events that passed event filters " << nEvts_EvtFilterPass << " events that pass event filters and GJets selection " << nEvts_GJetsPass << endl;
}




void SuperClusterSkimmer::AddSample(TFile* file){

	

}
