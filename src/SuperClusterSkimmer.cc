#include "SuperClusterSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"

#include <TFile.h>
#include <TH2D.h>
//make cluster param histograms
void SuperClusterSkimmer::Skim(){

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

	
	//set iso cuts
	if(_isocuts) _prod->SetIsoCut();
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
	int scidx;
	vector<JetPoint> bhRhs;
	double BHclusterPass = 0;
	double BHclusterFail = 0;
	double BHPass = 0;
	double BHFail = 0;
	int nscran = 0;
	vector<pair<int,int>> icoords, icoords_nocenter;
	vector<double> Es, Es_nocenter;	
	//genpt of photons whose mom is ~40 (ISR)
	//genpt of photons whose mom is ~50 (meson decay)
	int nspikes = 0;
	int nBH = 0;
	for(int e = _evti; e < _evtj; e++){
		_base->GetEntry(e);
		_obs.at("evt") = (double)e;
		_obs.at("evt_wt") = _weight;
		_obs.at("MET") = _base->Met_pt;
		cout << "evt " << e << " ntuple event " << _base->Evt_event << endl;//" base is nullptr? " << (_base == nullptr) << endl;
        	if(_BHFilter != notApplied){
        	        if(_BHFilter == applied){
        	                //apply beam halo filter - other noise filters needed for full Run2 recommendations
        	                if(!_base->Flag_globalSuperTightHalo2016Filter) continue;
        	        }
        	        else{
        	                //inversely apply beam halo filter - other noise filters needed for full Run2 recommendations
        	                if(_base->Flag_globalSuperTightHalo2016Filter) continue;
        	        }
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
			int ret = RunClustering(scs[s], bhc_sc, true, s); //fully remove PU clusters
			if(ret < 0){
				continue;
			}
			rhs.clear();
			bhc_sc.GetJets(rhs); 
	
			//make BH filter
			/*
			_clusterSize = MakeBHFilterCluster(bhc_sc,bhRhs);
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
			*/
		//if(_clusterSize > 3) cout << "~BH " << _BHcluster << " cluster size " << _clusterSize << " true BH filter " << _base->Flag_globalSuperTightHalo2016Filter << endl;
		
			map<string,double> mapobs;
			//get id_idx of procCat that matches sample - still 1/sample but with correct labels now
			int id_idx = -999;
			//skip "total" procCat for always separated hists (id_idx == 1)
			

			nscran++;
			//add CMS benchmark variables - R9, Sietaieta, Siphiiphi, Smajor, Sminor
			//add CMS benchmark variable - isolation information
			//need to find associated photon
			//SuperCluster_photonIndx
			int np = _base->SuperCluster_PhotonIndx->at(scidx);
			//get jet points (rhs with IDs) for CNN grid
			bhc_sc.GetJetPoints(rh_pts);
			
			//int label = GetTrainingLabel(scidx, k, gmm);
			int label = GetTrainingLabel(scidx, bhc_sc);
			//make CNN training grid
			MakeCNNInputGrid(rh_pts, mapobs);
			FillBranches(bhc_sc);
			
			mapobs["event"] = e;
			mapobs["event_weight"] = 1.;
			mapobs["object"] = scidx;
			mapobs["label"] = label;
			
			vector<double> ovalues; //discriminator output value, pass-by-ref
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
			if(label == 2) nBH++;
			if(label == 3) nspikes++;

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
				if(label == 1){ //phys bkg
					for(int ee = 0; ee < Es.size(); ee++)
						ENeighbors_physBkg->Fill(icoords[ee].first,icoords[ee].second,Es[ee]);
					for(int ee = 0; ee < Es_nocenter.size(); ee++)
						ENeighborsSkipCenter_physBkg->Fill(icoords_nocenter[ee].first,icoords_nocenter[ee].second,Es_nocenter[ee]);
						
				}
				else if(label == 2){ //BH
					for(int ee = 0; ee < Es.size(); ee++)
						ENeighbors_BH->Fill(icoords[ee].first,icoords[ee].second,Es[ee]);
					for(int ee = 0; ee < Es_nocenter.size(); ee++)
						ENeighborsSkipCenter_BH->Fill(icoords_nocenter[ee].first,icoords_nocenter[ee].second,Es_nocenter[ee]);
	
				}
				else if(label == 3){ //spike
					for(int ee = 0; ee < Es.size(); ee++)
						ENeighbors_spikes->Fill(icoords[ee].first,icoords[ee].second,Es[ee]);
					for(int ee = 0; ee < Es_nocenter.size(); ee++)
						ENeighborsSkipCenter_spikes->Fill(icoords_nocenter[ee].first,icoords_nocenter[ee].second,Es_nocenter[ee]);
				}
				else {}
	
			}

		//cout << "n obs " << mapobs.size() << " " << _inputs.size() << endl;			
			
			
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
	cout << "Total events that passed " << BHPass << " failed BH Filter " << BHFail << " \% events that passed cluster selection " << BHclusterPass/BHFail << "\% and \% of events that failed cluster selection given (for both) that the event failed the BH filter " << BHclusterFail/BHFail << "\%" << endl;
	cout << "total # spikes " << nspikes << " _nSpike_hist " << _nSpike_hist << " total # of BH " << nBH << " total # of bh hist " << _nBH_hist << endl;
}




void SuperClusterSkimmer::AddSample(TFile* file){

	

}
