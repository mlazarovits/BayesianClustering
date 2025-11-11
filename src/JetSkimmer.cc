#include "JetSkimmer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"
#include <TFile.h>
#include <time.h>
//#include <TH1D.h>
#include <TH2D.h>


//make cluster param histograms
//if specified, skim from events i to j
void JetSkimmer::Skim(){
	if(_jsonfile != "" && _applyLumiMask && _data){
		_jsonfile = "config/json/"+_jsonfile;
		cout << "Applying lumi mask " << _jsonfile << endl;
		_jsonTool.BuildMap(_jsonfile);
	}
	cout << "Writing skim to: " << _oname << endl;
	cout << "Using clustering strategy mixture model with pre-clustered AK4 jets (time calculated using MM components + naive methods)" << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	
	MakeTimeRecoCatHists();
	//create data smear matrix - smear in eta/phi
	double dphi = 2*acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	_smearMat = Matrix(3,3);
	_smearMat.SetEntry(deta*deta,0,0);
	_smearMat.SetEntry(dphi*dphi,1,1);
	_smearMat.SetEntry(0.,2,2); 

	//double alpha = 0.1;
	//double emAlpha = 0.5;
	//double thresh = 1.;
	
	map<string, Matrix> params;

	double jetSelEff = 0;
	double totEvt = 0;
	int totJet = 0;        
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}

	cout << setprecision(10) << "weight " << _weight << endl; 


	//set NN model + features - can move to .C for more flexibility
	SetNNModel("config/json/small3CNN_EMultr_2017and2018.json");
	//SetNNFeatures();
	_nnfeatures = {"EMultr"};

	double metThresh = 0.4;
	double geoAvgJets;
	double phogev = _gev;//1./30.;
	_prod->PrintPreselection();
	for(int i = _evti; i < _evtj; i+=_skip){
		_base->GetEntry(i);
		//apply lumi mask
		if(_applyLumiMask){
			if(!_jsonTool.IsGood(_base->Evt_run, _base->Evt_luminosityBlock) && _jsonfile != "" && _data){
				cout << "Skipping event " << i << " in run " << _base->Evt_run << " and lumi section " << _base->Evt_luminosityBlock << " due to lumi mask." << endl;
				continue;
			}
		}
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		vector<Jet> phos, bhc_phos;
		_prod->GetTruePhotons(phos, i, phogev);
		for(int p = 0; p < phos.size(); p++){
			Jet bhc_pho;
			int ret = RunClustering(phos[p], bhc_pho, true);
			if(ret == -1) continue;
			bhc_phos.push_back(bhc_pho);
		}
		if(i % (_skip) == 0) cout << "evt: " << i << " of " << _nEvts;
		//calc max time for photons
		//cout << "lead photon pt " << _phos[0].pt() << " max time " << CalcMaxTime(_phos[0]) << endl;
		//if(_phos.size() > 1) cout << "sublead photon pt " << _phos[1].pt() << " max time " << CalcMaxTime(_phos[1]) << endl;

		if(phos.size() > 1) FillTruePhotonHists(phos);
		totEvt++;	
	

		////fill true jet histograms
		vector<Jet> jets, bhc_jets;
		_prod->GetTrueJets(jets, i, _gev);
		if(jets.size() < 1){ cout << endl; continue; }
		
		//run BHC clustering
		for(int j = 0; j < jets.size(); j++){	
			Jet bhc_jet;
			int ret = RunClustering(jets[j], bhc_jet, true); //fully remove PU subclusters
			if(ret == -1) continue;
			bhc_jets.push_back(bhc_jet);
		}
		FillTrueJetHists(bhc_jets);	
		totJet += (int)jets.size();
		if(_data){
			//cut on ratio of MET pt to geo average of 2 leading jets
			geoAvgJets = jets[0].pt();
			if(jets.size() > 1){
				geoAvgJets *= jets[1].pt();
				geoAvgJets = sqrt(geoAvgJets);
			}
			if(_base->Met_pt/geoAvgJets > metThresh){ cout << endl; continue; }
		}

	
		if(i % (_skip) == 0) cout << " with " << jets.size() << " jets to cluster and " << phos.size() << " photons";
		cout << endl;
		for(int i = 0; i < trCats.size(); i++){
			vector<Jet> injets, inphos;
			if(i == 2){ //mmavg
				injets = bhc_jets;
				inphos = bhc_phos;
			}
			else{
				injets = jets;
				inphos = phos;
			}
			FillPVTimeHists(injets, inphos, i);
		}
		
		jetSelEff++;

	}


	WriteHists(ofile);
	cout << "Total number of events ran over: " << totEvt << " events that had at least two jets that passed selection: " << jetSelEff << " fraction: " << jetSelEff/totEvt << endl;
	cout << "Total jets ran over: " << totJet << endl;
	cout << "Wrote skim to: " << _oname << endl;
}








