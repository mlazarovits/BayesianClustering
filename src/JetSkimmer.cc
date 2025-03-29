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
	SetNNModel("json/small8CNN_EMultr.json");
	//SetNNFeatures();
	_nnfeatures = {"EMultr"};

	double metThresh = 0.4;
	double geoAvgJets;
	double phogev = _gev;//1./30.;
	_prod->PrintPreselection();
	for(int i = _evti; i < _evtj; i+=_skip){
		//do data MET selection
		_base->GetEntry(i);
		if(_BHFilter != notApplied){
                        if(_BHFilter == applied){
                                //apply beam halo filter - other noise filters needed for full Run2 recommendations
                                if(!_base->Flag_globalSuperTightHalo2016Filter){ cout << "BH Filter flagged - skipping" << endl; continue;}
                        }
                        else{
                                //inversely apply beam halo filter - other noise filters needed for full Run2 recommendations
                                if(_base->Flag_globalSuperTightHalo2016Filter) continue;
                        }
                }
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		_prod->GetTruePhotons(_phos, i, phogev);
		_scprod->GetTrueSuperClusters(_SCs, i, phogev);
		if(i % (_skip) == 0) cout << "evt: " << i << " of " << _nEvts;
		//calc max time for photons
		//cout << "lead photon pt " << _phos[0].pt() << " max time " << CalcMaxTime(_phos[0]) << endl;
		//if(_phos.size() > 1) cout << "sublead photon pt " << _phos[1].pt() << " max time " << CalcMaxTime(_phos[1]) << endl;

		if(_phos.size() > 1) FillTruePhotonHists(_phos);
		totEvt++;	
	

		////fill true jet histograms
		vector<Jet> jets;
		_prod->GetTrueJets(jets, i, _gev);
		if(jets.size() < 1){ cout << endl; continue; }
		FillTrueJetHists(jets);	
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

	
		if(i % (_skip) == 0) cout << " with " << jets.size() << " jets to cluster and " << _phos.size() << " photons";
		cout << endl;
		for(int i = 0; i < trCats.size(); i++)	
			//make sure time smearing doesn't happen here when it's turned off by the flag
			FillPVTimeHists(jets, i);
		
		jetSelEff++;

	}


	WriteHists(ofile);
	cout << "Total number of events ran over: " << totEvt << " events that had at least two jets that passed selection: " << jetSelEff << " fraction: " << jetSelEff/totEvt << endl;
	cout << "Total jets ran over: " << totJet << endl;
	cout << "Wrote skim to: " << _oname << endl;
}








