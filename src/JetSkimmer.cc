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
	Matrix smear = Matrix(3,3);
	double dphi = 2*acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(0.,2,2); 
	//for time smearing (energy dependent)
	double tres_c = 0.2;
	double tres_n = sqrt(1 - tres_c*tres_c)*_gev;	

	double alpha = 0.1;
	double emAlpha = 0.5;
	double thresh = 1.;
	
	map<string, Matrix> params;
	vector<node*> trees;

	double jetSelEff = 0;
	double totEvt = 0;
        
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}

	cout << setprecision(10) << "weight " << _weight << endl; 

	double metThresh = 0.4;
	double geoAvgJets;
	double phogev = 1./30.;
	_prod->SetTimeSmear(_timesmear);
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
		if(i % (_skip) == 0) cout << "evt: " << i << " of " << _nEvts;

		FillTruePhotonHists(_phos);
		totEvt++;	
	

		////fill true jet histograms
		vector<Jet> jets;
		_prod->GetTrueJets(jets, i, _gev);
		FillTrueJetHists(jets);	
		if(jets.size() < 1){ cout << endl; continue; }
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
			FillPVTimeHists(jets, i, smear, emAlpha, alpha, tres_c, tres_n);
		
		jetSelEff++;

	}
	//get bin errors
	//get maps for LH ratio
	double err = -999;
	double bin = -999;
	Matrix E_tNeg12(5,5);
	Matrix E_t0(5,5);
	Matrix E_rh(5,5);
	for(int i = 0; i < 5; i++){
		for(int j = 0; j < 5; j++){
			//t ~ -12
			err = sqrt(trCats[0].procCats[1].hists2D[0][37]->GetBinError(i+1,j+1)/trCats[0].procCats[1].hists2D[0][37]->GetBinContent(i+1,j+1));
			trCats[0].procCats[1].hists2D[0][43]->Fill(i-2,j-2,err);
			bin = trCats[0].procCats[1].hists2D[0][37]->GetBinContent(i+1,j+1);
			E_tNeg12.SetEntry(bin,i,j);

			//t ~ -5
			err = sqrt(trCats[0].procCats[1].hists2D[0][38]->GetBinError(i+1,j+1)/trCats[0].procCats[1].hists2D[0][38]->GetBinContent(i+1,j+1));
			trCats[0].procCats[1].hists2D[0][44]->Fill(i-2,j-2,err);
			
			//t ~ 0
			err = sqrt(trCats[0].procCats[1].hists2D[0][39]->GetBinError(i+1,j+1)/trCats[0].procCats[1].hists2D[0][39]->GetBinContent(i+1,j+1));
			trCats[0].procCats[1].hists2D[0][45]->Fill(i-2,j-2,err);
			bin = trCats[0].procCats[1].hists2D[0][39]->GetBinContent(i+1,j+1);
			E_t0.SetEntry(bin,i,j);
		}
	}
//cout << "unnorm neg12" << endl; E_tNeg12.Print(); cout << "t0" << endl; E_t0.Print();
	//normalize timed maps
	double tneg12prod,t0prod, rh_norm;
	double tneg12norm = 0;
	double t0norm = 0;
	for(int i = 0; i < 5; i++){
		for(int j = 0; j < 5; j++){
			tneg12norm += E_tNeg12.at(i,j);	
			t0norm += E_t0.at(i,j);	
		}
	}		
			
	E_t0.mult(E_t0,1./t0norm);
	E_tNeg12.mult(E_tNeg12,1./tneg12norm);

//cout << "neg12" << endl; E_tNeg12.Print(); cout << "t0" << endl; E_t0.Print();

	//get jets for PV calc
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
		vector<Jet> jets;
		_prod->GetTrueJets(jets, i, _gev);
		if(jets.size() < 2) continue;
		pair<Jet,Jet> hardjets;
		int pair = FindJetPair(jets, hardjets);
		std::vector<std::pair<int, int>> icoords;
		vector<double> neighborEs;
		if(pair == -999) continue; //jets not did pass selection	
		vector<JetPoint> rhs1 = hardjets.first.GetJetPoints();
		vector<JetPoint> rhs2 = hardjets.second.GetJetPoints();

		for(int r = 0; r < rhs1.size(); r++){
			GetNeighborE(rhs1, r, icoords, neighborEs);
			for(int rr = 0; rr < neighborEs.size(); rr++){
				//cout << "idx " << rr << " E " << neighborEs[rr] << " icoord " << icoords[rr].first << ", " << icoords[rr].second << " for entry " << icoords[rr].first+2 << " , " << icoords[rr].second+2 << endl;
				E_rh.SetEntry(neighborEs[rr],icoords[rr].first+2,icoords[rr].second+2);
			}
			//E_rh.Print();
			rh_norm = 0;
			for(int i = 0; i < 5; i++){
				for(int j = 0; j < 5; j++){
					rh_norm += E_rh.at(i,j);	
				}
			}		
			E_rh.mult(E_rh,1./rh_norm);
			tneg12prod = E_tNeg12.FrobProd(E_rh);
			t0prod = E_t0.FrobProd(E_rh);
			//cout << "tneg12prod " << tneg12prod << " t0prod " << t0prod << endl;	
			
			trCats[0].procCats[1].hists1D[0][66]->Fill(tneg12prod/t0prod);
			E_rh.reset();
		}
		for(int r = 0; r < rhs2.size(); r++){
			GetNeighborE(rhs2, r, icoords, neighborEs);
			for(int rr = 0; rr < neighborEs.size(); rr++){
				//cout << "idx " << rr << " E " << neighborEs[rr] << " icoord " << icoords[rr].first << ", " << icoords[rr].second << " for entry " << icoords[rr].first+2 << " , " << icoords[rr].second+2 << endl;
				E_rh.SetEntry(neighborEs[rr],icoords[rr].first+2,icoords[rr].second+2);
			}
			//E_rh.Print();
			rh_norm = 0;
			for(int i = 0; i < 5; i++){
				for(int j = 0; j < 5; j++){
					rh_norm += E_rh.at(i,j);	
				}
			}		
			E_rh.mult(E_rh,1./rh_norm);
			tneg12prod = E_tNeg12.FrobProd(E_rh);
			t0prod = E_t0.FrobProd(E_rh);
			//cout << "tneg12prod " << tneg12prod << " t0prod " << t0prod << endl;	
			trCats[0].procCats[1].hists1D[0][66]->Fill(tneg12prod/t0prod);
			E_rh.reset();


		}
		
	}

	WriteHists(ofile);
	cout << "Total number of events ran over: " << totEvt << " events that had at least two jets that passed selection: " << jetSelEff << " fraction: " << jetSelEff/totEvt << endl;
	cout << "Wrote skim to: " << _oname << endl;
}








