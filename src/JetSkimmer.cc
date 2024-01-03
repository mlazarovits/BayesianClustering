#include "JetSkimmer.hh"
#include "Clusterizer.hh"
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
	cout << "Using clustering strategy";
	if(_strategy == NlnN)
		cout << " NlnN (Delauney)" << endl;
	else if(_strategy == N2)
		cout << " N2 (naive)" << endl;
	else if(_strategy == MM)
		cout << " mixture model with pre-clustered AK4 jets (time calculated using MM components + naive methods)" << endl;
	else
		cout << " undefined. Please use SetStrategy(i) with i == 0 (NlnN), 1 (N2), 2 (MM)" << endl;
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	//set differences in samples (ie GMSB, data) here
	
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
	vector<Jet> rhs;

	double jetSelEff = 0;
	double totEvt = 0;
        
	//for computational time
	vector<double> x_nrhs, y_time;
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	


	int SKIP = 1;
	for(int i = _evti; i < _evtj; i++){
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		_prod->GetTruePhotons(_phos, i);
		if(i % (SKIP) == 0) cout << "evt: " << i << " of " << _nEvts;
		_prod->GetRecHits(rhs, i);
		x_nrhs.push_back((double)rhs.size());
		for(int r = 0; r < rhs.size(); r++){
			rhTime->Fill(rhs[r].t());
		}
		
		totEvt++;	
	
		//fill true jet histograms
		vector<Jet> jets;
		_prod->GetTrueJets(jets, i);
		if(jets.size() < 1){ cout << endl; continue; }
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;

	
		if(i % (SKIP) == 0) cout << " with " << jets.size() << " jets to cluster and " << _phos.size() << " photons";
		FillTrueJetHists(jets);
		for(int i = 0; i < trCats.size(); i++)	
			//make sure time smearing doesn't happen here when it's turned off by the flag
			//FillPVTimeHists(jets, i, smear, emAlpha, alpha, tres_c, tres_n);
			FillPVTimeHists(jets, i); //turn off time and spatial smearing
		
		jetSelEff++;
		

		//fill mm only jet hists - done above
		if(_strategy == MM){
			cout << endl;
			continue;	
		}

		clock_t t;
		BayesCluster* algo = new BayesCluster(rhs);
		if(_smear) algo->SetDataSmear(smear);
		//algo->SetTimeResSmear(tres_c, tres_n);
		algo->SetThresh(thresh);
		algo->SetAlpha(alpha);
		algo->SetSubclusterAlpha(emAlpha);
		algo->SetVerbosity(0);
		//run clustering
		//delauney NlnN version
		if(_strategy == NlnN){
			if(i % SKIP == 0) cout << " with " << rhs.size() << " rhs" << endl;
			//start clock
			t = clock();
			trees = algo->NlnNCluster();
		}
		//N^2 version
		else if(_strategy == N2){
			if(i % SKIP == 0) cout << " with " << rhs.size() << " rhs" << endl;
			//start clock
			t = clock();
			trees = algo->N2Cluster();
		}
		t = clock() - t;
		y_time.push_back((double)t/CLOCKS_PER_SEC);
		
		FillPredJetHists(trees);
		_phos.clear();
		
	}

	//do computational time graph
	nrhs_comptime = new TGraph(_nEvts, &x_nrhs[0], &y_time[0]);
	WriteHists(ofile);

	cout << "Total number of events ran over: " << totEvt << " events that had at least two jets that passed selection: " << jetSelEff << " fraction: " << jetSelEff/totEvt << endl;
	cout << "Wrote skim to: " << _oname << endl;
}








