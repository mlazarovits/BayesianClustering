#include "JetSkimmer.hh"
#include "JetProducer.hh"
#include "Clusterizer.hh"
#include "BayesCluster.hh"
#include "Matrix.hh"
#include <TFile.h>
#include <time.h>
//#include <TH1D.h>
#include <TH2D.h>

JetSkimmer::JetSkimmer(){ 
	_evti = 0;
	_evtj = 0;
	_mmonly = false;
};



JetSkimmer::~JetSkimmer(){ 
}


JetSkimmer::JetSkimmer(TFile* file) : BaseSkimmer(file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	_prod = new JetProducer(file);
	_base = _prod->GetBase();
		//	_base = new ReducedBase(tree);
	_nEvts = _base->fChain->GetEntries();
	_evti = 0;
	_evtj = _nEvts;
	_oname = "plots/jet_skims_"+_cms_label+".root";
	_mmonly = false;


	//true jet hists
	hists1D.push_back(nTrueJets);
	
	//true jets pv times
	hists1D.push_back(PVtime_median);	
	hists1D.push_back(PVtime_eAvg);	
	//time differences for back-to-back jets for time resolution
	hists1D.push_back(PVtimeDiff_median);	
	hists1D.push_back(PVtimeDiff_eAvg);	

	//mm only jets
	hists1D.push_back(PVtime_mmAvg);	
	hists1D.push_back(PVtimeDiff_mmAvg);	

	//predicted jets pv times - from BHC
	hists1D.push_back(nClusters);
	hists1D.push_back(rhTime); 
	hists1D.push_back(comptime);
	hists1D.push_back(PVtime_median_pred);	
	hists1D.push_back(PVtime_eAvg_pred);	
	hists1D.push_back(PVtime_mmAvg_pred);	
	//time differences for back-to-back jets for time resolution
	hists1D.push_back(PVtimeDiff_median_pred);	
	hists1D.push_back(PVtimeDiff_eAvg_pred);	
	hists1D.push_back(PVtimeDiff_mmAvg_pred);	
	hists2D.push_back(e_nRhs);
	graphs.push_back(nrhs_comptime);


}



//make cluster param histograms
//if specified, skim from events i to j
void JetSkimmer::Skim(){
	cout << "Writing skim to: " << _oname << endl;
	if(!_mmonly){
		cout << "Using clustering strategy";
		if(_strategy == NlnN)
			cout << " NlnN (Delauney)" << endl;
		else if(_strategy == N2)
			cout << " N2 (naive)" << endl;
		else
			cout << " undefined. Please use SetStrategy(i) with i == 0 (NlnN), 1 (N2)" << endl;
	}
	else{
		cout << "Using pre-clustered AK4 jets (time calculated using MM components + naive methods)" << endl;
	}
	TFile* ofile = new TFile(_oname.c_str(),"RECREATE");
	
	//create data smear matrix - smear in eta/phi
	Matrix smear = Matrix(3,3);
	double dphi = acos(-1)/360.; //1 degree in radians
	double deta = dphi;//-log( tan(1./2) ); //pseudorap of 1 degree
	//diagonal matrix
	smear.SetEntry(deta*deta,0,0);
	smear.SetEntry(dphi*dphi,1,1);
	smear.SetEntry(1.,2,2); //no smear in time	
	
	double alpha = 0.1;
	double emAlpha = 0.5;
	double thresh = 1.;
	
	map<string, Matrix> params;
	vector<node*> trees;
	vector<Jet> rhs;
        
	//for computational time
	vector<double> x_nrhs, y_time;
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nEvts;
	}
	int SKIP = 1;	
	for(int i = _evti; i < _evtj; i++){
		_base->GetEntry(i);
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		if(i % (SKIP) == 0) cout << "evt: " << i << " of " << _nEvts;
		_prod->GetRecHits(rhs, i);
		x_nrhs.push_back((double)rhs.size());
		

		for(int r = 0; r < rhs.size(); r++){
			rhTime->Fill(rhs[r].t());
		}
		
		//fill true jet histograms
		vector<Jet> jets;
		GetTrueJets(jets);
		if(jets.size() < 1) continue;
		//cout << "\33[2K\r"<< "evt: " << i << " of " << _nEvts << " with " << rhs.size() << " rhs" << flush;
		if(i % (SKIP) == 0) cout << " with " << jets.size() << " jets to cluster";
		FillTrueJetHists(jets);
		
		if(_mmonly){
			//fill mm only jet hists
			FillMMOnlyJetHists(jets);
			cout << endl;
			continue;	
		}

		if(i % SKIP == 0) cout << " with " << rhs.size() << " rhs" << endl;

		clock_t t;
		BayesCluster* algo = new BayesCluster(rhs);
		algo->SetDataSmear(smear);
		algo->SetTimeResSmear(0.2, 0.3*_gev);
		algo->SetThresh(thresh);
		algo->SetAlpha(alpha);
		algo->SetSubclusterAlpha(emAlpha);
		algo->SetVerbosity(0);
		
		//run clustering
		//delauney NlnN version
		if(_strategy == NlnN){
			//start clock
			t = clock();
			trees = algo->NlnNCluster();
		}
		//N^2 version
		else if(_strategy == N2){
			//start clock
			t = clock();
			trees = algo->N2Cluster();
		}
		t = clock() - t;
		y_time.push_back((double)t/CLOCKS_PER_SEC);
		
		FillPredJetHists(trees);

		
	}

	//do computational time graph
	nrhs_comptime = new TGraph(_nEvts, &x_nrhs[0], &y_time[0]);

	WriteHists(ofile);
}



void JetSkimmer::GetTrueJets(vector<Jet>& jets){
	jets.clear();
	double px, py, pz;
	int nrhs, rhid;
	int njets = _base->Jet_energy->size();
	Jet jet;
	Point vtx(3);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);
	bool jetid;
	for(int j = 0; j < njets; j++){
		nrhs = _base->Jet_drRhIds->at(j).size();
		px = _base->Jet_pt->at(j)*cos( _base->Jet_phi->at(j) );
		py = _base->Jet_pt->at(j)*sin( _base->Jet_phi->at(j) );
		pz = _base->Jet_pt->at(j)*sinh( _base->Jet_eta->at(j) );
		
		//rough preselection
		if(nrhs < 2) continue;
		if(_base->Jet_pt->at(j) < 30.) continue;
		if(fabs(_base->Jet_eta->at(j)) > 1.5) continue;
		//create jet id (based on tight 2017 Run II recommendations)
		//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
		jetid = (_base->Jet_neEmEF->at(j) < 0.9) && (_base->Jet_neHEF->at(j) < 0.9) && (_base->Jet_nConstituents->at(j) > 1)
			&& (_base->Jet_chHEF->at(j) > 0.0) && (_base->Jet_chHM->at(j) > 0);
		if(!jetid) continue;  		


		jet = Jet(px, py, pz, _base->Jet_energy->at(j));
		jet.SetVertex(vtx);
		//add rec hits
		for(int r = 0; r < nrhs; r++){
			rhid = _base->Jet_drRhIds->at(j).at(r);
			for(int rr = 0; rr < _base->ECALRecHit_ID->size(); rr++){
				if(_base->ECALRecHit_ID->at(rr) != rhid) continue;
				JetPoint rh = JetPoint(_base->ECALRecHit_rhx->at(rr), _base->ECALRecHit_rhy->at(rr), 
					_base->ECALRecHit_rhz->at(rr), _base->ECALRecHit_time->at(rr) + _base->ECALRecHit_TOF->at(rr));
				rh.SetEnergy(_base->ECALRecHit_energy->at(rr));
				rh.SetEta(_base->ECALRecHit_eta->at(rr));
				rh.SetPhi(_base->ECALRecHit_phi->at(rr));
				rh.SetWeight(_base->ECALRecHit_energy->at(rr)*_gev);
				jet.AddRecHit(rh);
				jet.SetUserIdx(j);
				break;
			}
		}
		jets.push_back(jet);
			
	}
}





