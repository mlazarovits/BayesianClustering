#include "JetSimProducer.hh"
#include "TSystem.h"
#include "RandomSample.hh"

JetSimProducer::JetSimProducer(){
	_gev = 1;
	//_minobjeta = 1.5;
	_recoptmin = 0;
	_minrhE = 0.5;
	_recoEmin = 0;
}

JetSimProducer::~JetSimProducer(){
	delete _base;
}

JetSimProducer::JetSimProducer(TFile* file){
	if(gSystem->AccessPathName(file->GetName())){ cout << "Error: file " << file->GetName() << " doesn't exist." << endl; return; }
	TTree* tree = (TTree*)file->Get("tree/llpgtree");
	_base = new ReducedBaseSim(tree);
	_nEvts = _base->fChain->GetEntries();
	//default to 1 GeV = 1 entry -> gev = 1
	_gev = 1;
	_recoptmin = 0;
	_minrhE = 0.5;
	_recoEmin = 0;
	//_minobjeta = 1.5;
	
}

void JetSimProducer::GetRecHits(vector<Jet>& rhs, int evt){
	double t, E, eta, phi;
	rhs.clear();
	double timecorr, drh, dpv, calibfactor;	

	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	int nRHs = (int)_base->ECALRecHit_energy->size();
	
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	//make weights - E/e_avg
	vector<double> ws;
	for(int r = 0; r < nRHs; r++){
		//not sure if below is needed as long as there is a clear and consistent time frame definition
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;
		
		//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
		if(_base->ECALRecHit_energy->at(r) < _minrhE) continue;
		JetPoint rh(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r),
		        _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r));
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetWeight(_base->ECALRecHit_energy->at(r)*_gev);

		Jet j(rh, vtx);
		rhs.push_back(j);
	}	
}

void JetSimProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	double t, E, eta, phi;
	rhs.clear();
	double timecorr, drh, dpv, calibfactor;	

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nRHs = (int)_base->ECALRecHit_energy->size();

	
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	//make weights - E/e_avg
	vector<double> ws;
	for(int r = 0; r < nRHs; r++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;
		
		//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
		if(_base->ECALRecHit_energy->at(r) < _minrhE) continue;
	
		JetPoint rh(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r),
		        _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r));
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetWeight(_base->ECALRecHit_energy->at(r)*_gev);

		rhs.push_back(rh);
	}
}	


void JetSimProducer::GetGenJets(vector<Jet>& genjets, int evt){
	double eta, phi, px, py, pz, pt;
	genjets.clear();

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nJets = (int)_base->Jet_genEnergy->size();
	
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	//make weights - E/e_avg
	vector<double> ws;
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		pt = _base->Jet_genPt->at(j);
		phi = _base->Jet_genPhi->at(j);
		eta = _base->Jet_genEta->at(j);


		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
		Jet jet(px, py,
		        pz, _base->Jet_genEnergy->at(j));
		
		jet.SetVertex(vtx);
		genjets.push_back(jet);
	}	
}

void JetSimProducer::GetRecoJets(vector<Jet>& recojets, int evt){
	double eta, phi, px, py, pz, pt;
	recojets.clear();

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nJets = (int)_base->Jet_energy->size();
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);

	//get rh info
	int nrhs = (int)_base->nRHs;
	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;
	vector<unsigned int> rhs;
	//make weights - E/e_avg
	vector<double> ws;
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		rhs.clear();

		pt = _base->Jet_pt->at(j);
		phi = _base->Jet_phi->at(j);
		eta = _base->Jet_eta->at(j);
		if(pt < _recoptmin) continue;
		if(_base->Jet_energy->at(j) < _recoEmin) continue;

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
		Jet jet(px, py,
		        pz, _base->Jet_energy->at(j));
		
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);

		int rhidx;
		rhs = _base->Jet_RhIDs->at(j);
		unsigned int rhid;
		double totE = 0;
		//add rhs
		if(rhs.size() < 1) continue;
		for(int r = 0; r < rhs.size(); r++){
			rhid = _base->Jet_RhIDs->at(j).at(r);
		        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;
				totE += _base->ECALRecHit_energy->at(rhidx);
				if(fabs(_base->ECALRecHit_time->at(rhidx)) > 20) continue;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), _base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx));
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*_gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				jet.AddRecHit(rh);
			}
		}		
		//cout << "jet " << j << " has " << jet.GetNRecHits() << " rhs - energy " << _base->Jet_energy->at(j) << " tot rh e " << totE << " ratio " << totE/_base->Jet_energy->at(j) << endl;
		//put cut on min n rhs (ie 2)
		recojets.push_back(jet);
	}
	//sort by pt
	if(recojets.size() > 0) SortJets(recojets);
}

void JetSimProducer::GetPrimaryVertex(BayesPoint& vtx, int evt){
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);

}

void JetSimProducer::SortJets(vector<Jet>& jets){
	vector<Jet> low;
	vector<Jet> high;
	vector<Jet> same;

	int n = jets.size();
	if(n < 2) return;
	RandomSample rs(111);
	rs.SetRange(0,n);
	int idx = rs.SampleFlat();	
	Jet pivot = jets[idx];
	for(int i = 0; i < n; i++){
		if( pivot.pt() > jets[i].pt() )
			low.push_back(jets[i]);
		else if( pivot.pt() == jets[i].pt() )
			same.push_back(jets[i]);
		else
			high.push_back(jets[i]);
	}
	SortJets(low);
	SortJets(high);
	
	jets.clear();
	for(int i = 0; i < high.size(); i++) jets.push_back(high[i]);
	for(int i = 0; i < same.size(); i++) jets.push_back(same[i]);
	for(int i = 0; i < low.size(); i++) jets.push_back(low[i]);
}
