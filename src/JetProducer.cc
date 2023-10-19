#include "JetProducer.hh"

#include "Clusterizer.hh"
#include "Matrix.hh"
#include <TFile.h>
//#include <TH1D.h>
#include <TH2D.h>



JetProducer::JetProducer(){ 
}



JetProducer::~JetProducer(){ 
	_file->Close();
	delete _base;
	delete _file;
}


JetProducer::JetProducer(TFile* file) : BaseProducer(file){
	//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
	//tof = (d_rh-d_pv)/c
	//in ntuplizer, stored as rh time

	//grab rec hit values
	//x, y, z, time (adjusted), energy, phi, eta
	//_file = file;
	//cout << "b" << endl;
	//cout << "file name: " << file->GetName() << endl;
	//TTree* tree = (TTree*)file->Get("tree/llpgtree");
	//cout << "c" << endl;
	//_base = new ReducedBase(tree);
	//cout << "d" << endl;
	//_nEvts = _base->fChain->GetEntries();

}

void JetProducer::GetRecHits(vector<vector<JetPoint>>& rhs){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned long id;
	int nRHs;
	rhs.clear();
	for(int i = 0; i < _nEvts; i++){
		_base->GetEntry(i);
		//TODO: switch to _base->nRHs when that's in the ntuples
		nRHs = (int)_base->ECALRecHit_ID->size();
		rhs.push_back({});
		for(int r = 0; r < nRHs; r++){
			//add tof = d_pv to time to get correct RH time
			//t = rh_time - d_rh/c + d_pv/c = tadj + TOF
			rh = JetPoint(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
			
			rh.SetEnergy(_base->ECALRecHit_energy->at(r));
			rh.SetEta(_base->ECALRecHit_eta->at(r));
			rh.SetPhi(_base->ECALRecHit_phi->at(r));
			rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
	
			rhs[i].push_back(rh);
		}
	}
}


void JetProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	JetPoint rh;
	double x, y, z, t, E, eta, phi;
	unsigned int rhId;
	int nRHs;
	rhs.clear();
	

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nJets = (int)_base->Jet_energy->size();
	//cout << "# jets in evt " << evt << ": " << nJets << endl;

	//int nJetRHs;
	//int nJetRHs_total = 0;
	//for(int j = 0; j < nJets; j++){
	//	nJetRHs = (int)_base->Jet_drRhIds->at(j).size();
	//	//cout << "jet #" << j << " has " << nJetRHs << " rhs - eta = " << fabs(_base->Jet_eta->at(j)) << endl;
	//	if(fabs(_base->Jet_eta->at(j)) > 1.479) continue;
	//	nJetRHs_total += _base->Jet_drRhIds->at(j).size(); 
	//}
//		cout << nJetRHs_total << " total rhs in all jets" << endl;
	//actually get rhs for clustering
	nRHs = (int)_base->ECALRecHit_ID->size();
	for(int r = 0; r < nRHs; r++){
		//clean out photon ids - if rh id is in phoIds continue
		rhId = _base->ECALRecHit_ID->at(r);
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		rh = JetPoint(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
	
		rhs.push_back(rh);
	}	
}




void JetProducer::GetRecHits(vector<Jet>& jets, int evt){
	JetPoint rh;
	Jet j;
	double t, E, eta, phi;
	unsigned int rhId;
	int nRHs;
	jets.clear();

	if(evt > _nEvts) return;

	_base->GetEntry(evt);

	Point vtx = Point(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);

	//actually get rhs for clustering
	nRHs = (int)_base->ECALRecHit_ID->size();
	for(int r = 0; r < nRHs; r++){
		//clean out photon ids - if rh id is in phoIds continue
		rhId = _base->ECALRecHit_ID->at(r);
		//add tof = d_pv to time to get correct RH time
		//t = rh_time - d_rh/c + d_pv/c
		rh = JetPoint(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r), _base->ECALRecHit_rhz->at(r), _base->ECALRecHit_time->at(r)+_base->ECALRecHit_TOF->at(r));
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));

		j = Jet(rh, vtx);
		jets.push_back(j);
	}	
}


void JetProducer::GetPrimaryVertex(Point& vtx, int evt){
	//reset to empty 3-dim point	
	vtx = Point(3);

	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	vtx.SetValue(_base->PV_x, 0);
	vtx.SetValue(_base->PV_y, 1);
	vtx.SetValue(_base->PV_z, 2);

}


void JetProducer::GetTrueJets(vector<Jet>& jets, int evt){
	double px, py, pz, pt, phi, eta;
	jets.clear();

	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nJets = (int)_base->Jet_energy->size();
	//cout << "# jets in evt " << evt << ": " << nJets << endl;
	int nJetRHs;
	int nJetRHs_total = 0;
	//for dR matched rec hits
	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;
	int rhidx;
	for(int j = 0; j < nJets; j++){
		pt = _base->Jet_pt->at(j);
		phi = _base->Jet_phi->at(j);
		eta = _base->Jet_eta->at(j);

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet = Jet(px, py, pz, _base->Jet_energy->at(j));	
		
		//set rec hits in jet
		vector<unsigned int> rhs = _base->Jet_drRhIds->at(j);
		vector<unsigned int> egidxs = _base->Jet_egIndxs->at(j);
		for(int r = 0; r < rhs.size(); r++){
			unsigned int rhid = rhs[r];
			rhit = std::find(rhids.begin(), rhids.end(), rhid);
			if(rhit != rhids.end()){
				rhidx = rhit - rhids.begin();
				JetPoint rh = JetPoint(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), 
					_base->ECALRecHit_rhz->at(rhidx), _base->ECALRecHit_time->at(rhidx));
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
				jet.AddRecHit(rh);
			}

		}
		if(jet.GetNConstituents() < 1) continue;
		jets.push_back(jet);
	}


}


