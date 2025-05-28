#include "JetSimProducer.hh"
#include "TSystem.h"
#include "RandomSample.hh"

JetSimProducer::JetSimProducer(){
	_gev = 1;
	//_minobjeta = 1.5;
	_minrhE = 0.5;
	_Emin_reco = 0;
	_ptmin_reco = 0;
	_Emin_gen = 0;
	_ptmin_gen = 0;
	_minNrhs = 1;
	_nConstsmin = 0;
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
	_minrhE = 0.5;
	_Emin_reco = 0;
	_ptmin_reco = 0;
	_Emin_gen = 0;
	_ptmin_gen = 0;
	_minNrhs = 1;
	_nConstsmin = 0;
}

void JetSimProducer::GetRecHits(vector<Jet>& rhs, int evt){
	double t, E, eta, phi, x, y, z;
	rhs.clear();
	double time, timecorr, drh;	

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
		//need to correct for geometric effects in detector (ie takes longer to get to more forward areas than central ones)
		/////TOF from 0 to rh location
		time = _base->ECALRecHit_time->at(r);
		x = _base->ECALRecHit_rhx->at(r);
		y = _base->ECALRecHit_rhy->at(r);
		z = _base->ECALRecHit_rhz->at(r);
		drh = sqrt(x*x + y*y + z*z)/_c;
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;
		
		timecorr = drh;
		time = time - timecorr;

		
		//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
		if(_base->ECALRecHit_energy->at(r) < _minrhE) continue;
		if(fabs(time) > 20) continue;
		JetPoint rh(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r),
		        _base->ECALRecHit_rhz->at(r), time);
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetWeight(_base->ECALRecHit_energy->at(r)*_gev);
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));
		Jet j(rh, vtx);
		rhs.push_back(j);
	}	
}

void JetSimProducer::GetRecHits(vector<JetPoint>& rhs, int evt){
	double t, E, eta, phi, x, y, z;
	rhs.clear();
	double time, timecorr, drh;	

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
		//need to correct for geometric effects in detector (ie takes longer to get to more forward areas than central ones)
		/////TOF from 0 to rh location
		time = _base->ECALRecHit_time->at(r);
		x = _base->ECALRecHit_rhx->at(r);
		y = _base->ECALRecHit_rhy->at(r);
		z = _base->ECALRecHit_rhz->at(r);
		drh = sqrt(x*x + y*y + z*z)/_c;
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;
		
		timecorr = drh;
		time = time - timecorr;
		
		//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
		if(_base->ECALRecHit_energy->at(r) < _minrhE) continue;
		if(fabs(time) > 20) continue;
	
		JetPoint rh(_base->ECALRecHit_rhx->at(r), _base->ECALRecHit_rhy->at(r),
		        _base->ECALRecHit_rhz->at(r), time);
		
		rh.SetEnergy(_base->ECALRecHit_energy->at(r));
		rh.SetEta(_base->ECALRecHit_eta->at(r));
		rh.SetPhi(_base->ECALRecHit_phi->at(r));
		rh.SetWeight(_base->ECALRecHit_energy->at(r)*_gev);
		rh.SetRecHitId(_base->ECALRecHit_ID->at(r));

		rhs.push_back(rh);
	}
}	

void JetSimProducer::GetGenJets(vector<Jet>& genAK4jets, vector<Jet>& genAK8jets, vector<Jet>& genAK15jets, int evt){
	double eta, phi, px, py, pz, pt;
	if(evt > _nEvts) return;
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	_base->GetEntry(evt);
	
	//get AK4 jets
	genAK4jets.clear();
	int nJets = _base->AK4Jet_genNJet;
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		pt = _base->AK4Jet_genPt->at(j);
		phi = _base->AK4Jet_genPhi->at(j);
		eta = _base->AK4Jet_genEta->at(j);
		if(pt < _ptmin_gen) continue;
		if(_base->AK4Jet_genEnergy->at(j) < _Emin_gen) continue;
		//multiplicity requirement
		if(_base->AK4Jet_genNConstituents->at(j) < _nConstsmin) continue;		
		//parton-matching requirement - also serves as lepton disambiguation
		//double mindr = 999;
		//double dr = 0;
		//int bestidx = -1;
		//for(int g = 0; g < _base->genpart_ngenpart; g++){
		//	dr = dR(_base->genpart_eta->at(g), _base->genpart_phi->at(g), eta, phi);
		//	if(dr < mindr){
		//		mindr = dr;
		//		bestidx = g;
		//	}
		//}
		////set some min dR
		//if(mindr > 0.1) continue;
		////set some E ratio window
		//double Eratio = _base->Jet_genEnergy->at(j)/_base->genpart_energy->at(bestidx);
		//if(Eratio > 1.5 || Eratio < 0.5) continue;
		//cout << "\ngen jet #" << j << " has best gen match (id) " << _base->genpart_id->at(bestidx) << " with dr " << mindr << " and jet/particle energy " << _base->Jet_genEnergy->at(j)/_base->genpart_energy->at(bestidx) << " with pt " << _base->Jet_genPt->at(j) << " and # constituents " << _base->Jet_genNConstituents->at(j) << endl;
		//int genid = fabs(_base->genpart_id->at(bestidx));
		//if(find(lepIds.begin(), lepIds.end(), genid) != lepIds.end()) continue; //can't match to gen lepton

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet(px, py,
		        pz, _base->AK4Jet_genEnergy->at(j));
		//check that mass is the same
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);
		//set # constituents
		genAK4jets.push_back(jet);
	}
	//sort by pt
	if(genAK4jets.size() > 0) SortJets(genAK4jets);
	
	//get AK8 jets
	genAK8jets.clear();
	nJets = _base->AK8Jet_genNJet;
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		pt = _base->AK8Jet_genPt->at(j);
		phi = _base->AK8Jet_genPhi->at(j);
		eta = _base->AK8Jet_genEta->at(j);
		if(pt < _ptmin_gen) continue;
		if(_base->AK8Jet_genEnergy->at(j) < _Emin_gen) continue;
		//multiplicity requirement
		if(_base->AK8Jet_genNConstituents->at(j) < _nConstsmin) continue;		
		//parton-matching requirement - also serves as lepton disambiguation
		//double mindr = 999;
		//double dr = 0;
		//int bestidx = -1;
		//for(int g = 0; g < _base->genpart_ngenpart; g++){
		//	dr = dR(_base->genpart_eta->at(g), _base->genpart_phi->at(g), eta, phi);
		//	if(dr < mindr){
		//		mindr = dr;
		//		bestidx = g;
		//	}
		//}
		////set some min dR
		//if(mindr > 0.1) continue;
		////set some E ratio window
		//double Eratio = _base->Jet_genEnergy->at(j)/_base->genpart_energy->at(bestidx);
		//if(Eratio > 1.5 || Eratio < 0.5) continue;
		//cout << "\ngen jet #" << j << " has best gen match (id) " << _base->genpart_id->at(bestidx) << " with dr " << mindr << " and jet/particle energy " << _base->Jet_genEnergy->at(j)/_base->genpart_energy->at(bestidx) << " with pt " << _base->Jet_genPt->at(j) << " and # constituents " << _base->Jet_genNConstituents->at(j) << endl;
		//int genid = fabs(_base->genpart_id->at(bestidx));
		//if(find(lepIds.begin(), lepIds.end(), genid) != lepIds.end()) continue; //can't match to gen lepton

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet(px, py,
		        pz, _base->AK8Jet_genEnergy->at(j));
		//check that mass is the same
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);
		//set # constituents
		genAK8jets.push_back(jet);
	}
	//sort by pt
	if(genAK8jets.size() > 0) SortJets(genAK8jets);
	
	//get AK15 jets
	genAK15jets.clear();
	nJets = _base->AK15Jet_genNJet;
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		pt = _base->AK15Jet_genPt->at(j);
		phi = _base->AK15Jet_genPhi->at(j);
		eta = _base->AK15Jet_genEta->at(j);
		if(pt < _ptmin_gen) continue;
		if(_base->AK15Jet_genEnergy->at(j) < _Emin_gen) continue;
		//multiplicity requirement
		if(_base->AK15Jet_genNConstituents->at(j) < _nConstsmin) continue;		
		//parton-matching requirement - also serves as lepton disambiguation
		//double mindr = 999;
		//double dr = 0;
		//int bestidx = -1;
		//for(int g = 0; g < _base->genpart_ngenpart; g++){
		//	dr = dR(_base->genpart_eta->at(g), _base->genpart_phi->at(g), eta, phi);
		//	if(dr < mindr){
		//		mindr = dr;
		//		bestidx = g;
		//	}
		//}
		////set some min dR
		//if(mindr > 0.1) continue;
		////set some E ratio window
		//double Eratio = _base->Jet_genEnergy->at(j)/_base->genpart_energy->at(bestidx);
		//if(Eratio > 1.5 || Eratio < 0.5) continue;
		//cout << "\ngen jet #" << j << " has best gen match (id) " << _base->genpart_id->at(bestidx) << " with dr " << mindr << " and jet/particle energy " << _base->Jet_genEnergy->at(j)/_base->genpart_energy->at(bestidx) << " with pt " << _base->Jet_genPt->at(j) << " and # constituents " << _base->Jet_genNConstituents->at(j) << endl;
		//int genid = fabs(_base->genpart_id->at(bestidx));
		//if(find(lepIds.begin(), lepIds.end(), genid) != lepIds.end()) continue; //can't match to gen lepton

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet(px, py,
		        pz, _base->AK15Jet_genEnergy->at(j));
		//check that mass is the same
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);
		//set # constituents
		genAK15jets.push_back(jet);
	}
	//sort by pt
	if(genAK15jets.size() > 0) SortJets(genAK15jets);
}


//get gen particles - only those that couldve originated a jet
//tops, bs, quarks (prompt and from Ws) 
void JetSimProducer::GetGenParticles(vector<Jet>& genparts, int evt){
	double eta, phi, px, py, pz, pt;
	genparts.clear();
	if(evt > _nEvts) return;

	_base->GetEntry(evt);
	int nParts = _base->genpart_ngenpart;
	
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	for(int p = 0; p < nParts; p++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		pt = _base->genpart_pt->at(p);
		phi = _base->genpart_phi->at(p);
		eta = _base->genpart_eta->at(p);
		//if(pt < _ptmin) continue;
		//if(_base->Jet_genEnergy->at(j) < _Emin) continue;
		//multiplicity requirement
		//if(_base->Jet_genConstituentIdxs->at(j) < _nConstmin) continue;		

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet part(px, py,
		        pz, _base->genpart_energy->at(p));
		part.SetVertex(vtx);
		part.SetUserIdx(p); //use this to get the id from the ntuple
		//set # constituents
		genparts.push_back(part);
	}
}

void JetSimProducer::GetRecoJets(vector<Jet>& recoAK4jets, vector<Jet>& recoAK8jets, vector<Jet>& recoAK15jets, int evt){
	double t, E, eta, phi, x, y, z, time, timecorr, drh, pt, px, py, pz;	
	if(evt > _nEvts) return;
	_base->GetEntry(evt);
	BayesPoint vtx = BayesPoint(3);
	vtx.SetValue(_base->PV_x,0);
	vtx.SetValue(_base->PV_y,1);
	vtx.SetValue(_base->PV_z,2);
	//get rh info
	int nrhs = (int)_base->nRHs;
	vector<unsigned int> rhids = *_base->ECALRecHit_ID;
	vector<unsigned int>::iterator rhit;
	vector<unsigned int> rhs;

	//get AK4 jets	
	recoAK4jets.clear();
	int nJets = (int)_base->AK4Jet_energy->size();
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		rhs.clear();

		pt = _base->AK4Jet_pt->at(j);
		phi = _base->AK4Jet_phi->at(j);
		eta = _base->AK4Jet_eta->at(j);
		if(pt < _ptmin_reco) continue;
		if(_base->AK4Jet_energy->at(j) < _Emin_reco) continue;

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet(px, py,
		        pz, _base->AK4Jet_energy->at(j));
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);

		int rhidx;
		rhs = _base->AK4Jet_RhIDs->at(j);
		unsigned int rhid;
		double totE = 0;
		//add rhs
		for(int r = 0; r < rhs.size(); r++){
			rhid = _base->AK4Jet_RhIDs->at(j).at(r);
		        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;
				totE += _base->ECALRecHit_energy->at(rhidx);
				
				//need to correct for geometric effects in detector (ie takes longer to get to more forward areas than central ones)
				/////TOF from 0 to rh location
				time = _base->ECALRecHit_time->at(r);
				x = _base->ECALRecHit_rhx->at(r);
				y = _base->ECALRecHit_rhy->at(r);
				z = _base->ECALRecHit_rhz->at(r);
				drh = sqrt(x*x + y*y + z*z)/_c;
				/////TOF from PV to rh location
				///dpv = _base->ECALRecHit_pvTOF->at(r); 
				///timecorr = drh - dpv;
				
				timecorr = drh;
				time = time - timecorr;
				if(fabs(time) > 20) continue;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), _base->ECALRecHit_rhz->at(rhidx), time);
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*_gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				jet.AddRecHit(rh);
			}
		}		
		if(jet.GetNRecHits() < _minNrhs) continue;
		jet.CalculateCenter();
		jet.CalculateCovariance();
		//cout << "jet " << j << " has " << jet.GetNRecHits() << " rhs - energy " << _base->AK4Jet_energy->at(j) << " tot rh e " << totE << " ratio " << totE/_base->AK4Jet_energy->at(j) << endl;
		//put cut on min n rhs (ie 2)
		vector<JetPoint> jet_rhs = jet.GetJetPoints();
		//cout << "jet #" << recoAK4jets.size() << " has";
		//for(auto rh : jet_rhs) cout << " rh time " << rh.time() << " energy " << rh.E() << endl;
		recoAK4jets.push_back(jet);
	}
	//sort by pt
	if(recoAK4jets.size() > 0) SortJets(recoAK4jets);
	
	//get AK8 jets	
	recoAK8jets.clear();
	nJets = (int)_base->AK8Jet_energy->size();
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		rhs.clear();

		pt = _base->AK8Jet_pt->at(j);
		phi = _base->AK8Jet_phi->at(j);
		eta = _base->AK8Jet_eta->at(j);
		if(pt < _ptmin_reco) continue;
		if(_base->AK8Jet_energy->at(j) < _Emin_reco) continue;

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet(px, py,
		        pz, _base->AK8Jet_energy->at(j));
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);

		int rhidx;
		rhs = _base->AK8Jet_RhIDs->at(j);
		unsigned int rhid;
		double totE = 0;
		//add rhs
		for(int r = 0; r < rhs.size(); r++){
			rhid = _base->AK8Jet_RhIDs->at(j).at(r);
		        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;
				totE += _base->ECALRecHit_energy->at(rhidx);
				
				//need to correct for geometric effects in detector (ie takes longer to get to more forward areas than central ones)
				/////TOF from 0 to rh location
				time = _base->ECALRecHit_time->at(r);
				x = _base->ECALRecHit_rhx->at(r);
				y = _base->ECALRecHit_rhy->at(r);
				z = _base->ECALRecHit_rhz->at(r);
				drh = sqrt(x*x + y*y + z*z)/_c;
				/////TOF from PV to rh location
				///dpv = _base->ECALRecHit_pvTOF->at(r); 
				///timecorr = drh - dpv;
				
				timecorr = drh;
				time = time - timecorr;
				if(fabs(time) > 20) continue;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), _base->ECALRecHit_rhz->at(rhidx), time);
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*_gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				jet.AddRecHit(rh);
			}
		}		
		if(jet.GetNRecHits() < _minNrhs) continue;
		jet.CalculateCenter();
		jet.CalculateCovariance();
		//cout << "jet " << j << " has " << jet.GetNRecHits() << " rhs - energy " << _base->AK8Jet_energy->at(j) << " tot rh e " << totE << " ratio " << totE/_base->AK8Jet_energy->at(j) << endl;
		//put cut on min n rhs (ie 2)
		vector<JetPoint> jet_rhs = jet.GetJetPoints();
		//cout << "jet #" << recoAK8jets.size() << " has";
		//for(auto rh : jet_rhs) cout << " rh time " << rh.time() << " energy " << rh.E() << endl;
		recoAK8jets.push_back(jet);
	}
	//sort by pt
	if(recoAK8jets.size() > 0) SortJets(recoAK8jets);
	
	//get AK15 jets	
	recoAK15jets.clear();
	nJets = (int)_base->AK15Jet_energy->size();
	for(int j = 0; j < nJets; j++){
		/////TOF from 0 to rh location
		///drh = _base->ECALRecHit_0TOF->at(r);
		/////TOF from PV to rh location
		///dpv = _base->ECALRecHit_pvTOF->at(r); 
		///timecorr = drh - dpv;

		rhs.clear();

		pt = _base->AK15Jet_pt->at(j);
		phi = _base->AK15Jet_phi->at(j);
		eta = _base->AK15Jet_eta->at(j);
		if(pt < _ptmin_reco) continue;
		if(_base->AK15Jet_energy->at(j) < _Emin_reco) continue;

		px = pt*cos(phi);
		py = pt*sin(phi);
		pz = pt*sinh(eta);

		Jet jet(px, py,
		        pz, _base->AK15Jet_energy->at(j));
		jet.SetVertex(vtx);
		jet.SetUserIdx(j);

		int rhidx;
		rhs = _base->AK15Jet_RhIDs->at(j);
		unsigned int rhid;
		double totE = 0;
		//add rhs
		for(int r = 0; r < rhs.size(); r++){
			rhid = _base->AK15Jet_RhIDs->at(j).at(r);
		        rhit = std::find(rhids.begin(), rhids.end(), rhid);
                        if(rhit != rhids.end()){
                                rhidx = rhit - rhids.begin();
				if(_base->ECALRecHit_energy->at(rhidx) < _minrhE) continue;
				totE += _base->ECALRecHit_energy->at(rhidx);
				
				//need to correct for geometric effects in detector (ie takes longer to get to more forward areas than central ones)
				/////TOF from 0 to rh location
				time = _base->ECALRecHit_time->at(r);
				x = _base->ECALRecHit_rhx->at(r);
				y = _base->ECALRecHit_rhy->at(r);
				z = _base->ECALRecHit_rhz->at(r);
				drh = sqrt(x*x + y*y + z*z)/_c;
				/////TOF from PV to rh location
				///dpv = _base->ECALRecHit_pvTOF->at(r); 
				///timecorr = drh - dpv;
				
				timecorr = drh;
				time = time - timecorr;
				if(fabs(time) > 20) continue;

				//t_meas = t_raw + TOF_0^rh - TOF_pv^rh
				JetPoint rh(_base->ECALRecHit_rhx->at(rhidx), _base->ECALRecHit_rhy->at(rhidx), _base->ECALRecHit_rhz->at(rhidx), time);
				rh.SetEnergy(_base->ECALRecHit_energy->at(rhidx));
                                rh.SetEta(_base->ECALRecHit_eta->at(rhidx));
                                rh.SetPhi(_base->ECALRecHit_phi->at(rhidx));
                                rh.SetWeight(_base->ECALRecHit_energy->at(rhidx)*_gev);
                                rh.SetRecHitId(_base->ECALRecHit_ID->at(rhidx));
				jet.AddRecHit(rh);
			}
		}		
		if(jet.GetNRecHits() < _minNrhs) continue;
		jet.CalculateCenter();
		jet.CalculateCovariance();
		//cout << "jet " << j << " has " << jet.GetNRecHits() << " rhs - energy " << _base->AK15Jet_energy->at(j) << " tot rh e " << totE << " ratio " << totE/_base->AK15Jet_energy->at(j) << endl;
		//put cut on min n rhs (ie 2)
		vector<JetPoint> jet_rhs = jet.GetJetPoints();
		//cout << "jet #" << recoAK15jets.size() << " has";
		//for(auto rh : jet_rhs) cout << " rh time " << rh.time() << " energy " << rh.E() << endl;
		recoAK15jets.push_back(jet);
	}
	//sort by pt
	if(recoAK15jets.size() > 0) SortJets(recoAK15jets);
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
