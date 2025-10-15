#include "ClusterAnalyzer.hh"

ClusterAnalyzer::ClusterAnalyzer(){
	_algo = nullptr;
	_gev = 1;
	_radius = 1.29;
	_detCenter = BayesPoint({0., 0., 0.});
	_PV = BayesPoint({0., 0., 0.});
	SetupDetIDsEB();
	_verb = -1;
}


//add rechit to list of rechits to be clustered
//can be run in a for-loop, looping over all rechits
void ClusterAnalyzer::AddRecHit(double rhx, double rhy, double rhz, double rhE, double rht, int rhId, bool invalidTime){
	double dx = rhx - _detCenter.at(0);
	double dy = rhy - _detCenter.at(1);
	double dz = rhz - _detCenter.at(2);
	double d_rh = sqrt(dx*dx + dy*dy + dz*dz)/_SOL;
	rht += d_rh;	

	JetPoint rh(rhx, rhy, rhz, rht);
	rh.SetEnergy(rhE);
	rh.SetWeight(rhE*_gev);
	rh.SetRecHitId(rhId);
	if(invalidTime) rh.SetInvalidTime();
	Jet jet(rh, _PV);
	//cout << "jet eta " << jet.eta() << " phi " << jet.phi() << " time " << jet.t() << " w " << rh.GetWeight() << " invalid time? " << invalidTime << " " << rh.InvalidTime() << endl;	
	_rhs.push_back(jet);
}

void ClusterAnalyzer::ClearRecHitList(){
	_rhs.clear();
}


//should be run after all rechits for clustering have been added
ClusterObj ClusterAnalyzer::RunClustering(){
	_algo = new BayesCluster(_rhs);	
	//hard coding parameters that won't change
	double cell = acos(-1)/180;
	//time resolution parameters - set by detector time resolution
	double tresCte = 0.1727;//times given in ns//0.133913 * 1e-9;
	double tresStoch = 0.5109;//1.60666 * 1e-9; 
	double tresNoise = 2.106;//0.00691415 * 1e-9;
	_algo->SetMeasErrParams(cell, tresCte, tresStoch*_gev, tresNoise*_gev); 
	//set time resolution smearing
	//if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n);
	double thresh = 0.1;
	_algo->SetThresh(thresh);
	double alpha = 1e-300;
	_algo->SetAlpha(alpha);
	double emAlpha = 1e-5;
	_algo->SetSubclusterAlpha(emAlpha);
	_algo->SetVerbosity(_verb);
	map<string, Matrix> prior_params;
	//beta
	prior_params["scale"] = Matrix(1e-3);
	//nu
	prior_params["dof"] = Matrix(3);
	//W
	Matrix W(3,3);
	W.SetEntry(0.013333,0,0);
	W.SetEntry(0.013333,1,1);
	W.SetEntry(33.33333,2,2);
	prior_params["scalemat"] = W;
	//m
	prior_params["mean"] = Matrix(3,1);
	_algo->SetPriorParameters(prior_params);


	//do hierarchical clustering for subcluster constraints
	vector<node*> trees = _algo->NlnNCluster();
	vector<ClusterObj> objs;
	_treesToObjs(trees, objs);
	sort(objs.begin(), objs.end(), Esort);

	objs[0].SetupDetIDsEB(_detIDmap, _invDetIDmap);
	return objs[0];

}
void ClusterAnalyzer::_iEtaiPhi(JetPoint rh, int& ieta, int& iphi){
	double eta = rh.eta();
	double phi = rh.phi();

	double pi = 4*atan(1);
	double deta = 2*pi / 360;
       	int eta_max_ring = 85;
	double eta_max = eta_max_ring * deta;
	if(eta > 1.479)
		ieta = -1; //endcap
	else{
		int ieta_mag = int(floor(fabs(eta) / deta + 0.5));
		if(ieta_mag < 1) ieta_mag = 1;
		if(ieta_mag > eta_max_ring) ieta_mag = eta_max_ring;

		ieta = (eta >= 0 ? +ieta_mag : -ieta_mag);
		if(eta == 0.0) ieta = 1;	
	}
	//put phi on [0,2pi)
	double ph = fmod(phi, 2*pi);
	if(ph < 0) ph += 2*pi;
	phi = ph;
	iphi = int(floor(phi / (2*pi) * 360)) + 1;
	if(iphi < 1) iphi = 1;
	if(iphi > 360) iphi = 360;

}


void ClusterAnalyzer::_setRhIds(Jet& jet){
	vector<Jet> rhs;
	jet.GetJets(rhs);
	int nmatch = 0;
	int ieta, iphi;
	for(int r = 0; r < rhs.size(); r++){
		JetPoint rh = rhs[r].GetJetPoints()[0];
		_iEtaiPhi(rh, ieta, iphi);
		//check if iphi in map
		pair<int, int> icoords = make_pair(ieta, iphi);
		//cout << "rh #" << r << " id from ieta iphi " << _invDetIDmap.at(icoords) << " ieta " << ieta << " iphi " << iphi << " eta " << rhs[r].eta() << " phi " << rhs[r].phi() << " time " << rhs[r].t() << " energy " << rhs[r].e() << endl;
		for(int rr = 0; rr < _rhs.size(); rr++){
			//placeholder logic for now
			if(rhs[r].eta() == _rhs[rr].eta() || rhs[r].phi() == _rhs[rr].phi() || rhs[r].t() == _rhs[rr].t()){
				//set ID for rh r - add function to Jet
				JetPoint rrh = _rhs[rr].GetJetPoints()[0];
				unsigned int id = rrh.rhId();
				//cout << "id from eta, phi, time match " << id << " _rhs eta " << _rhs[rr].eta() << " _rhs phi " << _rhs[rr].phi() << " time " << _rhs[rr].t() << " energy " << _rhs[rr].e() << endl;
				jet.SetRecHitId(r, id);
				//TODO - this strategy is implemented (see above) but is off from the id provided from the placeholder logic by 10 
				//figure out how to make (eta, phi) -> (ieta, iphi) -> rechit ID from functions in KUCMSTimeCaliFiles/KUCMS_TimeCalibration.hh
				//pass (eta, phi) -> rh id map to RunClustering to assign rh ids
				//then pass rh id -> res map (as done already)
				nmatch++;
				break;
			}
		}
	}
	//cout << "found " << nmatch << " matches for " << rhs.size() << " rhs and " << _rhs.size() << " _rhs" << endl;
}

void ClusterAnalyzer::_treesToObjs(vector<node*>& trees, vector<ClusterObj>& objs){
	objs.clear();
	double x, y, z, eta, phi, t, theta, px, py, pz;
	int njets_tot = 0;
	for(int i = 0; i < trees.size(); i++){
		if(trees[i] == nullptr) continue;
		//get points from tree
		PointCollection* pc = trees[i]->points;
		//at least 2 points (rhs)
		if(pc->GetNPoints() < 2) continue;
		Jet predJet(trees[i]->model, _PV, _gev, _radius);
		//add rh ids
		_setRhIds(predJet);
		//add Jet to jets	
		objs.push_back(ClusterObj(predJet));	
	}
}
