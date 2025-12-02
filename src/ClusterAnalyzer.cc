#include "ClusterAnalyzer.hh"

ClusterAnalyzer::ClusterAnalyzer(){
	_gev = 1;
	_radius = 129;
	_detCenter = BayesPoint({0., 0., 0.});
	_PV = BayesPoint({0., 0., 0.});
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
	//cout << "adding rh with e " << rhE << " x " << rhx << " y " << rhy << " z " << rhz << " t " << rht << " time from 0 to rh " << d_rh << " eta " << rh.eta() << " phi " << rh.phi() << " id " << rh.rhId() << endl;
	if(invalidTime) rh.SetInvalidTime();
	Jet jet(rh, _PV);
	//cout << "jet eta " << jet.eta() << " phi " << jet.phi() << " time " << jet.t() << " w " << rh.GetWeight() << " invalid time? " << invalidTime << " " << rh.InvalidTime() << endl;	
	_rhs.push_back(jet);
}

void ClusterAnalyzer::ClearRecHitList(){
	_rhs.clear();
}


//should be run after all rechits for clustering have been added
int ClusterAnalyzer::RunClustering(ClusterObj& retobj, bool pho){
	//cout << "ClusterAnalyzer::RunClustering - start" << endl;
	std::unique_ptr<BayesCluster> algo = std::make_unique<BayesCluster>(_rhs);
	//hard coding parameters that won't change
	double cell = acos(-1)/180;
	//time resolution parameters - set by detector time resolution
	double tresCte = 0.1727;//times given in ns//0.133913 * 1e-9;
	double tresStoch = 0.5109;//1.60666 * 1e-9; 
	double tresNoise = 2.106;//0.00691415 * 1e-9;
	algo->SetMeasErrParams(cell, tresCte, tresStoch*_gev, tresNoise*_gev); 
	//set time resolution smearing
	//if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n);
	double thresh = 0.1;
	algo->SetThresh(thresh);
	double alpha = 1e-300;
	algo->SetAlpha(alpha);
	double emAlpha = 1e-5;
	algo->SetSubclusterAlpha(emAlpha);
	algo->SetVerbosity(_verb);
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
	algo->SetPriorParameters(prior_params);


	//do hierarchical clustering for subcluster constraints
	vector<std::shared_ptr<BaseTree::node>> trees;
	algo->NlnNCluster(trees);
	vector<ClusterObj> objs;
	_treesToObjs(trees, objs);
	//safety for no trees found
	if(objs.size() < 1){
		if(_verb > -1) cout << "No BHC clusters found. Returning initial ClusterObj." << endl;
	cout << "ClusterAnalyzer::RunClustering - end" << endl;
		return -1;
	}
	sort(objs.begin(), objs.end(), Esort);
	//returns lead (energy) cluster
	retobj = objs[0];
	if(_verb > -1 && _detIDmap.size() == 0){
		cout << "Warning: detIDmap not set for this ClusterAnalyzer. This map will be empty for the resulting ClusterObj." << endl;
	}
	retobj.SetupDetIDs(_detIDmap);
	retobj.SetCNNModel(_detbkgmodel.get());
	retobj.SetDNNModel(_photonidmodel.get());
	//cout << "ClusterAnalyzer::RunClustering - end" << endl;
	return 0;
}

//does not run BHC and creates a clusterobj from given rhs (this clusterobj's jet will not have any subcluster info)
int ClusterAnalyzer::NoClusterRhs(ClusterObj& retobj, bool pho){
	Jet jet = Jet(_rhs);
	retobj = ClusterObj(jet);
	if(_verb > -1 && _detIDmap.size() == 0){
		cout << "Warning: detIDmap not set for this ClusterAnalyzer. This map will be empty for the resulting ClusterObj." << endl;
	}
	retobj.SetupDetIDs(_detIDmap);
	retobj.SetCNNModel(_detbkgmodel.get());
	retobj.SetDNNModel(_photonidmodel.get());
	return 0;
}


void ClusterAnalyzer::SetCNNModel(string model){
        if(_verb > -1) cout << "Using model " << model << " for instrumental background." << endl;
        _detbkgmodel = std::make_unique<fdeep::model>(fdeep::load_model(model));
}

void ClusterAnalyzer::SetDNNModel(string model){
        if(_verb > -1) cout << "Using model " << model << " for photon ID." << endl;
        _photonidmodel = std::make_unique<fdeep::model>(fdeep::load_model(model));
}





void ClusterAnalyzer::_treesToObjs(const vector<std::shared_ptr<BaseTree::node>>& trees, vector<ClusterObj>& objs){
	cout << "ClusterAnalyzer::_treesToObjs - start" << endl;
	objs.clear();
	double x, y, z, eta, phi, t, theta, px, py, pz;
	int njets_tot = 0;
	for(int i = 0; i < trees.size(); i++){
		if(trees[i] == nullptr) continue;
		//at least 2 points (rhs)
		if(trees[i]->points->GetNPoints() < 2) continue;
		Jet predJet(trees[i].get(), _PV, _gev, _radius);
		//add Jet to jets	
		objs.push_back(ClusterObj(predJet));
	}
	cout << "ClusterAnalyzer::_treesToObjs - end" << endl;
}
