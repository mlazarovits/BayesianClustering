#include "MergeTree.hh"
#include "VarEMCluster.hh"

//BHC with varEM
//assuming Dirichlet Process Model (sets priors)
node* MergeTree::CalculateMerge(node *l, node* r){
	//get points from l and points from r
	//get number of points in merged tree
	double n = l->points->GetNPoints() + r->points->GetNPoints();		
	double d, pi, phidist;
	if(_constraint){
		//calculate constraint probability rho
		//rho = 1 -> merge likely, rho = 0 -> no merge likely
		double rho = DistanceConstraint(l, r)/2.; //normalize by dividing by 2 for double counting (-pi/2 to 0 put into 0 to pi/2)
		
		//don't want to change the probability of a non-merge -> no factor on d_l*d_r term
		d = rho*_alpha*tgamma(n) + l->d*r->d;
		pi = rho*_alpha*tgamma(n)/d;
	}
	else{
		d = _alpha*tgamma(n) + l->d*r->d;
		pi = _alpha*tgamma(n)/d;
	}
	PointCollection* points = new PointCollection();
	points->AddPoints(*l->points);
	points->AddPoints(*r->points);
	struct node* x = (struct node*)malloc(sizeof *x);
	x->points = points;
	x->d = d;
	x->l = l;
	x->r = r;
	x->mirror = nullptr;
	double nndist = 1e300;
	//find nndist for x (should be O(n) operation)
	for(int i = 0; i < (int)_clusters.size(); i++){
		if(_clusters[i] == nullptr) continue;
		if(_euclidean_2d(l, r) < nndist) x->nndist = _euclidean_2d(l, r);
	}
	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from exp(Evidence()) = exp(ELBO) \approx exp(log(LH)) from Variational EM algorithm
	double p_dk_h1 = exp(Evidence(x));
	//marginal prob of t_k = null + alterantive hypo (separate trees)
	double p_dk_tk = pi*p_dk_h1 + ((l->d*r->d)/d)*l->prob_tk*r->prob_tk;
		
	double rk = -1;
	if(p_dk_h1 == 0 && p_dk_tk == 0) rk = 0; //the ELBO can get so negative s.t. sometimes exp(Evidence(x)) = 0 which results in rk = nan
	else rk = pi*p_dk_h1/p_dk_tk;
	if(std::isnan(rk)){
		cout << "rk " << rk << " pi " << pi << " p_dk_h1 " << p_dk_h1 << " p_dk_tk " << p_dk_tk << " d_l " << l->d << " d_r " << r->d << " d " << d << " p(D_l | T_l) " << l->prob_tk << " p(D_r | T_r) }" << r->prob_tk << endl;
		cout << "evidence is 0? " << (p_dk_h1 == 0) << " p_dk_tk == 0? " << (p_dk_tk == 0) << endl;
	}
//		cout << " rk " << rk << " pi " << pi << " p_dk_h1 " << p_dk_h1 << " p_dk_tk " << p_dk_tk << " d_l " << l->d << " d_r " << r->d << " d " << d << " p(D_l | T_l) " << l->prob_tk << " p(D_r | T_r) }" << r->prob_tk << endl;//" points " << endl; x->model->GetData()->Print(); cout << endl;
		
	//if total weight of tree is below threshold, break into separate points (ie dont merge, ie low posterior)
	//removing subclusters whose weight (ie norm) is below threshold is done within the GMM, but is not done at the BHC level
	//can put a requirement on predicted jets that # pts >= 2
	if(x->points->Sumw() < _thresh) rk = 0;

	x->val = rk;
	x->prob_tk = p_dk_tk;
	
	return x;
}


double MergeTree::CalculateMerge(int i, int j){
	node* n1 = Get(i);
	node* n2 = Get(j);
	if(n1 == nullptr) cout << "cluster for " << i << " null" << endl;
	if(n2 == nullptr) cout << "cluster for " << j << " null" << endl;
	node* x = CalculateMerge(n1, n2);
	return x->val;	
}


//if node needs to be mirrored across 0-2pi boundary
//create new node according to below and add it to _clusters
//based on Dnn2piCylinder::_CreateNecessaryMirrorPoints from FastJet
//----------------------------------------------------------------------
///// For each plane point specified in the vector plane_indices,
///// establish whether there is a need to create a mirror point
///// according to the following criteria:
/////
///// . phi < pi
///// . mirror does not already exist
///// . phi < NearestNeighbourDistance 
/////   (if this is not true then there is no way that its mirror point
/////   could have a nearer neighbour).
/////
///// If conditions all hold, then create the mirror point, insert it
///// into the _DNN structure, adjusting any nearest neighbours, and
///// return the list of plane points whose nearest neighbours have
///// changed (this will include the new neighbours that have just been
///// added)
void MergeTree::CreateMirrorNode(node* x){
	double phi = x->points->mean().at(1);
	double twopi = 2*acos(-1);
	//require absense of mirror
	if(x->mirror != nullptr) return;


	// check that we are sufficiently close to the border --
	//// i.e. closer than nearest neighbour distance.
	//// we actually sqrt the distance this time so no need to square phi :)
	double nndist = x->nndist;
	if(_verb > 1) cout << "checking to see if mirror point is necessary for point " << x->idx << " with nndist " << nndist << " and phi " << phi << endl;
	if( phi >= nndist && (twopi - phi) >= nndist) return;
	//if( fabs(phi-(0.1)) >= nndist && ((twopi+0.1) - phi) >= nndist) return;
	if(_verb > 1) cout << "creating mirror point for point " << x->idx << "  with phi center: " << phi << " with nndist: " << nndist << endl;

	//copy node x into node y so it has all the same info
	node* y = (node*)malloc(sizeof *y);
	y->points = new PointCollection(*x->points);
	y->l = new node(*x->l);
	y->r = new node(*x->r);	
	y->val = x->val;
	y->d = x->d;
	y->model = x->model;
	y->prob_tk = x->prob_tk;
	y->nndist = x->nndist;
	//map points across 0-2pi boundary
	_remap_phi(*y->points);

	//set x's mirror node to y and vice versa
	x->mirror = y;
	y->mirror = x;

	
	//add mirror point to list of clusters
	Insert(y);

	if(_verb > 1) cout << "created mirror point " << y->idx << " with mean phi: " <<  y->points->mean().at(1) << " for point " << x->idx << " with mean phi: " << x->points->mean().at(1) << endl;
}

void MergeTree::_remap_phi(PointCollection& points){
	double pi = acos(-1);
	double twopi = 2*pi;
	double phi = points.mean().at(1);
	double shift;
	if (phi < pi) { shift = twopi ;} else {shift = -twopi;}
	BayesPoint transl = BayesPoint({0., -shift, 0.});
	points.Translate(transl);
}
