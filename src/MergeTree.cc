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
	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from exp(Evidence()) = exp(ELBO) \approx exp(log(LH)) from Variational EM algorithm
	double p_dk_h1 = exp(Evidence(x));
	//marginal prob of t_k = null + alterantive hypo (separate trees)
	double p_dk_tk = pi*p_dk_h1 + ((l->d*r->d)/d)*l->prob_tk*r->prob_tk;	
	double rk = pi*p_dk_h1/p_dk_tk;
//cout << "pi: " << pi << " p_dk_h1: " << p_dk_h1 << " l->prob_tk: " << l->prob_tk << " l->d: " << l->d << " r->prob_tk: " << r->prob_tk << " r->d: " << r->d << " d: " << d << " p_dk_tk: " << p_dk_tk << " rk: " << rk << " points " << endl;
//points->Print();
//cout << " with val: " << rk << endl;	
	x->val = rk;
	x->prob_tk = p_dk_tk;
	
//	if(isnan(rk)) cout << "rk is nan - pi " << pi << " p_dk_h1: " << p_dk_h1 << " p_dk_tk: " << p_dk_tk << endl;
	//cout << "MergeTree::CalculateMerge - merging " << l->name << " + " << r->name << " val: " << rk << endl;	
//	x->name = "("+l->name + "+"  + r->name+")";
//cout << "MergeTree::Merge - end" << endl;

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

void MergeTree::Merge(int i, int j){
	node* n1 = Get(i);
	node* n2 = Get(j);
	Merge(n1, n2);	
}


node* MergeTree::Merge(node* l, node* r){
	node* x = CalculateMerge(l, r);
	x->l = l;
	x->r = r;
	Remove(l);
	Remove(r);
	Insert(x);
	return x;
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
	double nndist = x->nn3dist;
	if(phi >= nndist && (twopi - phi) >= nndist) return;
	if(_verb > 1) cout << "creating mirror point for point with phi center: " << phi << " with nndist: " << nndist << endl;

	//copy node x into node y so it has all the same info
	node* y = new node(*x);

	//map points across 0-2pi boundary
	_remap_phi(*y->points);

	//set x's mirror node to y and vice versa
	x->mirror = y;
	y->mirror = x;
	
	//add mirror point to list of clusters
	Insert(y);


}

void MergeTree::_remap_phi(PointCollection& points){
	double pi = acos(-1);
	double twopi = 2*pi;
	double phi = points.mean().at(1);
	double shift;
	if (phi < pi) { shift = twopi ;} else {shift = -twopi;}
	Point transl = Point({0., -shift, 0.});
	points.Translate(transl);
}
