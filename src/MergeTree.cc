#include "MergeTree.hh"
#include "VarEMCluster.hh"
//BHC with varEM
//assuming Dirichlet Process Model (sets priors)
node* MergeTree::CalculateMerge(node *l, node* r){
	//get points from l and points from r
	//get number of points in merged tree
	//double n = l->points->GetNPoints() + r->points->GetNPoints();	
	double n = l->points->Sumw() + r->points->Sumw();	
	double d =  _alpha*tgamma(n) + (double)l->d*(double)r->d;
	double pi = _alpha*tgamma(n)/(double)d;

	double pi_a = log(_alpha) + lgamma(n);
	double pi_b = (double)log(l->d) + (double)log(r->d);
	double pi_max = fmax(pi_a,pi_b);
	double pi_stable = exp(pi_a - pi_max)/(exp(pi_a - pi_max) + exp(pi_b - pi_max));

	//d = a*gam(n) + dl*dr = a*gam(n)/dl*dr + 1
	//d - 1 = a*gam(n)/dl*dr
	//double dmin1 = log(_alpha) + lgamma(n) - log(l->d) - log(r->d);
	cpp_bin_float_100 dmin1 = log(_alpha) + lgamma(cpp_bin_float_100(n)) - log(l->d) - log(r->d);
	cpp_bin_float_100 d_stable = exp(dmin1) + 1; 
	cpp_bin_float_100 d_100 = cpp_bin_float_100(_alpha)*tgamma(cpp_bin_float_100(n)) + cpp_bin_float_100(l->d)*cpp_bin_float_100(r->d);
	cpp_bin_float_100 pi_100 = cpp_bin_float_100(_alpha)*tgamma(cpp_bin_float_100(n))/d_100;
	
	PointCollection* points = new PointCollection();
	points->AddPoints(*l->points);
	points->AddPoints(*r->points);
	struct node* x = (struct node*)malloc(sizeof *x);
	x->points = points;
	x->d = d_100;
	x->l = l;
	x->r = r;
	x->mirror = nullptr;
	x->ismirror = l->ismirror || r->ismirror;
	cout << "calcmerge ismirror " << x->ismirror << " l ismirror " << l->ismirror << " r ismirror " << r->ismirror << endl;
	double nndist = 1e300;
	//find nndist for x (should be O(n) operation)
	for(int i = 0; i < (int)_clusters.size(); i++){
		if(_clusters[i] == nullptr) continue;
		if(_euclidean_2d(l, r) < nndist) x->nndist = _euclidean_2d(l, r);
	}
	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from exp(Evidence()) = exp(ELBO) \approx exp(log(LH)) from Variational EM algorithm
	double elbo = Evidence(x);
	double p_dk_h1 = exp(elbo);
	//marginal prob of t_k = null + alterantive hypo (separate trees)
	double p_dk_tk = pi_stable*p_dk_h1 + (double)((l->d*r->d)/d)*(double)l->prob_tk*(double)r->prob_tk;
	//p_dk_tk = pi*p_dk_h1 + (dl*dr/d)*p_dr_tr*p_dl_tl 
	//double p_dk_tkmin1 = log(pi_stable) + elbo - log(l->d) - log(r->d) + log(d) - log(l->prob_tk) - log(r->prob_tk);
	double p_dk_tkmin1 = log(_alpha) + lgamma(n) + elbo - (double)log(l->d) - (double)log(r->d) - (double)log(l->prob_tk) - (double)log(r->prob_tk);
	double p_dk_tk_stable = exp(p_dk_tkmin1) + 1; 
	
	//deal with numerical instability - rk = exp(A)/(exp(A) + exp(B))
	//exp(A) = pi*p_dk_h1 = pi*exp(ELBO(x)) = exp(log(pi))exp(ELBO(x)) = exp(log(pi) + ELBO(x))
	//exp(B) = (1-pi)*p_dk_tk = (dr*dl/d)*p_dk_tk = (dr*dl/d)*p_dr_tr*p_dl_tl = exp(log(dr*dl/d) + log(p_dr_tr) + log(p_dl_tl))
	//A = log(pi) + ELBO(x)
	cpp_bin_float_100 a = log(pi_stable) + elbo;
	//B = log(dr*dl/d) + log(p_dr_tr) + log(p_dl_tl)
	cpp_bin_float_100 b = log(((l->d*r->d)/d_stable)) + log(l->prob_tk) + log(r->prob_tk);
	//find m = max(A,B)
	cpp_bin_float_100 m = fmax(a,b);
	//rewrite rk as rk = exp(A)/(exp(A) + exp(B)) = exp(A - m)/(exp(A - m) + exp(B - m))
	cpp_bin_float_100 rk_stable = exp(a - m)/(exp(a - m) + exp(b - m));

	cpp_bin_float_100 elbo_100 = cpp_bin_float_100(elbo);
	cpp_bin_float_100 p_dk_h1_100 = exp(elbo_100);
	cpp_bin_float_100 p_dl = cpp_bin_float_100(l->prob_tk);

	
	cpp_bin_float_100 p_dk_tk_100 = pi_100*p_dk_h1_100 + ((cpp_bin_float_100(l->d)*cpp_bin_float_100(r->d))/d_100)*cpp_bin_float_100(l->prob_tk)*cpp_bin_float_100(r->prob_tk);
	cpp_bin_float_100 rk_100 = pi_100*p_dk_h1_100/p_dk_tk_100;
	double rk_d100 = (double)rk_100;
	
	double rk = (double)rk_stable;//pi*p_dk_h1/p_dk_tk;
	//if(p_dk_h1 == 0 && p_dk_tk == 0) rk = 0; //the ELBO can get so negative s.t. sometimes exp(Evidence(x)) = 0 which results in rk = nan
	//else rk = pi*p_dk_h1/p_dk_tk;
	if(std::isnan(rk)){
	cout << " rk " << rk << " rk stable " << rk_stable <<  " " << (double)rk_stable << " " << (double)rk_100 << " p(D_l | T_l) " << l->prob_tk << " p(D_r | T_r) }" << r->prob_tk << " log(p_dl) " << log(l->prob_tk) << " log(p_dr) " << log(r->prob_tk) << " n " << n << " npts " << l->points->GetNPoints() + r->points->GetNPoints() << endl;
cout << "d_100 " << d_100 << " doublecast d_100 " << (double)d_100 << " pi_100 " << pi_100 << " p_dk_h1_100 " << p_dk_h1_100 << " p_dk_tk_100 " << p_dk_tk_100 << " rk_100 " << rk_100 << " rk_d100 " << rk_d100 << " p_dl_100 " << p_dl << " log(p_dl_100) " << log(p_dl) << endl;
		//cout << "evidence is 0? " << (p_dk_h1 == 0) << " p_dk_tk == 0? " << (p_dk_tk == 0) << endl;
	//cout << "rk " << rk << " gamma(n) " << tgamma(n) << " lgamma(n) " << lgamma(n) << " pi " << pi <<  " logpi " << log(pi) << " pi_stable " << pi_stable << " d " << d << " logd " << log(d) << " dr " << r->d << " logdr " << log(r->d) << " dl " << l->d << " logdl " << (l->d) << " n " << n << " dmin1 " << dmin1 << " d_stable " << d_stable << " p_dk_tkmin1 " << p_dk_tkmin1 << " p_dk_tk_stable " << p_dk_tk_stable << endl;
	}
	//if(_verb > 1){	
	//cout << " rk " << pi*p_dk_h1/p_dk_tk << " rk stable " << rk_stable << " a " << a << " b " << b << " m " << m << " pi " << pi << " p_dk_h1 " << p_dk_h1 << " p_dk_tk " << p_dk_tk << " ((l->d*r->d)/d)*p(D_l | T_l)*p(D_r | T_r) " << ((l->d*r->d)/d)*l->prob_tk*r->prob_tk << " pi*p_dk_h1 " << pi*p_dk_h1 << " p(D_l | T_l) " << l->prob_tk << " p(D_r | T_r) }" << r->prob_tk << " (l->d*r->d)/d " << (l->d*r->d)/d << endl;
	//	cout << " points " << endl; x->model->GetData()->Print(); cout << endl;
	//}	
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
