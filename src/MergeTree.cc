#include "MergeTree.hh"
#include "VarEMCluster.hh"
//BHC with varEM
//assuming Dirichlet Process Model (sets priors)
node* MergeTree::CalculateMerge(node *l, node* r){
	//get points from l and points from r
	//get number of points in merged tree
	//double n = l->points->GetNPoints() + r->points->GetNPoints();	
	double n = l->points->Sumw() + r->points->Sumw();	

	double pi_a = log(_alpha) + lgamma(n);
	double pi_b = (double)log(l->d) + (double)log(r->d);
	double pi_max = fmax(pi_a,pi_b);
	double pi_stable = exp(pi_a - pi_max)/(exp(pi_a - pi_max) + exp(pi_b - pi_max));

	//d = alpha*gam(n) + dl*dr 
	//  = dr*dl*(alpha*gam(n)/(dl*dr) + 1)
	//  = dr*dl*(exp( log(alpha*gam(n)/(dl*dr)) ) + 1)
	//  = dr*dl*(exp( log(alpha) + lgam(n) - log(dl) - log(dr) ) + 1)
	//  = exp( log(dr) + log(dl) )*(exp( log(alpha) + lgam(n) - log(dl) - log(dr) ) + 1)
	//cpp_bin_float_100 d_100 = cpp_bin_float_100(_alpha)*tgamma(cpp_bin_float_100(n)) + cpp_bin_float_100(l->d)*cpp_bin_float_100(r->d);
	cpp_bin_float_100 d_100 = (exp( log(_alpha) + lgamma(n) - log(l->d) - log(r->d) ) + 1)*exp( log(l->d) + log(r->d)); 
	
	PointCollection* points = new PointCollection();
	points->AddPoints(*l->points);
	points->AddPoints(*r->points);
	struct node* x = (struct node*)malloc(sizeof *x);
	x->points = points;
//if(l->points->GetNPoints() + r->points->GetNPoints() != x->points->GetNPoints()){
//	cout << "mismatched pts" << endl;
//	cout << "x pts " << x->points->GetNPoints() << " l pts " << l->points->GetNPoints() << " r pts " << r->points->GetNPoints() << endl;
//	cout << "l pts" << endl; l->points->Print();
//	cout << "r pts" << endl; r->points->Print();
//	cout << "x pts" << endl; x->points->Print();
//}
	x->d = d_100;
	x->l = l;
	x->r = r;
	x->mirror = nullptr;
	x->ismirror = false;//l->ismirror || r->ismirror;
	x->idx = -1; //hasn't been added to merge tree yet
	//cout << "calcmerge ismirror " << x->ismirror << " l ismirror " << l->ismirror << " r ismirror " << r->ismirror << endl;
	double nndist = 1e300;
	//find nndist for x (should be O(n) operation)
	for(int i = 0; i < (int)_clusters.size(); i++){
		if(_clusters[i] == nullptr) continue;
		if(_euclidean_2d(l, r) < nndist) x->nndist = _euclidean_2d(l, r);
	}
	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from exp(Evidence()) = exp(ELBO) \approx exp(log(LH)) from Variational EM algorithm
	double elbo = Evidence(x);
	cpp_bin_float_100 elbo_100 = cpp_bin_float_100(elbo);
	//cpp_bin_float_100 p_dk_h1_100 = exp(elbo_100);
	cpp_bin_float_100 p_dk_h1 = exp(elbo_100);
	//marginal prob of t_k = null + alterantive hypo (separate trees) - need to save for future recursions
	//p_dk_tk = pi*p_dk_h1 + dr*dl*p_dl*p_dr/d
	// = pi*p_dk_h1(1 + dr*dl*p_dl*p_dr/(d*pi*p_dk_h1))
	// = exp(log(alpha) + 
	cpp_bin_float_100 p_dk_tk_100 = pi_stable*p_dk_h1 + ((l->d*r->d)/d_100)*l->prob_tk*r->prob_tk;


	//deal with numerical instability - rk = pi*p(Dk|H1)/p(Dk|Tk) = exp(A)/(exp(A) + exp(B))
        //exp(A) = pi*p_dk_h1 = pi*exp(ELBO(x)) = (alpha*gamma(n)/dk)*exp(ELBO) = exp(log(alpha*gamma(n)/dk))exp(ELBO(x)) = exp(log(alpha) + log(gamma(n) + ELBO(x) - log(dk))
        //exp(B) = p(Dr|Tr)*p(Dl|Tl)*dl*dr/dk = exp(log(p(Dr|Tr)*p(Dl|Tl)*dl*dr/dk)) = exp(log(p(Dr|Tr)) + log(p(Dl|Tl)) + log(dl*dr) - log(dk))
        //log(dk) appears in all terms in rk, so it can be factored out and cancelled (this is the recursively defined normalization term)
        //can rewrite exp(A) = exp(log(alpha) + log(gamma(n) + ELBO(x))
        //and rewrite exp(B) = exp(log(p(Dr|Tr)) + log(p(Dl|Tl)) + log(dl*dr))
        //where A = log(alpha) + log(gamma(n) + ELBO(x)
        cpp_bin_float_100 a = log(_alpha) + lgamma(n) + elbo_100;
        // and B = log(dr*dl/d) + log(p_dr_tr) + log(p_dl_tl)
        cpp_bin_float_100 b = log(l->d*r->d) + log(l->prob_tk) + log(r->prob_tk);
        //find m = max(A,B)
        cpp_bin_float_100 m = fmax(a,b);
        //rewrite rk as rk = exp(A)/(exp(A) + exp(B)) = exp(A - m)/(exp(A - m) + exp(B - m))
        cpp_bin_float_100 rk_stable = exp(a - m)/(exp(a - m) + exp(b - m));	


	
	//cpp_bin_float_100 p_dk_tk_100 = pi_100*p_dk_h1_100 + ((cpp_bin_float_100(l->d)*cpp_bin_float_100(r->d))/d_100)*cpp_bin_float_100(l->prob_tk)*cpp_bin_float_100(r->prob_tk);
	//cpp_bin_float_100 rk_100 = pi_100*p_dk_h1_100/p_dk_tk_100;
	//double rk_d100 = (double)rk_100;
	
	double rk = (double)rk_stable;//pi*p_dk_h1/p_dk_tk;
	//if total weight of tree is below threshold, break into separate points (ie dont merge, ie low posterior)
	//removing subclusters whose weight (ie norm) is below threshold is done within the GMM, but is not done at the BHC level
	//can put a requirement on predicted jets that # pts >= 2
	//if(x->points->Sumw() < _thresh) rk = 0;
	if(std::isnan(rk)){
        cout << std::setprecision(10) << "rk " << rk << " elbo " << elbo << " gamma(n) " << tgamma(n) << " lgamma(n) " << lgamma(n) << " pi_stable " << pi_stable << " d " << d_100 <<  " a " << a << " b " << b << " m " << m <<  " with # subclusters " << x->model->GetNClusters() << " n " << n << " npts " << l->points->GetNPoints() + r->points->GetNPoints() << endl;
                //cout << "evidence is 0? " << (p_dk_h1 == 0) << " p_dk_tk == 0? " << (p_dk_tk == 0) << endl;
        }
        //cout << std::setprecision(10) << "rk " << rk << " elbo " << elbo << " gamma(n) " << tgamma(n) << " lgamma(n) " << lgamma(n) << " dl " << l->d << " dr " << r->d << " log(dl*dr) " << log(l->d*r->d) << " p(Dl|Tl) " << l->prob_tk << " log(p(Dl|Tl)) " << log(l->prob_tk)  <<" p(Dr|Tr) " << r->prob_tk <<   " log(p(Dr|Tr)) " << log(r->prob_tk)  <<  " a " << a << " b " << b << " m " << m <<  " with # subclusters " << x->model->GetNClusters() << " n " << n << " npts " << l->points->GetNPoints() + r->points->GetNPoints() << " p_dk_tk " << p_dk_tk << " p_dk_tk100 " << p_dk_tk_100 << endl;
	double loga = elbo + log(_alpha) + lgamma(n);
	double logb = (double)log(l->prob_tk) + (double)log(r->prob_tk) + (double)log(l->d) + (double)log(r->d);
        //cout << std::setprecision(10) << "rk " << rk << " elbo " << elbo << " lgam " << lgamma(n) << " log(alpha) " << log(_alpha) << endl;
	//cout << "log(p_dl) " << log(l->prob_tk) << " log(p_dr) " << log(r->prob_tk) << " log(d_l) " << log(l->d) << " log(d_r) "  << log(r->d) << endl;
	//cout << "log(a) " << loga << " log(b) " << logb << endl;
	//cout << "log(a) - log(b) " << loga - logb << endl;
	//cout << " 3d distance from centroid " << endl;


if(_verb > 1)cout << "p_dk_tk_100 " << p_dk_tk_100 << " pi_stable " << pi_stable << " p_dk_h1 " << p_dk_h1 << " dl " << l->d << " dr " << r->d << " d " << d_100 << " p(dl) " << l->prob_tk << " p(dr) " << r->prob_tk << " n " << n << endl;
      x->val = rk;
      x->log_h1_prior = loga;
if(_verb > 1)cout << "log h1 prior = elbo " << elbo << " + log(alpha) " << log(_alpha) << " + lgam(n) " << lgamma(n) << endl;
      x->log_didj = logb;
if(_verb > 1)cout << "log didj = log(p_dl) " << log(l->prob_tk) << " " << (double)log(l->prob_tk) << " + log(p_dr) " << log(r->prob_tk) << " " << (double)log(r->prob_tk) << " + log(dl) " << log(l->d) << " " << (double)log(l->d) << " + log(dr) " << log(r->d) << " " << (double)log(r->d) << endl;
	x->prob_tk = p_dk_tk_100;
	

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
