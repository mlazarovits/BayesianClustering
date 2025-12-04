#include "MergeTree.hh"
#include "VarEMCluster.hh"

double MergeTree::_log_sum_exp(vector<double>& terms){
	double m = *std::max_element(terms.begin(), terms.end());
	double sum = 0.0;
	for(double t : terms){
		sum += exp(t - m);
	}
	return m + log(sum);
}

//BHC with varEM
//assuming Dirichlet Process Model (sets priors)
std::shared_ptr<BaseTree::node> MergeTree::CalculateMerge(node *l, node* r){
	//cout << "MergeTree::CalculateMerge - start" << endl;


	clock_t t;
	t = clock();
	//get points from l and points from r
	//get number of points in merged tree
	double n = l->points->Sumw() + r->points->Sumw();	
	//d = alpha*gam(n) + dl*dr 
	//  = dr*dl*(alpha*gam(n)/(dl*dr) + 1)
	//  = dr*dl*(exp( log(alpha*gam(n)/(dl*dr)) ) + 1)
	//  = dr*dl*(exp( log(alpha) + lgam(n) - log(dl) - log(dr) ) + 1)
	//  = exp( log(dr) + log(dl) )*(exp( log(alpha) + lgam(n) - log(dl) - log(dr) ) + 1)
	vector<double> logd_terms = {log(_alpha) + lgamma(n), l->log_d + r->log_d};
	double logd_LSE = _log_sum_exp(logd_terms);
	double log_pi = log(_alpha) + lgamma(n) - logd_LSE;
	
	auto x = make_unique<node>();
	//x->points = std::make_unique<PointCollection>();
	//if l is completely mirrored in r (ie all of l's points are mirrored in r) or vice versa, do NOT do merge
	//this merge has already been done (and executed), leading to the current parent nodes
	//merges should only be considered IF l contains unique info NOT in r (or vice versa)
	PointCollection uniquepts;
	uniquepts.GetUniquePoints({*l->points, *r->points});
	if(false) cout << "unique pts - # pts " << uniquepts.GetNPoints() << " sumw " << uniquepts.Sumw() << endl;
	x->points = std::make_unique<PointCollection>(uniquepts);
	

	//x->points->AddPoints(*l->points);
	//x->points->AddPoints(*r->points);
 if(false)cout << "set pts " << x->points->GetNPoints() << " from l pts " << l->points->GetNPoints() << " and r pts " << r->points->GetNPoints() << endl;
 if(false)cout << "total sumw " << x->points->Sumw() << " l sumw " << l->points->Sumw() << " r sumw " << r->points->Sumw() << endl;
	x->l = l;
	x->r = r;

	//x->d = d_100;
	x->log_d = logd_LSE;
	double nndist = 1e300;
	x->nndist = _euclidean_2d(l, r);
	//find nndist for x (should be O(n) operation)
	//for(int i = 0; i < (int)_clusters.size(); i++){
	//	if(_clusters[i] == nullptr) continue;
	//	if(_euclidean_2d(l, r) < nndist) x->nndist = _euclidean_2d(l, r);
	//}
	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from exp(Evidence()) = exp(ELBO) \approx exp(log(LH)) from Variational EM algorithm
	


if(false)cout << "doing evidence" << endl;
	double elbo;
	if(uniquepts.GetNPoints() < 1){ //all points are mirrored between both nodes - do not merge
		elbo = -1e308;
	}
	else if(uniquepts.Sumw() == l->points->Sumw() || uniquepts.Sumw() == r->points->Sumw()){ //all points in l are mirrored in r or vice versa - do not merge (redudant info)
		elbo = -1e308;
	}
	else{
		elbo = Evidence(x.get());
	}
if(false) cout << "did evidence" << endl;

	//marginal prob of t_k = null + alterantive hypo (separate trees) - need to save for future recursions
	//p_dk_tk = pi*p_dk_h1 + dr*dl*p_dl*p_dr/d
	// = pi*p_dk_h1(1 + dr*dl*p_dl*p_dr/(d*pi*p_dk_h1))


	//calculate log_prob_dk_tk
	// = log(exp(log(a)) + exp(log(b)))
	// for log(a) = log(pi_k * exp(elbo)) = log(pi_k)*elbo = elbo*(log(alpha) + lgamma(n) - log(dk))
	//            log(dk) calculated above with LSE
	// and log(b) = log(p_dl_tl*p_dr_tr*dl*dr/dk) = log(p_dl_tl) + log(p_dr_tr) + log(dl) + log(dr) - log(dk)
	//	      log(dl), log(dr), log(p_dl_tl), log(p_dr_tr) from previous step
	//	      log(dk) calculated above

	vector<double> logpk_dk_terms = {elbo + log(_alpha) + lgamma(n) - logd_LSE, l->log_prob_tk + r->log_prob_tk + l->log_d + r->log_d - logd_LSE};
	double log_p_dk_tk_LSE = _log_sum_exp(logpk_dk_terms);

	double log_rk = log_pi + elbo - log_p_dk_tk_LSE;

	double rk;
	if(log_rk == -1e308)
		rk = -1e308;
	else
		rk = exp(log_rk);
	//if total weight of tree is below threshold, break into separate points (ie dont merge, ie low posterior)
	//removing subclusters whose weight (ie norm) is below threshold is done within the GMM, but is not done at the BHC level
	//can put a requirement on predicted jets that # pts >= 2
	if(std::isnan(rk)){
        cout << std::setprecision(10) << "log_rk " << log_rk << " rk " << rk << " elbo " << elbo << " lgamma(n) " << lgamma(n) << " log_pi " << log_pi << " logd_LSE " << logd_LSE << " with # subclusters " << x->model->GetNClusters() << " n " << n << " npts " << l->points->GetNPoints() + r->points->GetNPoints() << endl;
                //cout << "evidence is 0? " << (p_dk_h1 == 0) << " p_dk_tk == 0? " << (p_dk_tk == 0) << endl;
        }
	double loga = elbo + log(_alpha) + lgamma(n);
	double logb = l->log_prob_tk + r->log_prob_tk + l->log_d + r->log_d;
      x->val = rk;
      x->log_val = log_rk; 
      x->log_h1_prior = loga;
      x->log_didj = logb;
      x->log_prob_tk = log_p_dk_tk_LSE;
	t = clock() - t;
	_total_calcmerge_time += (double)t/CLOCKS_PER_SEC;
	_n_calcmerge_calls++;
if(_verb > 1)cout << "log didj = log(p_dl) " << l->log_prob_tk << " + log(p_dr) " << r->log_prob_tk << " + log(dl) " << l->log_d  << " + log(dr) " << r->log_d << endl;
if(_verb > 1) cout << "merge val = " << loga - logb << " for n pts " << x->points->GetNPoints() << endl;
if(_verb > 1) cout << "calculatemerge pts " << endl;
if(_verb > 1) x->points->Print();
	//cout << "MergeTree::CalculateMerge - end" << endl;
	return x;
}

/*
double MergeTree::CalculateMerge(int i, int j){
	node* n1 = Get(i);
	node* n2 = Get(j);
	if(n1 == nullptr) cout << "cluster for " << i << " null" << endl;
	if(n2 == nullptr) cout << "cluster for " << j << " null" << endl;
	auto x = CalculateMerge(n1, n2);
	return x->val;	
}
*/

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
void MergeTree::CreateMirrorNode(std::shared_ptr<node> x){
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
	std::shared_ptr<node> y = std::make_shared<node>();
	y->points = std::make_unique<PointCollection>(*x->points);
	y->l = x->l;
	y->r = x->r;	
	y->val = x->val;
	y->log_val = x->log_val;
	y->log_d = x->log_d;
	y->model = std::make_unique<GaussianMixture>(*x->model);
	y->log_prob_tk = x->log_prob_tk;
	y->nndist = x->nndist;
	//map points across 0-2pi boundary
	_remap_phi(*y->points);

	//set x's mirror node to y and vice versa
	x->mirror = y.get();
	y->mirror = x.get();

     	
	//add mirror point to list of clusters
	//Insert(y);

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
