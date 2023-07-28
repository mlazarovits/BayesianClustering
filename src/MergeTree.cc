#include "MergeTree.hh"



//assuming Dirichlet Process Model (sets priors)
double MergeTree::CalculateMerge(node *l, node* r){
	//get points from l and points from r
	//get number of points in merged tree
	double n = l->points->GetNPoints() + r->points->GetNPoints();		
	double d = _alpha*tgamma(n) + l->d*r->d;
	double pi = _alpha*tgamma(n)/d;
	PointCollection* points;
	points->add(*l->points);
	points->add(*r->points);

	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from ConjugateEvidence();
	double p_dk_h1 = _model->ConjugateEvidence(*points);
	//marginal prob of t_k = null + alterantive hypo (separate trees)
	double p_dk_tk = pi*p_dk_h1 + ((l->d*r->d)/d)*l->prob_tk*r->prob_tk;	

	double r_k = pi*p_dk_h1/p_dk_tk;
		
	return r_k;
}
