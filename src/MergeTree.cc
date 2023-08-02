#include "MergeTree.hh"



//assuming Dirichlet Process Model (sets priors)
node* MergeTree::CalculateMerge(node *l, node* r){
	//get points from l and points from r
	//get number of points in merged tree
	double n = l->points->GetNPoints() + r->points->GetNPoints();		
	double d = _alpha*tgamma(n) + l->d*r->d;
	double pi = _alpha*tgamma(n)/d;
	PointCollection* points = new PointCollection();
	points->AddPoints(*l->points);
	points->AddPoints(*r->points);

	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from ConjugateEvidence();
	double p_dk_h1 = _model->ConjugateEvidence(*points);
	//marginal prob of t_k = null + alterantive hypo (separate trees)
	double p_dk_tk = pi*p_dk_h1 + ((l->d*r->d)/d)*l->prob_tk*r->prob_tk;	

	double rk = pi*p_dk_h1/p_dk_tk;
	
	struct node* x = (struct node*)malloc(sizeof *x);
	x->points = points;
	x->val = rk;
	x->d = d;
	x->prob_tk = p_dk_tk;
	x->l = l;
	x->r = r;

	return x;
}
