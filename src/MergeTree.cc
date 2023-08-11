#include "MergeTree.hh"
#include "VarEMCluster.hh"

//BHC with varEM
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

	struct node* x = (struct node*)malloc(sizeof *x);
	x->points = points;
	x->d = d;
	x->l = l;
	x->r = r;
	if(x->r->color != -1)
		x->color = x->r->color;
	else{
		x->color = _c;
	}
	//null hypothesis - all points in one cluster
	//calculate p(dk | null) from Evidence() = ELBO from Variational EM algorithm
	double p_dk_h1 = Evidence(x);
	//marginal prob of t_k = null + alterantive hypo (separate trees)
	double p_dk_tk = pi*p_dk_h1 + ((l->d*r->d)/d)*l->prob_tk*r->prob_tk;	

	double rk = pi*p_dk_h1/p_dk_tk;
	
	x->val = rk;
	x->prob_tk = p_dk_tk;
	
//	if(isnan(rk)) cout << "rk is nan - pi " << pi << " p_dk_h1: " << p_dk_h1 << " p_dk_tk: " << p_dk_tk << endl;
	//cout << "MergeTree::CalculateMerge - merging " << l->name << " + " << r->name << " val: " << rk << endl;	
//	x->name = "("+l->name + "+"  + r->name+")";

	return x;
}

/*
BHC only
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
	
//	if(isnan(rk)) cout << "rk is nan - pi " << pi << " p_dk_h1: " << p_dk_h1 << " p_dk_tk: " << p_dk_tk << endl;
	//cout << "MergeTree::CalculateMerge - merging " << l->name << " + " << r->name << " val: " << rk << endl;	
//	x->name = "("+l->name + "+"  + r->name+")";

	return x;
}
*/
