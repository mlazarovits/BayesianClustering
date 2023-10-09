#include "MergeTree.hh"
#include "VarEMCluster.hh"

//BHC with varEM
//assuming Dirichlet Process Model (sets priors)
node* MergeTree::CalculateMerge(node *l, node* r){
	//get points from l and points from r
	//get number of points in merged tree
	double n = l->points->GetNPoints() + r->points->GetNPoints();		
	double d, pi;
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
	
	x->val = rk;
	x->prob_tk = p_dk_tk;
	
//	if(isnan(rk)) cout << "rk is nan - pi " << pi << " p_dk_h1: " << p_dk_h1 << " p_dk_tk: " << p_dk_tk << endl;
	//cout << "MergeTree::CalculateMerge - merging " << l->name << " + " << r->name << " val: " << rk << endl;	
//	x->name = "("+l->name + "+"  + r->name+")";

	return x;
}


node* MergeTree::Merge(node* l, node* r){
	node* x = CalculateMerge(l, r);
	x->l = l;
	x->r = r;
	cout << "Merge - l pts: " << x->l->points->GetNPoints() << " r pts: " << x->r->points->GetNPoints() << endl;
	Remove(l);
	Remove(r);
	Insert(x);
	return x;
}
