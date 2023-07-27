#include "MergeTree.hh"



//assuming Dirichlet Process Model (sets priors)
double MergeTree::CalculateMerge(node *l, node* r, BasePDFMixture* model){
	//get points from l and points from r
	//get number of points in merged tree
	double n = l->points->GetNPoints() + r->points->GetNPoints();		
	double d = _alpha*tgamma(n) + l->d*r->d;
	double pi = _alpha*tgamma(n)/d;

	//null hypothesis - all points in one cluster
	//double nullhypo = _model->ConjugateEvidence();

	
	//calculate p_lr
	//return
}
