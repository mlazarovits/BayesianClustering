#include <math.h>
#include "BayesHierCluster.hh"
#include "TriangularPDF.hh"
#include "UniformPDF.hh"
#include "NodeStack.hh"
#include "MergeTree.hh"

BayesHierCluster::BayesHierCluster(){  _verb = 0; _mergeTree = nullptr; _thresh = 0;  _constraint_min = -999; _constraint_max = -999; _constrain = false; _wraparound = false;}
		


BayesHierCluster::BayesHierCluster(double alpha){
	_mergeTree = new MergeTree(alpha);
	_verb = 0;
	_thresh = 1;
	_mergeTree->SetThresh(_thresh);
	_mergeTree->SetVerbosity(_verb);
	_constrain = false;
	_wraparound = false;
	 _constraint_min = -999; _constraint_max = -999;
}


void BayesHierCluster::AddData(PointCollection* pc){
	_mergeTree->AddData(pc);
}



//uses a triangular distribution to constraint certain clusterings
double BayesHierCluster::DistanceConstraint(node* i, node* j){
	double c = (_constraint_a + _constraint_b)/2.;
	double pi = acos(-1);
	double d;
	//phi
	double cent1 = i->points->Centroid(1);
	double cent2 = j->points->Centroid(1);
	TriangularPDF* tri = new TriangularPDF(_constraint_a,_constraint_b,c);
	
	//transform deltaPhi to be on [0,pi], wrapped s.t. 0 is close to 2pi (-3 close to 3)
	d = fabs(cent1 - cent2);
	if(_wraparound){
		if(d > pi) d = 2*pi - d; 
	}	

	double phi = 0;
	double theta = 0;
	if(d >= _constraint_a && d <= _constraint_b) phi = cos(d);
	
	//eta to theta
	cent1 = i->points->Centroid(0); 
	cent2 = j->points->Centroid(0); 

	cent1 = 2*atan(exp(-cent1));
	cent2 = 2*atan(exp(-cent2));
	//don't need to wrap eta -> only goes from 0 to pi in theta
	d = fabs(cent1 - cent2);	
	if(d >= _constraint_a && d <= _constraint_b) theta = cos(d);

	return theta*phi;//tri->Prob(d)/tri->Prob(c);

}

vector<node*> BayesHierCluster::Cluster(){
	int n;
	node* di; node* dj;
	double rk;
	double rk_max;
	//update cluster history - set layer with nodes + max rk for cut

//	m_nclusters = (int)m_pts.size(); //needs to be updated

	int it = 0;

	//while(_mergeTree->GetNPoints() > 1){
	int m_npts = _mergeTree->GetNClusters();
	//construct rk list (NodeStack)
	//start algorithm with all combinations - O(n^2)
	for(int i = 0; i < m_npts; i++){
		for(int j = i; j < m_npts; j++){
			if(i == j) continue;
			//get subtrees i, j	
			di = _mergeTree->Get(i);
			dj = _mergeTree->Get(j);
//		cout << "i: " << i << " npts: " << di->points->GetNPoints() << " j: " << j << " npts: " << dj->points->GetNPoints() << endl;
			//calculate posterior for potential merge
			node *x = _mergeTree->CalculateMerge(di, dj);
//		cout << "potential merge x has " << x->points->GetNPoints() << " points" << endl; x->points->Print();	
			//modify probability of merge (rk) by distance constraint if specified
			if(_constrain){
				//cout << "distance constraining" << endl;
				x->val *= DistanceConstraint(di,dj);		
		
			}	
		
		//	//insert into search tree
			_list.insert(x);	
	
		}
	}
	vector<node*> nodes;

	//loop over possible merges
	while(_mergeTree->GetNClusters() > 1){
		if(_verb > 1) cout << "---------- iteration: " << it << " ----------" << endl;
		nodes = _mergeTree->GetClusters();
		for(int i = 0; i < (int)nodes.size(); i++){
			int kmax = nodes[i]->model->GetNClusters();
		if(_verb > 1){	cout << "cluster " << i << " has " << kmax << " subclusters and " << nodes[i]->model->GetData()->GetNPoints() << " points - rk: " << nodes[i]->val << endl; }

		}
		//get max rk as top of sorted list - quicksort search tree (list) - get top value (pop)
		_list.sort();
	//	if(_verb > 0){
	//	cout << "pre pop" << endl;
	//	_list.Print();}
		//remove all combinations containing one subtree from list
		node* max = _list.fullpop();
		if(max == nullptr) break;
		//cout << "max merge rk: " << max->val << endl;

		//if rk < 0.5: cut tree
		double maxval = 0.5;
		if(max->val < maxval){
			if(_verb > 0) cout << "reached min rk = " << max->val << " <  " << maxval << " - final iteration: " << it <<  " - " << _mergeTree->GetNClusters() << " clusters" << endl;
			break;
		}
		//cout << "post pop" << endl;
		//_list.Print();
		//merge corresponding subtrees in merge tree: merge = x (node)
		if(_verb > 4){
		cout << "merging clusters" << endl;
		max->points->Print();	
		cout << "removing cluster - left" << endl;
		max->l->points->Print();
		cout << "removing cluster - right" << endl;
		max->r->points->Print();
		
			cout << "left - left subtree val " << max->l->l->val << endl;
			if(max->l->l->val != -1) max->l->l->points->Print();
			cout << "left - right subtree val " << max->l->r->val << endl;
			if(max->l->r->val != -1) max->l->r->points->Print();
			cout << "right - left subtree val " << max->r->l->val << endl;
			if(max->r->l->val != -1) max->r->l->points->Print();
			cout << "right - right subtree val " << max->r->r->val << endl;
			if(max->r->r->val != -1) max->r->r->points->Print();
		}
		_mergeTree->Insert(max);
		_mergeTree->Remove(max->l);
		_mergeTree->Remove(max->r);
		
	//	cout << "n clusters in merge tree: " << _mergeTree->GetNClusters() << endl;
	//	cout << "max merge has: " << max->points->GetNPoints() << " points" << endl;
	//	max->points->Print();
		
		//for all nodes in merge tree: check against newly formed cluster
		NodeStack _list1;
		for(int i = 0; i < _mergeTree->GetNClusters(); i++){
			di = _mergeTree->Get(i);
			if(di == max) continue;
			//cout << "calculating merge for max and " << endl; di->points->Print();
			node* x = _mergeTree->CalculateMerge(di, max);
			if(isnan(x->val)){  return _mergeTree->GetClusters(); }
			//modify probability of merge (rk) by distance constraint if specified
			if(_constrain){
				x->val *= DistanceConstraint(di,dj);		
			}	
			
			_list1.insert(x);
		}
		//cout << "list1" << endl;
		_list1.sort();
		//_list1.Print();
		_list.merge(_list1);
		//cout << "merging" << endl;
		//_list.Print();
		it++;
		//cout << "\n" << endl;
	}
	
	if(_verb > 0) cout << "---------- final iteration: " << it <<  " - " << _mergeTree->GetNClusters() << " clusters ----------" << endl;
	nodes = _mergeTree->GetClusters();
	for(int i = 0; i < (int)nodes.size(); i++){
			int kmax = nodes[i]->model->GetNClusters();
			if(_verb > 0) cout << "cluster " << i << " has " << kmax << " subclusters - rk: " << nodes[i]->val << " and " << nodes[i]->points->GetNPoints() << " points " << endl;	

	}

	return _mergeTree->GetClusters();
}







