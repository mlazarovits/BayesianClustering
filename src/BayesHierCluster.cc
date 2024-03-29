#include <math.h>
#include "BayesHierCluster.hh"
#include "TriangularPDF.hh"
#include "UniformPDF.hh"
#include "NodeStack.hh"
#include "MergeTree.hh"

BayesHierCluster::BayesHierCluster(){
	_verb = 0;
	_mergeTree = nullptr;
	_thresh = 0;
	_npts = 0;
}
		


BayesHierCluster::BayesHierCluster(double alpha){
	_mergeTree = new MergeTree(alpha);
	_verb = 0;
	_thresh = 1;
	_mergeTree->SetThresh(_thresh);
	_mergeTree->SetVerbosity(_verb);
	_npts = 0;
}


void BayesHierCluster::AddData(PointCollection* pc){
	_mergeTree->AddData(pc);
	_npts = _mergeTree->GetNClusters();
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
		
		//	//insert into search tree
			_list.insert(x);	
	
		}
	}
	vector<node*> nodes;

	//loop over possible merges
	while(_mergeTree->GetNClusters() > 1){
		if(_verb > 1) cout << "---------- iteration: " << it << " ----------" << endl;
		if(_verb > 1) _list.Print();
		//_mergeTree->GetClusters(nodes);
		//for(int i = 0; i < (int)nodes.size(); i++){
			//int kmax = nodes[i]->model->GetNClusters();
		//if(_verb > 1){	cout << "cluster " << i << " has " << kmax << " subclusters and " << nodes[i]->model->GetData()->GetNPoints() << " points - rk: " << nodes[i]->val << endl; }

		//}
		//get max rk as top of sorted list - quicksort search tree (list) - get top value (pop)
		_list.sort();
		if(_verb > 5){
		cout << "pre pop" << endl;
		_list.Print();}
		//remove all combinations containing one subtree from list
		node* max = _list.fullpop();
		if(max == nullptr) break;
		if(_verb > 2) cout << "max merge rk: " << max->val << " " << endl;

		//if rk < 0.5: cut tree
		double maxval = 0.5;
		if(max->val < maxval){
			if(_verb > 0) cout << "reached min rk = " << max->val << " <  " << maxval << " - final iteration: " << it <<  " - " << _mergeTree->GetNClusters() << " clusters" << endl;
			break;
		}


		if(_verb > 5){
		cout << "post pop" << endl;
		_list.Print();
		}
		//merge corresponding subtrees in merge tree: merge = x (node)
		if(_verb > 5){
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
	
		if(_verb > 4){	
		cout << "n clusters in merge tree: " << _mergeTree->GetNClusters() << endl;
		cout << "max merge has: " << max->points->GetNPoints() << " points" << endl;
		max->points->Print();
		}	
		
		if(max->points->GetNPoints() == _npts){
			cout << "All points clustered. Exiting clustering." << endl;
			break;
		} 

		//for all nodes in merge tree: check against newly formed cluster
		NodeStack _list1;
		for(int i = 0; i < _mergeTree->GetNClusters(); i++){
			di = _mergeTree->Get(i);
			if(di == max) continue;
			if(di == nullptr) continue;
			if(_verb > 4){ cout << "calculating merge for max and " << endl; di->points->Print();}
			node* x = _mergeTree->CalculateMerge(di, max);
			if(isnan(x->val)){  return _mergeTree->GetClusters(); }
			//modify probability of merge (rk) by distance constraint if specified
			
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







