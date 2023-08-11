#include "BayesHierCluster.hh"
#include "BaseTree.hh"
#include "NodeList.hh"
#include "MergeTree.hh"


using node = BaseTree::node;
BayesHierCluster::BayesHierCluster(){  }


BayesHierCluster::BayesHierCluster(double alpha){
	m_mergeTree = MergeTree(alpha);
}

void BayesHierCluster::AddData(PointCollection* pc){
	m_mergeTree.AddData(pc);
}

void BayesHierCluster::SetAlpha(double a){
	m_alpha = a;
	m_mergeTree.SetAlpha(m_alpha);
}


vector<node*> BayesHierCluster::Cluster(){
	cout << "BHC::Cluster" << endl;
	int n;
	node* di; node* dj;
	double rk;
	double rk_max;
	//update cluster history - set layer with nodes + max rk for cut

//	m_nclusters = (int)m_pts.size(); //needs to be updated

	int it = 0;

	//while(m_mergeTree.GetNPoints() > 1){
	int m_npts = m_mergeTree.GetNClusters();
	//construct rk list (NodeList)
	//start algorithm with all combinations - O(n^2)
	for(int i = 0; i < m_npts; i++){
		for(int j = i; j < m_npts; j++){
			if(i == j) continue;
			
			//get subtrees i, j	
			di = m_mergeTree.Get(i);
			dj = m_mergeTree.Get(j);
			//calculate posterior for potential merge
			node *x = m_mergeTree.CalculateMerge(di, dj);
		//	//insert into search tree
			_list.insert(x);	
		}
	}
cout << "Made first r_k list" << endl;
	vector<node*> nodes;

	//loop over possible merges
	while(m_mergeTree.GetNClusters() > 1){
		cout << "---------- iteration: " << it << " ----------" << endl;
		nodes = m_mergeTree.GetClusters();
		for(int i = 0; i < (int)nodes.size(); i++){
			int kmax = nodes[i]->model->GetNClusters();
			cout << "cluster " << i << " has " << kmax << " subclusters and " << nodes[i]->model->GetData()->GetNPoints() << " points - rk: " << nodes[i]->val << endl;
			map<string, Matrix> params;
		//	for(int k = 0; k < kmax; k++){
		//	}

		}


		//get max rk as top of sorted list - quicksort search tree (list) - get top value (pop)
		_list.sort();
		//cout << "pre pop" << endl;
		//_list.Print();
		//remove all combinations containing one subtree from list
		node* max = _list.fullpop();
		//cout << "max merge rk: " << max->val << endl;
		//if rk < 0.5: cut tree
		if(max->val < 0.5){
			cout << "reached min rk = 0.5 - final iteration: " << it <<  " - " << m_mergeTree.GetNClusters() << " clusters" << endl;
			break;
		}
		//cout << "post pop" << endl;
		//_list.Print();
		//merge corresponding subtrees in merge tree: merge = x (node)
		m_mergeTree.Insert(max);
		m_mergeTree.Remove(max->l);
		m_mergeTree.Remove(max->r);
		//cout << "n clusters in merge tree: " << m_mergeTree.GetNClusters() << endl;
		
		//for all nodes in merge tree: check against newly formed cluster
		NodeList _list1;
		for(int i = 0; i < m_mergeTree.GetNClusters(); i++){
			di = m_mergeTree.Get(i);
			if(di == max) continue;
			node* x = m_mergeTree.CalculateMerge(di, max);
			if(isnan(x->val)){  return m_mergeTree.GetClusters(); }
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
	
	cout << "---------- final iteration: " << it <<  " - " << m_mergeTree.GetNClusters() << " clusters ----------" << endl;
	nodes = m_mergeTree.GetClusters();
	for(int i = 0; i < (int)nodes.size(); i++){
			int kmax = nodes[i]->model->GetNClusters();
			cout << "cluster " << i << " has " << kmax << " subclusters - rk: " << nodes[i]->val << endl;
			map<string, Matrix> params;
		//	for(int k = 0; k < kmax; k++){
		//		params = nodes[i]->model->GetParameters(k);
		//		cout << "mean " << k << endl;
		//		params["mean"].Print();
		//		cout << "cov " << k << endl;
		//		params["cov"].Print();
		//	}

	}

	return m_mergeTree.GetClusters();
}







