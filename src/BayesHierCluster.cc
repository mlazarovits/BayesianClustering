#include "BayesHierCluster.hh"
#include "BaseTree.hh"
#include "NodeList.hh"
#include "MergeTree.hh"


using node = BaseTree::node;
BayesHierCluster::BayesHierCluster(){ m_nclusters = 999; }


BayesHierCluster::BayesHierCluster(BasePDF* model){
	m_pdf = model;
	m_mergeTree = MergeTree(m_pdf);
}

void BayesHierCluster::AddData(PointCollection* pc){
	for(int i = 0; i < pc->GetNPoints(); i++){
		m_pts.push_back(new PointCollection(pc->at(i)));
	}
	m_nclusters = (int)m_pts.size();
	m_mergeTree.AddData(pc);
}

void BayesHierCluster::SetAlpha(double a){
	m_alpha = a;
	m_mergeTree.SetAlpha(m_alpha);
}

vector<PointCollection*> BayesHierCluster::GetClusters(double rk){
	return m_clusterHist.at(rk);
}

vector<PointCollection*> BayesHierCluster::GetClusters(int d){
	return m_clusterHist.at(d);
}


void BayesHierCluster::Cluster(){
	int n;
	node* di; node* dj;
	double rk;
	double rk_max;
	double rk_max_old = -999;
	node* merge_node;
	pair<int, int> idxs;
	//update cluster history - set layer with nodes + max rk for cut

	m_nclusters = (int)m_pts.size(); //needs to be updated


	//construct rk list (NodeList)
	//start algorithm with all combinations - O(n^2)
	for(int i = 0; i < m_nclusters; i++){
		for(int j = i; j < m_nclusters; j++){
			if(i == j) continue;

			//get subtrees i, j	
			di = m_mergeTree.Get(i);
			dj = m_mergeTree.Get(j);
			//calculate posterior for potential merge
			node *x = m_mergeTree.CalculateMerge(di, dj);
			cout << "i: " << i << " j: " << j << " r_ij: " << x->val << " points" << endl;
		//	//insert into search tree
			_list.insert(x);	
		}
	}
	//get max rk as top of sorted list - quicksort search tree (list) - get top value (pop)
	_list.sort();
cout << "pre pop" << endl;
_list.Print();
	//remove all combinations containing one subtree from list
	node* max = _list.fullpop();
cout << "max merge rk: " << max->val << endl;
cout << "post pop" << endl;
_list.Print();
	//merge corresponding subtrees in merge tree: merge = x (node)
	m_mergeTree.Insert(max);
	m_mergeTree.Remove(max->l);
	m_mergeTree.Remove(max->r);
		
	//add depth to history
//	m_clusterHist.AddLayer(m_pts, rk_max);
	
	//for all nodes in merge tree:
	NodeList _list1;
	for(int i = 0; i < m_mergeTree.GetNPoints(); i++){
		di = m_mergeTree.Get(i);
		if(di == max) continue;
		node* x = m_mergeTree.CalculateMerge(di, max);
		_list1.insert(x);
	}
cout << "list1" << endl;
	_list1.sort();
	_list1.Print();
	_list.merge(_list1);
	cout << "merging" << endl;
	_list.Print();
	return;

/*
	while(m_nclusters > 1){
		rk_max = -999;
		//update current rk value for merge
		rk_max_old = rk_max;

		
		//find max posterior in search tree, get idxs of nodes to merge
		pair<int, int> idxs_merge = m_searchTree.search(rk_max);
		int merge_i = idxs_merge.first;
		int merge_j = idxs_merge.second;
		di = m_mergeTree.Get(idxs_merge.first);
		dj = m_mergeTree.Get(idxs_merge.second);
		//merge di + dj = dk, remove di + dj, insert dk
		merge_node = m_mergeTree.Merge(di,dj);
		//remove all impossible merges from search tree (ie merges involving either di or dj)
		m_searchTree.del(merge_i);
		m_searchTree.del(merge_j);
		
		//remove all impossible merges from m_pts (ie merges involving either di or dj)
		for(int i = 0; i < m_nclusters; i++){
			if(i == merge_i || i == merge_j)
				m_pts.erase(m_pts.begin(),m_pts.begin()+i);
		}
		m_pts.push_back(merge_node->points);
		
		//update cluster count
		m_nclusters--; //should equal m_pts.size
	}

*/
}







