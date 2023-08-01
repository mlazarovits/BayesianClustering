#include "BayesHierCluster.hh"
#include "BaseTree.hh"
#include "NodeList.hh"
#include "MergeTree.hh"


using node = BaseTree::node;
BayesHierCluster::BayesHierCluster(){ m_nclusters = 999; }


BayesHierCluster::BayesHierCluster(const PointCollection* pc){
	for(int i = 0; i < pc->GetNPoints(); i++){
		m_pts.push_back(new PointCollection(pc->at(i)));
	}
	m_nclusters = (int)m_pts.size();
}

void BayesHierCluster::SetModel(BasePDF* pdf){
	m_pdf = pdf;
}

void BayesHierCluster::SetAlpha(double a){
	m_alpha = a;
}

void BayesHierCluster::Init(){
	PointCollection* pts = new PointCollection();
	for(int i = 0; i < (int)m_pts.size(); i++) pts->AddPoints(*m_pts[i]);
	m_mergeTree = MergeTree(pts);
	m_mergeTree.SetAlpha(m_alpha);
	//make sure prior within model is set
	m_mergeTree.SetModel(m_pdf);
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
		for(int j = 0; j < m_nclusters; j++){
			if( i > j ) break; //don't double count
			if(i == j) continue;

			//get subtrees i, j	
			di = m_mergeTree.Get(i);
			dj = m_mergeTree.Get(j);

			//calculate posterior for potential merge
			node *x = m_mergeTree.CalculateMerge(di, dj);

			//insert into search tree
			_list.insert(x);	
			
		}
	}
	//get max rk as top of sorted list - quicksort search tree (list) - get top value (pop)
	//merge corresponding subtrees in merge tree: merge = x (node)
	//remove subtrees from list
	//remove all combinations containing one subtree from list
	//for all nodes in merge tree:
	//if node == x: skip
	//add rk to list (insertion sort)
	//repeat from first step
 

/*
	while(m_nclusters > 1){
		rk_max = -999;
		//update current rk value for merge
		rk_max_old = rk_max;

		//add depth to history
		m_clusterHist.AddLayer(m_pts, rk_max);
		
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







