#include "BayesHierCluster.hh"
#include "BaseTree.hh"

using node = BaseTree::node;


BayesHierCluster::BayesHierCluster(){ m_nclusters = 999; }


BayesHierCluster::BayesHierCluster(PointCollection* pc){
	m_pts = pc;
	m_nclusters = m_pts->GetNPoints();
}

void BayesHierCluster::SetClusterPDF(BasePDF* pdf){
	m_pdf = pdf;
}

void BayesHierCluster::SetAlpha(double a){
	m_alpha = a;
}

void BayesHierCluster::Init(){
	m_mergeTree = MergeTree(m_pts);
	m_mergeTree.SetAlpha(m_alpha);
	//make sure prior within model is set
	m_mergeTree.SetModel(m_pdf);
}

vector<PointCollection> BayesHierCluster::GetClusters(double rk){
	return m_clusterHist.at(rk);
}

vector<PointCollection> BayesHierCluster::GetClusters(int d){
	return m_clusterHist.at(d);
}


void BayesHierCluster::Cluster(){
	int n;
	node* di; node* dj;
	double rk;
	double rk_max = -999;
	//update cluster history - set layer with nodes + max rk for cut


	m_nclusters = m_pts->GetNPoints(); //needs to be updated
	while(m_nclusters > 1){
		for(int i = 0; i < m_nclusters; i++){
			for(int j = 0; j < m_nclusters; j++){
				if( i > j ) break; //don't double count
				if(i == j) continue;
				//also check if merge is already in search tree because skip if so
			
				//get subtrees i, j	
				di = m_mergeTree.Get(i);
				dj = m_mergeTree.Get(j);

				//calculate posterior for potential merge
				rk = m_mergeTree.CalculateMerge(di, dj);

				//update maximum posterior value
				if(rk > rk_max) rk_max = rk;

				//insert into search tree
				m_searchTree.insert(rk, std::make_pair(i,j));	
				
			}
		} 

		//find max posterior in search tree, get idxs of nodes to merge
		pair<int, int> idxs_merge = m_searchTree.search(rk_max);
		int merge_i = idxs_merge.first;
		int merge_j = idxs_merge.second;
		di = m_mergeTree.Get(idxs_merge.first);
		dj = m_mergeTree.Get(idxs_merge.second);
		//merge di + dj = dk, remove di + dj, insert dk
		m_mergeTree.Merge(di,dj);
		//remove all impossible merges from search tree (ie merges involving either di or dj)
		m_searchTree.del(merge_i);
		m_searchTree.del(merge_j);
		//update cluster count
		m_nclusters--;
	}


}







