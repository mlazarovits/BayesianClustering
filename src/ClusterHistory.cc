#include "ClusterHistory.hh"
#include <stdexcept>

using std::logic_error;

ClusterHistory::ClusterHistory(){ };

ClusterHistory::~ClusterHistory(){
	m_tree.clear();
	m_rks.clear();
}




void ClusterHistory::AddLayer(vector<node*> pcs, double rk){
	m_tree.push_back(pcs);
	m_rks.push_back(rk);
}

vector<node*> ClusterHistory::at(int i){
	return m_tree[i];
}

node* ClusterHistory::at(int i, int j){
	return m_tree[i][j];
}

//return clusters at depth given by posterior value
vector<node*> ClusterHistory::at(double rk){
	//define
	vector<double>::iterator it = m_rks.begin();
	int idx;
	//find closest vector element to given rk (rk must be <= it)
	while(rk > *it){
		it++;		
	}
	idx = it - m_rks.begin();
	return m_tree[idx];

}



int ClusterHistory::GetNLayers(){
	return (int)m_tree.size();
}

int ClusterHistory::GetNClusters(int layer){
	return (int)m_tree[layer].size();
}
int ClusterHistory::GetNClusters(){
	int cnt = 0;
	for(int i = 0; i < (int)m_tree.size(); i++)
		cnt += (int)m_tree[i].size();
	return cnt;
}


