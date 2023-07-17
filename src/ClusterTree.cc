#include "ClusterTree.hh"
#include <stdexcept>

using std::logic_error;

ClusterTree::ClusterTree(){ };

ClusterTree::~ClusterTree(){
	m_tree.clear();
}

void ClusterTree::AddLayer(vector<PointCollection> pcs){
	m_tree.push_back(pcs);
}

vector<PointCollection> ClusterTree::at(int i){
	return m_tree[i];
}

PointCollection ClusterTree::at(int i, int j){
	return m_tree[i][j];
}

//return clusters at depth given by posterior value
vector<PointCollection> ClusterTree::at(double rk){
	auto it = std::find(m_rks.begin(), m_rks.end(), rk);
	
	return m_tree[it - m_rks.begin()];
}



int ClusterTree::GetNLayers(){
	return (int)m_tree.size();
}

int ClusterTree::GetNClusters(){
	int cnt = 0;
	for(int i = 0; i < (int)m_tree.size(); i++)
		cnt += (int)m_tree[i].size();
	return cnt;
}


