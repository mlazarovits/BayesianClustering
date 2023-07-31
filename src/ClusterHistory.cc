#include "ClusterHistory.hh"
#include <stdexcept>

using std::logic_error;

ClusterHistory::ClusterHistory(){ };

ClusterHistory::~ClusterHistory(){
	m_tree.clear();
	m_rks.clear();
}




void ClusterHistory::AddLayer(vector<PointCollection*> pcs, double rk){
	m_tree.push_back(pcs);
	m_rks.push_back(rk);
}

vector<PointCollection*> ClusterHistory::at(int i){
	return m_tree[i];
}

PointCollection* ClusterHistory::at(int i, int j){
	return m_tree[i][j];
}

//return clusters at depth given by posterior value
vector<PointCollection*> ClusterHistory::at(double rk){
	auto it = std::find(m_rks.begin(), m_rks.end(), rk);
	
	return m_tree[it - m_rks.begin()];
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


