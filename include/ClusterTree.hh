#ifndef CLUSTERTREE_HH
#define CLUSTERTREE_HH

#include "Jet.hh"

using std::vector;

//this class keeps track of the cluster (merge) history
class ClusterTree{
	public:
		ClusterTree();
		virtual ~ClusterTree();
	
	private:
		//m_tree[i] = vector<jet> jets at depth i -> m_tree[i][j] = jet j at depth i
		vector<vector<Jet>> m_tree;
















};
#endif
