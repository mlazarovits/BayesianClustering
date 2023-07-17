#ifndef JETTREE_HH
#define JETTREE_HH

#include "Jet.hh"

using std::vector;

//this class keeps track of the cluster (merge) history
class ClusterTree{
	public:
		ClusterTree();
		virtual ~ClusterTree();

		void AddLayer(vector<PointCollection> pcs);	
	
		//clusters at depth i
		vector<PointCollection> at(int i);

		//clusters at rk value
		vector<PointCollection> at(double rk);

		//cluster j at depth i
		PointCollection at(int i, int j);

		int GetNLayers();
		int GetNClusters();

		//removal functions?


	private:
		//m_tree[i] = vector<cluster> clusters at depth i -> m_tree[i][j] = cluster j at depth i
		vector<vector<PointCollection>> m_tree;
		//posterior values, index matched to point collections
		vector<double> m_rks;















};
#endif
