#ifndef ClusterHistory_HH
#define ClusterHistory_HH

#include "SearchTree.hh"
#include "Jet.hh"
#include <vector>


using std::vector;
using node = SearchTree::node;

//this class keeps track of the cluster (merge) history
class ClusterHistory{
	public:
		ClusterHistory();
		virtual ~ClusterHistory();


		//init tree stack with leaf nodes (connected to null, external termination nodes)
		void Init();

		void AddLayer(vector<PointCollection*> pcs, double rk);	
	
		//clusters at depth i
		vector<PointCollection*> at(int i);

		//clusters at rk value
		vector<PointCollection*> at(double rk);

		//cluster j at depth i
		PointCollection* at(int i, int j);

		int GetNLayers();
		int GetNClusters();
		int GetNClusters(int layer);


		void clear(){
			m_tree.clear();
			m_rks.clear();
		}

	private:
		//m_tree[i] = vector<cluster> clusters at depth i -> m_tree[i][j] = cluster j at depth i
		vector<vector<PointCollection*>> m_tree;
		//posterior values, index matched to point collections
		vector<double> m_rks;











};
#endif