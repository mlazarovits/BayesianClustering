#ifndef BAYESHIERCLUSTER_HH
#define BAYESHIERCLUSTER_HH

#include "PointCollection.hh"
#include "BasePDF.hh"
#include "NodeList.hh"
#include "MergeTree.hh"
#include "BaseTree.hh"

using node = BaseTree::node;
using std::vector;
class BayesHierCluster{
	public:
		BayesHierCluster();
		BayesHierCluster(double alpha);
		virtual ~BayesHierCluster(){ };
	
		void AddData(PointCollection* pc);
	

		
		//cluster trees (points in pc)
		//might need arg for recursion
		//void Cluster(const vector<PointCollection>& pc);
		vector<node*> Cluster();

		void SetAlpha(double a);
	private:
		//int m_nclusters;

		vector<PointCollection*> m_pts;
		vector<node*> _clusters;

		//calculates and tracks merges
		MergeTree m_mergeTree;

		//tracks posterior values
		NodeList _list;

		//dirichlet param
		double m_alpha;


};
#endif




