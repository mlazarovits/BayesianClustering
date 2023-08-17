#ifndef BAYESHIERCLUSTER_HH
#define BAYESHIERCLUSTER_HH

#include "PointCollection.hh"
#include "BasePDF.hh"
#include "NodeStack.hh"
#include "MergeTree.hh"
//#include "BaseTree.hh"

//using node = BaseTree::node;
using std::vector;
class BayesHierCluster{
	public:
		BayesHierCluster();
		BayesHierCluster(double alpha);
		virtual ~BayesHierCluster(){ };
	
		void AddData(PointCollection* pc);
	
		void SetVerbosity(int verb){ _verb = verb; }
		
		void SetThresh(double t){ _thresh = t; _mergeTree->SetThresh(_thresh); }
		
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
		MergeTree *_mergeTree;

		//tracks posterior values
		NodeStack _list;

		//dirichlet param
		double _alpha;

		//threshold for varEM
		double _thresh;
		//verbosity
		int _verb;

};
#endif




