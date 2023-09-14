#ifndef BAYESHIERCLUSTER_HH
#define BAYESHIERCLUSTER_HH

#include "PointCollection.hh"
#include "BasePDF.hh"
#include "NodeStack.hh"
#include "MergeTree.hh"

using std::vector;
class BayesHierCluster{
	public:
		BayesHierCluster();
		BayesHierCluster(double alpha);
		virtual ~BayesHierCluster(){ };
	
		void AddData(PointCollection* pc);
	
		void SetVerbosity(int verb){ _verb = verb; }
		
		void SetThresh(double t){ _thresh = t; _mergeTree->SetThresh(_thresh); }
		void SetDataSmear(const Matrix& cov){ _mergeTree->SetDataSmear(cov); }		
		void SetPriorParameters(map<string, Matrix> params){ _mergeTree->SetPriorParameters(params); }		

		//cluster trees (points in pc)
		//might need arg for recursion
		//void Cluster(const vector<PointCollection>& pc);
		vector<node*> Cluster();

		void SetAlpha(double a);

		void SetDistanceConstraint(int d, double thresh, double a = 0, double b = 1){ _constraint_d = d; _constraint_thresh = thresh; _constraint_a = a; _constraint_b = b; }
		double DistanceConstraint(node* i, node* j);


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

		//hierarchical cluster distance constraints
		//dimension along which to constrain
		int _constraint_d;
		//cutoff threshold (clusters farther away than thresh will not be clustered)
		double _constraint_thresh;
		double _constraint_a, _constraint_b;
};
#endif




