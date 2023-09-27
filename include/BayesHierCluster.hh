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
	
		void SetVerbosity(int verb){ _verb = verb; _mergeTree->SetVerbosity(_verb);}
		
		void SetThresh(double t){ _thresh = t; _mergeTree->SetThresh(_thresh); }
		void SetDataSmear(const Matrix& cov){ _mergeTree->SetDataSmear(cov); }		
		void SetPriorParameters(map<string, Matrix> params){ _mergeTree->SetPriorParameters(params); }		
		//cluster trees (points in pc)
		//might need arg for recursion
		//void Cluster(const vector<PointCollection>& pc);
		vector<node*> Cluster();

		void SetAlpha(double a){_mergeTree->SetAlpha(a);}
		void SetSubclusterAlpha(double a){ _mergeTree->SetSubclusterAlpha(a); }

		void SetDistanceConstraint(double a = 0, double b = 1){ _constrain = true; _constraint_a = a; _constraint_b = b; }
		void SetPhiWraparound(bool phi){ _wraparound = phi; _mergeTree->SetPhiWraparound(_wraparound); }
		double DistanceConstraint(node* i, node* j);


	private:
		//int m_nclusters;

		vector<PointCollection*> m_pts;
		vector<node*> _clusters;

		//calculates and tracks merges
		MergeTree *_mergeTree;

		//tracks posterior values
		NodeStack _list;


		//threshold for varEM
		double _thresh;
		//verbosity
		int _verb;

		//hierarchical cluster distance constraints
		//dimension along which to constrain
		//int _constraint_d;
		//transform distance to this interval
		double _constraint_min, _constraint_max;
		//cutoff interval (clusters outside of [a,b] will not be clustered)
		double _constraint_a, _constraint_b;
		bool _constrain;
		//if phi is used for distance constraint
		bool _wraparound;
};
#endif




