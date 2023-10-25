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

		void SetDistanceConstraint(double a = 0, double b = 1){ _mergeTree->SetDistanceConstraint(a, b); }


	private:
		//int m_nclusters;
		int _npts;
		vector<node*> _clusters;

		//calculates and tracks merges
		MergeTree *_mergeTree;

		//tracks posterior values
		NodeStack _list;


		//threshold for varEM
		double _thresh;
		//verbosity
		int _verb;

};
#endif




