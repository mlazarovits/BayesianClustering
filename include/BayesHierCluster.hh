#ifndef BAYESHIERCLUSTER_HH
#define BAYESHIERCLUSTER_HH

#include "PointCollection.hh"
#include "BasePDF.hh"
#include "ClusterHistory.hh"
#include "NodeList.hh"
#include "MergeTree.hh"

using std::vector;
class BayesHierCluster{
	public:
		BayesHierCluster();
		BayesHierCluster(const PointCollection* pc);
		virtual ~BayesHierCluster(){ };
		
		void Init();
		void SetClusterPDF(BasePDF* pdf);

		//get clusters at specified rk cut-off
		vector<PointCollection*> GetClusters(double rk);	

		//get clusters at specified depth
		vector<PointCollection*> GetClusters(int d);	
		
		//cluster trees (points in pc)
		//might need arg for recursion
		//void Cluster(const vector<PointCollection>& pc);
		void Cluster();

		void SetAlpha(double a);
	private:
		int m_nclusters;

		vector<PointCollection*> m_pts;
		BasePDF* m_pdf; //p(x | theta)

		//clustering history
		ClusterHistory m_clusterHist;
		
		//calculates and tracks merges
		MergeTree m_mergeTree;

		//tracks posterior values
		NodeList _list;

		//dirichlet param
		double m_alpha;


};
#endif




