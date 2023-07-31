#ifndef BAYESHIERCLUSTER_HH
#define BAYESHIERCLUSTER_HH

#include "PointCollection.hh"
#include "BasePDF.hh"
#include "ClusterHistory.hh"
#include "SearchTree.hh"
#include "MergeTree.hh"

using std::vector;
class BayesHierCluster{
	public:
		BayesHierCluster();
		BayesHierCluster(PointCollection* pc);
		virtual ~BayesHierCluster(){ };
		
		void Init();
		void SetClusterPDF(BasePDF* pdf);

		//get clusters at specified rk cut-off
		vector<PointCollection> GetClusters(double rk);	

		//get clusters at specified depth
		vector<PointCollection> GetClusters(int d);	
		
		//cluster trees (points in pc)
		//might need arg for recursion
		//void Cluster(const vector<PointCollection>& pc);
		void Cluster();

		void SetAlpha(double a);
	private:
		int m_nclusters;

		PointCollection* m_pts;
		BasePDF* m_pdf; //p(x | theta)

		//clustering history
		ClusterHistory m_clusterHist;
		
		//calculates and tracks merges
		MergeTree m_mergeTree;

		//tracks posterior values
		SearchTree m_searchTree;

		//dirichlet param
		double m_alpha;

		//weighting params: pi (multinomial)
		vector<double> m_pi;

		//d_k for recursive prior computation 
		vector<double> m_d;


};
#endif




