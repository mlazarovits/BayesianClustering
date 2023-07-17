#ifndef BAYESHIERCLUSTER_HH
#define BAYESHIERCLUSTER_HH

#include "PointCollection.hh"
#include "BasePDF.hh"
#include "ClusterTree.hh"

using std::vector;
class BayesHierCluster{
	public:
		BayesHierCluster();
		BayesHierCluster(PointCollection pc);
		virtual ~BayesHierCluster(){ };
		
		void SetClusterPDF(BasePDF* pdf);

		//get clusters at specified rk cut-off
		vector<PointCollection> GetClusters(double rk);	

		//get clusters at specified depth
		vector<PointCollection> GetClusters(int d);	
		
		//cluster trees (points in pc)
		void Cluster(const vector<PointCollection>& pc);

	private:
		PointCollection m_pts;
		BasePDF* m_pdf; //p(x | theta)

		//clustering history
		ClusterTree m_clusterTree;

		//dirichlet param
		double m_alpha;

		//weighting params: pi (multinomial)
		vector<double> m_pi;

		//d_k for recursive prior computation 
		vector<double> m_d;


};
#endif




