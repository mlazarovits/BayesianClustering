#include "BayesHierCluster.hh"



BayesHierCluster::BayesHierCluster(){ };


BayesHierCluster::BayesHierCluster(PointCollection pc){
	m_pts = pc;


}

BayesHierCluster::SetClusterPDF(BasePDF* pdf){
	m_pdf = pdf;
}

BayesHierCluster::SetAlpha(double a){
	m_alpha = a;
}

BayesHierCluster::Init(){
	//init values: 
	//P(D_i | T_i) = p(D_i | H_1^i) based on PDF given
	//d_i = alpha
	//pi_i = 1

}

vector<PointCollection> BayesHierCluster::GetClusters(double rk){
	return m_clusterTree.at(rk);
}

vector<PointCollection> BayesHierCluster::GetClusters(int d){
	return m_clusterTree.at(d);
}


void Cluster(const vector<PointCollection>& pc){
	//find pairs of (sub)trees
	//for(int i = 0; i < pc.size(); i++){
	//	for(int j = 0; j < pc.size(); j++){
	//		if(i > j) continue;
	//



}







