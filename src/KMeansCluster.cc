#include "KMeansCluster.hh"
#include "RandomSample.hh"


void KMeansCluster::Initialize(unsigned long long seed){
	if(m_k == 0){
		cout << "Error: number of clusters hasn't been set yet (via ctor or SetNClusters())." << endl;
		return;
	}
	//initialize means to be random points in data
	PointCollection initpts = m_data->SelectPoints(m_k,seed);		
	//create separate matrix for each mean
	for(int k = 0; k < m_k; k++){
		m_means.push_back(Matrix(initpts.at(k)));
		//only keeps track of which cluster assigned to (only 1 col - not multiple for multiple clusters)
		m_counts.push_back(0);
	}
	for(int n = 0; n < m_n; n++) 
		m_assigns.push_back(0);
}
void KMeansCluster::Initialize(const PointCollection& pc){
	if(m_k == 0){
		cout << "Error: number of clusters hasn't been set yet (via ctor or SetNClusters())." << endl;
		return;
	}
	//initialize from given points
	//create separate matrix for each mean
	for(int k = 0; k < m_k; k++){
		m_means.push_back(Matrix(pc.at(k)));
		//only keeps track of which cluster assigned to (only 1 col - not multiple for multiple clusters)
		m_counts.push_back(0);
	}
	for(int n = 0; n < m_n; n++) 
		m_assigns.push_back(0);
}

//E-step: estimate assignments
void KMeansCluster::Estimate(){
	double dmin, dist;
	int kmin;
	//reset stats
	m_nchg = 0;
	for(int k = 0; k < m_k; k++) m_counts[k] = 0;
	for(int n = 0; n < m_n; n++){
		dmin = 999.;
		for(int k = 0; k < m_k; k++){
			dist = 0.;
			//calculate euclidean distance
			//for(int d = 0; d < m_dim; d++) dist += pow(m_data->at(n).Value(d) - m_means[k].at(d,0),0.5);
			for(int d = 0; d < m_dim; d++) dist += pow(m_data->at(n).Value(d) - m_means[k].at(d,0),2);
			dist = sqrt(dist);
			//track mean and minimum distance to mean per point
			if(dist < dmin){ dmin = dist; kmin = k; }
		}
		//track number of points that change assignments
		if(kmin != m_assigns[n]) m_nchg++;
		//update assignment
		m_assigns[n] = kmin;
		//keep track of number of points assigned to cluster kmin
		m_counts[kmin]++;
	}
}

//M-step: update means
void KMeansCluster::Update(){
	for(int k = 0; k < m_k; k++){
		m_means[k].clear();
		m_means[k].InitEmpty();
	}
	double val;
	//sum over data points in cluster they've been assigned to
	for(int n = 0; n < m_n; n++){
		for(int d = 0; d < m_dim; d++){
			val = m_means[m_assigns[n]].at(d,0);
			m_means[m_assigns[n]].SetEntry(val+m_data->at(n).Value(d),d,0);
		}
	}
	//normalize to get mean for each cluster
	for(int k = 0; k < m_k; k++){
		if(m_counts[k] > 0){
			for(int d = 0; d < m_dim; d++){
				val = m_means[k].at(d,0);
				m_means[k].SetEntry(val/m_counts[k],d,0);
			}
		}
	}
}
