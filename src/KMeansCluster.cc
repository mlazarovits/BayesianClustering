#include "KMeansCluster.hh"
#include "RandomSample.hh"


void KMeansCluster::Initialize(unsigned long long seed){
	if(m_k == 0){
		cout << "Error: number of clusters hasn't been set yet (via ctor or SetNClusters())." << endl;
		return;
	}
	//initialize means to be random points in data
	PointCollection initpts = m_data->SelectPoints(m_k,seed);		
	if(initpts.GetNPoints() < m_k && initpts.GetNPoints() != 0) m_k = initpts.GetNPoints();
	//cout << "KMeans starting points: " << initpts.GetNPoints() << " m_k: " << m_k << endl;
	//initpts.Print();
	//create separate matrix for each mean
	for(int k = 0; k < m_k; k++){
		m_means.push_back(Matrix(initpts.at(k)));
		//only keeps track of which cluster assigned to (only 1 col - not multiple for multiple clusters)
		m_counts.push_back(0);
	}
	for(int n = 0; n < m_n; n++) 
		m_assigns.push_back(-1);
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
		m_assigns.push_back(-1);
}


//kmeans++ initialization
void KMeansCluster::Initialize_pp(unsigned long long seed){
	//randomly select first point
	BayesPoint c = m_data->SelectPoints(1, seed).at(0);
	PointCollection* means = new PointCollection();
	PointCollection* pts = new PointCollection;
	means->AddPoint(c);

	for(int i = 0; i < m_n; i++){
		if(m_data->at(i) == c) continue;
		pts->AddPoint(m_data->at(i));
	}

	RandomSample rs;
	int nit = 0;
	while(means->GetNPoints() < m_k){
		//calculate distance of all other points from c
		vector<double> dist;
		double d;
		for(int i = 0; i < pts->GetNPoints(); i++){
			d = 0;
			for(int j = 0; j < m_dim; j++) d += pow(pts->at(i).Value(j) - c.Value(j),2);
			d *= pts->at(i).w();
			dist.push_back(d);
		}

		int idx = rs.SampleCategorical(dist);
		means->AddPoint(pts->at(idx));
		c = pts->at(idx);
		pts->Remove(idx);

	}
	//cout << "centers" << endl;
	double N = 0;
	for(int k = 0; k < m_k; k++){ 
	//	cout << "k: " << k << endl; means->at(k).Print(); 
		m_means.push_back(Matrix(means->at(k)));
		m_counts.push_back(0);
	}
	for(int n = 0; n < m_n; n++) N += m_data->at(n).w();
 
//cout << "total norm: " << N << endl;
	
	for(int n = 0; n < m_n; n++) 
		m_assigns.push_back(-1);
}

//E-step: estimate assignments
void KMeansCluster::Estimate(){
	double dmin, dist;
	int kmin = -1;
	//reset stats
	m_nchg = 0;
	for(int k = 0; k < m_k; k++) m_counts[k] = 0;
	for(int n = 0; n < m_n; n++){
		dmin = 1e300;
		for(int k = 0; k < m_k; k++){
			dist = 0.;
			//calculate euclidean distance squared as optimization metric
			//for(int d = 0; d < m_dim; d++) dist += pow(m_data->at(n).Value(d) - m_means[k].at(d,0),0.5);
			for(int d = 0; d < m_dim; d++) dist += pow(m_data->at(n).Value(d) - m_means[k].at(d,0),2);
			dist = sqrt(dist);
			//dist *= m_data->at(n).w();
			//track mean and minimum distance to mean per point
			//considering saving all dists and sorting to get min
			if(dist < dmin){ dmin = dist; kmin = k; }
			//cout << "k: " << k << " current kmin: " << kmin << " dist: " << dist << " current dmin: " << dmin << endl;
		}
		//if no cluster assigned - assign to first cluster
		if(kmin == -1) kmin = 0;
		//track number of points that change assignments
		if(kmin != m_assigns[n]) m_nchg++; //+= m_data->at(n).w();//m_nchg++; 
		//update assignment
		m_assigns[n] = kmin;
		//keep track of number of points assigned to cluster kmin
		//unweighted data -> w = 1
		m_counts[kmin] += m_data->at(n).w();
	}
}

//M-step: update means
void KMeansCluster::Update(){
	for(int k = 0; k < m_k; k++){
		m_means[k].clear();
		m_means[k].InitEmpty();
	}
	double val, w;
	double test = 0;
	//	cout << "m_n " << m_n << " m_k " << m_k << " d " << m_dim << " means " << m_means.size() << " assigns " << m_assigns.size() << " counts " << m_counts.size() << " data " << m_data->GetNPoints() << " data dim " << m_data->Dim() << endl;
	//sum over data points in cluster they've been assigned to
	for(int n = 0; n < m_n; n++){
		w = m_data->at(n).w();
		for(int d = 0; d < m_dim; d++){
			val = m_means[m_assigns[n]].at(d,0);
			//unweighted data -> w = 1
			//cout << "setting entry for n " << n << " d " << d << " assignment " << m_assigns[n] << " means size  " << m_means[m_assigns[n]].GetDims()[0] << " " << m_means[m_assigns[n]].GetDims()[1] << endl;
			m_means[m_assigns[n]].SetEntry(val + w*m_data->at(n).Value(d),d,0);
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


void KMeansCluster::GetAssignments(vector<PointCollection*>& pcs){
	pcs.clear();
	//push back k point collections of m_counts[k] points
	for(int k = 0; k < m_k; k++){
		pcs.push_back(new PointCollection());
	}
	int kass;
	for(int n = 0; n < m_n; n++){
		kass = m_assigns[n];
		pcs[kass]->AddPoint(m_data->at(n));	
	}


}
