#ifndef KMEANSCLUSTER_HH
#define KMEANSCLUSTER_HH

#include "BaseCluster.hh"
#include "Matrix.hh"



class KMeansCluster : public BaseCluster{
	//inheriting constructors
	using BaseCluster::BaseCluster;
	public:
		//KMeansCluster();
		virtual ~KMeansCluster(){ m_means.clear(); };

		void Initialize();

		//E-step
		void Estimate();
		//M-step
		void Update();
		//for convergence
		double EvalLogL(){ return double(m_nchg); }

		void GetMeans(vector<Matrix>& m){ m.clear(); m = m_means; }
		void GetAssignments(vector<int>& a){ a.clear(); a = m_assigns; }	
		int GetAssignment(int n){ if(n < m_data->GetNPoints()) return m_assigns[n]; else return -999.; }
	

	private:
		//parameters - means (1 per cluster)
		vector<Matrix> m_means;
		//assignment values
		vector<int> m_assigns;
		//number of points assigned to each cluster
		vector<int> m_counts;
		//number of points change assignment
		int m_nchg;
};

#endif
