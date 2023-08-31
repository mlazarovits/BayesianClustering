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

		//init from randomly selected points
		void Initialize(unsigned long long seed = 123);
		//init means from given points
		void Initialize(const PointCollection& pc);
		void Initialize_pp(unsigned long long seed = 123);

		//E-step
		void Estimate();
		//M-step
		void Update();
		//for convergence
		double EvalLogL(){ return double(m_nchg); }

		void GetMeans(vector<Matrix>& m){ m.clear(); m = m_means; }
		map<string, vector<Matrix>> GetModel(){ map<string,vector<Matrix>> ret;
			ret["means"] = m_means;
			return ret;
		 }
		//get assignments for each point
		void GetAssignments(vector<int>& a){ a.clear(); a = m_assigns; }
		//get number of assignments per cluster	
		int GetAssignment(int n){ if(n < m_data->GetNPoints()) return m_assigns[n]; else return -999.; }
		//get points separated into point collections for each cluster
		void GetAssignments(vector<PointCollection*>& pcs);
	
		void GetCounts(vector<double>& c){ c.clear(); c = m_counts; }
	private:
		//parameters - means (1 per cluster)
		vector<Matrix> m_means;
		//assignment values
		vector<int> m_assigns;
		//number of points assigned to each cluster
		vector<double> m_counts;
		//number of points change assignment
		double m_nchg;
};

#endif
