#ifndef BASECLUSTER_HH
#define BASECLUSTER_HH

#include "PointCollection.hh"
#include "Matrix.hh"


class BaseCluster{
	public:
		BaseCluster(){ m_k = 0; m_dim = 0; m_n = 0;}
		BaseCluster(int k){ m_k = k; m_dim = 0; m_n = 0;}
		BaseCluster(PointCollection* pc){ m_data = pc; m_dim = m_data->Dim(); m_k = 0; m_n = m_data->GetNPoints();}
		virtual ~BaseCluster(){ };

		virtual void Initialize(unsigned long long seed) = 0;
		//E-step
		virtual void Estimate() = 0;
		//M-step
		virtual void Update() = 0;
		//for convergence
		virtual double EvalLogL() = 0;

		void SetNClusters(int k){ m_k = k; }
		int GetNClusters(){ return m_k; }
		void AddData(PointCollection* pc){ m_data = pc; m_dim = m_data->Dim(); }
		PointCollection* GetData(){ return m_data; }
		
		PointCollection* m_data;
		//number of clusters
		int m_k;
		//dimension
		int m_dim;
		//number of points
		int m_n;
		//posterior matrix of n data points for k clusters
		Matrix m_post;

};
#endif
