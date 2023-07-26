#ifndef BASECLUSTER_HH
#define BASECLUSTER_HH

#include "PointCollection.hh"
#include "Matrix.hh"
#include "BasePDFMixture.hh"

class BaseCluster{
	public:
		BaseCluster(){ m_k = 0; m_dim = 0; m_n = 0;}
		BaseCluster(PointCollection* pc, int k){ m_k = k; m_data = pc; m_dim = m_data->Dim(); m_n = m_data->GetNPoints(); }
		BaseCluster(BasePDFMixture* pdf, int k){ m_pdfmix = pdf; m_k = k; m_data = m_pdfmix->GetData(); m_dim = m_data->Dim(); m_n = m_data->GetNPoints();}
		virtual ~BaseCluster(){ };

		//E-step
		virtual void Estimate() = 0;
		//M-step
		virtual void Update() = 0;
		//for convergence
		virtual double EvalLogL() = 0;

		void SetData(PointCollection* pc){ m_data = pc; m_dim = m_data->Dim(); m_k = 0; m_n = m_data->GetNPoints();}
		
		double Cluster(BasePDFMixture* pdfmix){
			//E step
			Estimate();
			//M step
			Update();
			
			//Check for convergence
			double newLogL = EvalLogL();
			
			pdfmix = GetModel();
			return newLogL;
		};

		void GetParameters(){ };

		void SetNClusters(int k){ m_k = k; }
		int GetNClusters(){ return m_k; }
		PointCollection* GetData(){ return m_data; }

		BasePDFMixture* GetModel(){ return m_pdfmix; }

		//model
		BasePDFMixture* m_pdfmix;
		//data
		PointCollection* m_data;
		//number of clusters
		int m_k;
		//dimension
		int m_dim;
		//number of points
		int m_n;
		//posterior matrix of n data points for k clusters
//		Matrix m_post;

		//1 mixing param for each cluster k
		vector<double> m_coeffs;
		//normalizations for each cluster, N_k = sum_n(gamma(z_nk)) (k entries)
		vector<double> m_norms;
};
#endif
