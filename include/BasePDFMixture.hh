#ifndef BASEPDFMIXTURE_HH
#define BASEPDFMIXTURE_HH

#include "PointCollection.hh"
#include "Matrix.hh"


class BasePDFMixture{
	public:
		BasePDFMixture(){ };
		
		virtual ~BasePDFMixture(){ };
		
		virtual void Initialize(unsigned long long seed = 111) = 0;
		//E-step
		virtual void CalculatePosterior() = 0;
		//M-step
		virtual void UpdateParameters() = 0;
		//some stuff for BHC

		void AddData(const PointCollection& pc){
			m_x = pc; 
			m_n = pc.GetNPoints();
			m_dim = pc.Dim();
		
		}
		
		PointCollection GetData() const{
			return m_x;
		}
	
		Matrix GetPosterior() const{
			return m_post;
		}
		
		int GetNClusters() const{
			return m_k;
		}

		//data to fit
		PointCollection m_x;
		//TODO: need to initialize dim from data
		int m_dim;
		//TODO: number of data points - also need to init from data
		int m_n;
		//number of clusters k needs to be user specified
		int m_k;
		//posterior matrix of n data pts for k clusters
		Matrix m_post;		
};










#endif
