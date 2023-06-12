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
		//log likelihood
		virtual double EvalLogL() = 0;

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
		
		//gaussian for one data point
		double Gaus(const Point& x, const Matrix& mu, const Matrix& cov){
			double det = cov.det();
			int dim = x.Dim();
			if(dim != mu.GetDims()[0]){
				cout << "Error: x and mu length do not match. " << dim << " " << mu.GetDims()[0] << endl;
				return 0;
			}
			if(mu.GetDims()[0] != cov.GetDims()[0]){
				cout << "Error: covariance and mu dimensions do not match. " << cov.GetDims()[0] << " " << mu.GetDims()[0] << endl;
				return 0;
			}
			if(!cov.square()){
				cout << "Error: non-square covariance matrix." << endl;
				return 0;
			}
			//construct x - mu
			Matrix x_mat = Matrix(x.Value());
			Matrix x_min_mu;
			x_min_mu.mult(mu,-1.);
			x_min_mu.add(x_mat);
			
			//transpose x - mu
			Matrix x_min_muT;
			x_min_muT.transpose(x_min_mu);
			Matrix cov_inv;
			cov_inv.invert(cov);
			
			double coeff = 1/(pow(det,0.5)*pow(2*acos(-1),0.5*dim));
			//should only be 1 element matrix
			//muT*cov*mu = 1xd * dxd * dx1
			Matrix mat_expon = Matrix(1,1);
			Matrix cov_mu = Matrix(m_dim,1);
			
			cov_mu.mult(cov_inv,x_min_mu);
			mat_expon.mult(x_min_muT,cov_mu);
			double expon = mat_expon.at(0,0);
			return coeff*exp(-0.5*expon);
			
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
