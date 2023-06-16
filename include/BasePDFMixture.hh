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

		//fill vectors with estimated parameters
		virtual void GetGausParameters(vector<Matrix>& mus, vector<Matrix>& covs) = 0;
		virtual void GetDirichletParameters(vector<double>& alphas) = 0;
		virtual void GetMixingCoeffs(vector<double>& pis) = 0;

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
		double Gaus(const Point& x, const Matrix& mu, const Matrix& cov);


		//ln of Dirichlet coefficient (eq. B.23)
		double Dir_lnC(double alpha, int k);
		//ln of Dirichlet coefficient (eq. B.23)
		double Dir_lnC(vector<double> alphas);

		//log of Wishart entropy normalization (eq. B.79)
		double Wish_lnB(Matrix W, double nu);

		//Wishart entropy (eq. B.82)
		double Wish_H(Matrix W, double nu);

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
