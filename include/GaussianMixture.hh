#ifndef GAUSSIANMixture_HH
#define GAUSSIANMixture_HH
#include "BasePDFMixture.hh"
#include "Matrix.hh"
#include "PointCollection.hh"
#include <vector>
using std::vector;


class GaussianMixture : public BasePDFMixture{
	public:
		GaussianMixture();
		GaussianMixture(int k);
		virtual ~GaussianMixture(){ };
	
		void Initialize(unsigned long long seed = 111);
		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();
		double Gaus(const Point& x, const Matrix& mu, const Matrix& cov);
	private:
		//parameters (mu, sigma, and mixing coeffs) for k clusters
		//d-length vector (d x 1 matrix) for each cluster k
		vector<Matrix> m_mus;
		//d x d matrix for each cluster k
		vector<Matrix> m_covs;
		//1 mixing param for each cluster k
		vector<double> m_coeffs;
		//posterior matrix of n data pts for k clusters, post_nk = gamma(z_nk)
		Matrix m_post;		
		//normalizations for each cluster, N_k = sum_n(gamma(z_nk)) (k entries)
		vector<double> m_norms;


};


#endif
