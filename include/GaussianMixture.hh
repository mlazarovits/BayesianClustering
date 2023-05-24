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
	
		void Initialize();
		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters();
		//eval
		void EvalLogL();
		double Gaus(const Point& x, const Matrix& mu, const Matrix& cov);
	private:
		//data to fit
		PointCollection m_x;
		//TODO: need to initialize dim from data
		int m_dim;
		//TODO: number of data points - also need to init from data
		int m_n;
		//number of clusters k needs to be user specified
		int m_k;
		//parameters (mu, sigma, and mixing coeffs) for k clusters
		//d-length vector (d x 1 matrix) for each cluster k
		vector<Matrix> m_mus;
		//d x d matrix for each cluster k
		vector<Matrix> m_covs;
		//1 mixing param for each cluster k
		vector<double> m_coeffs;
		//posterior matrix of n data pts for k clusters
		Matrix m_post;		



};


#endif
