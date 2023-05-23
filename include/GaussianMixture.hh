#ifndef GAUSSIANMixture_HH
#define GAUSSIANMixture_HH
#include "BaseMixture.hh"
#include "Matrix.hh"
#include <vector>
using std::vector;


class GaussianMixture : public BasePDFMixture{
	public:
		GaussianMixture();
		virtual ~GaussianMixture(){ };
	
		void Initialize();
		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters(BaseMixture* post);
	private:
		//TODO: need to initialize dim from data
		int m_dim;
		//TODO: number of data points - also need to init from data
		int m_n;
		//number of clusters k needs to be user specified
		int m_k;
		//parameters (mu, sigma, and mixing coeffs) for k clusters
		vector<double> m_mus;
		vector<Matrix> m_covs;
		vector<double> m_coeffs;
		//posterior matrix of n data pts for k clusters
		Matrix m_post;		






#endif
