#ifndef GaussianMixture_HH
#define GaussianMixture_HH


#include "BasePDFMixture.hh"
#include "Gaussian.hh"

class GaussianMixture : public BasePDFMixture{
	//inheriting all constructors from base class
	//using BasePDFMixture::BasePDFMixture;

	public:
		GaussianMixture(){ m_k = 0; m_seed = 123; }
		GaussianMixture(int k);
		virtual ~GaussianMixture(){ };

		void SetSeed(unsigned long long seed){m_seed = seed;}
		void InitParameters();
		//for EM algorithm
		void CalculatePosterior(Matrix& post);
		void UpdateParameters(const Matrix& post);
		//returns mu, cov, and mixing coeffs
		map<string, vector<Matrix>> GetParameters(); 
		
		//for variational EM algorithm
		void CalculateVariationalPosterior(Matrix& post);
		void UpdateVariationalParameters(const Matrix& post);
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart)
		map<string, vector<Matrix>> GetPriorParameters(); 

		double EvalLogL();
		double EvalVariationalLogL();

	private:
		//model
		vector<Gaussian*> m_gaus;

		//parameters (mu, sigma, and mixing coeffs) for k clusters
		//d-length vector (d x 1 matrix) for each cluster k
		vector<Matrix> m_mus;
		//d x d matrix for each cluster k
		vector<Matrix> m_covs;
	
		//parameters on conjugate prior
		//k scaling parameters for NW (normal wishart) distribution
		vector<double> m_betas;
		//k dx1 means for normal component of NW
		vector<Matrix> m_means;
		//k dxd matrices - covariance matrix of normal distribution that occurs in the construction of a Wishart distribution
		//G_i = (g1_i, ..., gp_i)T ~ N_p(0, W)
		//S = GG_T ~ W_p(V,n) for n degrees of freedom
		vector<Matrix> m_Ws;
		//k degrees of freedom for NW
		vector<double> m_nus;


		unsigned long long m_seed;
};
#endif
