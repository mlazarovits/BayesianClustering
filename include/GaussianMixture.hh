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
	
		void InitParameters(unsigned long long seed = 111);
		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();
		
		//fill vectors with parameters
		//returns mu, cov, and mixing coeffs for cluster k
		map<string, Matrix> GetParameters(int k); 
 

			
		//for variational EM algorithm
		void InitPriorParameters(unsigned long long seed = 111);
		void InitParameters(const PointCollection& pc);
		void CalculateVariationalPosterior();
		void CalculateExpectations();
		void CalculateRStatistics();

		void UpdateVariationalParameters();
		void UpdatePriorParameters();
		double EvalVariationalLogL();
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		map<string, Matrix> GetPriorParameters(int k); 

		

	private:
		/*
		//parameters (mu, sigma, and mixing coeffs) for k clusters
		//d-length vector (d x 1 matrix) for each cluster k
		vector<Matrix> m_mus;
		//d x d matrix for each cluster k
		vector<Matrix> m_covs;


		//variational stuff
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
		*/

		//initial parameters
		double m_beta0;
		double m_nu0;
		Matrix m_W0;
		Matrix m_W0inv;
		Matrix m_mean0;
		Matrix m_meanBeta0;

		/*
		//responsbility statistics
		//k weighted means - dx1 each
		vector<Matrix> m_xbars;
		//k weighted covariance matrices - dxd each
		vector<Matrix> m_Ss;
		*/

		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;

};


#endif
