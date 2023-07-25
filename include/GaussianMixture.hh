#ifndef GaussianMixture_HH
#define GaussianMixture_HH


#include "BasePDFMixture.hh"
#include "Gaussian.hh"

class GaussianMixture : public BasePDFMixture{
	public:
		GaussianMixture(){ m_k = 0; m_seed = 123; }
		GaussianMixture(int k);
		virtual ~GaussianMixture(){ };

		void SetSeed(unsigned long long seed = 123){m_seed = seed;}
		void InitParameters();
		void InitPriorParameters();
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
		void CalculateExpectations();

		double EvalLogL();
		double EvalVariationalLogL(const Matrix& post);

	protected:
		void CalculateRStatistics(const Matrix& post);
		void UpdatePriorParameters();

	private:
		//model
		vector<Gaussian*> m_gaus;

		//parameters (mu, sigma, and mixing coeffs) for k clusters
		//d-length vector (d x 1 matrix) for each cluster k
		//vector<Matrix> m_mus;
		//k weighted means - dx1 each
		vector<Matrix> m_xbars;
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
		//initial parameters
		double m_alpha0;
		double m_beta0;
		double m_nu0;
		Matrix m_W0;
		Matrix m_mean0;
		Matrix m_meanBeta0;
		Matrix m_W0inv;

		unsigned long long m_seed;
		
		//responsbility statistics
		//k normalization factors
		vector<double> m_norms;
		//k weighted covariance matrices - dxd each
		vector<Matrix> m_Ss;

		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;
};
#endif

