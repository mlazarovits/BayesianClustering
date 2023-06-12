#ifndef VARGAUSSIANMIXTURE_HH
#define VARGAUSSIANMIXTURE_HH

#include "BasePDFMixture.hh"

class VarGaussianMixture : public BasePDFMixture{

	public:
		VarGaussianMixture();
		VarGaussianMixture(int k);
		virtual ~VarGaussianMixture(){ };

		void Initialize(unsigned long long seed = 111);

		void InitAlpha(double a){
			m_alpha0 = a;
		}


		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();
	
		//UPDATE - fill vectors with estimated parameters
		void GetParameters(vector<Matrix>& mus, vector<Matrix>& covs);

	protected:
		//pre-E step (don't have a million for loops)
		void CalculateExpectations();


	private:
		//hyperparameters
		//number of clusters - can change
		int m_k;
		//number of data points - doesn't change
		int m_n;
		//initial values of parameters
		double m_beta0, m_nu0, m_alpha0;
		Matrix m_W0inv, m_mean0;
		Matrix m_meanBeta0;	
	
		//parameters
		//k concentrations parameters for each Dirichlet prior on mixing coeff
		vector<double> m_alphas;
		//k scaling parameters for NW (normal wishart) distribution
		vector<double> m_betas;
		//k dx1 means for normal component of NW
		vector<Matrix> m_means;
		//k dxd matrices - covariance matrix of normal distribution that occurs in the construction of a Wishart distribution
		//G_i = (g1_i, ..., gp_i)T ~ N_p(0, W)
		//S = GG_T ~ W_p(V,n) for n degrees of freedom
		vector<Matrix> m_W;
		//k degrees of freedom for NW
		vector<double> m_nus;

		//responsbility statistics
		//k normalization factors
		vector<double> m_norms;
		//k weighted means - dx1 each
		vector<Matrix> m_mus;
		//k weighted covariance matrices - dxd each
		vector<Matrix> m_covs;

		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;



};


#endif





