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
	
		void GetGausParameters(vector<Matrix>& mus, vector<Matrix>& covs){
			mus.clear();
			covs.clear();
			mus = m_xbars;
			covs = m_Ss;
		};
		void GetDirichletParameters(vector<double>& alphas){
			alphas.clear();
			alphas = m_alphas;
		};
		//E[pi_k] = (alpha_k + N_k)/(K*alpha_0 + N)
		void GetMixingCoeffs(vector<double>& pis){
			pis.clear();
			for(int k = 0; k < m_k; k++){
				pis.push_back((m_alphas[k] + m_norms[k])/(m_k*m_alpha0 + m_n));
			}
		};


		//get responsibilities for point n
		void GetPointPosterior(int n, vector<double>& r){
			r.clear();
			for(int k = 0; k < m_k; k++){
				r.push_back(m_post.at(n,k));
			}
		};


	protected:
		//pre-E step (don't have a million for loops)
		void CalculateExpectations();
		//pre-M step
		void CalculateRStatistics();
	private:
		//hyperparameters
		//initial values of parameters for k clusters
		double m_beta0, m_nu0, m_alpha0;
		Matrix m_W0, m_W0inv, m_mean0, m_meanBeta0;	
	
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
		vector<Matrix> m_Ws;
		//k degrees of freedom for NW
		vector<double> m_nus;

		//responsbility statistics
		//k normalization factors
		vector<double> m_norms;
		//k weighted means - dx1 each
		vector<Matrix> m_xbars;
		//k weighted covariance matrices - dxd each
		vector<Matrix> m_Ss;

		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;



};


#endif





