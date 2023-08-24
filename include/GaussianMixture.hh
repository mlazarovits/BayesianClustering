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

		//estimates data points as Gaussians with mean = pt and covariance = set here
		void SetDataSmear(const Matrix& cov){ _data_cov = cov; }
		
		//fill vectors with parameters
		//returns mu, cov, and mixing coeffs for cluster k
		map<string, Matrix> GetParameters(int k); 
 
		void SetPriorParameters(map<string, Matrix> params){
			m_beta0 = params["scale"].at(0,0);
			m_nu0 = params["dof"].at(0,0);
			m_W0 = params["scalemat"];
			m_mean0 = params["mean"];

		}
			
		//for variational EM algorithm
		void InitPriorParameters(unsigned long long seed = 111);
		void CalculateVariationalPosterior();
		void CalculateExpectations();
		void CalculateRStatistics();

		void UpdateVariationalParameters();
		void UpdatePriorParameters();
		double EvalVariationalLogL();
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		map<string, Matrix> GetPriorParameters(int k); 

		

	private:
		//initial parameters
		double m_beta0;
		double m_nu0;
		Matrix m_W0;
		Matrix m_W0inv;
		Matrix m_mean0;
		Matrix m_meanBeta0;

		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;

		//data smear
		Matrix _data_cov;

};


#endif
