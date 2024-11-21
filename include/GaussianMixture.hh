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
		virtual ~GaussianMixture(){
			m_Elam.clear();
			m_Epi.clear();
		};
	
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

 
		void SetPriorParameters(map<string, Matrix> params){
			m_beta0 = params["scale"].at(0,0);
			m_nu0 = params["dof"].at(0,0);
			m_W0 = params["scalemat"];
			m_mean0 = params["mean"];
		//	if(_verb > 0){
		//		cout << "GaussianMixture::SetPriorParameters - setting prior parameters" << endl;
		//		cout << "beta0 = " << m_beta0 << " nu0 = " << m_nu0 << endl;
		//		cout << "W0 = " << endl;
		//		m_W0.Print();
		//		cout << "m0 = " << endl;
		//		m_mean0.Print();
		//	}

		}


		//shift data + learned model parameters
		void Shift(const BayesPoint& pt){
			m_data->Translate(pt);

			//only need to shift Gaussian means + prior mean
			for(int k = 0; k < m_k; k++){
				Matrix mean = m_model[k]->GetParameter("mean");
				PointCollection m = mean.MatToPoints();
				m.Translate(pt);
				m_model[k]->SetParameter("mean",m);
			}
			PointCollection m = m_mean0.MatToPoints();
			m.Translate(pt);
			m_mean0.PointsToMat(m);


		}

			
		//for variational EM algorithm
		void InitPriorParameters(unsigned long long seed = 111);
		void CalculateVariationalPosterior();
		void CalculateExpectations();
		void CalculateRStatistics();

		void UpdateVariationalParameters();
		void UpdatePriorParameters();
		double EvalVariationalLogL();
		double EvalEntropyTerms();
		double EvalNLLTerms();
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		map<string, Matrix> GetPriorParameters(int k) const; 
		map<string, Matrix> GetOnlyPriorParameters(int k); 



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


		

};


#endif
