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
		double GetPi(int k){ return (m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_data->Sumw()); }
		double GetCoeff(int k){ return m_coeffs[k]; }

		void SetJetPriorParameters(map<string, Matrix>& params){ }
		void SetJetParameters(map<string, Matrix>& params){ }
		void GetJetParameters(map<string, Matrix>& params){ }
 
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
			m_meanBeta0 = Matrix(m_dim, 1);
			m_meanBeta0.mult(m_mean0, m_beta0);

		}


		//shift learned model parameters
		void ShiftParameters(const BayesPoint& pt){
			//only need to shift Gaussian means + prior mean
			Matrix mean, priormean;
			Matrix shift;
			shift.PointToMat(pt);
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				mean.minus(shift);
				m_model[k]->SetParameter("mean",mean);
				
				//translate posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				priormean.minus(shift);
				m_model[k]->GetPrior()->SetParameter("mean",priormean);
			}
		}
		
		//scale learned model parameters
		void ScaleParameters(Matrix sc){
			//scale means + matrices
			Matrix mean, priormean, cov, W;
			Matrix scT;
			scT.transpose(sc);
			for(int k = 0; k < m_k; k++){
				//scale posterior mean in likelihood 
				mean = m_model[k]->GetParameter("mean");
				mean.mult(sc,mean);	
				m_model[k]->SetParameter("mean",mean);
				
				//scale posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				priormean.mult(sc,priormean);
				m_model[k]->GetPrior()->SetParameter("mean",priormean);

				//scale posterior cov in likelihood
				//var(AX) = Avar(X)A^T 
				cov = m_model[k]->GetParameter("cov");
				cov.mult(sc,cov); //Avar(X)
				cov.mult(cov,scT); //Avar(X)A^T
				m_model[k]->SetParameter("cov",cov);

				//scale posterior W in prior distribution
				W = m_model[k]->GetPrior()->GetParameter("scalemat");
				W.mult(sc,W); //Avar(X)
				W.mult(W,scT); //Avar(X)A^T
				m_model[k]->GetPrior()->SetParameter("scalemat",W);
			}

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
