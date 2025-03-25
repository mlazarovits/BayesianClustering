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
		map<string, Matrix> GetLikelihoodParameters(int k); 

		//void SetJetPriorParameters(map<string, Matrix>& params){ }
		//void SetJetParameters(map<string, Matrix>& params){ }
		//void GetJetParameters(map<string, Matrix>& params){ }
 
		void SetPriorParameters(map<string, Matrix> params){
			m_beta0 = params["scale"].at(0,0);
			m_nu0 = params["dof"].at(0,0);
			m_W0 = params["scalemat"];
			m_W0inv.invert(m_W0);
			m_mean0 = params["mean"];
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


				//shift data statistics
				_xbar[k].minus(shift);
			}
		}
		
		void PutPhi02pi_params(){
			//put means on [0,2pi]
			Matrix mean, priormean;
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				PointCollection mean_pt = mean.MatToPoints();
				mean_pt.Put02pi(1);	
				m_model[k]->SetParameter("mean",Matrix(mean_pt));
				
				//translate posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				PointCollection priormean_pt = priormean.MatToPoints();
				priormean_pt.Put02pi(1);	
				m_model[k]->SetParameter("priormean",Matrix(priormean_pt));


				//shift data statistics
				PointCollection xbar_pt = _xbar[k].MatToPoints();
				xbar_pt.Put02pi(1);
				_xbar[k] = Matrix(xbar_pt);	
			}
		}
		
		//scale learned model parameters
		void ScaleParameters(Matrix sc){
			//scale means + matrices
			Matrix mean, priormean, cov, W;
			Matrix scT, scInv, scTinv;
			scT.transpose(sc);
			scInv.invert(sc);
			scTinv.transpose(scInv);
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
				W.mult(scInv,W); //Avar(X)
				W.mult(W,scTinv); //Avar(X)A^T
				m_model[k]->GetPrior()->SetParameter("scalemat",W);
		
				//scale data parameters - Sbar
				_Sbar[k].mult(sc,_Sbar[k]); //Avar(X)
				_Sbar[k].mult(_Sbar[k],scT); //Avar(X)A^T

			}

		}

			
		//for variational EM algorithm
		void InitPriorParameters(unsigned long long seed = 111);
		void CalculateVariationalPosterior();
		void CalculateExpectations();
		void CalculateRStatistics();

		void UpdateVariationalParameters();
		void UpdatePosteriorParameters();
		double EvalVariationalLogL();
		double EvalEntropyTerms();
		double EvalNLLTerms();
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		map<string, Matrix> GetLHPosteriorParameters(int k) const; 
		map<string, Matrix> GetOnlyPosteriorParameters(int k);
		map<string, Matrix> GetDataStatistics(int k) const; 


		


	private:
		//initial parameters - prior is normal-wishart over mean and precision parameters of likelihood
		double m_beta0;
		double m_nu0;
		Matrix m_W0;
		Matrix m_W0inv;
		Matrix m_mean0;
		Matrix m_meanBeta0;

		//data statistics
		vector<Matrix> _xbar;
		vector<Matrix> _Sbar;
		
		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;


		

};


#endif
