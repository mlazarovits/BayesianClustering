#ifndef HIERGAUSSIANMixture_HH
#define HIERGAUSSIANMixture_HH
#include "BasePDFMixture.hh"
#include "Matrix.hh"
#include "PointCollection.hh"
#include "Gaussian.hh"
#include <vector>
using std::vector;

//A class that fits a hierarchical Gaussian Mixture model
//with k subclusters and a prior on the location parameter of these subclusters
//such that the joint likelihood of the model is
//p(X, Z, pi, mu, lam, m) = p(pi)*p(Z | pi)*p(X | Z, mu, lam)*p(mu, lam)*p(m_j, lam_j)
//with prior parameters and functional forms,
//p(X, Z, pi, mu, lam, m) = Dir(pi | alpha0)*Categorical(Z | pi)*N(X | Z, mu, lam)*N(mu | lam, m_j, beta0)*Wishart(lam | W0, nu0)*p(m_j | lam_j, m_j0, betaj_0)*Wishart(lam_j | W_j0, nu_j0)
//where subscript 0 parameters are prior hyperparameters to be set and subsequently updated with the EM algorithm
//and subscript j parameters describe the overall "jet" that the MM is a part of
class HierGaussianMixture : public BasePDFMixture{
	public:
		HierGaussianMixture();
		HierGaussianMixture(int k);
		virtual ~HierGaussianMixture(){
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

			m_meanBeta0 = Matrix(m_dim, 1);
			m_meanBeta0.mult(m_mean0, m_beta0);

			m_W0inv = Matrix(m_dim,m_dim);
			m_W0inv.invert(m_W0);

		}
		void SetJetPriorParameters(map<string, Matrix>& params){	
			m_betaj0 = params["scale"].at(0,0);
			m_nuj0 = params["dof"].at(0,0);
			m_Wj0 = params["scalemat"];
			m_meanj0 = params["mean"];
			
			m_Wj0inv.invert(m_Wj0);
			m_meanBetaj0 = Matrix(m_dim, 1);
			m_meanBetaj0.mult(m_meanj0,m_betaj0);
		//	if(_verb > 0){
		//		cout << "HierGaussianMixture::SetPriorParameters - setting prior parameters" << endl;
		//		cout << "beta0 = " << m_beta0 << " nu0 = " << m_nu0 << endl;
		//		cout << "W0 = " << endl;
		//		m_W0.Print();
		//		cout << "m0 = " << endl;
		//		m_mean0.Print();
		//	}
			m_modelj->SetPrior(new NormalWishart(m_Wj0, m_meanj0, m_nuj0, m_betaj0));
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
		void CalculateJStatistics();

		void UpdateVariationalParameters();
		void UpdatePriorParameters();
		double EvalVariationalLogL();
		double EvalEntropyTerms();
		double EvalNLLTerms();
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		map<string, Matrix> GetPriorParameters(int k) const; 
		map<string, Matrix> GetOnlyPriorParameters(int k); 
		void GetJetParameters(map<string, Matrix>& params){
			params.clear();
			params["mean"] = m_modelj->GetParameter("mean");
			params["cov"] = m_modelj->GetParameter("cov");
		}
		void SetJetParameters(map<string, Matrix>& params){
			m_modelj->SetParameter("mean",params["mean"]);
			m_modelj->SetParameter("cov",params["cov"]);
		}
		


	private:
		//initial parameters for MM
		double m_beta0;
		double m_nu0;
		Matrix m_W0;
		Matrix m_W0inv;
		Matrix m_mean0;
		Matrix m_meanBeta0;

		//initial parameters for jet
		double m_betaj0;
		double m_nuj0;
		Matrix m_Wj0;
		Matrix m_Wj0inv;
		Matrix m_meanj0;
		Matrix m_meanBetaj0;

		Gaussian* m_modelj;

		//expectation values - used in ELBO + calculated in an earlier step
		//E_lam = E[ln|lambda_k|] (eq. 10.65)
		//E_pi = E[ln(pi_k)] (eq. 10.66)
		vector<double> m_Elam, m_Epi;
		double m_Elamj;

		

};


#endif
