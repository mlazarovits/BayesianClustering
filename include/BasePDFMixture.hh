#ifndef BasePDFMixture_HH
#define BasePDFMixture_HH

#include "Matrix.hh"
#include "BasePDF.hh"
#include <vector>
#include <map>
#include <string>
using std::vector;
using std::map;
using std::string;

class BasePDFMixture : public BasePDF{
	public:
		BasePDFMixture(){ m_k = 0; m_n = 0; }
		BasePDFMixture(int k){ m_k = k; 
			for(int k = 0; k < m_k; k++){
				m_weights.push_back(0.);
				m_alphas.push_back(0.);
			}
			m_n = 0;
		}

		virtual ~BasePDFMixture(){ m_weights.clear(); m_alphas.clear(); m_model.clear(); }

		double Prob(const Point& x);

		void SetData(PointCollection* data){m_data = data; m_n = m_data->GetNPoints();}
		PointCollection* GetData(){ return m_data; }


		//for EM algorithm
		virtual void CalculatePosterior(Matrix& post) = 0;
		virtual void UpdateParameters(const Matrix& post) = 0;
		//returns mu, cov, and mixing coeffs
		virtual map<string, vector<Matrix>> GetParameters() = 0; 
		
		//for variational EM algorithm
		virtual void CalculateVariationalPosterior(Matrix& post) = 0;
		virtual void UpdateVariationalParameters(const Matrix& post) = 0;
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart)
		virtual map<string, vector<Matrix>> GetPriorParameters() = 0; 

		void GetMixingCoeffs(vector<double>& weights){ weights.clear(); weights = m_weights; }	
		void GetDirichletParams(vector<double>& alphas){ alphas.clear(); alphas = m_alphas; }

		virtual double EvalLogL() = 0;
	
		PointCollection* m_data;
		//number of data points
		int m_n;
		//number of components
		int m_k;
		//probabilistic model
		vector<BasePDF*> m_model;
		//mixture weights (probabilities sum to 1, multinomal dist.)
		vector<double> m_weights;
		//dirichlet (prior) parameters (on weights)
		vector<double> m_alphas;
		//normalization on posterior
		vector<double> m_norms;



};
#endif
