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

class BasePDFMixture{
	public:
		BasePDFMixture(){ m_k = 0; m_n = 0; }
		BasePDFMixture(int k){ m_k = k; 
			for(int k = 0; k < m_k; k++){
				m_weights.push_back(0.);
				m_alphas.push_back(0.);
			}
			m_n = 0;
		}

		virtual void InitParameters(unsigned long long seed = 123) = 0;
		virtual ~BasePDFMixture(){ m_weights.clear(); m_alphas.clear(); m_model.clear(); }

		double Prob(const Point& x);

		void SetData(PointCollection* data){m_data = data; m_n = m_data->GetNPoints(); m_dim = m_data->Dim(); }
		PointCollection* GetData(){ return m_data; }


		//for EM algorithm
		virtual void CalculatePosterior() = 0;
		virtual void UpdateParameters() = 0;
		//returns mu, cov, and mixing coeffs
		virtual map<string, vector<Matrix>> GetParameters() = 0; 
		
		//for variational EM algorithm
		virtual void InitPriorParameters(unsigned long long seed = 123) = 0;
		virtual void CalculateVariationalPosterior() = 0;
		virtual void UpdateVariationalParameters() = 0;
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart)
		virtual map<string, vector<Matrix>> GetPriorParameters() = 0; 

		void GetMixingCoeffs(vector<double>& weights){ weights.clear(); weights = m_weights; }	
		void GetDirichletParams(vector<double>& alphas){ alphas.clear(); alphas = m_alphas; }

		virtual double EvalLogL() = 0;
		virtual double EvalVariationalLogL() = 0;

		int Dim(){ return m_dim; }
	
		Matrix GetPosterior() const{
			return m_post;
		}

		int GetNClusters(){ return m_k; }		

		int GetMaxPointAssignment(int n){
			int max_post = 0;
			int max_k = 0;
			for(int k = 0; k < m_k; k++)
				if(m_post.at(n,k) > max_post){
					max_post = m_post.at(n,k);
					max_k = k;
				}
			return max_k;
		}

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
		Matrix m_post;
	
		int m_dim;



};
#endif
