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
		BasePDFMixture(){ m_k = 0; m_n = 0; m_alpha0 = 0.; _verb = 0; _smear = false;}
		BasePDFMixture(int k){ 
			m_k = k; 
			for(int k = 0; k < m_k; k++){
				m_coeffs.push_back(0.);
				m_norms.push_back(1.);
				m_norms_unwt.push_back(1.);
				m_alphas.push_back(0.);
			}
			m_n = 0;
			//alpha > 0
			//choose the same value for all alpha_0k by symmetry (see Bishop eq. 10.39)
			m_alpha0 = 0.1; _verb = 0; _smear = false;
		}

		//virtual void InitParameters(unsigned long long seed = 123) = 0;
		virtual ~BasePDFMixture(){ m_coeffs.clear(); m_alphas.clear(); m_model.clear(); }

		void SetVerbosity(int v){ _verb = v; }
		double Prob(const Point& x);
		double Prob(const PointCollection& x);

		void SetData(PointCollection* data){
			m_data = data; m_n = m_data->GetNPoints(); m_dim = m_data->Dim(); 
			if(data->GetNPoints() < m_k){
				//remove extra models
				for(int i = 0; i < m_k - data->GetNPoints(); i++) RemoveModel(i);
				m_k = data->GetNPoints();
			}
		}
		PointCollection* GetData(){ return m_data; }

		//estimates data points as Gaussians with mean = pt and covariance = set here
		void SetDataSmear(const Matrix& cov){ _data_cov = cov; _smear = true; }

		//for EM algorithm
		virtual void CalculatePosterior() = 0;
		virtual void UpdateParameters() = 0;
		//returns mu, cov, and mixing coeffs for cluster k
		virtual map<string, Matrix> GetParameters(int k) = 0; 
		
		//for variational EM algorithm
		virtual void SetPriorParameters(map<string, Matrix> params) = 0;
		virtual void InitPriorParameters(unsigned long long seed = 123) = 0;
		virtual void CalculateVariationalPosterior() = 0;
		virtual void UpdateVariationalParameters() = 0;
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		virtual map<string, Matrix> GetPriorParameters(int k) = 0; 

		void GetMixingCoeffs(vector<double>& coeffs){ coeffs.clear(); coeffs = m_coeffs; }	
		void GetDirichletParameters(vector<double>& alphas){ alphas.clear(); alphas = m_alphas; }
		void SetDirichletParameter(double alpha){ m_alpha0 = alpha; }
		void SetAlpha(double alpha){  m_alpha0 = alpha; }

		BasePDF* GetModel(int k){ return m_model[k]; }
		void RemoveModel(int k){
			//erase model
			m_model.erase(m_model.begin()+k); 
			//erase corresponding dirichlet param
			m_alphas.erase(m_alphas.begin()+k);
			//erase corresponding number of associated points (N_k)
			m_norms.erase(m_norms.begin()+k);
			//don't need to update m_coeffs - only used for normal EM algorithm
			//update number of clusters
			m_k--; }

		//removes components that do not contribute to overall likelihood
		void UpdateMixture(double thresh){
		//if Dirichlet parameter (m_alpha) is below some threshold, remove cluster
			for(int k = 0; k < m_k; k++){
				if(m_alphas[k] < thresh){
					if(_verb > 1) cout << "Removing cluster " << k << " with alpha: " << m_alphas[k] << endl;
					//remove model + update number of clusters
					RemoveModel(k);
					k--; //make sure to check following model
				}
			}
			//cout if all clusters are removed
			if(m_k < 1) cout << "Error: all clusters have been removed. Update threshold accordingly." << endl;

		}


		virtual double EvalLogL() = 0;
		virtual double EvalVariationalLogL() = 0;

		int Dim(){ return m_dim; }
	
		Matrix GetPosterior() const{
			return m_post;
		}

		int GetNClusters(){ return m_k; }		

		int GetMaxPointAssignment(int n){
			if(m_post.empty()) return -1;
			int max_post = 0;
			int max_k = 0;
			for(int k = 0; k < m_k; k++)
				if(m_post.at(n,k) > max_post){
					max_post = m_post.at(n,k);
					max_k = k;
				}
			return max_k;
		}
		
		//ie gets subcluster average energy (within subcluster) - returns average weight per subcluster
		void GetAvgWeights(vector<double>& v){
			v.clear(); for(int k = 0; k < m_k; k++){ 
			v.push_back(m_norms[k]/m_norms_unwt[k]); }
		}  
		void GetAvgVarWeights(vector<double>& v){
			for(int k = 0; k < m_k; k++) v.push_back(m_norms[k]/m_norms_unwt[k]);

			sortByVarCoeff(v);

		}  


		void sortByVarCoeff(vector<double>& v){
			if((int)v.size() != m_k) return;
			vector<pair<int,double>> mweight;
			for(int k = 0; k < m_k; k++) mweight.push_back(std::make_pair(k,(m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_n)));
			//sort by mixing coeff weight - insertion sort
			int i, j;
			pair<int, double> t;
			for(i = 2; i < m_k; i++){
				t = mweight[i]; j = i;
				while(mweight[j-1].second > t.second){
					mweight[j] = mweight[j-1]; j--;
				}
				mweight[j] = t;
			}

			v.clear();
			vector<double> w;
			int idx; 
			for(int k = 0; k < m_k; k++){
				idx = mweight[k].first; 
				w.push_back(v[idx]);
			}
			v = w;
		}


		void GetNormsUnwt(vector<double>& v){
			for(int k = 0; k < m_k; k++) v.push_back(m_norms_unwt[k]);
			sortByVarCoeff(v);
		}

		PointCollection* m_data;
		//number of data points
		int m_n;
		//number of components
		int m_k;
		//probabilistic model
		vector<BasePDF*> m_model;
		//mixture coeffs (probabilities sum to 1, multinomal dist.)
		vector<double> m_coeffs;
		//dirichlet (prior) parameters (on coeffs)
		vector<double> m_alphas;
		//initial alpha
		double m_alpha0;
		//normalization on posterior
		vector<double> m_norms;
		vector<double> m_norms_unwt;
		Matrix m_post;
	
		int m_dim;

		int _verb;
		//data smear
		Matrix _data_cov;
		bool _smear;

};
#endif
