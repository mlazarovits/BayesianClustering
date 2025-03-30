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
		BasePDFMixture(){
			m_k = 0; m_n = 0; m_alpha0 = 0.; _verb = 0; _smear = false; m_post.SetDims(m_n, m_k);
			_cell = 0;
			_tresCte = 0;
			_tresStoch = 0;
			_tresNoise = 0;
		}
		BasePDFMixture(int k){ 
			m_k = k; 
			for(int k = 0; k < m_k; k++){
				m_coeffs.push_back(0.);
				m_norms.push_back(1.);
				m_alphas.push_back(0.);
				m_model.push_back(nullptr);
			
			}
			m_n = 0;
			//alpha > 0
			//choose the same value for all alpha_0k by symmetry (see Bishop eq. 10.39)
			m_alpha0 = 0.1; _verb = 0; _smear = false;
			m_post.SetDims(m_n, m_k);
			_cell = acos(-1)/180; //default is CMS ECAL cell size
			//default is EGamma res - https://github.com/jking79/GammaResTool/blob/main/macros/ecal_config/caliSmearConfig.txt
			//times should already be given in ns
			_tresCte = 0.133913;
			_tresStoch = 1.60666; 
			_tresNoise = 0.00691415;
		}

		//virtual void InitParameters(unsigned long long seed = 123) = 0;
		virtual ~BasePDFMixture(){ m_coeffs.clear(); m_alphas.clear();
			for(auto& pointer : m_model)
				delete pointer;
		}

		void SetVerbosity(int v){ _verb = v; }
		double Prob(const BayesPoint& x);
		double Prob(const PointCollection& x);

		double _cell, _tresCte, _tresStoch, _tresNoise;
		void SetMeasErrParams(double spatial, double tresCte, double tresStoch, double tresNoise){
			_cell = spatial;
			_tresCte = tresCte;
			_tresStoch = tresStoch;
			_tresNoise = tresNoise;
		}
		
		void SetData(PointCollection* data){
			if(_verb > 1) cout << "Using tres_cte = " << _tresCte << " ns, tres_stoch = " << _tresStoch << " ns and tres_noise = " << _tresNoise << endl;
			m_data = data; 
			m_n = m_data->GetNPoints(); 
			m_dim = m_data->Dim(); 
			m_post.SetDims(m_n, m_k);
			if(data->GetNPoints() < m_k){
				//remove extra models
				for(int i = 0; i < m_k - data->GetNPoints(); i++) RemoveModel(i);
				m_k = data->GetNPoints();
			}
			//set up lambda_n
			double tresSq;
			for(int n = 0; n < m_n; n++){
				Matrix lamStar(m_dim, m_dim);
				//need to make sure tResRate = tResRate_true*gev for units to match
				tresSq = _tresCte*_tresCte + _tresStoch*_tresStoch/(m_data->at(n).w()) + _tresNoise*_tresNoise/(m_data->at(n).w()*m_data->at(n).w());
				tresSq /= 2;
				lamStar.SetEntry(1/(_cell*_cell),0,0);
				lamStar.SetEntry(1/(_cell*_cell),1,1);
				lamStar.SetEntry(1/(tresSq),2,2);
				if(_verb > 3){ cout << "point " << n << " has weight " << m_data->at(n).w() << " and sigma_t " << sqrt(tresSq) << " ns " << endl; m_data->at(n).Print(); cout << "lamstar" << endl;lamStar.Print();}
				_lamStar.push_back(lamStar); //r = 1 for all k on init
				
			}
			for(int k = 0; k < m_k; k++){
				m_model[k]->SetDim(m_dim);
				m_model[k]->GetPrior()->SetDim(m_dim);
			}
		}
		PointCollection* GetData(){ return m_data; }

		//estimates data points as Gaussians with mean = pt and covariance = set here
		void SetDataSmear(const Matrix& cov){ _data_cov = cov; _smear = true; }
	
		//for EM algorithm
		virtual void CalculatePosterior() = 0;
		virtual void UpdateParameters() = 0;
		//returns mu, cov, and mixing coeffs for cluster k
		virtual map<string, Matrix> GetLikelihoodParameters(int k) = 0; 
		
		//for variational EM algorithm
		virtual void SetPriorParameters(map<string, Matrix> params) = 0;
		//virtual void SetJetPriorParameters(map<string, Matrix>& params) = 0;
		//virtual void SetJetParameters(map<string, Matrix>& params) = 0;
		//virtual void GetJetParameters(map<string, Matrix>& params) = 0;
		virtual void InitPriorParameters(unsigned long long seed = 123) = 0;
		virtual void CalculateVariationalPosterior() = 0;
		virtual void UpdateVariationalParameters() = 0;
		//returns params on priors (alpha, W, nu, m, beta - dirichlet + normalWishart) for cluster k
		virtual map<string, Matrix> GetLHPosteriorParameters(int k) const = 0; 
		virtual map<string, Matrix> GetDataStatistics(int k) const = 0; 
		double GetPi(int k){ return (m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_data->Sumw()); }
		double GetCoeff(int k){ return m_coeffs[k]; }

		void GetMixingCoeffs(vector<double>& coeffs){ coeffs.clear(); coeffs = m_coeffs; }	
		void GetDirichletParameters(vector<double>& alphas){ alphas.clear(); alphas = m_alphas; }
		void SetDirichletParameter(double alpha){ m_alpha0 = alpha; }
		void SetAlpha(double alpha){  m_alpha0 = alpha; }

		BasePDF* GetModel(int k){ return m_model[k]; }
		void RemoveModel(int j){
			//if this is removing the last cluster
			//break and do resetting of posterior in UpdateModel()
			if(m_k - 1 == 0) return;
			//erase model
			//m_model.erase(m_model.begin()+j); 
			//for(int i = 0; i < m_model.size(); i++){
			for(auto& pointer : m_model){
				//avoid dangling pointers
				if(pointer == m_model[j]){
					delete pointer;
					pointer = nullptr;
					break;
				}
			}			
			m_model.erase(std::remove(m_model.begin(), m_model.end(), nullptr), m_model.end());

			//erase corresponding dirichlet param
			m_alphas.erase(m_alphas.begin()+j);
			//erase corresponding number of associated points (N_k)
			m_norms.erase(m_norms.begin()+j);
			//don't need to update m_coeffs - only used for normal EM algorithm
			//update dimensions of posterior
			Matrix newpost = Matrix(m_n,m_k-1);
			int l = 0;
			for(int k = 0; k < m_k; k++){
				if(k == j){ l++; continue; }
				for(int n = 0; n < m_n; n++){
					newpost.SetEntry(m_post.at(n,k),n,k-l);
				}
			}
		
			m_post = newpost;
			//update number of clusters
			m_k--;
		 }

		//removes components that do not contribute to overall likelihood
		void UpdateMixture(double thresh){
		//if Dirichlet parameter (m_alpha) is below some threshold, remove cluster
		//if(m_k > 1){ for(int k = 0; k < m_k; k++) cout << "cluster " << k << " has " << m_norms[k] + m_alpha0 << " points - norm " << m_norms[k] << endl; m_post.Print(); if(m_n < 3) m_data->Print(); }
	//	cout << "points" << endl; m_data->Print();
	//	cout << "weights" << endl; for(int n = 0; n < m_n; n++) cout << m_data->at(n).w() << endl;
			for(int k = 0; k < m_k; k++){
				//alpha_k = norms_k + alpha0 -> may need to remove before all parameters have been updated
				//if(m_norms[k] + m_alpha0 < thresh){
				if(m_norms[k] < thresh){
					if(_verb > 1) 
						cout << "Removing cluster " << k << " with alpha " << m_alphas[k] << " and norm " << m_norms[k] << endl;
					//remove model + update number of clusters
					RemoveModel(k);
					//if the above call removes all clusters
					//break and continue to set single gaussian
					if(m_k == 1) break;
					k--; //make sure to check following model
				}
			}
			//if all clusters removed -> set model to single gaussian
			if(m_k < 1) m_k = 1;//cout << "Error: all clusters for " << m_n << " points have been removed. Update threshold accordingly." << endl; 
			//if single gaussian -> all points need to be under this gaussian
			if(m_k == 1){
				for(int n = 0; n < m_n; n++){
					if(m_post.at(n,0) < m_data->at(n).w()){
						m_post.SetEntry(m_data->at(n).w(), n, 0);
					}
				}
			}
			//update effective counts + corresponding parameters - the corresponding 
			//entry r_nk is the point weight
			UpdateVariationalParameters();

		}


		virtual double EvalLogL() = 0;
		virtual double EvalVariationalLogL() = 0;
		virtual double EvalEntropyTerms() = 0;
		virtual double EvalNLLTerms() = 0;

		int Dim(){ return m_dim; }
	
		Matrix GetPosterior() const{
			return m_post;
		}

		int GetNClusters() const{ return m_k; }		

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
		

		void GetNorms(vector<double>& v){
			v.clear();
			for(int k = 0; k < m_k; k++) v.push_back(m_norms[k]);
		}
	
	
		//sorts smallest first - ascending order
		void SortIdxs(vector<int>& idxs){
			vector<pair<int,double>> mweight;
			for(int k = 0; k < m_k; k++) mweight.push_back(std::make_pair(k,(m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_n)));
			//sort by mixing coeff weight - insertion sort
			int i, j;
			pair<int, double> t;
	//		if(m_k < 3){
				if(m_k < 2){ idxs = {0};
		//		if(mweight[0].second < mweight[1].second){
	//				idxs = {1, 0};
	//			}
	//			else{ idxs = {0, 1}; }
				return;
			}
			for(i = 1; i < m_k; i++){
				t = mweight[i]; j = i;
				while(mweight[j-1].second > t.second){
					mweight[j] = mweight[j-1]; j--;
				}
				mweight[j] = t;
			}

			
			int idx;	
			for(int k = 0; k < m_k; k++){
				idx = mweight[k].first; 
				idxs.push_back(idx);
			}

		}

		//shift data 
		void ShiftData(const BayesPoint& pt){
			m_data->Translate(pt.at(0),0);
			//m_data->Translate(pt.at(1),1);
			m_data->CircularTranslate(pt.at(1),1);
			m_data->Translate(pt.at(2),2);
			
		}
		virtual void PutPhi02pi_params() = 0;

		void PutPhi02pi(){
			m_data->Put02pi(1);
			PutPhi02pi_params();
		}
		//scale data 
		void ScaleData(Matrix sc){
			Matrix x;
			x.PointsToMat(*m_data);
			//Matrix scaled_data(x.GetDims()[0], x.GetDims()[1]);
			x.mult(sc,x);
			//get weights from original data
			vector<double> ws;
			m_data->GetWeights(ws);
			m_data = new PointCollection(x.MatToPoints());
			m_data->SetWeights(ws);

			//also scale smear - if datacov is set
			if(_smear){
				_data_cov.mult(sc,_data_cov);
				Matrix scT;
				scT.transpose(sc);
				_data_cov.mult(_data_cov,scT);
			}
			//measurement err ~ lambda is a precision matrix
			//therefore, it needs to be inverted to be in correct units (ie cell^2 instead of cell^-2 or ns^2 instead of ns^-2)
			else{ //else scale meas err lamba*_n
				for(int n = 0; n < m_n; n++){
					//if(n == 0){ cout << "scale mat " << endl; sc.Print(); cout << "data pt #" << n << " w " << m_data->at(n).w() << "pre scale" << endl; _lamStar[n].Print();}

					Matrix sigStar(m_dim, m_dim);
					sigStar.invert(_lamStar[n]);
					
					sigStar.mult(sc,sigStar);
					Matrix scT;
					scT.transpose(sc);
					sigStar.mult(sigStar,scT);
					_lamStar[n].invert(sigStar);
					//if(n == 0){ cout << "post scale" << endl; _lamStar[n].Print();}
				}

			}
		}

		//shift learned parameters
		virtual void ShiftParameters(const BayesPoint& pt) = 0;
		//scale learned parameters 
		virtual void ScaleParameters(Matrix sc) = 0;

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
		Matrix m_post;

	
		int m_dim;

		int _verb;
		//data smear
		Matrix _data_cov;
		bool _smear;
		vector<Matrix> _lamStar; //lamStar_n - measurement error

};
#endif
