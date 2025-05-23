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
	
		void InitParameters(map<string, Matrix> priors = {}, unsigned long long seed = 111);
		//E-step
		void CalculatePosterior();
		//M-step
		void UpdateParameters();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();

		//fill vectors with parameters
		//returns mu, cov, and mixing coeffs for cluster k
		map<string, Matrix> GetLikelihoodParameters(int k); 

		/* 
		void SetPriorParameters(map<string, Matrix> params){
			m_beta0 = params["scale"].at(0,0);
			m_nu0 = params["dof"].at(0,0);
			m_W0 = params["scalemat"];
			m_W0inv.invert(m_W0);
			m_mean0 = params["mean"];
			m_meanBeta0 = Matrix(m_dim, 1);
			m_meanBeta0.mult(m_mean0, m_beta0);
			if(_verb > 6){
				cout << "Prior Parameters" << endl;
				cout << "beta0" << endl;
				params["scale"].Print();
				cout << "mean0" << endl;
				m_mean0.Print();
				cout << "nu0" << endl;
				params["dof"].Print();
				cout << "W0" << endl;
				m_W0.Print();
			} 
		}
		*/

		//shift learned model parameters
		void ShiftParameters(const BayesPoint& pt){
			//only need to shift Gaussian means + prior mean
			Matrix mean, priormean;
			Matrix shift;
			shift.PointToMat(pt);
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				PointCollection meanpts = mean.MatToPoints();
				meanpts.Translate(pt.at(0),0);
				meanpts.CircularTranslate(pt.at(1),1);
				meanpts.Translate(pt.at(2),2);
				m_model[k]->SetParameter("mean",Matrix(meanpts));
				
				//translate posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				meanpts = priormean.MatToPoints();
				meanpts.Translate(pt.at(0),0);
				meanpts.CircularTranslate(pt.at(1),1);
				meanpts.Translate(pt.at(2),2);
				m_model[k]->GetPrior()->SetParameter("mean",Matrix(meanpts));


				//shift data statistics
				meanpts = _xbar[k].MatToPoints();
				meanpts.Translate(pt.at(0),0);
				meanpts.CircularTranslate(pt.at(1),1);
				meanpts.Translate(pt.at(2),2);
				_xbar[k] = Matrix(meanpts);
			}
		}
		
		void ThetaToEta_params(){
			Matrix mean, priormean;
			double theta;
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				theta = mean.at(0,0);
				mean.SetEntry(-log(tan(theta/2)),0,0);
				m_model[k]->SetParameter("mean",mean);
				
				//posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				theta = priormean.at(0,0);	
				priormean.SetEntry(-log(tan(theta/2)),0,0);
				m_model[k]->GetPrior()->SetParameter("mean",priormean);


				//data statistics
				theta = _xbar[k].at(0,0);	
				_xbar[k].SetEntry(-log(tan(theta/2)),0,0);
			}
		}
		
		void UnprojectTheta_params(){
			Matrix mean, priormean;
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				PointCollection mean_pt = mean.MatToPoints();
				//cout << "mean" << endl; mean_pt.Print();
				mean_pt.PlaneToAngleProject(0);	
				//cout << "unproj mean" << endl; mean_pt.Print();
				m_model[k]->SetParameter("mean",Matrix(mean_pt));
				
				//translate posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				PointCollection priormean_pt = priormean.MatToPoints();
				priormean_pt.PlaneToAngleProject(0);	
				m_model[k]->GetPrior()->SetParameter("mean",Matrix(priormean_pt));


				//shift data statistics
				PointCollection xbar_pt = _xbar[k].MatToPoints();
				xbar_pt.PlaneToAngleProject(0);	
				_xbar[k] = Matrix(xbar_pt);	
			}
		}
	
		void UnprojectPhi_params(){
			Matrix mean, priormean;
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				PointCollection mean_pt = mean.MatToPoints();
				//cout << "mean" << endl; mean_pt.Print();
				mean_pt.PlaneToAngleProject(1);	
				//cout << "unproj mean" << endl; mean_pt.Print();
				m_model[k]->SetParameter("mean",Matrix(mean_pt));
				
				//translate posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				PointCollection priormean_pt = priormean.MatToPoints();
				priormean_pt.PlaneToAngleProject(1);	
				m_model[k]->GetPrior()->SetParameter("mean",Matrix(priormean_pt));


				//shift data statistics
				PointCollection xbar_pt = _xbar[k].MatToPoints();
				xbar_pt.PlaneToAngleProject(1);	
				_xbar[k] = Matrix(xbar_pt);	
			}
		}

	
		void PutPhi02pi_params(){
			//put means on [0,2pi]
			Matrix mean, priormean;
			for(int k = 0; k < m_k; k++){
				mean = m_model[k]->GetParameter("mean");
				PointCollection mean_pt = mean.MatToPoints();
				//cout << "mean" << endl; mean_pt.Print();
				mean_pt.Put02pi(1);	
				//cout << "02pi mean" << endl; mean_pt.Print();
				m_model[k]->SetParameter("mean",Matrix(mean_pt));
				
				//translate posterior mean in prior distribution
				priormean = m_model[k]->GetPrior()->GetParameter("mean");
				PointCollection priormean_pt = priormean.MatToPoints();
				priormean_pt.Put02pi(1);	
				m_model[k]->SetParameter("mean",Matrix(priormean_pt));


				//shift data statistics
				PointCollection xbar_pt = _xbar[k].MatToPoints();
				xbar_pt.Put02pi(1);
				_xbar[k] = Matrix(xbar_pt);	
			}
		}
	
		//scale prior parameters
		void ScalePriorParameters(Matrix sc){
			//scale means + matrices
			Matrix mean, priormean, cov, W;
			Matrix scT, scInv, scTinv;
			scT.transpose(sc);
			scInv.invert(sc);
			scTinv.transpose(scInv);
			
			m_W0.mult(sc,m_W0); //Avar(X)
			m_W0.mult(m_W0,scT); //Avar(X)A^T
			m_W0inv.invert(m_W0);
			
			m_mean0.mult(sc,m_mean0);	
			m_meanBeta0.mult(m_mean0, m_beta0);
				/*
				cout << "Scaled prior Parameters" << endl;
				cout << "beta0" << endl;
				cout << m_beta0 << endl;
				cout << "mean0" << endl;
				m_mean0.Print();
				cout << "nu0" << endl;
				cout << m_nu0 << endl;
				cout << "W0" << endl;
				m_W0.Print();
			*/
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

				//scale r stat mean
				_xbar[k].mult(sc,_xbar[k]);
	

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


		
		void RemoveModel(int j){
			//if only one cluster, don't remove
			if(m_k > 1 && _xbar.size() > 1){
				_xbar.erase(_xbar.begin()+j);
				_Sbar.erase(_Sbar.begin()+j);
				m_Elam.erase(m_Elam.begin()+j);
				m_Epi.erase(m_Epi.begin()+j);
			}
			BaseRemoveModel(j);
		}


	protected:	
		//for variational EM algorithm
		//this needs to be called BEFORE init parameters because InitParameters can remove subclusters if the k-means algo
		//finds nothing assigned to the subcluster BUT it relies on the priors still
		//bc the UpdateMixture method also updates the posterior parameters (ie the dims of m_post, etc)
		void InitPriorParameters(map<string, Matrix> params = {});


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
		
		void InitParameters(unsigned long long seed){ };

		

};


#endif
