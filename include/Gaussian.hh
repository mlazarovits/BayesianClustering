#ifndef GAUSSIAN_HH
#define GAUSSIAN_HH
#include "BasePDF.hh"
#include "Matrix.hh"

class Gaussian : public BasePDF{
	public:
		Gaussian();
		Gaussian(int d);
		Gaussian(Point mu, Matrix cov);
		Gaussian(Matrix mu, Matrix cov);
		virtual ~Gaussian(){ };


		void InitParameters();
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		map<string, vector<Matrix>> GetParameters();
		void UpdateParameters(){ m_mu = m_params["mean"]; m_cov = m_params["cov"]; }	
/*	
		void SetParameters(map<string, Matrix> params){ 
			if(params.find("mean") == params.end()) cout << "Specify Gaussian mean with 'mean'." << endl;
			if(params.find("cov") == params.end()) cout << "Specify Gaussian covariance with 'cov'." << endl;
			if(!params["cov"].square()){ cout << "Error: non-square covariance." << endl; return;}
			if(!params["cov"].symmetric()){ cout << "Error: non-symmetric covariance." << endl; return; }
			m_dim = params["mean"].Dim(0);
			if(m_dim != params["cov"].Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }
			m_mu = params["mean"];
			m_cov = params["cov"];
			m_params = params;	
		}	
*/
		double Prob(const Point& x);


	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
