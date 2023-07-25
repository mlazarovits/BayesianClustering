#ifndef GAUSSIAN_HH
#define GAUSSIAN_HH
#include "BasePDF.hh"
#include "Matrix.hh"

class Gaussian : public BasePDF{
	public:
		Gaussian();
		Gaussian(Point mu, Matrix cov);
		Gaussian(Matrix mu, Matrix cov);
		virtual ~Gaussian(){ };


		void InitParameters();
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		map<string, vector<Matrix>> GetParameters();
		
		void SetParameters(map<string, Matrix> params){ 
			if(params.find("mean") == params.end()) cout << "Specify Gaussian mean with 'mean'." << endl;
			if(params.find("cov") == params.end()) cout << "Specify Gaussian covariance with 'cov'." << endl;
			m_mu = params["mean"];
			m_cov = params["cov"];
			m_params = params;	
		}	
		double Prob(const Point& x);

	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
