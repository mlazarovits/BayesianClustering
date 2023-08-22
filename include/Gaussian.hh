#ifndef GAUSSIAN_HH
#define GAUSSIAN_HH
#include "BasePDF.hh"
#include "Matrix.hh"
#include "NormalWishart.hh"

class Gaussian : public BasePDF{
	public:
		Gaussian();
		Gaussian(int d);
		Gaussian(Point mu, Matrix cov);
		Gaussian(Matrix mu, Matrix cov);
		virtual ~Gaussian(){ };


		void InitParameters(unsigned long long seed = 123);
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		void UpdateParameters(){ m_mu = m_params["mean"]; m_cov = m_params["cov"]; }	
		double Prob(const Point& x);
		double Prob(const PointCollection& x);
		NormalWishart* Posterior();

		Gaussian* mult(Gaussian* p1); 


	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
