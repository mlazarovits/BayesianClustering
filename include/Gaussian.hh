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
		
		double Prob(const Point& x);

	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
