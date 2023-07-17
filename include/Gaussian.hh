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
		
		double Prob(const Point& x);

	private:
		int m_dim;
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
