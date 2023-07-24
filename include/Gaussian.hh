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


		void SetParameters(const Matrix& mu, const Matrix& cov){ m_mu = mu; m_cov = cov; }		
		void GetParameters(Matrix& mu, Matrix& cov){ mu = m_mu; cov = m_cov; }		
		
		double Prob(const Point& x);

	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
