#ifndef Wishart_HH
#define Wishart_HH

#include "BasePDF.hh"
#include "Matrix.hh"

class Wishart : public BasePDF{
	public:
		Wishart(){ m_dim = 0; }
		Wishart(const Matrix& W, double nu);
		virtual ~Wishart(){ };

		void SetParameters(Matrix W, double nu);
		void GetParameters(Matrix& W, double& nu){ W = m_W; nu = m_nu; }

		double Prob(const Matrix& x);
		double Prob(const Point& x);

		double lnB();
		double H(); //entropy
	private:
		Matrix m_W;
		double m_nu;


};
#endif
