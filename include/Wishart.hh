#ifndef Wishart_HH
#define Wishart_HH

#include "BasePDF.hh"
#include "Matrix.hh"

class Wishart : public BasePDF{
	public:
		Wishart(){ m_dim = 0; }
		Wishart(const Matrix& W, double nu);
		virtual ~Wishart(){ };

		void SetParameters(Matrix W, double nu){ m_W = W; m_nu = nu; }
		void InitParameters();
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		map<string, vector<Matrix>> GetParameters();

		double Prob(const Matrix& x);
		double Prob(const Point& x);

		double lnB();
		double H(); //entropy
	private:
		Matrix m_W;
		double m_nu;


};
#endif