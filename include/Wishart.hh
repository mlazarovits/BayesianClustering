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
		void InitParameters(unsigned long long seed = 123);
		
		BasePDF* mult(BasePDF* p1){ return nullptr; }
		void UpdateParameters(){ 
			m_nu = m_params["nu"].at(0,0);
			m_W = m_params["W"];
		}	
		BasePDF* Posterior(){ return nullptr; }
		double Prob(const Matrix& x);
		double Prob(const BayesPoint& x);
		double Prob(const PointCollection& x){ return -1.; }

		double lnB();
		double H(); //entropy

		void SetDim(int d){
			if(m_dim != 0) return;
			BaseSetDim(d);
			m_W = Matrix(d,d);
			m_nu = 0;
		}

	private:
		Matrix m_W;
		double m_nu;


};
#endif
