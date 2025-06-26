#ifndef UNIFORMPDF_HH
#define UNIFORMPDF_HH

#include "BasePDF.hh"

class UniformPDF : public BasePDF{
	public:
		UniformPDF();
		UniformPDF(int d);
		UniformPDF(double a, double b){ _a = a; _b = b; 
			m_params["a"] = Matrix(BayesPoint(_a));
			m_params["b"] = Matrix(BayesPoint(_b));
		}
		virtual ~UniformPDF(){ };

		double Prob(double x);
		double Prob(const BayesPoint& x);
		double Prob(const PointCollection& x);

		void InitParameters(unsigned long long seed = 123);

		void UpdateParameters(){ _a = m_params["a"].at(0,0); _b = m_params["b"].at(0,0);
			m_params["a"] = Matrix(BayesPoint(_a));
			m_params["b"] = Matrix(BayesPoint(_b));
		 }	

		void SetDim(int d){
			if(m_dim != 0) return;
			BaseSetDim(d);
		}

	private:
		double _a, _b;


};

#endif
