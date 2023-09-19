#ifndef UNIFORMPDF_HH
#define UNIFORMPDF_HH

#include "BasePDF.hh"

class UniformPDF : public BasePDF{
	public:
		UniformPDF();
		UniformPDF(int d);
		UniformPDF(double a, double b){ _a = a; _b = b; 
			m_params["a"] = Matrix(Point(_a));
			m_params["b"] = Matrix(Point(_b));
		}
		virtual ~UniformPDF(){ };

		double Prob(double x);
		double Prob(const Point& x);
		double Prob(const PointCollection& x);

		void InitParameters(unsigned long long seed = 123);

		void UpdateParameters(){ _a = m_params["a"].at(0,0); _b = m_params["b"].at(0,0);
			m_params["a"] = Matrix(Point(_a));
			m_params["b"] = Matrix(Point(_b));
		 }	

	private:
		double _a, _b;


};

#endif
