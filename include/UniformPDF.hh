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


	private:
		double _a, _b;


};

#endif
