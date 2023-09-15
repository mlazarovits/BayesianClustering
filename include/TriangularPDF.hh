#ifndef TriangularPDF_HH
#define TriangularPDF_HH

#include "BasePDF.hh"
#include "Point.hh"


class TriangularPDF : public BasePDF{
	public:
		TriangularPDF();
		TriangularPDF(int d);
		TriangularPDF(double a, double b, double c){
			if(b >= a || c > a){ cout << "Error: c: a <= c <= b and b: a < b" << endl; }
			_a = a; _b = b; _c = c;
		}
		virtual ~TriangularPDF(){ };

		double Prob(const Point& x);
		double Prob(double x);
		double Prob(const PointCollection& x);

		void InitParameters(unsigned long long seed = 123);
		void UpdateParameters(){ _a = m_params["a"].at(0,0); _b = m_params["b"].at(0,0); _c = m_params["c"].at(0,0);
			m_params["a"] = Matrix(Point(_a));
			m_params["b"] = Matrix(Point(_b));
			m_params["c"] = Matrix(Point(_c));
		 }	

	private:
		double _a, _b, _c;

};
#endif
