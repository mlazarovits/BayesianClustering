#ifndef TriangularPDF_HH
#define TriangularPDF_HH

#include "BasePDF.hh"
#include "Point.hh"

class TriangularPDF : public BasePDF{
	public:
		TriangularPDF();
		TriangularPDF(int d);
		TriangularPDF(double a, double b, double c){
			if(b <= a || c < a || c > b){ cout << "Error: c: a <= c <= b and b: a < b" << endl; cout << "a: " << a << " b: " << b << " c: " << c << endl;}
			_a = a; _b = b; _c = c;
			m_params["a"] = Matrix(BayesPoint(_a));
			m_params["b"] = Matrix(BayesPoint(_b));
		}
		virtual ~TriangularPDF(){ };
		std::unique_ptr<BasePDF> clone() const override{
			return std::unique_ptr<BasePDF>(new TriangularPDF(*this));
		}

		double Prob(const BayesPoint& x) override;
		double Prob(double x);
		double Prob(const PointCollection& x) override;

		void InitParameters(unsigned long long seed = 123) override;
		void UpdateParameters() override{ _a = m_params["a"].at(0,0); _b = m_params["b"].at(0,0); _c = m_params["c"].at(0,0);
			m_params["a"] = Matrix(BayesPoint(_a));
			m_params["b"] = Matrix(BayesPoint(_b));
			m_params["c"] = Matrix(BayesPoint(_c));
		 }	

		void SetDim(int d) override{
			if(m_dim != 0) return;
			BaseSetDim(d);
		}

	private:
		double _a, _b, _c;

};
#endif
