#include "TriangularPDF.hh"
#include "RandomSample.hh"

TriangularPDF::TriangularPDF(){ _a = 0; _b = 0; _c = 0; m_dim = 1;
	m_params["a"] = Matrix(Point(_a));
	m_params["b"] = Matrix(Point(_b));
	m_params["c"] = Matrix(Point(_c));
}

TriangularPDF::TriangularPDF(int d){ 
	if(d != 1){ cout << "Error: triangular PDF is only defined for univariate case." << endl; }
	_a = 0; _b = 0; _c = 0;
	m_dim = 1;

	m_params["a"] = Matrix(Point(_a));
	m_params["b"] = Matrix(Point(_b));
	m_params["c"] = Matrix(Point(_c));
}

//x is nonzero within a < x < b
double TriangularPDF::Prob(double x){
	if(x == _c) return 2./(_b - _a);
	else if(x >= _a && x < _c) return (2.*(x - _a))/( (_b - _a)*(_c - _a) );
	else if(x >= _c && x < _b) return (2.*(_b - x))/( (_b - _a)*(_b - _c) );
	else return 0;  
}

double TriangularPDF::Prob(const Point& x){
	if(x.Dim() != m_dim){ cout << "Error: point dimension is " << x.Dim() << ". Should be d = 1." << endl; return -999; }

	return Prob(x.at(0));
}

double TriangularPDF::Prob(const PointCollection& x){
	if(x.Dim() != m_dim){ cout << "Error: point dimension is " << x.Dim() << ". Should be d = 1." << endl; return -999; }

	double ret = 1;
	for(int i = 0; i < x.GetNPoints(); i++) ret *= Prob(x.at(i));

	return ret;
}

void TriangularPDF::InitParameters(unsigned long long seed){
	_a = 0;
	_b = 1;
	_c = 0;

	m_params["a"] = Matrix(Point(_a));
	m_params["b"] = Matrix(Point(_b));
	m_params["c"] = Matrix(Point(_c));

}

