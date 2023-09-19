#include "UniformPDF.hh"


UniformPDF::UniformPDF(){ _a = 0; _b = 0; m_dim = 0; 
	m_params["a"] = Matrix(Point(_a));
	m_params["b"] = Matrix(Point(_b));


}

UniformPDF::UniformPDF(int d){
	if(d != 1){ cout << "Error: uniform PDF is only defined for univariate case." << endl; }
	_a = 0; _b = 1; m_dim = d;
	m_params["a"] = Matrix(Point(_a));
	m_params["b"] = Matrix(Point(_b));
}


double UniformPDF::Prob(double x){
	if(x < _a || x > _b) return 0;
	else return 1./(_b - _a);

}

double UniformPDF::Prob(const Point& x){
	if(x.Dim() != m_dim){ cout << "Error: point dimension is " << x.Dim() << ". Should be d = 1." << endl; return -999; }
	
	return Prob(x.at(0));
}

double UniformPDF::Prob(const PointCollection& x){
	if(x.Dim() != m_dim){ cout << "Error: point dimension is " << x.Dim() << ". Should be d = 1." << endl; return -999; }

	double ret = 1;
	for(int i = 0; i < x.GetNPoints(); i++) ret *= Prob(x.at(i));

	return ret;
}

void UniformPDF::InitParameters(unsigned long long seed){
	_a = 0;
	_b = 1;

	m_params["a"] = Matrix(Point(_a));
	m_params["b"] = Matrix(Point(_b));

}
