#ifndef BASEPDF_HH
#define BASEPDF_HH

#include "Point.hh"
#include "Matrix.hh"
#include <string>
using std::string;


class BasePDF{
	public:
		BasePDF(){ };
		virtual ~BasePDF(){ delete m_prior; }	

		virtual double Prob(const Point& x) = 0;
		virtual void SetPrior(BasePDF* pdf){ m_prior = pdf; }

		virtual void InitParameters() = 0;
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		virtual map<string, vector<Matrix>> GetParameters() = 0;
		
		int m_dim;
		BasePDF* m_prior;








};

#endif
