#ifndef BASEPDF_HH
#define BASEPDF_HH

#include "Point.hh"

class BasePDF{
	public:
		BasePDF(){ };
		virtual ~BasePDF(){ delete m_prior; }	

		virtual double Prob(Point x) = 0;
		virtual void SetPrior(BasePDF* pdf){ m_prior = pdf; }

		


	private:
		BasePDF* m_prior;








};

#endif