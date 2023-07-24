#include "BasePDFMixture.hh"

double Prob(const Point& x){
	double ret;
	for(int i = 0; i < k; i++)
		ret += m_weights[i]*m_model[i]->Prob(x);
	return ret;
}
