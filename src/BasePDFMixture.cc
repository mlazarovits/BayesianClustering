#include "BasePDFMixture.hh"

double BasePDFMixture::Prob(const Point& x){
	double ret;
	for(int i = 0; i < m_k; i++)
		ret += m_weights[i]*m_model[i]->Prob(x);
	return ret;
}
