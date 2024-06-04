#include "BasePDFMixture.hh"

double BasePDFMixture::Prob(const BayesPoint& x){
	double ret = 0.;
	for(int i = 0; i < m_k; i++)
		ret += m_coeffs[i]*m_model[i]->Prob(x);
	return ret;
}

double BasePDFMixture::Prob(const PointCollection& x){
	double ret = 1.;
	int n = x.GetNPoints();
	for(int i = 0; i < n; i++)
		ret *= Prob(x.at(i));	
	return ret;
}
