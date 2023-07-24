#include <boost/math/special_functions/gamma.hpp>
#include "Dirichlet.hh"

using boost::math::tgamma;
using boost::math::lgamma;

double Dirichlet::Prob(const Point& x){
	//x needs to be k-dimensional
	if(x.Dim() != m_dim){
		cout << "Error: x dimensions"  << x.Dim() << " and PDF dimension " << (int)m_alphas.size() << " incompatible." << endl;
		return -1;
	}

	//make sure sum over x dim is one
	double norm = 0;
	for(int i = 0; i < m_dim; i++){
		norm += x.at(i);
		if(m_alphas[i] < 0){
			cout << "Error: alpha #" << i << " is less than zero: " << m_alphas[0] << endl;
			return -1;
		}
	}
	if(norm != 1){
		cout << "Error: x is not on standard K - 1 simplex (sum over k != 1): " << norm << endl;
		return -1;
	}

	double num = 1.;
	double denom = 0.;
	double pdf = 1.;
	for(int k = 0; k < m_dim; k++){
		num *= tgamma(m_alphas[k]);
		denom += m_alphas[k];
		pdf *= pow(x.at(k), m_alphas[k] - 1);
	}
	denom = tgamma(denom);
	//normalizing constant for PDF now
	norm = 0;
	norm = num/denom;

	return pdf/norm;

}



double Dirichlet::lnC(){
	double alpha_hat = 0;
	double norm = 0;
	
	for(int i = 0; i < m_dim; i++){
		alpha_hat += m_alphas[i];
		norm += lgamma(m_alphas[i]);
	}
	return lgamma(alpha_hat) - norm;

}



void Dirichlet::InitParameters(){
	//init alphas to 1
	for(int i = 0; i < m_dim; i++){
		m_alphas.push_back(1.);
	}
}


map<string, vector<Matrix>> Dirichlet::GetParameters(){
	map<string, vector<Matrix>> ret;
	ret["alphas"] = {Matrix(m_alphas)};
	return ret;

}
