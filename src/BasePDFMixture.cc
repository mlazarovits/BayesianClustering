#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include "BasePDFMixture.hh"

using boost::math::tgamma;
using boost::math::digamma;

//gaussian for one data point
double BasePDFMixture::Gaus(const Point& x, const Matrix& mu, const Matrix& cov){
	double det = cov.det();
	int dim = x.Dim();
	if(dim != mu.GetDims()[0]){
		cout << "Error: x and mu length do not match. " << dim << " " << mu.GetDims()[0] << endl;
		return 0;
	}
	if(mu.GetDims()[0] != cov.GetDims()[0]){
		cout << "Error: covariance and mu dimensions do not match. " << cov.GetDims()[0] << " " << mu.GetDims()[0] << endl;
		return 0;
	}
	if(!cov.square()){
		cout << "Error: non-square covariance matrix." << endl;
		return 0;
	}
	//construct x - mu
	Matrix x_mat = Matrix(x.Value());
	Matrix x_min_mu;
	x_min_mu.mult(mu,-1.);
	x_min_mu.add(x_mat);
	
	//transpose x - mu
	Matrix x_min_muT;
	x_min_muT.transpose(x_min_mu);
	Matrix cov_inv;
	cov_inv.invert(cov);
	
	double coeff = 1/(pow(det,0.5)*pow(2*acos(-1),0.5*dim));
	//should only be 1 element matrix
	//muT*cov*mu = 1xd * dxd * dx1
	Matrix mat_expon = Matrix(1,1);
	Matrix cov_mu = Matrix(m_dim,1);
	
	cov_mu.mult(cov_inv,x_min_mu);
	mat_expon.mult(x_min_muT,cov_mu);
	double expon = mat_expon.at(0,0);
	return coeff*exp(-0.5*expon);
	
};



double BasePDFMixture::Dir_lnC(double alpha, int k){
	double alpha_hat = k*alpha;
	
	double norm = lgamma(alpha)*k;
	return lgamma(alpha_hat) - norm;

};

//returns ln(C(alpha))
double BasePDFMixture::Dir_lnC(vector<double> alphas){
	double alpha_hat = 0;
	double norm = 0;
	int k = (int)alphas.size();
	
	for(int i = 0; i < k; i++){
		alpha_hat += alphas[i];
		norm += lgamma(alphas[i]);
	}
	return lgamma(alpha_hat) - norm;

};


//log of B from Wishart distribution
double BasePDFMixture::Wish_lnB(Matrix W, double nu){
	int d = W.GetDims()[0];
	double det = W.det();
	double pi = acos(-1);
	double lgam = 0;
	for(int i = 0; i < d; i++)
		lgam += lgamma( (nu + 1 - i)/2.);
	

	return -(nu/2.)*log(det) - (nu*d/2.)*log(2) - (d*(d-1)/4.)*log(pi) - lgam; 

};


double BasePDFMixture::Wish_H(Matrix W, double nu){
	int d = W.GetDims()[0];
	double E_lam = 0;
	double det = W.det();
	double lnB;
	for(int i = 1; i < d+1; i++)
		E_lam += digamma( (nu + 1 - i)/2. );
	E_lam += d*log(2);
	E_lam += log(det);
	lnB = Wish_lnB(W,nu);

	return -lnB - ((nu - d - 1)/2.)*E_lam + (nu*d)/2.;
};


