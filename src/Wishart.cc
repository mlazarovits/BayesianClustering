#include <boost/math/special_functions/digamma.hpp>
#include "Wishart.hh"
#include "RandomSample.hh"

using boost::math::digamma;

Wishart::Wishart(const Matrix& W, double nu){
	if(!W.square()){ cout << "Error: non-square scale matrix." << endl; return;}
	m_dim = W.Dim(0);

	if(nu < m_dim - 1){ cout << "Error: degrees of freedom " << nu << " less than dimensions: " << m_dim << endl; return; } 

	m_W = W;
	m_nu = nu;

}

//Wishart is only for matrices
double Wishart::Prob(const Point& x){
	return -1;
}

double Wishart::Prob(const Matrix& x){
	if(!x.square()){ cout << "Error: non-square random variable matrix." << endl; return -1;}

	double det_x = x.det();
	double det_w = m_W.det();
	double pi = acos(-1);

	double gam_p = multidim_gam(m_nu/2.);

	Matrix W_inv_mult_x;
	W_inv_mult_x.invert(m_W);
	W_inv_mult_x.mult(W_inv_mult_x,x);
	
	double tr = W_inv_mult_x.trace();

	return (1./(pow(2,m_dim*m_nu/2.)*pow(det_w,m_nu/2.)*gam_p))*pow(det_x, (m_nu - m_dim - 1)/2.)*exp(-0.5*tr);

}


double Wishart::lnB(){
	double det = m_W.det();
	double pi = acos(-1);
	double lgam = 0;
	for(int i = 0; i < m_dim; i++)
		lgam += lgamma( (m_nu + 1 - i)/2.);
	

	return -(m_nu/2.)*log(det) - (m_nu*m_dim/2.)*log(2) - (m_dim*(m_dim-1)/4.)*log(pi) - lgam; 
}

double Wishart::H(){
	double E_lam = 0;
	double det = m_W.det();
	for(int i = 1; i < m_dim+1; i++)
		E_lam += digamma( (m_nu + 1 - i)/2. );
	E_lam += m_dim*log(2);
	E_lam += log(det);

	return -lnB() - ((m_nu - m_dim - 1)/2.)*E_lam + (m_nu*m_dim)/2.;
}



void Wishart::InitParameters(unsigned long long seed){
	//init W to identity
	m_W.InitIdentity();
	//nu > d - 1 (degrees of freedom)
	RandomSample rs;
	rs.SetRange(m_dim - 1, m_dim+2);
	m_nu = rs.SampleFlat();

}

