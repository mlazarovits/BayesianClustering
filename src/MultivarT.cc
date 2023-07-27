#include <boost/math/special_functions/gamma.hpp>
#include "MultivarT.hh"
using boost::math::tgamma;

MultivarT::MultivarT(int d) : BasePDF(d){
	m_params["mean"] = {};
	m_params["scale_mat"] = {};
	m_params["dof"] = {};
}

MultivarT::MultivarT(Matrix mean, Matrix scale, double dof){
	if(!scale.square()){ cout << "Error: non-square covariance." << endl; return;}
	if(m_dim == 0) m_dim = mean.Dim(0);
	if(m_dim != scale.Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }
	m_mean = mean;
	m_scale = scale;
	m_dof = dof;
	m_params["mean"] = m_mean;
	m_params["scale"] = m_scale;
	vector<double> dofs; dofs.push_back(dof);	
	m_params["dof"] = Matrix(dofs);

}


void MultivarT::InitParameters(){
	m_mean.InitEmpty();
	m_scale.InitIdentity();
	m_dof = 1.1;
}


double MultivarT::Prob(const Point& x){
	//undefined case
	if(m_dof <= 1) return -999;

	double pi = acos(-1);
	double det = m_scale.det();
	double prefactor = (tgamma(m_dof + m_dim)/2.)/( tgamma(m_dof/2.)*pow(m_dof,m_dim/2.)*pow(pi,m_dim/2.)*pow(det,0.5) );

	Matrix mat = Matrix(x);
	mat.minus(m_mean);
	Matrix matT = Matrix(1,m_dim);
	matT.transpose(mat); 
	Matrix scale_inv = Matrix(m_dim, m_dim);
	scale_inv.invert(m_scale);
	
	Matrix factor = Matrix(1,1);
	matT.mult(matT,scale_inv);
	factor.mult(matT,mat);

	return prefactor*(1 + (1./m_dof)*factor.at(0,0));
}

