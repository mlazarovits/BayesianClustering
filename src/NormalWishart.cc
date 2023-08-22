#include "NormalWishart.hh"

NormalWishart::NormalWishart(){ 
	m_params["mean"] = Matrix(); m_params["scalemat"] = Matrix();
	m_dim = 0; m_params["scale"] = Matrix(1,1); m_params["dof"] = Matrix(1,1); 
	m_prior = nullptr;
}

NormalWishart::NormalWishart(int d) : BasePDF(d){
	m_params["mean"] = Matrix(d,d); m_params["scalemat"] = Matrix(d,1);
	m_dim = d; m_params["scale"] = Matrix(1,1); m_params["dof"] = Matrix(1,1); 
	m_prior = nullptr;

}

NormalWishart::NormalWishart(Matrix scalemat, Matrix mean, double dof, double scale){
	if(!scalemat.square()){ cout << "Error: non-square covariance." << endl; return;}
	if(!scalemat.symmetric()){ cout << "Error: non-symmetric covariance." << endl; return; }
	if(m_dim == 0) m_dim = mean.Dim(0);
	if(m_dim != scalemat.Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }
	if(scale <= 0){ cout << "Error: scale is less than or equal to zero." << endl; return; }
	if(dof <= m_dim - 1){ cout << "Error: dof is less than or equal to dims." << endl; return; }

	m_mean = mean;
	m_scalemat = scalemat;
	m_dof = dof;
	m_scale = scale;
	m_params["mean"] = m_mean;
	m_params["scalemat"] = m_scalemat;	
	m_params["dof"] = Matrix(dof);
	m_params["scale"] = Matrix(scale);

}

//not defined
double NormalWishart::Prob(const Point& x){ return -999; }

double NormalWishart::Prob(const Matrix& mu, const Matrix& cov){
	double pi = acos(-1);
	double p_gam = multidim_gam(m_dof/2);
	double det = m_scalemat.det();
	Matrix cov_inv = Matrix(m_dim, m_dim);
	Matrix scale_cov_inv = Matrix(m_dim, m_dim);
	cov_inv.invert(cov);
	scale_cov_inv.mult(cov_inv,m_scale);
	double scale_cov_inv_det = scale_cov_inv.det();
	double prefactor = pow(2,m_dof*m_dim/2.)*pow(det,m_dof/2.)*p_gam*pow(2*pi,m_dim/2.)*pow(scale_cov_inv_det,0.5);

	Matrix mat = Matrix(1,1);
	Matrix x = mu;
	x.minus(m_mean);
	Matrix xT = Matrix(1, m_dim);
	xT.transpose(x);

	xT.mult(xT,scale_cov_inv);
	mat.mult(xT,x);

	double tr;
	Matrix cov_scale = Matrix(m_dim, m_dim);
	cov_scale.mult(cov_inv,m_scalemat);
	tr = cov_scale.trace();

	double expon = 0.5*mat.at(0,0) - 0.5*tr;

	return (1./prefactor)*pow(det,(m_dof - m_dim - 1.)/2.)*exp(expon);

}



void NormalWishart::InitParameters(unsigned long long seed){ 
	m_mean.InitEmpty();
	m_scalemat.InitIdentity();
	m_scale = 1.;
	m_dof = m_dim;
};

void NormalWishart::UpdateParameters(){ 
	m_mean = m_params["mean"];
	m_scalemat = m_params["scale_mat"];
	m_dof = m_params["dof"].at(0,0);
	m_scale = m_params["scale"].at(0,0);
};
