#include "Gaussian.hh"
#include "Matrix.hh"
#include "MultivarT.hh"


Gaussian::Gaussian(){ 
	m_params["mean"] = {}; m_params["cov"] = {};
	m_dim = 0; m_params["mean"] = {}; m_params["cov"] = {}; 
	m_prior = nullptr;
}

Gaussian::Gaussian(int d) : BasePDF(d){ 
	m_params["mean"] = {}; m_params["cov"] = {};
	m_prior = nullptr;
}

Gaussian::Gaussian(Point mu, Matrix cov){
	m_params["mean"] = {}; m_params["cov"] = {};
	if(!cov.square()){ cout << "Error: non-square covariance." << endl; return;}
	if(!cov.symmetric()){ cout << "Error: non-symmetric covariance." << endl; return; }
	m_dim = mu.Dim();
	if(mu.Dim() != cov.Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }

	if(m_dim == 0){ m_dim = mu.Dim(); }

	m_cov = cov;
	PointCollection pc;
	pc += mu;
	m_mu.PointsToMat(pc);


	m_params["mean"] = m_mu;
	m_params["cov"] = m_cov;	
}

Gaussian::Gaussian(Matrix mu, Matrix cov){
	if(!cov.square()){ cout << "Error: non-square covariance." << endl; return;}
	if(!cov.symmetric()){ cout << "Error: non-symmetric covariance." << endl; return; }
	if(m_dim == 0) m_dim = mu.Dim(0);
	if(m_dim != cov.Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }
	m_mu = mu;
	m_cov = cov;
	m_params["mean"] = m_mu;
	m_params["cov"] = m_cov;	
}

double Gaussian::Prob(const Point& x){
	if(m_dim != x.Dim()){ cout << "Point dimensions " << x.Dim() << " incompatible. " << m_dim << endl; return -999;}

	double det = m_cov.det();
	Matrix x_min_mu = Matrix(x);
	x_min_mu.minus(m_mu);
	//transpose x - mu
	Matrix x_min_muT;
	x_min_muT.transpose(x_min_mu);
	Matrix cov_inv;
	cov_inv.invert(m_cov);

	double coeff = 1./(pow(det,0.5)*pow(2*acos(-1),0.5*m_dim));
	//should only be 1 element matrix
	//muT*cov*mu = 1xd * dxd * dx1
	Matrix mat_expon = Matrix(1,1);
	Matrix cov_mu = Matrix(m_dim,1);

	cov_mu.mult(cov_inv,x_min_mu);

	mat_expon.mult(x_min_muT,cov_mu);

	double expon = mat_expon.at(0,0);
	return coeff*exp(-0.5*expon);


}

void Gaussian::InitParameters(){
	//init cov to identity and mean to 0
	m_cov.InitIdentity();
	m_mu.InitEmpty();
}


double Gaussian::ConjugateEvidence(const Point& x){
//assuming conjugate prior - for multidim gaussian with unknown mean + variance (precision) is normal inverse wishart (normal wishart)
	if(m_prior == nullptr){
		return -999;
	}

	//prior should be normal inverse wishart
	Matrix mean = m_prior->GetParameter("mean");
	Matrix scalemat = m_prior->GetParameter("scalemat");
	//check that these parameters exist
	if(mean.empty() || scalemat.empty()){ cout << "Error: incorrect prior specified. Needs to be Normal inverse-Wishart." << endl; return -1.; }


	double dof = m_prior->GetParameter("dof").at(0,0);
	double scale = m_prior->GetParameter("scale").at(0,0);

	scalemat.mult(scalemat,(scale + 1)/scale*(dof - m_dim + 1));

	//make distribution parameters
	MultivarT* dist = new MultivarT(mean,scalemat,dof - m_dim + 1); 

	return dist->Prob(x);
}
