#include "Gaussian.hh"
#include "Matrix.hh"

Gaussian::Gaussian(){ 
	m_params["mean"] = {}; m_params["cov"] = {};
	m_dim = 0; m_params["mean"] = {}; m_params["cov"] = {}; }

Gaussian::Gaussian(int d) : BasePDF(d){ 
	m_params["mean"] = {}; m_params["cov"] = {};
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

cout << "gaus Prob" << endl;
cout << "mu dim: " << m_mu.Dim(0) << " cov dim: " << m_cov.Dim(0) << endl;


	double det = m_cov.det();
	double denom = sqrt(det)*pow(2*acos(-1),1/m_dim);

	Matrix x_min_mu = Matrix(x);
	x_min_mu.minus(m_mu);
	
	Matrix x_min_muT;
	x_min_muT.transpose(x_min_mu);
	Matrix cov_inv;
	cov_inv.invert(m_cov);

	x_min_muT.mult(x_min_muT, cov_inv);
	Matrix mat_expon = Matrix(1,1);
	mat_expon.mult(x_min_muT, x_min_mu);

	double expon = -0.5*mat_expon.at(0,0);

	return exp(expon)/denom;	

}

void Gaussian::InitParameters(){
	//init cov to identity and mean to 0
	m_cov.InitIdentity();
	m_mu.InitEmpty();
}



map<string, vector<Matrix>> Gaussian::GetParameters(){
	map<string, vector<Matrix>> ret;
	ret["mu"] = {m_mu};
	ret["cov"] = {m_cov};
	return ret;
}
