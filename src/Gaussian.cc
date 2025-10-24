#include "Gaussian.hh"
#include "Matrix.hh"


Gaussian::Gaussian(){ 
	m_dim = 0; m_params["mean"] = Matrix(); m_params["cov"] = Matrix(); 
	m_prior = nullptr;
}

Gaussian::Gaussian(int d) : BasePDF(d){ 
	m_params["mean"] = Matrix(d,1); m_params["cov"] = Matrix(d,d);
	m_prior = nullptr;
}

Gaussian::Gaussian(BayesPoint mu, Matrix cov){
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
	m_prefactor = 1./(pow(m_cov.det(),0.5)*pow(2*acos(-1),0.5*m_dim));
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
	m_prefactor = 1./(pow(m_cov.det(),0.5)*pow(2*acos(-1),0.5*m_dim));
}

double Gaussian::Prob(const BayesPoint& x){
	if(m_dim != x.Dim()){ cout << "Point dimensions " << x.Dim() << " incompatible. " << m_dim << endl; return -999;}

	Matrix x_min_mu = Matrix(x);
	x_min_mu.minus(m_mu);
	//transpose x - mu
	Matrix x_min_muT;
	x_min_muT.transpose(x_min_mu);
	Matrix cov_inv;
	cov_inv.invert(m_cov);

	//should only be 1 element matrix
	//muT*cov*mu = 1xd * dxd * dx1
	Matrix mat_expon = Matrix(1,1);
	Matrix cov_mu = Matrix(m_dim,1);

	cov_mu.mult(cov_inv,x_min_mu);

	mat_expon.mult(x_min_muT,cov_mu);

	double expon = mat_expon.at(0,0);
	return m_prefactor*exp(-0.5*expon);


}

double Gaussian::Prob(const PointCollection& x){
	if(m_dim != x.Dim()){ cout << "Point dimensions " << x.Dim() << " incompatible. " << m_dim << endl; return -999;}
	int n = x.GetNPoints();
	double prefactor = pow(m_prefactor, n);
	Matrix m = Matrix(m_dim,1);
	m.mean(x);
	Matrix scatter = Matrix(m_dim, m_dim);
	scatter.scatter(x);
	//from Bishop 4.196
	Matrix cov_scatter = Matrix(m_dim, m_dim);
	Matrix cov_inv = Matrix(m_dim, m_dim);
	cov_inv.invert(m_cov);
	cov_scatter.mult(cov_inv, scatter);
	
	double trace = cov_scatter.trace();
	prefactor *= exp(-n/2*trace);
	
	m.mult(m,sqrt(n));

	//N/2(mu - mean)*covinv*(mu - mean);
	Gaussian* prob = new Gaussian(m,m_cov);
	Matrix p = Matrix(m_dim, 1);
	p.mult(m_mu,sqrt(n));
	BayesPoint pt = BayesPoint(m_dim);
	for(int i = 0; i < m_dim; i++) pt.SetValue(p.at(i,0),i);

	prob->SetPrefactor(prefactor);
	return prob->Prob(pt);


}

void Gaussian::InitParameters(unsigned long long seed){
	//init cov to identity and mean to 0
	m_cov.InitIdentity();
	m_mu.InitEmpty();
}

Gaussian* Gaussian::mult(Gaussian* p1){ 
	//check that p1 is a Gaussian by looking at parameters
	Matrix mu2; p1->GetParameter("mean", mu2);
	Matrix mu1 = m_mu;

	Matrix cov2; p1->GetParameter("cov", cov2);
	cov2.invert(cov2);
	Matrix cov1 = Matrix(m_dim, m_dim);
	cov1.invert(m_cov);
	

	Gaussian* ret = new Gaussian(m_dim);

	//form from: https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf (Ch. 8.1.8)
	Matrix new_cov = Matrix(m_dim, m_dim);
	new_cov.add(cov1,cov2);
	
	Matrix new_mu1 = Matrix(m_dim, 1);
	Matrix new_mu2 = Matrix(m_dim, 1);
	Matrix new_mu = Matrix(m_dim, 1);
	
	new_mu1.mult(cov1,mu1);
	new_mu2.mult(cov2,mu2);

	new_mu.add(new_mu1, new_mu2);
	new_mu.mult(new_cov,new_mu);

	//use added covariances in this distribution
	Gaussian* gaus_coeff = new Gaussian(mu2, new_cov);
	//then invert added covariances for overall distribution
	new_cov.invert(new_cov); 

	BayesPoint x; mu1.MatToPoint(x);
	double coeff = gaus_coeff->Prob(x);	

	ret->SetParameter("mean",new_mu);
	ret->SetParameter("cov",new_cov);
	ret->SetPrefactor(coeff);

	return ret;	
};
