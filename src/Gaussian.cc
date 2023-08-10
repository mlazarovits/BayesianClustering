#include "Gaussian.hh"
#include "Matrix.hh"
#include "MultivarT.hh"


Gaussian::Gaussian(){ 
	m_dim = 0; m_params["mean"] = Matrix(); m_params["cov"] = Matrix(); 
	m_prior = nullptr;
}

Gaussian::Gaussian(int d) : BasePDF(d){ 
	m_params["mean"] = Matrix(d,1); m_params["cov"] = Matrix(d,d);
	m_prior = nullptr;
}

Gaussian::Gaussian(Point mu, Matrix cov){
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

double Gaussian::Prob(const Point& x){
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
	Point pt = Point(m_dim);
	for(int i = 0; i < m_dim; i++) pt.SetValue(p.at(i,0),i);

	prob->SetPrefactor(prefactor);
	return prob->Prob(pt);


}

void Gaussian::InitParameters(unsigned long long seed){
	//init cov to identity and mean to 0
	m_cov.InitIdentity();
	m_mu.InitEmpty();
}


double Gaussian::ConjugateEvidence(const Point& x){
	PointCollection pc;
	pc += x;
	return ConjugateEvidence(pc);

/*
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
cout << "Gaussian::ConjugateEvidence - dof: " << dof - m_dim + 1 << " scale matrix ";
scalemat.Print(); cout << "mean "; mean.Print();
	//make distribution parameters
	MultivarT* dist = new MultivarT(mean,scalemat,dof - m_dim + 1); 

	return dist->Prob(x);
*/
}

double Gaussian::ConjugateEvidence(const PointCollection& x){
//cout << "Gaussian::ConjugateEvidence" << endl;
	//assuming conjugate prior - for multidim gaussian with unknown mean + variance (precision) is normal inverse wishart (normal wishart)
	if(m_prior == nullptr){
		return -999;
	}
	int n = x.GetNPoints();
	Matrix xbar = Matrix(m_dim, 1);
	xbar.mean(x);
//cout << "xbar" << endl;
//xbar.Print();

	//prior should be normal inverse wishart
	Matrix m0 = m_prior->GetParameter("mean");
	Matrix S0 = m_prior->GetParameter("scalemat");
	//check that these parameters exist
	double nu0 = m_prior->GetParameter("dof").at(0,0);
	double beta0 = m_prior->GetParameter("scale").at(0,0);
	if(m0.empty() || S0.empty()){ cout << "Error: incorrect prior specified. Needs to be Normal inverse-Wishart." << endl; return -1.; }

//cout << "S0" << endl;
//S0.Print();
//cout << "\n" << endl;

	//from Bishop
	//using posterior parameters in eqns (4.209) - (4.214)
	Matrix Sn = Matrix(m_dim, m_dim);
	Matrix mn = Matrix(m_dim, 1);
	double beta_n = beta0 + n;
	double nu_n = nu0 + n;
//cout << "beta0 " << beta0 << " beta_n: " << beta_n << " nu0: " << nu0 << " nu_n: " << nu_n << endl;
	//posterior mean of normal part of NIW
	mn.mult(m0,beta0/beta_n);
//cout << "m0 scaled by beta0/betan" << endl;
//mn.Print();
//cout << "\n" << endl;
	Matrix xbar_scaled;
	xbar_scaled.mult(xbar, n/beta_n);
//cout << "xbar scaled by n/betan" << endl;
//xbar_scaled.Print();
//cout << "\n" << endl;
	mn.add(mn,xbar_scaled);
//cout << "mn" << endl;
//mn.Print();
//cout << "\n" << endl;

	//posterior scale of wishart part of NIW
	Matrix S = Matrix(m_dim, m_dim);
	Matrix xmat, xmatT, xmat_xmatT;
	for(int i = 0; i < n; i++){
		xmat.PointToMat(x.at(i));
		xmat.minus(xbar);
		xmatT.transpose(xmat);
		xmat_xmatT.mult(xmat,xmatT);
		S.add(S,xmat_xmatT);
	}
//cout << "S" << endl;
//S.Print();
//cout << "\n" << endl;	

	Sn.add(S0,S);

	Matrix m0T;
	Matrix m0_m0T = Matrix(m_dim, m_dim);
	m0T.transpose(m0);
	m0_m0T.mult(m0,m0T);
//cout << "m0m0T" << endl;
//m0_m0T.Print();
//cout << "\n" << endl;
	m0_m0T.mult(m0_m0T,beta0);
//cout << "m0m0T*beta0" << endl;
//m0_m0T.Print();
//cout << "\n" << endl;
	Sn.add(Sn,m0_m0T);

	Matrix mnT = Matrix(1, m_dim);
	Matrix mn_mnT = Matrix(m_dim, m_dim);
	mnT.transpose(mn);
	mn_mnT.mult(mn,mnT);
//cout << "mnmnT" << endl;
//mn_mnT.Print();
//cout << "\n" << endl;
	mn_mnT.mult(mn_mnT,beta_n);
//cout << "mnmnT*betan" << endl;
//mn_mnT.Print();
//cout << "\n" << endl;
	Sn.add(Sn,mn_mnT);
//cout << "Sn" << endl;
//Sn.Print();
//cout << "\n" << endl;
	//create evidence from Bishop (5.29)
	double pi = acos(-1);
	double ret = pow(pi,-n*m_dim/2.);
	double det0 = S0.det();
	double detn = Sn.det();
	//cout << "detn: " << detn << " nu_n/2: " << nu_n/2. << endl;
	//cout << "ret 1: " << ret << endl;
	ret *= pow(beta0/beta_n,m_dim/2.);
	//cout << "ret 2: " << ret << endl;
	ret *= pow(det0,nu0/2)/pow(detn,nu_n/2);
	//cout << "pow(det0,nu0/2) " << pow(det0,nu0/2) << endl;
	//cout << "pow(detn,nu_n/2) " << pow(detn,nu_n/2) << endl;
	//cout << "ret 3: " << ret << endl;

	ret *= multidim_gam(nu_n/2);
//cout << "gam_d(nu_n/2) " << multidim_gam(nu_n/2) << endl;
//	cout << "ret 5: " << ret << endl;
//cout << "gam_d(nu_0/2) " << multidim_gam(nu0/2) << endl;
	ret /= multidim_gam(nu0/2);	

//cout << "Gaussian::ConjugateEvidence - end: " << ret << endl;
	return ret;	
}


Gaussian* Gaussian::mult(Gaussian* p1){ 
	//check that p1 is a Gaussian by looking at parameters
	Matrix mu2 = p1->GetParameter("mean");
	Matrix mu1 = m_mu;

	Matrix cov2 = p1->GetParameter("cov");
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

	Point x = mu1.MatToPoints().at(0);
	double coeff = gaus_coeff->Prob(x);	

	ret->SetParameter("mean",new_mu);
	ret->SetParameter("cov",new_cov);
	ret->SetPrefactor(coeff);

	return ret;	
};
