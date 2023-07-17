#include "Gaussian.hh"
#include "Matrix.hh"

Gaussian::Gaussian(){ }

Gaussian::Gaussian(Point mu, Matrix cov){
	if(!cov.square()){ cout << "Error: non-square covariance." << endl; return;}
	if(!cov.symmetric()){ cout << "Error: non-symmetric covariance." << endl; return; }
	m_dim = mu.Dim();
	if(mu.Dim() != cov.Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }

	m_cov = cov;
	PointCollection pc;
	pc += mu;
	m_mu.PointsToMat(pc);
	
}

Gaussian::Gaussian(Matrix mu, Matrix cov){
	if(!cov.square()){ cout << "Error: non-square covariance." << endl; return;}
	if(!cov.symmetric()){ cout << "Error: non-symmetric covariance." << endl; return; }
	m_dim = mu.Dim(0);
	if(m_dim != cov.Dim(0)){ cout << "Error: mean and covariance dimensions not compatible." << endl; return; }
	m_mu = mu;
	m_cov = cov;
}

double Gaussian::Prob(const Point& x){
	if(m_dim != x.Dim()){ cout << "Point dimensions incompatible." << endl; return -999;}

	Point pt = Point(m_dim);

	double det = m_cov.det();
	double denom = sqrt(det)*pow(2*acos(-1),1/m_dim);

	Matrix x_min_mu = Matrix(x);
	x_min_mu.minus(m_mu);

	Matrix x_min_muT;
	x_min_muT.transpose(x_min_mu);

	Matrix cov_inv;
	cov_inv.invert(m_cov);

	x_min_muT.mult(x_min_muT, cov_inv);
	x_min_muT.mult(x_min_muT, x_min_mu);

	double expon = -0.5*x_min_muT.at(0,0);

	return exp(expon)/denom;	

}


