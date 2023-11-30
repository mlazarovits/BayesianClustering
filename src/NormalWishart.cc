#include "NormalWishart.hh"
#include "Gaussian.hh"
#include "Wishart.hh"

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

double NormalWishart::Prob(const Matrix& mu, const Matrix& precision){
	Matrix gaus_cov = Matrix(m_dim, m_dim);
	gaus_cov.mult(precision,m_scale);
	gaus_cov.invert(gaus_cov);
	
	Gaussian* gaus = new Gaussian(m_mean, gaus_cov);
	Wishart* wish = new Wishart(m_scalemat, m_dof);
	Point x = Point(m_dim);
	for(int i = 0; i < m_dim; i++) x.SetValue(mu.at(i,0), i);

	return gaus->Prob(x)*wish->Prob(precision);





}



void NormalWishart::InitParameters(unsigned long long seed){ 
	m_mean.InitEmpty();
	m_scalemat.InitIdentity();
	m_scale = 1.;
	m_dof = m_dim;
};

void NormalWishart::UpdateParameters(){ 
	m_mean = m_params["mean"];
	m_scalemat = m_params["scalemat"];
	m_dof = m_params["dof"].at(0,0);
	m_scale = m_params["scale"].at(0,0);
};
