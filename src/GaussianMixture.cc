#include "GaussianMixture.hh"
#include "RandomSample.hh"

#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

//k = # of clusters (cols)
//n = # of data pts (rows)
GaussianMixture::GaussianMixture(){ 
	m_k = 0;
	m_n = 0;
}

GaussianMixture::GaussianMixture(int k){
	m_k = k;
	m_n = 0;
}

//don't forget to include weights (eventually)
void GaussianMixture::Initialize(unsigned long long seed){
	cout << "Gaussian Mixture Model with " << m_k << " clusters for " << m_n << " " << m_dim << "-dimensional points." << endl;
	//randomly initialize mean, covariance + mixing coeff.
	RandomSample randy(seed);
	//TODO: should use data to set the range on the possible parameter values
	randy.SetRange(0.,1.);
	double mu_lower = m_x.min()-0.1;
	double mu_upper = m_x.max()+0.1;
	double coeff_norm = 0;
	for(int k = 0; k < m_k; k++){
		m_mus.push_back(Matrix(m_dim,1));
		m_mus[k].InitRandom(mu_lower,mu_upper,seed+k);
	//	cout << "k: " << k << endl;
	//	cout << "mu" << endl;
	//	m_mus[k].Print();
		m_covs.push_back(Matrix(m_dim,m_dim));
		m_covs[k].InitIdentity();
//cout << "cov" << endl;
//m_covs[k].Print();
		m_coeffs.push_back(randy.SampleFlat());
		//make sure sum_k m_coeffs[k] = 1
		coeff_norm += m_coeffs[k];
		m_norms.push_back(0.);
	}
	//make sure sum_k m_coeffs[k] = 1
	for(int k = 0; k < m_k; k++) m_coeffs[k] /= coeff_norm;
	m_post.SetDims(m_n, m_k);
}


//E-step
void GaussianMixture::CalculatePosterior(){
	//each kth vector is a cluster of n data points
	double val;
	Matrix gaus = Matrix(m_n,m_k);
	//these are NOT N_k, they are the normalizations for the posterior gamma(z_nk)
	//summed over k for each n (n entries)
	vector<double> post_norms;
	//calculate norms for each data pt (summed over k)
	for(int	n = 0; n < m_n; n++){
		double norm = 0.;
		for(int k = 0; k < m_k; k++){
		//fill posterior with values according to Bishop eq. (9.23)	
		//gamma(z_nk) = pi_k*N(x_n | mu_k, sigma_k)/sum_{j=1}^K(pi_j*N(x_n | mu_j, sigma_j))
//			cout << "n: " << n << " k: " << k << endl;
			double g = Gaus(m_x.at(n),m_mus[k],m_covs[k]);
//		cout << "n: " << n << " k: " << k << " coef: " << m_coeffs[k] << " gaus: " << g << endl;
			gaus.SetEntry(g,n,k);
			norm += m_coeffs[k]*gaus.at(n,k);
		}
		post_norms.push_back(norm);
	}
	//calculate posterior for each data pt n in each cluster k
	for(int	n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			val = m_coeffs[k]*gaus.at(n,k)/post_norms[n];
		//	cout << "k: " << k << " n: " << n << " coeff: " << m_coeffs[k] << " gaus: " << gaus.at(n,k) << " norm: " << norms[n] << " val: " << val << endl;
			m_post.SetEntry(val,n,k);		
		}
	}
	//m_post.Print();
}



//M-step
//equations were derived from maximizing the posterior calculated in the E-step
void GaussianMixture::UpdateParameters(){
	//re-calculate normalization (overwrites)
	m_norms.clear();
	//this is for N_k - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms.push_back(0.);
		for(int n = 0; n < m_n; n++){
			m_norms[k] += m_post.at(n,k);
		}
	}
//cout << "new means" << endl;
	//set new means 	
	for(int k = 0; k < m_k; k++){
		m_coeffs[k] = m_norms[k]/m_n;
		//clear and overwrite m_mu[k]
		m_mus[k].clear();
		m_mus[k].InitEmpty();
		//check above does what we want (clear the entries and reset an empty matrix)
		for(int n = 0; n < m_n; n++){
			//add data pt x_n,
			Matrix x = Matrix(m_x.at(n).Value());
			//weighted by posterior value gamma(z_nk),
			x.mult(x,m_post.at(n,k));
			//to new mu for cluster k
			m_mus[k].add(x);
		}
		//normalized by sum of posteriors for cluster k
		m_mus[k].mult(m_mus[k],1./m_norms[k]);
		//cout << "k: " << k << endl;
		//m_mus[k].Print();
	}



//cout << "new covs" << endl;

//sigma_k = 1/N_k sum_n(gamma(z_nk)*(x_n - mu_k)*(x_n - mu_k)T) for mu_k = mu^new_k
	for(int k = 0; k < m_k; k++){
		//create (x_n - mu)*(x_n - mu)T matrices for each data pt
		Matrix new_cov = Matrix(m_dim,m_dim);
		for(int n = 0; n < m_n; n++){
			Matrix cov_k = Matrix(m_dim, m_dim);

			//construct x - mu
			Matrix x_mat = Matrix(m_x.at(n).Value());
			Matrix x_min_mu;
			x_min_mu.mult(m_mus[k],-1.);
			x_min_mu.add(x_mat);
	//	cout << "x - mu" << endl;
//		x_min_mu.Print();

	
			//transpose x - mu
			Matrix x_min_muT;
			x_min_muT.transpose(x_min_mu);
	//	cout << "(x - mu)T" << endl;
	//	x_min_muT.Print();
			//(x_n - mu_k)*(x_n - mu_k)T
			cov_k.mult(x_min_mu,x_min_muT);
	//	cout << "(x_n - mu_k)*(x_n - mu_k)T" << endl;
	//	cov_k.Print();
			//weighting by posterior gamma(z_nk)
			cov_k.mult(cov_k,m_post.at(n,k));
	//	cout << "gam(z_nk)*(x_n - mu_k)*(x_n - mu_k)T" << endl;
	//	cov_k.Print();
	//	cout << "gam(z_nk) = " << m_post.at(n,k) << endl;	
			//sum over n
			new_cov.add(cov_k);
	//		cout << "running sum" << endl;
	//		new_cov.Print();
	//	cout << "\n" << endl;
		}	
		//normalize by N_k
		new_cov.mult(new_cov,1./m_norms[k]);
		//overwrites m_covs[k]
		m_covs[k] = new_cov;
//	cout << "k: " << k << endl;	
//	m_covs[k].Print();
//	cout << "new" << endl;
//	new_cov.Print();
		
	}

}


//want to maximize
//this is used to test for algorithm convergence
//ln[p(X | mu, sigma, pi) = sum_n( ln[sum_k(pi_k*Gaus(x_n | mu_k, sigma_k))] )
double GaussianMixture::EvalLogL(){
	double L;
	double sum_k;
	for(int n = 0; n < m_n; n++){
		sum_k = 0;
		for(int k = 0; k < m_k; k++){
			sum_k += m_coeffs[k]*Gaus(m_x.at(n),m_mus[k],m_covs[k]);
			
		}
		L += log(sum_k);
	}
	return L;
}






//consider doing this for whole class or moving Gaus function to more general class
//gaussian for one data point
double GaussianMixture::Gaus(const Point& x, const Matrix& mu, const Matrix& cov){
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
//	cout << "	x - mu:" << endl;
//	x_min_mu.Print();

	//transpose x - mu
	Matrix x_min_muT;
	x_min_muT.transpose(x_min_mu);
	Matrix cov_inv;
//	cout << "	cov:" << endl;
//	cov.Print();
	cov_inv.invert(cov);
//	cout << "	cov_inv:" << endl;
//	cov_inv.Print();

	double coeff = 1/(pow(det,0.5)*pow(2*acos(-1),0.5*dim));
	//should only be 1 element matrix
	//muT*cov*mu = 1xd * dxd * dx1
	Matrix mat_expon = Matrix(1,1);
	Matrix cov_mu = Matrix(m_dim,1);

	cov_mu.mult(cov_inv,x_min_mu);
//	cout << "	cov_inv*(x-mu):" << endl;
//	cov_mu.Print();

	mat_expon.mult(x_min_muT,cov_mu);
//	cout << "	(x-mu)T*cov_inv*(x-mu):" << endl;
//	mat_expon.Print();


	double expon = mat_expon.at(0,0);
//	cout << "	x:";
//	x.Print();
//	cout << "	expon: " << expon << endl;
	return coeff*exp(-0.5*expon);

}



//fill vectors with params for each cluster
void GaussianMixture::GetParameters(vector<Matrix>& mus, vector<Matrix>& covs){
	mus.clear();
	covs.clear();
	mus = m_mus;
	covs = m_covs;
}

