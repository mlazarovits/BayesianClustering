#include "GaussianMixture.hh"
#include "Gaussian.hh"
#include "KMeansCluster.hh"


GaussianMixture::GaussianMixture(int k) : BasePDFMixture(k){
	for(int i = 0; i < k; k++){
		m_model.push_back(new Gaussian());
	}

}

void GaussianMixture::InitParameters(unsigned long long seed){
	m_n = m_data->GetNPoints();
	//run k-means clustering on data to convergence to get means
	//get means for each cluster 
	KMeansCluster kmc = KMeansCluster(m_data);
	kmc.SetNClusters(m_k);
	kmc.Initialize();

	//use number of points that change assignment at E-step to track convergence
	int nit = 0;
	int nchg = 999;
	while(nchg > 0){
		kmc.Estimate();
		//check for convergence with number of points that 
		//change assignment
		nchg = (int)kmc.EvalLogL();
		
		kmc.Update();
		nit++;
	}
	
	kmc.GetMeans(m_mus);

	//init covs to identity
	//randomly initialize mixing coeffs
	RandomSample randy(seed);
	randy.SetRange(0.,1.);
	double coeff_norm = 0;
	for(int k = 0; k < m_k; k++){
		m_covs.push_back(Matrix(m_dim,m_dim));
		m_covs[k].InitIdentity();
		m_weights.push_back(randy.SampleFlat());
		coeff_norm += m_weights[k];
	}
	//make sure sum_k m_weights[k] = 1
	for(int k = 0; k < m_k; k++) m_weights[k] /= coeff_norm;


}


//for normal EM algorithm - no priors
void GaussianMixture::CalculatePosterior(Matrix& post){
	//each kth vector is a cluster of n data points
	double val, g;
	Matrix mat = Matrix(m_n,m_k);
	//these are NOT N_k, they are the normalizations for the posterior gamma(z_nk)
	//summed over k for each n (n entries)
	vector<double> post_norms;
	Gaussian gaus;
	//calculate norms for each data pt (summed over k)
	for(int	n = 0; n < m_n; n++){
		double norm = 0.;
		for(int k = 0; k < m_k; k++){
		//fill posterior with values according to Bishop eq. (9.23)	
		//gamma(z_nk) = pi_k*N(x_n | mu_k, sigma_k)/sum_{j=1}^K(pi_j*N(x_n | mu_j, sigma_j))
			gaus = Gaussian(m_mus[k],m_covs[k]);
			g = gaus.Prob(m_data->at(n));
			mat.SetEntry(g,n,k);
			norm += m_weights[k]*mat.at(n,k);
		}
		post_norms.push_back(norm);
	}
	//calculate posterior for each data pt n in each cluster k
	for(int	n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			val = m_weights[k]*mat.at(n,k)/post_norms[n];
			post.SetEntry(val,n,k);		
		}
	}

}

//for variational EM algorithm, assuming conjugate priors on mu, cov, and mixing coeffs
void GaussianMixture::CalculateVariationalPosterior(Matrix& post){


}

//for normal EM algorithm - no priors, updates mu, cov, and mixing coeffs
void GaussianMixture::UpdateParameters(const Matrix& post){
	//re-calculate normalization (overwrites)
	m_norms.clear();
	//this is for N_k - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms.push_back(0.);
		for(int n = 0; n < m_n; n++){
			m_norms[k] += post.at(n,k);
		}
	}
//cout << "new means" << endl;
	//set new means + mixing coeffs
	m_weights.clear();	
	for(int k = 0; k < m_k; k++){
		//overwrite and recalculate mixing coefficients
		m_weights.push_back(m_norms[k]/m_n);
		//clear and overwrite m_mu[k]
		m_mus[k].clear();
		m_mus[k].InitEmpty();
		for(int n = 0; n < m_n; n++){
			//add data pt x_n,
			Matrix x = Matrix(m_data->at(n).Value());
			//weighted by posterior value gamma(z_nk),
			x.mult(x,post.at(n,k));
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
			Matrix x_mat = Matrix(m_data->at(n).Value());
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
			cov_k.mult(cov_k,post.at(n,k));
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


//for variational EM algorithm, assuming conjugate priors on mu, cov (normalWishart), and mixing coeffs (Dirichlet)
//updates alpha (Dirichlet), W, nu (Wishart), beta, m (normal)
void GaussianMixture::UpdateVariationalParameters(const Matrix& post){



}

map<string, vector<Matrix>> GaussianMixture::GetParameters(){
	map<string, vector<Matrix>> params;

	params["mu"] = m_mus;
	params["cov"] = m_covs;
	params["pi"] = {Matrix(m_weights)};
	return params;
}

map<string, vector<Matrix>> GaussianMixture::GetPriorParameters(){
	map<string, vector<Matrix>> params;

	params["model:mu"] = m_mus;
	params["model:cov"] = m_covs;
	params["model:pi"] = {Matrix(m_weights)};
	params["dirPrior:alphas"] = {Matrix(m_alphas)};
	params["nwPrior:Ws"] = m_Ws;
	params["nwPrior:nus"] = {Matrix(m_nus)};
	params["nwPrior:ms"] = m_means;
	params["nwPrior:betas"] = {Matrix(m_betas)};	
	return params;
}



double GaussianMixture::EvalLogL(){
	double L;
	double sum_k;
	Gaussian gaus;
	for(int n = 0; n < m_n; n++){
		sum_k = 0;
		for(int k = 0; k < m_k; k++){
			//matrix mu = params[k][0];
			//matrix cov = params[k][1];
			gaus = Gaussian(m_mus[k], m_covs[k]);
			sum_k += m_weights[k]*gaus.Prob(m_data->at(n));
		}
		L += log(sum_k);
	}
	return L;

}
