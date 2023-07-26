#include "GaussianMixture.hh"
#include "RandomSample.hh"
#include "Gaussian.hh"
#include "KMeansCluster.hh"

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

GaussianMixture::GaussianMixture(int k) : BasePDFMixture(k){
	for(int i = 0; i < k; i++){
		m_model.push_back(new Gaussian());
	}
}

//don't forget to include weights (eventually)
void GaussianMixture::InitParameters(unsigned long long seed){
	cout << "Gaussian Mixture Model with " << m_k << " clusters for " << m_n << " " << m_dim << "-dimensional points." << endl;
	//randomly initialize mean, covariance + mixing coeff.
	RandomSample randy(seed);
	randy.SetRange(0.,1.);
	double coeff_norm = 0;
	for(int k = 0; k < m_k; k++){
		m_model[k]->SetDim(m_dim);

		m_covs.push_back(Matrix(m_dim,m_dim));
		m_covs[k].InitIdentity();
		m_model[k]->SetParameter("cov",m_covs[k]);		

		m_coeffs.push_back(randy.SampleFlat());
		//make sure sum_k m_coeffs[k] = 1
		coeff_norm += m_coeffs[k];
		m_norms.push_back(0.);
	}
	//init means
	KMeansCluster kmc = KMeansCluster(m_data, m_k);
	kmc.Initialize();
	//use number of points that change assignment at E-step to track convergence
	int nit = 0;
	int nchg = 999;
	while(nchg > 0){
		kmc.Estimate();
		kmc.Update();
		nchg = (int)kmc.EvalLogL();
		//check for convergence with number of points that 
		//change assignment
		nit++;
	}
	kmc.GetMeans(m_mus);
	for(int k = 0; k < m_k; k++)
	m_model[k]->SetParameter("mean",m_mus[k]);		

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
//			double g = Gaus(m_data->at(n),m_mus[k],m_covs[k]);
			double g = m_model[k]->Prob(m_data->at(n));	
			gaus.SetEntry(g,n,k);
			norm += m_coeffs[k]*gaus.at(n,k);
		}
		post_norms.push_back(norm);
	}
	//calculate posterior for each data pt n in each cluster k
	for(int	n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			val = m_coeffs[k]*gaus.at(n,k)/post_norms[n];
			m_post.SetEntry(val,n,k);		
		}
	}
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
	//set new means 	
	for(int k = 0; k < m_k; k++){
		m_coeffs[k] = m_norms[k]/m_n;
		//clear and overwrite m_mu[k]
		m_mus[k].clear();
		m_mus[k].InitEmpty();
		//check above does what we want (clear the entries and reset an empty matrix)
		for(int n = 0; n < m_n; n++){
			//add data pt x_n,
			Matrix x = Matrix(m_data->at(n).Value());
			//weighted by posterior value gamma(z_nk),
			x.mult(x,m_post.at(n,k));
			//to new mu for cluster k
			m_mus[k].add(x);
		}
		//normalized by sum of posteriors for cluster k
		m_mus[k].mult(m_mus[k],1./m_norms[k]);
		m_model[k]->SetParameter("mean",m_mus[k]);		
	}




//sigma_k = 1/N_k sum_n(gamma(z_nk)*(x_n - mu_k)*(x_n - mu_k)T) for mu_k = mu^new_k
	for(int k = 0; k < m_k; k++){
		//create (x_n - mu)*(x_n - mu)T matrices for each data pt
		Matrix new_cov = Matrix(m_dim,m_dim);
		for(int n = 0; n < m_n; n++){
			Matrix cov_k = Matrix(m_dim, m_dim);

			//construct x - mu
			Matrix x_min_mu = Matrix(m_data->at(n).Value());
			x_min_mu.minus(m_mus[k]);

	
			//transpose x - mu
			Matrix x_min_muT;
			x_min_muT.transpose(x_min_mu);
			//(x_n - mu_k)*(x_n - mu_k)T
			cov_k.mult(x_min_mu,x_min_muT);
			//weighting by posterior gamma(z_nk)
			cov_k.mult(cov_k,m_post.at(n,k));
			//sum over n
			new_cov.add(cov_k);
		}	
		//normalize by N_k
		new_cov.mult(new_cov,1./m_norms[k]);
		//overwrites m_covs[k]
		m_covs[k] = new_cov;
		m_model[k]->SetParameter("cov",m_covs[k]);		
		
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
			sum_k += m_coeffs[k]*m_model[k]->Prob(m_data->at(n));
			
		}
		L += log(sum_k);
	}
	return L;
}




map<string, vector<Matrix>> GaussianMixture::GetParameters(){ 
	map<string, vector<Matrix>> params;
	for(int k = 0; k < m_k; k++){
		params["means"].push_back(m_model[k]->GetParameter("mean"));
		params["covs"].push_back(m_model[k]->GetParameter("cov"));
		//params["means"].push_back(m_mus[k]);
		//params["covs"].push_back(m_covs[k]);
	}
		params["pis"] = {Matrix(m_weights)};
	return params;
}; 
/*
//fill vectors with params for each cluster
void GaussianMixture::GetParameters(vector<Matrix>& mus, vector<Matrix>& covs){
	mus.clear();
	covs.clear();
	mus = m_mus;
	covs = m_covs;
}
*/
