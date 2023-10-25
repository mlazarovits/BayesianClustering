#include "GaussianMixture.hh"
#include "RandomSample.hh"
#include "Gaussian.hh"
#include "KMeansCluster.hh"
#include "Dirichlet.hh"
#include "Wishart.hh"
#include <boost/math/special_functions/digamma.hpp>

#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;
using boost::math::digamma;

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

void GaussianMixture::InitParameters(unsigned long long seed){
//cout << "GaussianMixture::InitParameters alpha0 = " << m_alpha0 << endl;
//	cout << "Gaussian Mixture Model with " << m_k << " clusters for " << m_n << " " << m_dim << "-dimensional points." << endl;
	//cout << "InitParameters" << endl;
	//randomly initialize mean, covariance + mixing coeff.
	RandomSample randy(seed);
	randy.SetRange(0.,1.);
	double coeff_norm = 0;
	double norm_norm = 0;
	for(int k = 0; k < m_k; k++){
		//seed N_k to even posterior values (even probabilities for all clusters -> n*(1/kmax)) - make sure 0params are distinct for convergence
		m_norms[k] = double(m_n)/double(m_k);
		double w = 0;
		for(int n = 0; n < m_n; n++){ w += m_data->at(n).w(); }
		m_norms[k] *= w/double(m_n); 

		m_norms_unwt[k] = double(m_n)/double(m_k);
		m_model[k]->SetDim(m_dim);


		m_coeffs[k] = randy.SampleFlat();
		//make sure sum_k m_coeffs[k] = 1
		coeff_norm += m_coeffs[k];
		norm_norm += m_norms[k];
	}
	//make sure sum_k m_coeffs[k] = 1
	for(int k = 0; k < m_k; k++) m_coeffs[k] /= coeff_norm;
	
	//for(int k = 0; k < m_k; k++){ m_norms[k] *= double(m_n)/norm_norm;}
	//for(int k = 0; k < m_k; k++) cout << "k: " << k << " Nk: " << m_norms[k] << endl;
	//init means
	KMeansCluster kmc = KMeansCluster(m_data, m_k);
	kmc.Initialize(seed);
	//use number of points that change assignment at E-step to track convergence
	int nit = 0;
	int nchg = 999;
	while(nchg > 0 && nit < 50){
		kmc.Estimate();
		kmc.Update();
		nchg = (int)kmc.EvalLogL();
		//check for convergence with number of points that 
		//change assignment
		nit++;
	}
	vector<Matrix> mus;
	kmc.GetMeans(mus);
	for(int k = 0; k < m_k; k++)
		m_model[k]->SetParameter("mean",mus[k]);		
	

	//seed covariance matrix from identity 
	for(int k = 0; k < m_k; k++){
		Matrix S = Matrix(m_dim,m_dim);
		S.InitIdentity();
		m_model[k]->SetParameter("cov",S);
	}

	m_post.SetDims(m_n, m_k);
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++)
			m_post.SetEntry(1./double(m_k),n,k);

//	cout << "InitParameters - end" << endl;
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
			double g = m_model[k]->Prob(m_data->at(n));	
			gaus.SetEntry(g,n,k);
			norm += m_coeffs[k]*gaus.at(n,k);
		}
		post_norms.push_back(norm);
	}
	//calculate posterior for each data pt n in each cluster k
	for(int	n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			//weight by data weight
			val = m_data->at(n).w()*m_coeffs[k]*gaus.at(n,k)/post_norms[n];
			m_post.SetEntry(val,n,k);		
		}
	}
}



//M-step
//equations were derived from maximizing the posterior calculated in the E-step
void GaussianMixture::UpdateParameters(){
	//re-calculate normalization (overwrites)
	//this is for N_k - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms[k] = 0;
		m_norms_unwt[k] = 0;
		for(int n = 0; n < m_n; n++){
			m_norms[k] += m_post.at(n,k);
			//weighted in posterior calculation - need to unweight
			m_norms_unwt[k] += m_post.at(n,k)/m_data->at(n).w();
		}
	}
	//set new means + coeffs	
	Matrix mu = Matrix(m_dim, 1);
	m_coeffs.clear();

	for(int k = 0; k < m_k; k++){
		m_coeffs.push_back(m_norms[k]/m_n);
		//clear and overwrite m_mu[k]
		mu.clear();
		mu.InitEmpty();
		//m_mus[k].clear();
		//m_mus[k].initempty();
		//check above does what we want (clear the entries and reset an empty matrix)
		for(int n = 0; n < m_n; n++){
			//add data pt x_n,
			Matrix x = Matrix(m_data->at(n).Value());
			//weighted by posterior value gamma(z_nk),
			x.mult(x,m_post.at(n,k));
			//to new mu for cluster k
			mu.add(x);
			//m_mus[k].add(x);
		}
		//normalized by sum of posteriors for cluster k
		mu.mult(mu,1./m_norms[k]);
		//m_mus[k].mult(m_mus[k],1./m_norms[k]);
		m_model[k]->SetParameter("mean",mu);		
	}


	//sigma_k = 1/N_k sum_n(gamma(z_nk)*(x_n - mu_k)*(x_n - mu_k)T) for mu_k = mu^new_k
	for(int k = 0; k < m_k; k++){
		//create (x_n - mu)*(x_n - mu)T matrices for each data pt
		Matrix new_cov = Matrix(m_dim,m_dim);
		for(int n = 0; n < m_n; n++){
			Matrix cov_k = Matrix(m_dim, m_dim);
			mu = m_model[k]->GetParameter("mean");
			//construct x - mu
			Matrix x_min_mu = Matrix(m_data->at(n).Value());
			x_min_mu.minus(mu);
	
			//transpose x - mu
			Matrix x_min_muT;
			x_min_muT.transpose(x_min_mu);
			//(x_n - mu_k)*(x_n - mu_k)T
			cov_k.mult(x_min_mu,x_min_muT);
			//weighting by posterior gamma(z_nk)
			//sum over n
			new_cov.add(cov_k);
		}	
		//normalize by N_k
		new_cov.mult(new_cov,1./m_norms[k]);
		//overwrites m_covs[k]
		//m_covs[k] = new_cov;
		//m_model[k]->SetParameter("cov",m_covs[k]);		
		m_model[k]->SetParameter("cov",new_cov);		
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




map<string, Matrix> GaussianMixture::GetParameters(int k){ 
	map<string, Matrix> p;
	p["mean"] = m_model[k]->GetParameter("mean");
	p["cov"] = m_model[k]->GetParameter("cov");
	p["pi"] = Matrix(m_coeffs[k]);
	return p;
};



map<string, Matrix> GaussianMixture::GetPriorParameters(int k){ 
	map<string, Matrix> p;
	if(k >= m_k) return p;
	p["mean"] = m_model[k]->GetParameter("mean");
	p["cov"] = m_model[k]->GetParameter("cov");
	p["pi"] = Matrix((m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_data->Sumw()));
	p["scalemat"] = m_model[k]->GetPrior()->GetParameter("scalemat");
	p["m"] = m_model[k]->GetPrior()->GetParameter("mean");
	p["scale"] = m_model[k]->GetPrior()->GetParameter("scale");
	p["dof"] = m_model[k]->GetPrior()->GetParameter("dof");
	p["alpha"] = Matrix(m_alphas[k]);
	
	return p;
};


//variational stuff 
void GaussianMixture::InitPriorParameters(unsigned long long seed){
//cout << "INIT PRIOR PARAMS" << endl;
	if(m_dim == 0){
		cout << "VarGaussianMixture Initialize - Error: data has not been set." << endl;
		return;
	}

	//assuming conjugate prior - normal wishart (using precision matrix and corresponding update equations)
	for(int k = 0; k < m_k; k++){m_model[k]->SetDim(m_dim); m_model[k]->SetPrior(new NormalWishart(m_dim));}	

	//beta > 0
	m_beta0 = 1e-3;
	//cout << "beta0: " << m_beta0 << endl;
	//m > 0
	m_mean0 = Matrix(m_dim,1);
	//choose m_0 = 0 by symmetry (see Bishop eq. 10.40)
	m_mean0.InitEmpty();
	
	//nu > d - 1 (degrees of freedom)
	m_nu0 = m_dim;// - 1) + 1e-3;
	//cout << "nu0: " << m_nu0 << endl;
	
	m_meanBeta0 = Matrix(m_dim, 1);
	m_meanBeta0.mult(m_mean0, m_beta0);
	//W <- R in dxd space - is a covariance matrix
	m_W0 = Matrix(m_dim, m_dim);
	//this is the inverse prior covariance
	m_W0.InitIdentity();
	//least informative prior - nu0^-1*sigma0^-1
	m_W0.mult(m_W0,1./m_nu0);
	//need W0 inverse for parameter calculation
	m_W0inv = Matrix(m_dim,m_dim);
	m_W0inv.invert(m_W0);


	m_post.SetDims(m_n, m_k);
	
	m_Elam.clear();
	m_Epi.clear();
	for(int k = 0; k < m_k; k++){
		m_Elam.push_back(0.);
		m_Epi.push_back(0.);
	}

	//init parameters from alpha0
	//assuming a Dirichlet prior on the multinomial (categorical) assignment distribution (over latent variable z - sets pis)
	for(int k = 0; k < m_k; k++){
		m_alphas[k] = m_alpha0;
	}
	//to init prior parameters without calculating Rstats from posterior
	UpdatePriorParameters();

 
}






void GaussianMixture::CalculateExpectations(){
	//cout << "CALC EXPECTATIONS" << endl;
	//calculate alpha_hat
	double alpha_hat = 0.;
	//alpha_hat = sum_k alpha_k
	for(int k = 0; k < m_k; k++)
		alpha_hat += m_alphas[k];
	
	//calculate Elam (10.65) and Epi (10.66)
	double digam, dof;
	Matrix scalemat = Matrix(m_dim, m_dim);
	for(int k = 0; k < m_k; k++){
		scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
 		dof = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		digam = 0;
		for(int d = 1; d < m_dim+1; d++){
			digam += digamma((dof + 1 - d)/2.);
		}
	//	cout << "k: " << k << " digam: " << digam << " d*ln2: " << m_dim*log(2) << " lndet: " << log(scalemat.det()) << " det: " << scalemat.det() << endl; scalemat.Print();
		m_Elam[k] = digam + m_dim*log(2) + log(scalemat.det());
		m_Epi[k] = digamma(m_alphas[k]) - digamma(alpha_hat);
		if(isnan(m_Elam[k])){ cout << "NAN!!!!! k: " << k << " alpha: " << m_alphas[k] << " dof: " << dof << " Elam: " << m_Elam[k] << " Epi: " << m_Epi[k] << " detW[k]: " << scalemat.det() << " W[k]: " << endl;
		scalemat.Print(); cout << "W0" << endl; m_W0.Print();}
	}	
}





//E-step - calculate posterior
//E[z_nk] = r_nk
//(10.46) ln(rho_nk) = E[ln(pi_k)] + 1/2E[ln|lambda_k|] - D/2*ln(2pi) - 1/2E[(x_n - mu_k)Tlambda_k(x_n - mu_k)]
//(10.49) r_nk = rho_nk/sum_k rho_nk
//(10.64) ln(rho_nk) = psi(alpha_k) - psi(alpha_hat) + 1/2(sum^d_i psi( (nu_k + 1 - i) /2) + d*ln2 + ln|W_k| - D/2*ln(2pi) - 1/2( D*beta_k^inv + nu_k*(x_n - m_k)T*W_k*(x_n - m_k) )
void GaussianMixture::CalculateVariationalPosterior(){
//cout << "CALCULATE POSTERIOR - E STEP" << endl;
	//calculate necessary expectation values for E-step and ELBO
	CalculateExpectations();
	double E_mu_lam, post, norm, dof, scale;
	//vector<double> post_norms;
	vector<double> post_norms_adj; //norms for each pt (n vals)
	vector<double> post_k_vals; //find max value for normalization (k vals for each pt)
	vector<double> post_n_max; //max k val for nth pt (n vals)
	Matrix x_mat, x_min_m, x_min_mT;

	Matrix mean = Matrix(m_dim, 1);
	Matrix scalemat = Matrix(m_dim, m_dim);
	//cout << m_n << " points and " << m_k << " clusters " << endl;
	for(int n = 0; n < m_n; n++){
		norm = 0;
		for(int k = 0; k < m_k; k++){
			dof = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
			scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
			scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
			mean = m_model[k]->GetPrior()->GetParameter("mean");	
			//nu_k*(x_n - m_k)T*W_k*(x_n - m_k)
			x_min_m = Matrix(m_dim,1);
			x_mat = Matrix(m_data->at(n));
		//	if(k == 1){cout << "n: " << n << " x" << endl; x_mat.Print(); cout << "m[k]" << endl; mean.Print();}
			x_min_m.minus(x_mat,mean);		
			//if(n == 47){cout << "x - m[k]" << endl;x_min_m.Print();}
			x_min_mT = Matrix(1, m_dim);
			x_min_mT.transpose(x_min_m);
			//full term
			//if(n == 47){cout << "W[k]:" << endl; scalemat.Print();}
			Matrix transp_W = Matrix(1,m_dim);
			transp_W.mult(x_min_mT,scalemat);
			//if(n == 47){cout << "(x - m[k])T*W[k]" << endl; transp_W.Print();}
			Matrix full = Matrix(1,1);
			full.mult(transp_W,x_min_m);
			//if(n == 47){cout << "(x - m[k])T*W[k]*(x - m[k])" << endl; full.Print();}
			E_mu_lam = m_dim/scale + dof*full.at(0,0);	
			//gives ln(rho_nk)
			post = m_Epi[k] + 0.5*m_Elam[k] - (m_dim/2.)*log(2*acos(-1)) - 0.5*E_mu_lam;
		//	cout << "x - m[k]" << endl;
		//	x_min_m.Print();
		//	cout << "W[k]:" << endl;
		//	scalemat.Print();
		//	cout << "(x - m[k])T*W[k]*(x - m[k])" << endl; full.Print();
			post_k_vals.push_back(post);
		
			//post = exp(post);
			norm += post;
			//need to normalize
			m_post.SetEntry(post, n, k);
			//if(k == 1){ cout << std::setprecision(10) << "n: " << n << " k: " << k << " scale: " << scale << " dof: " << dof << " Elam: " << m_Elam[k] << " E_pi: " << m_Epi[k] << " E_mu_lam: " << E_mu_lam << " post: " << post << " mat post: " << m_post.at(n,k) << " full: " << full.at(0,0) << endl;}
		}
		if(post_k_vals.size() < 1) continue;
		post_n_max.push_back(*std::max_element(post_k_vals.begin(), post_k_vals.end()));
		post_norms_adj.push_back(0.);
		for(int k = 0; k < m_k; k++) post_norms_adj[n] += exp(post_k_vals[k] - post_n_max[n]);

		post_k_vals.clear();
		//post_norms.push_back(norm);
	}
	//cout << "posterior pre-norm" << endl;	
	//m_post.Print();
	//normalize
	for(int n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			//will lead to nan
	//	 cout << "Entry at n: " << n << " k: " << k << " is " << m_post.at(n,k) << " weight - " << m_data->at(n).w() << " point  " << endl; m_data->at(n).Print(); cout << "post_norms_adj: " << post_norms_adj[n] << " post_n_max: " << post_n_max[n] << " mu_k: " << endl; m_model[k]->GetPrior()->GetParameter("mean").Print();  
			//if(post_norms[n] == 0){ cout << "Entry at n: " << n << " k: " << k << " is " << m_post.at(n,k) << " weight - " << m_data->at(n).w() << " point  " << endl; m_data->at(n).Print(); cout << "m_k: " << endl; m_model[k]->GetPrior()->GetParameter("mean").Print(); } 
			//weight by data weight and adjusted by max ln(p_nk)
			m_post.SetEntry(m_data->at(n).w()*exp(m_post.at(n,k) - post_n_max[n])/post_norms_adj[n],n,k);
			
	
			//put in safeguard for computer precision for doubles (~1e\pm308)/rounding
			if(m_post.at(n,k) < 1e-308) m_post.SetEntry(0.,n,k);

			//if(m_post.at(n,k) > 0 && isinf(1/m_post.at(n,k))){ cout << "Entry at n: " << n << " k: " << k << " is " << m_post.at(n,k) << " weight - " << m_data->at(n).w() << " point  " << endl; m_data->at(n).Print(); cout << "post_norms_adj: " << post_norms_adj[n] << " post_n_max: " << post_n_max[n] << " mu_k: " << endl; m_model[k]->GetPrior()->GetParameter("mean").Print(); } 


			//m_post.SetEntry(m_data->at(n).w()*m_post.at(n,k)/post_norms[n],n,k);
			//uncomment here to check posterior values
		//	cout << "m post = " << m_post.at(n,k) << endl;
			//if(k == 1) cout << "k: " << k << " n: " << n << " post: " << m_post.at(n,k) << " norm: " << post_norms[n] << endl;
		}
	}

	//cout << "posterior normed" << endl;	
	//m_post.Print();
//cout << "\n" << endl;
};




void GaussianMixture::CalculateRStatistics(){
	//cout << "Calculate RStats" << endl;
	//responsibility statistics
	//this is for N_k (Bishop eq. 10.51) - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms[k] = 0;
		m_norms_unwt[k] = 0;
		for(int n = 0; n < m_n; n++){
			//if(k == 0) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
			//if(n == 3) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
			
			//weighted in posterior calculation - need to unweight
			m_norms[k] += m_post.at(n,k);
			m_norms_unwt[k] += m_post.at(n,k)/m_data->at(n).w();
			//cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
		}
		//cout << "k: " << k << " norm: " << m_norms[k] << endl;
		//if(m_norms[k] > 0 && isinf(1/m_norms[k])){ cout << "k: " << k << " m_norms: " << m_norms[k] << endl; m_post.Print();cout << " w1: " << m_data->at(0).w() << " w2: " << m_data->at(1).w() << endl; }
	}

	//this is for x_k (eq. 10.52) - k dx1 matrices
	for(int k = 0; k < m_k; k++){
		Matrix mu = Matrix(m_dim, 1);
		//to avoid nans, if there are no effective points in a cluster, the update equations dictate that its parameters are just the priors
		//however in calculating the r statistics, N_k = 0 can lead to nans
		if(m_norms[k] == 0){
			m_model[k]->SetParameter("mean",mu);		
			continue;	
		}
		for(int n = 0; n < m_n; n++){
			//add data pt x_n,
			Matrix x = Matrix(m_data->at(n));
			//if(n == 0){
			//cout << "n: " << n << " k: " << k << endl;
			//cout << "post: " << m_post.at(n,k) << endl;
			//}
			//cout << "x" << endl;
			//x.Print();
			//weighted by posterior value gamma(z_nk),
			//if(k == 1) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << " data weight: " << m_data->at(n).w() << endl;
			x.mult(x,m_post.at(n,k));
			//cout << "post*x" << endl;
			//x.Print();
			//add to new mu for cluster k
			mu.add(x);
		}
		//cout << "k: " << k << " sum_n post*x" << endl;
		//mu.Print();
		//normalized by sum of posteriors for cluster k
		mu.mult(mu,1./m_norms[k]);
		//cout << "k: " << k << " m_norm: " << m_norms[k] << endl;
		//cout << "1/N[k]*(sum_n post*x)" << endl;
		//mu.Print();
		m_model[k]->SetParameter("mean",mu);		
	}
	


//		cout << std::setprecision(10) << endl;	
	//this is for S_k (eq.10.53) - k dxd matrices
	for(int k = 0; k < m_k; k++){
		//create (x_n - mu)*(x_n - mu)T matrices for each data pt
		Matrix S = Matrix(m_dim,m_dim);
		Matrix mu = m_model[k]->GetParameter("mean");
		//to avoid nans, if there are no effective points in a cluster, the update equations dictate that its parameters are just the priors
		//however in calculating the r statistics, N_k = 0 can lead to nans
		if(m_norms[k] == 0){
			m_model[k]->SetParameter("cov",S);		
			continue;	
		}

//		cout << "k: " << k << " mu" << endl;
//		mu.Print();
		for(int n = 0; n < m_n; n++){
			//construct x - mu
			Matrix x_mat = Matrix(m_data->at(n));
			//cout << "n: " << n << " k: " << k << endl;
			//cout << "mu:" << endl;
			//mu.Print();
			//cout << "x:" << endl;
			//x_mat.Print();
			Matrix x_min_mu = Matrix(m_dim, 1);
			x_min_mu.minus(x_mat,mu);
			//cout << "x - mu" << endl;
			//x_min_mu.Print();	
			//transpose x - mu
			Matrix x_min_muT = Matrix(1, m_dim);
			x_min_muT.transpose(x_min_mu);
		
			Matrix S_k = Matrix(m_dim, m_dim);	
			//(x_n - mu_k)*(x_n - mu_k)T
			S_k.mult(x_min_mu,x_min_muT);
			//cout << "(x - mu)*(x - mu)T" << endl;	
			//S_k.Print();
			//weighting by posterior r_nk
			S_k.mult(S_k,m_post.at(n,k));
		//	if(n == 3) cout << "cov calc - post: " << m_post.at(n,k) << endl;
			//cout << "post*(x - mu)*(x - mu)T" << endl;	
			//S_k.Print();
			//sum over n
			S.add(S_k);
		}	
		//cout << "sum_n post*(x - mu)*(x - mu)T" << endl;	
		//S.Print();
		//normalize by N_k
		S.mult(S,1./m_norms[k]);
		//cout << "m_norm: " << m_norms[k] << endl;
		//cout << "k: " << k << " norm: " << m_norms[k] << " alpha: " << m_alphas[k] << " (1/N[k])*sum_n post*(x - mu)*(x - mu)T" << endl;	
		//if data smear is specified - provides lower bound on covariance -> regularization and provides nonzero covariance in single point case
		//N_k = 0 clusters do not get a smear because there are no points to smear
		if(_smear) S.add(_data_cov);
//		S.Print();
		m_model[k]->SetParameter("cov",S);
		//if(k == 1){ cout << "CalculateRStats - cov" << endl; m_model[k]->GetParameter("cov").Print();}
	}

}


//M-step
void GaussianMixture::UpdateVariationalParameters(){
	cout << "GaussianMixture::UpdateVariationalParameters - start" << endl;
	CalculateRStatistics();
	//can't remove N_k = 0 clusters because alpha0 keeps these clusters alive -> instead these parameters will be only priors
	UpdatePriorParameters();
	cout << "GaussianMixture::UpdateVariationalParameters - end" << endl;

}

//M-step
void GaussianMixture::UpdatePriorParameters(){
//cout << "UPDATE PARAMETERS - M STEP" << endl;

	//now update variational distribution parameters
	//updating based on first step (sub-0 params)
	for(int k = 0; k < m_k; k++){
		//to avoid nans, if there are no effective points in a cluster, the update equations dictate that its parameters are just the priors
		//however in calculating the r statistics, N_k = 0 can lead to nans
		if(m_norms[k] == 0){
			m_model[k]->GetPrior()->SetParameter("scale", Matrix(m_beta0));
			m_model[k]->GetPrior()->SetParameter("dof", Matrix(m_nu0));
			m_model[k]->GetPrior()->SetParameter("mean", m_mean0);
			m_model[k]->GetPrior()->SetParameter("scalemat", m_W0);
			continue;	
		}



		//already calculated
		Matrix mu = m_model[k]->GetParameter("mean");
		Matrix cov = m_model[k]->GetParameter("cov");
		
		//alphas - eq. 10.58 (all the same in vector)
		m_alphas[k] = m_alpha0 + m_norms[k];
		//cout << "k: " << k << " norm: " << m_norms[k] << " alpha: " << m_alphas[k] << endl;	
		//cout << "k: " << k << " alpha: " << m_alphas[k] << endl;	
	
		//betas - eq. 10.60
		double new_scale = m_beta0 + m_norms[k];
		m_model[k]->GetPrior()->SetParameter("scale", Matrix(new_scale));
		//cout << "k: " << k << " scale: " << m_model[k]->GetPrior()->GetParameter("scale").at(0,0) << endl;	
		//nus - eq. 10.63
		double new_dof = m_nu0 + m_norms[k];
		m_model[k]->GetPrior()->SetParameter("dof", Matrix(new_dof));
		//cout << "k: " << k << " dof: " << m_model[k]->GetPrior()->GetParameter("dof").at(0,0) << endl;	
		
		Matrix new_mean = Matrix(m_dim, 1);
		//means - eq. 10.61
		//N_k*bar{x}_k
		new_mean.mult(mu,m_norms[k]);
		//m_meanBeta0 + N_k*bar{x}_k
		new_mean.add(m_meanBeta0);
		//normalize by 1/beta
		new_mean.mult(new_mean,1./new_scale);
		m_model[k]->GetPrior()->SetParameter("mean", new_mean);
		//cout << "k: " << k << " mean: " << endl; 
		//m_model[k]->GetPrior()->GetParameter("mean").Print();	
		//cout << "k: " << k << " xbar: " << endl;
		//mu.Print();

		//Ws - eq. 10.62
		//caluclated for W inverse
		//beta0*N_k/(beta0 + N_k)*(bar{x}_k - m0)(bar{x}_k - m0)T
		Matrix new_scalemat = Matrix(m_dim, m_dim);
		//bar{x}_k - m0
		Matrix x_min_mean = Matrix(m_dim, 1);
		x_min_mean.minus(mu, m_mean0);
		//cout << "xbar" << endl;
		//mu.Print();
		Matrix x_min_meanT = Matrix(1, m_dim);
		x_min_meanT.transpose(x_min_mean);
		//(bar{x}_k - m0)(bar{x}_k - m0)T
		new_scalemat.mult(x_min_mean, x_min_meanT);
		//cout << "xxT" << endl;
		//new_scalemat.Print();
		double prefactor = m_beta0*m_norms[k]/(m_beta0 + m_norms[k]);
		new_scalemat.mult(new_scalemat, prefactor);
		//cout << "(b*N)/(b + N)*xxT" << endl;
		//new_scalemat.Print();
		//N_k*S_k
		Matrix scaledS = Matrix(m_dim, m_dim);
		scaledS.mult(cov,m_norms[k]);
		//add first two terms to last term
		new_scalemat.add(scaledS);
		//cout << "N*S + (b*N)/(b + N)*xxT" << endl;
		//new_scalemat.Print();
		//add W0inv
		new_scalemat.add(m_W0inv);
		//cout << "W-1 + N*S + (b*N)/(b + N)*xxT" << endl;
		//new_scalemat.Print();
		//invert (calculated for W_k inverse)
		new_scalemat.invert(new_scalemat);

		if(isnan(new_scalemat.at(0,0))) cout << "W IS NAN!!!!! for cluster " << k << " m_norms: " << m_norms[k] << endl;
		m_model[k]->GetPrior()->SetParameter("scalemat", new_scalemat);
		//cout << "k: " << k << " cov: " << endl;
		//cov.Print();
		//cout << "k: " << k << " scalemat: " << endl; 
		//m_model[k]->GetPrior()->GetParameter("scalemat").Print();	

	}
};





//calculates ELBO
//(10.70) ELBO = E[ln(p(X|Z,mu,lam))] + E[ln(p(Z|pi)] + E[ln(p(pi))] + E[ln(p(mu,lam))] - E[ln(q(Z))] - E[ln(q(pi))] - E[ln(q(mu,lam))]
double GaussianMixture::EvalVariationalLogL(){
	double E_p_all, E_p_Z, E_p_pi, E_p_muLam, E_q_Z, E_q_pi, E_q_muLam, E_lam;
	//E[ln p(X|Z,mu,lam)] = 0.5*sum_k( N_k*(ln~lam_k - m_dim/beta_k - nu_k*Tr(S_k*W_k) - nu_k*(mus_k - m_k)T*W_k*(mu_k - m_k) - D*log(2*pi) ))
	E_p_all = 0;
	//recalculate expectation values E_lam + E_pi with updated parameters
	CalculateExpectations();
	for(int k = 0; k < m_k; k++){
		double scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		double nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		Matrix mean = m_model[k]->GetPrior()->GetParameter("mean");
		Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		Matrix mu = m_model[k]->GetParameter("mean");
		Matrix cov = m_model[k]->GetParameter("cov");
		////cout << "k: " << k << " scale: " << scale << " dof: " << nu << " norm: " << m_norms[k] << " alpha: " << m_alphas[k] << endl;
	//	cout << "m" << endl;
	//	mean.Print();
	//	cout << "W" << endl;
	//	scalemat.Print();
	
	//	cout << "mean" << endl;
	//	m_model[k]->GetParameter("mean").Print();
	//	cout << "cov" << endl;
	//	cov.Print();
		
		//(x_n - m_k)
		//m_xbars[k].Print();
		Matrix xbar_min_m = Matrix(m_dim,1);
		xbar_min_m.minus(mu, mean);
		Matrix xbar_min_mT = Matrix(1, m_dim);
		xbar_min_mT.transpose(xbar_min_m);
		//(x_n - m_k)T*W_k*(x_n - m_k)
		Matrix xbarT_x_W = Matrix(1,m_dim);
		xbarT_x_W.mult(xbar_min_mT,scalemat);
		Matrix full = Matrix(1,1);
		full.mult(xbarT_x_W,xbar_min_m);
		//tr(s_k*w_k)
		//S_k*W_k = dxd matrix
		Matrix tmp_S_W = Matrix(m_dim, m_dim);
		tmp_S_W.mult(cov,scalemat);
		E_p_all += m_norms[k]*(m_Elam[k] - m_dim/scale - nu*tmp_S_W.trace()  - nu*full.at(0,0) - m_dim*log(2*acos(-1)));

	}
	E_p_all *= 0.5;
//	cout << "E_p_all: " << E_p_all << endl;

	//E[ln(p(Z|pi))] = sum_n sum_k r_nk*ln(m_Epi[k])
	E_p_Z = 0;
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++)
			E_p_Z += m_post.at(n,k)*m_Epi[k];
//	cout << "E_p_Z: " << E_p_Z << endl;
	
	//E[ln(p(pi))] = ln(C(alpha)) + (alpha_0 - 1)*sum_k m_Epi[k]
	E_p_pi = 0;
	for(int k = 0; k < m_k; k++)
		E_p_pi += m_Epi[k];
	E_p_pi *= (m_alpha0 - 1);
	vector<double> alpha0s(m_k, m_alpha0);
	Dirichlet* dir = new Dirichlet(alpha0s);
	E_p_pi += dir->lnC();
//	cout << "E_p_pi: " << E_p_pi << endl;

	//E[ln(p(mu,lambda))] = 1/2*sum_k(D*ln(beta0/2*pi) + m_Elam[k] - D*beta0/beta_k - beta0*nu_k(m_k - m0)T*W_kW*(m_k - m0))
	E_p_muLam = 0;
	double lam_sum = 0; 
	double half_sum = 0; 
	double tr_sum = 0;
	Matrix tmp, tmpT;
	for(int k = 0; k < m_k; k++){
		//(m_k - m_0)
		double scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		double nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		Matrix mean = m_model[k]->GetPrior()->GetParameter("mean");
		Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		
		Matrix mk_min_m0 = Matrix(m_dim,1);
		mk_min_m0.minus(mean,m_mean0);		
		Matrix mk_min_m0T = Matrix(1, m_dim);
		mk_min_m0T.transpose(mk_min_m0);
		//(m_k - m_0)T*W_k*(m_k - m_0)
		Matrix mT_x_W = Matrix(1,m_dim);
		mT_x_W.mult(mk_min_m0T,scalemat);
		Matrix full = Matrix(1,1);
		full.mult(mT_x_W,mk_min_m0);
		
		half_sum += (m_dim*log(m_beta0/(2*acos(-1))) + m_Elam[k] - m_dim*m_beta0/scale - m_beta0*nu*full.at(0,0));
		lam_sum += m_Elam[k];
		
		Matrix tr = Matrix(m_dim,m_dim);
		tr.mult(m_W0inv,scalemat);	
	
		tr_sum += nu*tr.trace();
	}
	E_p_muLam = 0.5*half_sum + ((m_nu0 - m_dim - 1)/2.)*lam_sum - 0.5*tr_sum;
	Wishart* wish = new Wishart(m_W0, m_nu0);
	E_p_muLam += m_k*wish->lnB();
//	cout << "E_p_muLam: " << E_p_muLam << endl;

	//E[ln(q(Z)]
	E_q_Z = 0;
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++){
			if(m_post.at(n,k) == 0) E_q_Z += 0;	
			else E_q_Z += m_post.at(n,k)*log(m_post.at(n,k));	
		}
//	cout << "E_q_Z: " <<  E_q_Z << endl;
	

	//E[ln(q(pi))]
	E_q_pi = 0;
	for(int k = 0; k < m_k; k++){
		E_q_pi += (m_alphas[k] - 1)*m_Epi[k];
	}
	Dirichlet* dir_k = new Dirichlet(m_alphas);
	E_q_pi += dir_k->lnC();
//	cout << "E_q_pi: " << E_q_pi << endl;
	
	//E[ln(q(mu, lam))]
	E_q_muLam = 0;
	double H;
	for(int k = 0; k < m_k; k++){
		Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		double nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		double scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		Wishart* wish_k = new Wishart(scalemat, nu);
		H = wish_k->H();
		E_q_muLam += 0.5*m_Elam[k] + m_dim/2.*log(scale/(2*acos(-1))) - m_dim/2. - H;

	}
//	cout << "E_q_muLam: " << E_q_muLam << endl;
	/*
	cout << "E_p_all: " << E_p_all << endl;
	cout << "E_p_Z: " << E_p_Z << endl;
	cout << "E_p_pi: " << E_p_pi << endl;
	cout << "E_p_muLam: " << E_p_muLam << endl;
	cout << "E_q_Z: " <<  E_q_Z << endl;
	cout << "E_q_pi: " << E_q_pi << endl;
	cout << "E_q_muLam: " <<  E_q_muLam << endl;
	cout << "m_post" << endl;
	*/
	return E_p_all+ E_p_Z + E_p_pi+ E_p_muLam - E_q_Z - E_q_pi - E_q_muLam;
	


}; //end ELBO

