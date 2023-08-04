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

//don't forget to include coeffs (eventually)
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


		m_coeffs[k] = randy.SampleFlat();
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
	Matrix mu = Matrix(m_dim, 1);
	for(int k = 0; k < m_k; k++){
		m_coeffs[k] = m_norms[k]/m_n;
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




vector<map<string, Matrix>> GaussianMixture::GetParameters(){ 
	vector<map<string, Matrix>> params;
	for(int k = 0; k < m_k; k++){
		map<string, Matrix> p;
		p["mean"] = m_model[k]->GetParameter("mean");
		p["cov"] = m_model[k]->GetParameter("cov");
		p["pi"] = Matrix(m_coeffs[k]);
		params.push_back(p);
		p.clear();
	}
	return params;
};



vector<map<string, Matrix>> GaussianMixture::GetPriorParameters(){ 
	vector<map<string, Matrix>> params;
	for(int k = 0; k < m_k; k++){
		map<string, Matrix> p;
		p["mean"] = m_model[k]->GetParameter("mean");
		p["cov"] = m_model[k]->GetParameter("cov");
		p["pi"] = Matrix((m_alphas[k] + m_norms[k])/(k*m_alpha0 + m_n));
		p["scalemat"] = m_model[k]->GetPrior()->GetParameter("scalemat");
		p["mean"] = m_model[k]->GetPrior()->GetParameter("mean");
		p["scale"] = m_model[k]->GetPrior()->GetParameter("scale");
		p["dof"] = m_model[k]->GetPrior()->GetParameter("dof");
		p["alpha"] = Matrix(m_alphas[k]);
		
		params.push_back(p);
		p.clear();
	}
	return params;
};


//variational stuff 
void GaussianMixture::InitPriorParameters(unsigned long long seed){
	if(m_dim == 0){
		cout << "VarGaussianMixture Initialize - Error: data has not been set." << endl;
		return;
	}

	//assuming conjugate prior - normal inverse wishart
	for(int k = 0; k < m_k; k++) m_model[k]->SetPrior(new NormalInvWishart());	


	//alpha > 0
		//choose the same value for all alpha_0k by symmetry (see Bishop eq. 10.39)
	m_alpha0 = 1;
	
	//beta > 0
	RandomSample rs;
	rs.SetRange(0.,1.);
	m_beta0 = rs.SampleFlat();
	//m > 0
	m_mean0 = Matrix(m_dim,1);
	//choose m_0 = 0 by symmetry (see Bishop eq. 10.40)
	m_mean0.InitEmpty();
	
	m_meanBeta0 = Matrix(m_dim, 1);
	m_meanBeta0.mult(m_mean0, m_beta0);
	//W <- R in dxd space - is a covariance matrix
	m_W0 = Matrix(m_dim, m_dim);
	m_W0.InitIdentity();
	//need W0 inverse for parameter calculation
	m_W0inv = Matrix(m_dim,m_dim);
	m_W0inv.invert(m_W0);

	//nu > d - 1 (degrees of freedom)
	rs.SetRange(m_dim - 1, m_dim+2);
	m_nu0 = rs.SampleFlat();
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
		m_alphas.push_back(m_alpha0);
	}

	rs.SetRange(0.,1.);
	for(int k = 0; k < m_k; k++){
	//	m_betas.push_back(k*rs.SampleFlat() + 0.1);//rs.SampleFlat());
		m_model[k]->GetPrior()->SetParameter("scale",Matrix(k*rs.SampleFlat() + 0.1));
	}
	
	Matrix mat;
	double mean_lower = m_data->min()-0.1;
	double mean_upper = m_data->max()+0.1;
	for(int k = 0; k < m_k; k++){
		mat = Matrix(m_dim, 1);
		mat.InitRandom(mean_lower, mean_upper, seed+k);
		m_means.push_back(mat);
		m_model[k]->GetPrior()->SetParameter("mean",m_means[k]);		
		mat.clear();
	}
	for(int k = 0; k < m_k; k++){
		mat = Matrix(m_dim, m_dim);
		mat.InitRandomSymPosDef(0.,1.,seed+k);//InitIdentity();
		//m_Ws.push_back(mat);
		m_model[k]->GetPrior()->SetParameter("scalemat",mat);		
		mat.clear();
	}
	rs.SetRange(m_dim,m_dim*1.5);
	//rs.SetRange(m_dim-1,m_dim+2.);
	for(int k = 0; k < m_k; k++){
		//m_nus.push_back(rs.SampleFlat());
		m_model[k]->GetPrior()->SetParameter("dof",Matrix(rs.SampleFlat()));	
	}

	//init responsibility statistics
	for(int k = 0; k < m_k; k++){
		//seed N_k to 1 - multiplicative factor in Update for Rstat quantities
		m_norms[k] = 1.;
		mat = Matrix(m_dim, m_dim);
		mat.InitIdentity();
		//m_Ss.push_back(Matrix(m_dim,m_dim));
		//m_Ss[k].InitIdentity();
		m_model[k]->SetParameter("cov",mat);
	}
	//to init xbars
	//m_xbars = get means
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
	vector<Matrix> xbars;
	kmc.GetMeans(xbars);
	for(int k = 0; k < m_k; k++) m_model[k]->SetParameter("mean",xbars[k]);
	//to init prior parameters without calculating Rstats from posterior
	UpdatePriorParameters();

}





void GaussianMixture::CalculateExpectations(){
	//calculate alpha_hat
	double alpha_hat = 0.;
	//alpha_hat = sum_k alpha_k
	for(int k = 0; k < m_k; k++)
		alpha_hat += m_alphas[k];
	
	//calculate Elam (10.65) and Epi (10.66)
	double digam, nu;
	Matrix scalemat = Matrix(m_dim, m_dim);
	for(int k = 0; k < m_k; k++){
		scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
 		nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		digam = 0;
		for(int d = 1; d < m_dim+1; d++)
			digam += digamma((nu + 1 - d)/2);
		m_Elam[k] = digam + m_dim*log(2) + log(scalemat.det());
		m_Epi[k] = digamma(m_alphas[k]) - digamma(alpha_hat);
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
	double E_mu_lam, post, norm, beta, nu;
	vector<double> post_norms;
	Matrix x_mat, x_min_m, x_min_mT;

	Matrix mean = Matrix(m_dim, 1);
	Matrix scalemat = Matrix(m_dim, m_dim);
	for(int n = 0; n < m_n; n++){
		norm = 0;
		for(int k = 0; k < m_k; k++){
			nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
			beta = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
			scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
			mean = m_model[k]->GetPrior()->GetParameter("mean");	
			nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);	
	
			//nu_k*(x_n - m_k)T*W_k*(x_n - m_k)
			x_min_m = Matrix(m_dim,1);
			x_mat = Matrix(m_data->at(n).Value());
			x_min_m.minus(x_mat,mean);		
			//x_min_m.minus(x_mat,m_means[k]);		
			x_min_mT = Matrix(1, m_dim);
			x_min_mT.transpose(x_min_m);
			//full term
			Matrix transp_W = Matrix(1,m_dim);
			transp_W.mult(x_min_mT,scalemat);
			Matrix full = Matrix(1,1);
			full.mult(transp_W,x_min_m);
			//E_mu_lam = m_dim/beta + nu*full.at(0,0);	
			E_mu_lam = m_dim/beta + nu*full.at(0,0);	
		//if(n == 0){ cout << "k: " << k << " beta: " << m_betas[k] << " nu: " << m_nus[k] << " full: " << full.at(0,0) << " x_n: " << endl; m_data->at(n).Print();}
			//gives ln(rho_nk)
			post = m_Epi[k] + 0.5*m_Elam[k] - (m_dim/2.)*log(2*acos(-1)) - 0.5*E_mu_lam;
		//	if(n == 0){ cout << "k: " << k << " n: " << n << " Epi: " << m_Epi[k] << " Elam: " << m_Elam[k] << " E_muLam: " << E_mu_lam << " beta: " << m_betas[k] << " nu: " << m_nus[k]  << " post: " << post << " exp(post): " << exp(post) << " W: " <<  endl;
		//	m_Ws[k].Print(); cout << "m: " << endl; m_means[k].Print(); cout << "x: " << endl; m_data->at(n).Print();
			//}	
			post = exp(post);
			norm += post;
			//need to normalize
			m_post.SetEntry(post, n, k);
		}
		post_norms.push_back(norm);
	}
//cout << "posterior pre-norm" << endl;	
	//m_post.Print();
	//normalize
	for(int n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			m_post.SetEntry(m_post.at(n,k)/post_norms[n],n,k);
		//uncomment here to check posterior values
		//	if(n == 0) cout << "k: " << k << " n: " << n << " post: " << m_post.at(n,k) << " norm: " << post_norms[n] << endl;
		}
	}
//cout << "posterior normed" << endl;	
//	m_post.Print();
//cout << "\n" << endl;
};




void GaussianMixture::CalculateRStatistics(){
	//responsibility statistics
	//this is for N_k (Bishop eq. 10.51) - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms[k] = 0;
		for(int n = 0; n < m_n; n++){
			if(k == 0 && n == 0 && isnan(m_post.at(n,k))) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
			m_norms[k] += m_post.at(n,k);
		}
//cout << "k: " << k << " norm: " << m_norms[k] << endl;
	}









	//this is for x_k (eq. 10.52) - k dx1 matrices
	Matrix mu = Matrix(m_dim, 1);
	for(int k = 0; k < m_k; k++){
//if(k == 0) cout << "k: " << k << " old x_bar: " << endl;
//if(k == 0) m_xbars[k].Print();
		//clear and overwrite m_xbar[k]
		mu.clear();
		mu.InitEmpty();
		//m_mus[k].clear();
		//m_mus[k].InitEmpty();
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
//if(k == 0) cout << "k: " << k << " new x_bar: " << endl;
//if(k == 0) m_xbars[k].Print();
	}
	


	//this is for S_k (eq.10.53) - k dxd matrices
	for(int k = 0; k < m_k; k++){
//if(k == 0)cout << "k: " << k << " old S: " << endl;
//if(k == 0)m_Ss[k].Print();
		//create (x_n - mu)*(x_n - mu)T matrices for each data pt
		Matrix new_S = Matrix(m_dim,m_dim);
		for(int n = 0; n < m_n; n++){
			Matrix S_k = Matrix(m_dim, m_dim);

			//construct x - mu
			Matrix x_mat = Matrix(m_data->at(n).Value());
			Matrix x_min_mu;
			x_min_mu.mult(m_mus[k],-1.);
			x_min_mu.add(x_mat);
			
			//transpose x - mu
			Matrix x_min_muT;
			x_min_muT.transpose(x_min_mu);
			
			//(x_n - mu_k)*(x_n - mu_k)T
			S_k.mult(x_min_mu,x_min_muT);
			
			//weighting by posterior r_nk
			S_k.mult(S_k,m_post.at(n,k));
			
			//sum over n
			new_S.add(S_k);
		}	
		//normalize by N_k
		new_S.mult(new_S,1./m_norms[k]);
		//overwrites m_Ss[k]
//		m_Ss[k] = new_S;
		m_model[k]->SetParameter("cov",new_S);
//if(k == 0) cout << "k: " << k << " new S: " << endl;
//if(k == 0) m_Ss[k].Print();
	}

}


//M-step
void GaussianMixture::UpdateVariationalParameters(){
	CalculateRStatistics();
	UpdatePriorParameters();

}

//M-step
void GaussianMixture::UpdatePriorParameters(){
//cout << "UPDATE PARAMETERS - M STEP" << endl;
	//now update variational distribution parameters
	//updating based on first step (sub-0 params)
//	m_betas.clear();
//	m_alphas.clear();
//	m_nus.clear();
	double prefactor;
	Matrix new_mean = Matrix(m_dim, 1);
	Matrix new_scalemat = Matrix(m_dim, m_dim);
	double new_dof = 0;
	double new_scale = 0;
	Matrix mu = Matrix(m_dim, 1);
	Matrix cov = Matrix(m_dim, m_dim);
	for(int k = 0; k < m_k; k++){
		//already calculated
		mu = m_model[k]->GetParameter("mean");
		cov = m_model[k]->GetParameter("cov");
		//alphas - eq. 10.58 (all the same in vector)
		m_alphas[k] = m_alpha0 + m_norms[k];
		//betas - eq. 10.60
//cout << "k: " << k << " old beta: " << m_betas[k] << endl;
		new_scale = m_beta0 + m_norms[k];	
		m_model[k]->GetPrior()->SetParameter("scale", new_scale);
		//m_betas[k] = m_beta0 + m_norms[k];
//cout << "k: " << k << " new beta: " << m_betas[k] << " norm: " << m_norms[k] << endl;	
		
		//means - eq. 10.61
//		m_means[k].clear();
//		m_means[k].InitEmpty();
		new_mean.clear(); new_mean.InitEmpty();
		//N_k*bar{x}_k
//		m_means[k].mult(mu,m_norms[k]);
		new_mean.mult(mu,m_norms[k]);
		//m_meanBeta0 + N_k*bar{x}_k
//		m_means[k].add(m_meanBeta0);
		new_mean.add(m_meanBeta0);
		//normalize by 1/beta
		//m_means[k].mult(m_means[k],1./m_betas[k]);
//		m_means[k].mult(m_means[k],1./new_scale);
		new_mean.mult(new_mean,1./new_scale);
		m_model[k]->GetPrior()->SetParameter("mean", new_mean);


		//Ws - eq. 10.62
//if(k == 0) cout << "k: " << k << " old W: " << endl;
//if(k == 0) m_Ws[k].Print();
		//caluclated for W inverse
		//beta0*N_k/(beta0 + N_k)*(bar{x}_k - m0)(bar{x}_k - m0)T
		new_scalemat.clear();
		new_scalemat.InitEmpty();		
	
		//m_Ws[k].clear();
		//m_Ws[k].InitEmpty();		
//if(k == 0) cout << "k: " << k << " 1 - clear W: " << endl;
//if(k == 0) m_Ws[k].Print();
		//bar{x}_k - m0
		Matrix x_min_mean = Matrix(m_dim, 1);
		x_min_mean.minus(x_min_mean,mu);
		Matrix x_min_meanT = Matrix(1, m_dim);
		x_min_meanT.transpose(x_min_mean);
//if(k == 0) cout << "k: " << k << " (x_n - m_0): " << endl;
//if(k == 0) x_min_mean.Print();
		//(bar{x}_k - m0)(bar{x}_k - m0)T
		new_scalemat.mult(x_min_mean, x_min_meanT);
//if(k == 0) cout << "k: " << k << " 2 - W = (x_n - m_0)(x_n - m_0)T: " << endl;
//if(k == 0) m_Ws[k].Print();
		prefactor = m_beta0*m_norms[k]/(m_beta0 + m_norms[k]);
//if(k == 0) cout << "beta0: " << m_beta0 << " norm: " << m_norms[k] << " prefactor: " << prefactor << endl;
		new_scalemat.mult(new_scalemat, prefactor);
//if(k == 0)cout << "k: " << k << " 3 - W = (beta0*N_k)/(beta0 + N_k)*(x_n - m_0)(x_n - m_0)T: " << endl;
//if(k == 0)m_Ws[k].Print();
		//N_k*S_k
		Matrix scaledS = Matrix(m_dim, m_dim);
		scaledS.mult(cov,m_norms[k]);
//if(k ==0 )cout << "k: " << k << " N_k*S_k: " << endl;
//if(k ==0 )scaledS.Print();
		//add first two terms to last term
		new_scalemat.add(scaledS);
//if(k == 0)cout << "k: " << k << " 4 - W = (beta0*N_k)/(beta0 + N_k)*(x_n - m_0)(x_n - m_0)T + N_k*S_k: " << endl;
//if(k == 0)m_Ws[k].Print();
		//add W0inv
		new_scalemat.add(m_W0inv);
//if(k == 0) cout << "W0inv" << endl;
//if(k == 0) m_W0inv.Print();
//if(k == 0) cout << "k: " << k << " 5 - W = (beta0*N_k)/(beta0 + N_k)*(x_n - m_0)(x_n - m_0)T + N_k*S_k + W0inv: " << endl;
//if(k == 0) m_Ws[k].Print();
		//Matrix testW;
		//testW.invert(m_Ws[k]);
		//invert (calculated for W_k inverse)
		new_scalemat.invert(new_scalemat);
		m_model[k]->GetPrior()->SetParameter("scalemat", new_scalemat);
//if(k == 0)cout << "k: " << k << " new W: " << endl;
//if(k == 0)m_Ws[k].Print();
//if(k == 0)cout << "k: " << k << " new W?: " << endl;
//if(k == 0)testW.Print();
//if(k == 0)cout << "k: " << k << " test for posdef" << endl;
//if(k == 0) Matrix testL = m_Ws[k].cholesky();
//if(k == 0) testL.Print();
	
	//cout << "old nu: " << m_nus[k] << endl;
	//	cout << "N_k: " << m_norms[k] << endl;
		//nus - eq. 10.63
		new_dof = m_nu0 + m_norms[k];
		m_model[k]->GetPrior()->SetParameter("dof", Matrix(new_dof));
	//cout << "new nu: " << m_nus[k] << endl;

	}
};





//calculates ELBO
//(10.70) ELBO = E[ln(p(X|Z,mu,lam))] + E[ln(p(Z|pi)] + E[ln(p(pi))] + E[ln(p(mu,lam))] - E[ln(q(Z))] - E[ln(q(pi))] - E[ln(q(mu,lam))]
double GaussianMixture::EvalVariationalLogL(){
//	cout << "EvalVarLogL begin" << endl;
	double E_p_all, E_p_Z, E_p_pi, E_p_muLam, E_q_Z, E_q_pi, E_q_muLam, E_lam;
	//E[ln p(X|Z,mu,lam)] = 0.5*sum_k( N_k*(ln~lam_k - m_dim/beta_k - nu_k*Tr(S_k*W_k) - nu_k*(mus_k - m_k)T*W_k*(mu_k - m_k) - D*log(2*pi) ))
	E_p_all = 0;
	double scale, nu;
	Matrix mean = Matrix(m_dim, 1);
	Matrix cov = Matrix(m_dim, m_dim);
	Matrix scalemat = Matrix(m_dim, m_dim);
	for(int k = 0; k < m_k; k++){
		scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		mean = m_model[k]->GetPrior()->GetParameter("mean");
		cov = m_model[k]->GetParameter("cov");
		scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		//(x_n - m_k)
	//	cout << "m_xbars" << endl;
		//m_xbars[k].Print();
		Matrix xbar_min_m = m_model[k]->GetParameter("mean");//m_mus[k];//Matrix(m_dim,1);
		xbar_min_m.minus(mean);
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
//	cout << "k: " << k << " norm: " << m_norms[k] << " Elam: " << m_Elam[k] << " m_dim: " << m_dim << " beta: " << m_betas[k] << " nu: " << m_nus[k] << " trace: " << tmp_S_W.trace() << " full: " << full.at(0,0) << endl;		
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
	vector<double> alpha0s(m_alpha0, m_k);
	Dirichlet* dir = new Dirichlet(alpha0s);
	E_p_pi += dir->lnC();
//	cout << "E_p_pi: " << E_p_pi << endl;

	//E[ln(p(mu,lambda))] = 1/2*sum_k(D*ln(beta0/2*pi) + m_Elam[k] - D*beta0/beta_k - beta0*nu_k(m_k - m0)T*W_kW*(m_k - m0))
	E_p_muLam = 0;
	Matrix tmp, tmpT;
	for(int k = 0; k < m_k; k++){
		//(m_k - m_0)
		mean = m_model[k]->GetPrior()->GetParameter("mean");
		scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		Matrix mk_min_m0 = Matrix(m_dim,1);
		mk_min_m0.minus(mean,m_mean0);		
		Matrix mk_min_m0T = Matrix(1, m_dim);
		mk_min_m0T.transpose(mk_min_m0);
		//(m_k - m_0)T*W_k*(m_k - m_0)
		Matrix mT_x_W = Matrix(1,m_dim);
		mT_x_W.mult(mk_min_m0T,scalemat);
		Matrix full = Matrix(1,1);
		full.mult(mT_x_W,mk_min_m0);
		
		E_p_muLam += 0.5*(m_dim*log(m_beta0/(2*acos(-1))) + m_Elam[k] - m_dim*m_beta0/scale - m_beta0*nu*full.at(0,0));
		E_p_muLam += (m_nu0 - m_dim - 1)/2.*m_Elam[k];
		
		Matrix tr = Matrix(m_dim,m_dim);
		tr.mult(m_W0inv,scalemat);	
	
		E_p_muLam += -0.5*nu*tr.trace();
	}
	Wishart* wish = new Wishart(m_W0, m_nu0);
	E_p_muLam += wish->lnB();
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
	for(int k = 0; k < 1; k++){
//	cout << "k: " << k << endl;
//	cout << "	Elam: " << m_Elam[k] << endl;
//	cout << "	beta: " << m_betas[k] << endl;
	scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
	Wishart* wish_k = new Wishart(scalemat, nu);
	H = wish_k->H();
//	cout << "	Wish_H: " << H << endl;
		E_q_muLam += 0.5*( m_Elam[k] + m_dim/2.*log(scale/(2*acos(-1))) - m_dim/2. );
		E_q_muLam -= H;
	}
//	cout << "E_q_muLam: " <<  E_q_muLam << endl;
/*
	cout << "E_p_all: " << E_p_all << endl;
	cout << "E_p_Z: " << E_p_Z << endl;
	cout << "E_p_pi: " << E_p_pi << endl;
	cout << "E_p_muLam: " << E_p_muLam << endl;
	cout << "E_q_Z: " <<  E_q_Z << endl;
	cout << "E_q_pi: " << E_q_pi << endl;
	cout << "E_q_muLam: " <<  E_q_muLam << endl;

*/



//	cout << "EvalVarLogL end" << endl;
	return E_p_all+ E_p_Z + E_p_pi+ E_p_muLam - E_q_Z - E_q_pi - E_q_muLam;
	


}; //end ELBO

