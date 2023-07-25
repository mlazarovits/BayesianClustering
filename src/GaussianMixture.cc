#include <boost/math/special_functions/digamma.hpp>
#include "GaussianMixture.hh"
#include "Gaussian.hh"
#include "KMeansCluster.hh"
#include "ClusterViz2D.hh"
#include "Dirichlet.hh"
#include "Wishart.hh"

using boost::math::digamma;

GaussianMixture::GaussianMixture(int k) : BasePDFMixture(k){
	for(int i = 0; i < k; i++){
		m_model.push_back(new Gaussian());
	}
	m_seed = 123;
}

void GaussianMixture::InitParameters(){
	m_n = m_data->GetNPoints();
	//run k-means clustering on data to convergence to get means
	//get means for each cluster 
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
	vector<int> cnts;
	kmc.GetCounts(cnts);
	kmc.GetMeans(m_xbars);

	//init covs to identity
	//randomly initialize mixing coeffs
	RandomSample randy(m_seed);
	randy.SetRange(0.,1.);
	double coeff_norm = 0;
	map<string, Matrix> params;
	for(int k = 0; k < m_k; k++){
		m_covs.push_back(Matrix(m_dim,m_dim));
		m_covs[k].InitIdentity();
		m_weights[k] = randy.SampleFlat();
		coeff_norm += m_weights[k];
		
		params["mean"] = m_xbars[k];
		params["cov"] = m_covs[k];	
	
		m_model[k]->SetParameters(params);
	}
	//make sure sum_k m_weights[k] = 1
	for(int k = 0; k < m_k; k++) m_weights[k] /= coeff_norm;
}


void GaussianMixture::InitPriorParameters(){
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
//	m_W0.InitRandomSymPosDef();
	//need W0 inverse for parameter calculation
	m_W0inv = Matrix(m_dim,m_dim);
	m_W0inv.invert(m_W0);

	//nu > d - 1 (degrees of freedom)
	rs.SetRange(m_dim - 1, m_dim+2);
	m_nu0 = rs.SampleFlat();
	//m_post.SetDims(m_n, m_k);
	

	m_Elam.clear();
	m_Epi.clear();
	for(int k = 0; k < m_k; k++){
		m_Elam.push_back(0.);
		m_Epi.push_back(0.);
	}

	//init parameters from alpha0
	for(int k = 0; k < m_k; k++){
		m_alphas.push_back(m_alpha0);
	}

	rs.SetRange(0.,1.);
	for(int k = 0; k < m_k; k++){
		m_betas.push_back(k*rs.SampleFlat() + 0.1);//rs.SampleFlat());
	}
	
	Matrix mat;
	double mean_lower = m_data->min()-0.1;
	double mean_upper = m_data->max()+0.1;
	for(int k = 0; k < m_k; k++){
		mat = Matrix(m_dim, 1);
		mat.InitEmpty();
		//mat.InitRandom(mean_lower, mean_upper, m_seed+k);
		m_means.push_back(mat);
		mat.clear();
	}
	for(int k = 0; k < m_k; k++){
		mat = Matrix(m_dim, m_dim);
		mat.InitRandomSymPosDef(0.,1.,m_seed+k);//InitIdentity();
		m_Ws.push_back(mat);
		mat.clear();
	}
	rs.SetRange(m_dim,m_dim*1.5);
	//rs.SetRange(m_dim-1,m_dim+2.);
	for(int k = 0; k < m_k; k++)
		m_nus.push_back(rs.SampleFlat());

	//init responsibility statistics
	for(int k = 0; k < m_k; k++){
		//seed N_k to 1 - multiplicative factor in Update for Rstat quantities
		m_norms.push_back(1.);
		m_Ss.push_back(Matrix(m_dim,m_dim));
		m_Ss[k].InitIdentity();
	}
	
	//seed means from kmeans
	InitParameters();
	UpdatePriorParameters();


}

//for normal EM algorithm - no priors
void GaussianMixture::CalculatePosterior(Matrix& post){
	//each kth vector is a cluster of n data points
	if(post.Dim(0) != m_n || post.Dim(1) != m_k){ cout << "GaussianMixture Error: post " << post.Dim(0) << " " << post.Dim(1) << " does not have the correct dimensions: " << m_n << " " << m_k << endl; return; }
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
			//gaus = Gaussian(m_xbars[k],m_covs[k]);
			g = m_model[k]->Prob(m_data->at(n));
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
//cout << "CALCULATE POSTERIOR - E STEP" << endl;
	//calculate necessary expectation values for E-step and ELBO
	CalculateExpectations();
	double E_mu_lam, post_val, norm;
	vector<double> post_norms;
	Matrix x_mat, x_min_m, x_min_mT;
	for(int n = 0; n < m_n; n++){
		norm = 0;
		for(int k = 0; k < m_k; k++){
			//nu_k*(x_n - m_k)T*W_k*(x_n - m_k)
			x_min_m = Matrix(m_dim,1);
			x_min_m.mult(m_means[k],-1.);
			x_mat = Matrix(m_data->at(n).Value());
			x_min_m.add(x_mat);		
			x_min_mT = Matrix(1, m_dim);
			x_min_mT.transpose(x_min_m);
			//full term
			Matrix transp_W = Matrix(1,m_dim);
			transp_W.mult(x_min_mT,m_Ws[k]);
			Matrix full = Matrix(1,1);
			full.mult(transp_W,x_min_m);
			E_mu_lam = m_dim/m_betas[k] + m_nus[k]*full.at(0,0);	
		//if(n == 0){ cout << "k: " << k << " beta: " << m_betas[k] << " nu: " << m_nus[k] << " full: " << full.at(0,0) << " x_n: " << endl; m_data->at(n).Print();}
			//gives ln(rho_nk)
			post_val = m_Epi[k] + 0.5*m_Elam[k] - (m_dim/2.)*log(2*acos(-1)) - 0.5*E_mu_lam;
		//	if(n == 0){ cout << "k: " << k << " n: " << n << " Epi: " << m_Epi[k] << " Elam: " << m_Elam[k] << " E_muLam: " << E_mu_lam << " beta: " << m_betas[k] << " nu: " << m_nus[k]  << " post: " << post << " exp(post): " << exp(post) << " W: " <<  endl;
		//	m_Ws[k].Print(); cout << "m: " << endl; m_means[k].Print(); cout << "x: " << endl; m_data->at(n).Print();
			//}	
			post_val = exp(post_val);
			norm += post_val;
			//need to normalize
			post.SetEntry(post_val, n, k);
		}
		post_norms.push_back(norm);
	}
//cout << "posterior pre-norm" << endl;	
	//m_post.Print();
	//normalize
	for(int n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			post.SetEntry(post.at(n,k)/post_norms[n],n,k);
		//uncomment here to check posterior values
		//	if(n == 0) cout << "k: " << k << " n: " << n << " post: " << m_post.at(n,k) << " norm: " << post_norms[n] << endl;
		}
	}
//cout << "posterior normed" << endl;	
//	m_post.Print();
//cout << "\n" << endl;


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
		//clear and overwrite m_xbar[k]
		m_xbars[k].clear();
		m_xbars[k].InitEmpty();
		for(int n = 0; n < m_n; n++){
			//add data pt x_n,
			Matrix x = Matrix(m_data->at(n).Value());
			//weighted by posterior value gamma(z_nk),
			x.mult(x,post.at(n,k));
			//to new mu for cluster k
			m_xbars[k].add(x);
		}
		//normalized by sum of posteriors for cluster k
		m_xbars[k].mult(m_xbars[k],1./m_norms[k]);
		m_model[k]->SetParameter("mean", m_xbars[k]); 
		//cout << "k: " << k << endl;
		//m_xbars[k].Print();
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
			x_min_mu.mult(m_xbars[k],-1.);
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
		m_model[k]->SetParameter("cov",m_covs[k]); 
//	cout << "k: " << k << endl;	
//	m_covs[k].Print();
//	cout << "new" << endl;
//	new_cov.Print();
		
	}

}


//for variational EM algorithm, assuming conjugate priors on mu, cov (normalWishart), and mixing coeffs (Dirichlet)
//updates alpha (Dirichlet), W, nu (Wishart), beta, m (normal)
//M-step
void GaussianMixture::UpdateVariationalParameters(const Matrix& post){
	CalculateRStatistics(post);
	UpdatePriorParameters();
}



void GaussianMixture::CalculateRStatistics(const Matrix& post){
	//responsibility statistics
	//this is for N_k (Bishop eq. 10.51) - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms[k] = 0;
		for(int n = 0; n < m_n; n++){
			//if(k == 0 && n == 0 && isnan(m_post.at(n,k))) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
			m_norms[k] += post.at(n,k);
		}
//cout << "k: " << k << " norm: " << m_norms[k] << endl;
	}

	//this is for x_k (eq. 10.52) - k dx1 matrices
	for(int k = 0; k < m_k; k++){
//if(k == 0) cout << "k: " << k << " old x_bar: " << endl;
//if(k == 0) m_xbars[k].Print();
		//clear and overwrite m_xbar[k]
		m_xbars[k].clear();
		m_xbars[k].InitEmpty();
		for(int n = 0; n < m_n; n++){
			Matrix x = Matrix(m_data->at(n).Value());
			//weighted by posterior value gamma(z_nk),
			x.mult(x,post.at(n,k));
			//to new mu for cluster k
			m_xbars[k].add(x);
		}
		//normalized by sum of posteriors for cluster k
		m_xbars[k].mult(m_xbars[k],1./m_norms[k]);
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
			x_min_mu.mult(m_xbars[k],-1.);
			x_min_mu.add(x_mat);
			
			//transpose x - mu
			Matrix x_min_muT;
			x_min_muT.transpose(x_min_mu);
			
			//(x_n - mu_k)*(x_n - mu_k)T
			S_k.mult(x_min_mu,x_min_muT);
			
			//weighting by posterior r_nk
			S_k.mult(S_k,post.at(n,k));
			
			//sum over n
			new_S.add(S_k);
		}	
		//normalize by N_k
		new_S.mult(new_S,1./m_norms[k]);
		//overwrites m_Ss[k]
		m_Ss[k] = new_S;
//if(k == 0) cout << "k: " << k << " new S: " << endl;
//if(k == 0) m_Ss[k].Print();
	}

}




void GaussianMixture::UpdatePriorParameters(){
//cout << "UPDATE PARAMETERS - M STEP" << endl;
	//now update variational distribution parameters
	//updating based on first step (sub-0 params)
//	m_betas.clear();
//	m_alphas.clear();
//	m_nus.clear();
	double prefactor;
	for(int k = 0; k < m_k; k++){
		//alphas - eq. 10.58 (all the same in vector)
		m_alphas[k] = m_alpha0 + m_norms[k];
		//betas - eq. 10.60
//cout << "k: " << k << " old beta: " << m_betas[k] << endl;	
		m_betas[k] = m_beta0 + m_norms[k];
//cout << "k: " << k << " new beta: " << m_betas[k] << " norm: " << m_norms[k] << endl;	
		//means - eq. 10.61
		m_means[k].clear();
		m_means[k].InitEmpty();
		//N_k*bar{x}_k
		m_means[k].mult(m_xbars[k],m_norms[k]);
		//m_meanBeta0 + N_k*bar{x}_k
		m_means[k].add(m_meanBeta0);
		//normalize by 1/beta
		m_means[k].mult(m_means[k],1./m_betas[k]);


		//Ws - eq. 10.62
//if(k == 0) cout << "k: " << k << " old W: " << endl;
//if(k == 0) m_Ws[k].Print();
		//caluclated for W inverse
		//beta0*N_k/(beta0 + N_k)*(bar{x}_k - m0)(bar{x}_k - m0)T
		m_Ws[k].clear();
		m_Ws[k].InitEmpty();		
//if(k == 0) cout << "k: " << k << " 1 - clear W: " << endl;
//if(k == 0) m_Ws[k].Print();
		//bar{x}_k - m0
		Matrix x_min_mean = Matrix(m_dim, 1);
		x_min_mean.mult(m_mean0,-1.);
		x_min_mean.add(x_min_mean,m_xbars[k]);
		Matrix x_min_meanT = Matrix(1, m_dim);
		x_min_meanT.transpose(x_min_mean);
//if(k == 0) cout << "k: " << k << " (x_n - m_0): " << endl;
//if(k == 0) x_min_mean.Print();
		//(bar{x}_k - m0)(bar{x}_k - m0)T
		m_Ws[k].mult(x_min_mean, x_min_meanT);
//if(k == 0) cout << "k: " << k << " 2 - W = (x_n - m_0)(x_n - m_0)T: " << endl;
//if(k == 0) m_Ws[k].Print();
		prefactor = m_beta0*m_norms[k]/(m_beta0 + m_norms[k]);
//if(k == 0) cout << "beta0: " << m_beta0 << " norm: " << m_norms[k] << " prefactor: " << prefactor << endl;
		m_Ws[k].mult(m_Ws[k], prefactor);
//if(k == 0)cout << "k: " << k << " 3 - W = (beta0*N_k)/(beta0 + N_k)*(x_n - m_0)(x_n - m_0)T: " << endl;
//if(k == 0)m_Ws[k].Print();
		//N_k*S_k
		Matrix scaledS = Matrix(m_dim, m_dim);
		scaledS.mult(m_Ss[k],m_norms[k]);
//if(k ==0 )cout << "k: " << k << " N_k*S_k: " << endl;
//if(k ==0 )scaledS.Print();
		//add first two terms to last term
		m_Ws[k].add(scaledS);
//if(k == 0)cout << "k: " << k << " 4 - W = (beta0*N_k)/(beta0 + N_k)*(x_n - m_0)(x_n - m_0)T + N_k*S_k: " << endl;
//if(k == 0)m_Ws[k].Print();
		//add W0inv
		m_Ws[k].add(m_W0inv);
//if(k == 0) cout << "W0inv" << endl;
//if(k == 0) m_W0inv.Print();
//if(k == 0) cout << "k: " << k << " 5 - W = (beta0*N_k)/(beta0 + N_k)*(x_n - m_0)(x_n - m_0)T + N_k*S_k + W0inv: " << endl;
//if(k == 0) m_Ws[k].Print();
		//Matrix testW;
		//testW.invert(m_Ws[k]);
		//invert (calculated for W_k inverse)
		m_Ws[k].invert(m_Ws[k]);
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
		m_nus[k] = m_nu0 + m_norms[k];
	//cout << "new nu: " << m_nus[k] << endl;

	}
};


map<string, vector<Matrix>> GaussianMixture::GetParameters(){
	map<string, vector<Matrix>> params;

	params["mean"] = m_xbars;
	params["cov"] = m_covs;
	params["pi"] = {Matrix(m_weights)};
	return params;
}

map<string, vector<Matrix>> GaussianMixture::GetPriorParameters(){
	map<string, vector<Matrix>> params;

	params["model:mu"] = m_xbars;
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
	double L, sum_k, prob;
	for(int n = 0; n < m_n; n++){
		sum_k = 0;
		for(int k = 0; k < m_k; k++){
			//matrix mu = params[k][0];
			//matrix cov = params[k][1];
			prob = m_model[k]->Prob(m_data->at(n));
			sum_k += m_weights[k]*prob;
		}
		L += log(sum_k);
	}
	return L;

}

//for use in ELBO + E-step
//posterior calculation uses these expectation values - need to calculate first
void GaussianMixture::CalculateExpectations(){
	//calculate alpha_hat
	double alpha_hat = 0.;
	//alpha_hat = sum_k alpha_k
	for(int k = 0; k < m_k; k++)
		alpha_hat += m_alphas[k];
	
	//calculate Elam (10.65) and Epi (10.66)
	double digam;
	for(int k = 0; k < m_k; k++){
 		digam = 0;
		for(int d = 1; d < m_dim+1; d++)
			digam += digamma((m_nus[k] + 1 - d)/2);
		m_Elam[k] = digam + m_dim*log(2) + log(m_Ws[k].det());
		m_Epi[k] = digamma(m_alphas[k]) - digamma(alpha_hat);
		
	}
}



double GaussianMixture::EvalVariationalLogL(const Matrix& post){
	double E_p_all, E_p_Z, E_p_pi, E_p_muLam, E_q_Z, E_q_pi, E_q_muLam, E_lam;
	//E[ln p(X|Z,mu,lam)] = 0.5*sum_k( N_k*(ln~lam_k - m_dim/beta_k - nu_k*Tr(S_k*W_k) - nu_k*(mus_k - m_k)T*W_k*(mu_k - m_k) - D*log(2*pi) ))
	E_p_all = 0;
	for(int k = 0; k < m_k; k++){
		//(x_n - m_k)
		Matrix xbar_min_m = Matrix(m_dim,1);
		xbar_min_m.mult(m_means[k],-1.);
		xbar_min_m.add(m_xbars[k]);		
		Matrix xbar_min_mT = Matrix(1, m_dim);
		xbar_min_mT.transpose(xbar_min_m);
		//(x_n - m_k)T*W_k*(x_n - m_k)
		Matrix xbarT_x_W = Matrix(1,m_dim);
		xbarT_x_W.mult(xbar_min_mT,m_Ws[k]);
		Matrix full = Matrix(1,1);
		full.mult(xbarT_x_W,xbar_min_m);
		//tr(S_k*W_k)
		//S_k*W_k = dxd matrix
		Matrix tmp_S_W = Matrix(m_dim, m_dim);
		tmp_S_W.mult(m_Ss[k],m_Ws[k]);
//	cout << "k: " << k << " norm: " << m_norms[k] << " Elam: " << m_Elam[k] << " m_dim: " << m_dim << " beta: " << m_betas[k] << " nu: " << m_nus[k] << " trace: " << tmp_S_W.trace() << " full: " << full.at(0,0) << endl;		
		E_p_all += m_norms[k]*(m_Elam[k] - m_dim/m_betas[k] - m_nus[k]*tmp_S_W.trace()  - m_nus[k]*full.at(0,0) - m_dim*log(2*acos(-1)));

	}
	E_p_all *= 0.5;

	//E[ln(p(Z|pi))] = sum_n sum_k r_nk*ln(m_Epi[k])
	E_p_Z = 0;
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++)
			E_p_Z += post.at(n,k)*m_Epi[k];
	
	//E[ln(p(pi))] = ln(C(alpha)) + (alpha_0 - 1)*sum_k m_Epi[k]
	E_p_pi = 0;
	for(int k = 0; k < m_k; k++)
		E_p_pi += m_Epi[k];
	E_p_pi *= (m_alpha0 - 1);
	vector<double> alpha0s(m_alpha0, m_k);
	Dirichlet* dir = new Dirichlet(alpha0s);
	E_p_pi += dir->lnC();


	//E[ln(p(mu,lambda))] = 1/2*sum_k(D*ln(beta0/2*pi) + m_Elam[k] - D*beta0/beta_k - beta0*nu_k(m_k - m0)T*W_kW*(m_k - m0))
	E_p_muLam = 0;
	Matrix tmp, tmpT;
	for(int k = 0; k < m_k; k++){
		//(m_k - m_0)
		Matrix mk_min_m0 = Matrix(m_dim,1);
		mk_min_m0.mult(m_mean0,-1.);
		mk_min_m0.add(m_means[k]);		
		Matrix mk_min_m0T = Matrix(1, m_dim);
		mk_min_m0T.transpose(mk_min_m0);
		//(m_k - m_0)T*W_k*(m_k - m_0)
		Matrix mT_x_W = Matrix(1,m_dim);
		mT_x_W.mult(mk_min_m0T,m_Ws[k]);
		Matrix full = Matrix(1,1);
		full.mult(mT_x_W,mk_min_m0);
		
		E_p_muLam += 0.5*(m_dim*log(m_beta0/(2*acos(-1))) + m_Elam[k] - m_dim*m_beta0/m_betas[k] - m_beta0*m_nus[k]*full.at(0,0));
		E_p_muLam += (m_nu0 - m_dim - 1)/2.*m_Elam[k];
		
		Matrix tr = Matrix(m_dim,m_dim);
		tr.mult(m_W0inv,m_Ws[k]);	
	
		E_p_muLam += -0.5*m_nus[k]*tr.trace();
	}
	Wishart* wish = new Wishart(m_W0, m_nu0);
	E_p_muLam += wish->lnB();

	//E[ln(q(Z)]
	E_q_Z = 0;
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++){
			if(post.at(n,k) == 0) E_q_Z += 0;	
			else E_q_Z += post.at(n,k)*log(post.at(n,k));	
		}
	//E[ln(q(pi))]
	E_q_pi = 0;
	for(int k = 0; k < m_k; k++){
		E_q_pi += (m_alphas[k] - 1)*m_Epi[k];
	}
	//TODO: TURN ON
	Dirichlet* dir_k = new Dirichlet(m_alphas);
	E_q_pi += dir_k->lnC();
	//E_q_pi += dir_norm;
	
	//E[ln(q(mu, lam))]
	E_q_muLam = 0;
	double H;
	for(int k = 0; k < 1; k++){
//	cout << "k: " << k << endl;
//	cout << "	Elam: " << m_Elam[k] << endl;
//	cout << "	beta: " << m_betas[k] << endl;
	//TODO: TURN ON
	Wishart* wish_k = new Wishart(m_Ws[k], m_nus[k]);
	H = wish_k->H();
//	cout << "	Wish_H: " << H << endl;
		E_q_muLam += 0.5*( m_Elam[k] + m_dim/2.*log(m_betas[k]/(2*acos(-1))) - m_dim/2. );
		E_q_muLam -= H;
	}
/*
	cout << "E_p_all: " << E_p_all << endl;
	cout << "E_p_Z: " << E_p_Z << endl;
	cout << "E_p_pi: " << E_p_pi << endl;
	cout << "E_p_muLam: " << E_p_muLam << endl;
	cout << "E_q_Z: " <<  E_q_Z << endl;
	cout << "E_q_pi: " << E_q_pi << endl;
	cout << "E_q_muLam: " <<  E_q_muLam << endl;

*/

	return E_p_all+ E_p_Z + E_p_pi+ E_p_muLam - E_q_Z - E_q_pi - E_q_muLam;
	






}
