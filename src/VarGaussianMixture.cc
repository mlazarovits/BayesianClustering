#include <boost/math/special_functions/digamma.hpp>
#include "VarGaussianMixture.hh"

using boost::math::digamma;

VarGaussianMixture::VarGaussianMixture(){
	m_k = 0;
	m_n = 0;
	m_dim = 0;
};


VarGaussianMixture::VarGaussianMixture(int k){
	m_k = k;
	m_n = 0;
	m_dim = 0;
};



void VarGaussianMixture::Initialize(unsigned long long seed){
	if(m_dim == 0){
		cout << "VarGaussianMixture Initialize - Error: data has not been set." << endl;
		return;
	}


	//alpha > 0
	for(int k = 0; k < m_k; k++){
		m_alphas0.push_back(1);
	}
	
	RandomSample rs;
	rs.SetRange(0.,20.);
	//beta > 0
	for(int k = 0; k < m_k; k++){
		m_betas0.push_back(rs.SampleFlat());
	}
	//m > 0
	for(int k = 0; k < m_k; k++){
		m_means0.push_back(Matrix(m_dim,1));
		m_means0[k].InitRandom();
		m_meanBetas0.push_back(Matrix(m_dim, 1));
		m_meanBetas0[k].mult(m_means0[k], m_betas0[k]);
	}
	//W <- R in dxd space - is a covariance matrix
	for(int k = 0; k < m_k; k++){
		m_Ws0.push_back(Matrix(m_dim, m_dim));
		m_Ws0[k].InitRandomSymPosDef();
		//need W0 inverse for parameter calculation
		m_Ws0inv.push_back(Matrix(m_dim,m_dim));
		m_Ws0inv.invert(m_Ws0);
	}

	//nu > d - 1 (degrees of freedom)
	rs.SetRange(m_dim - 1, 20.);
	for(int k = 0; k < m_k; k++){
		m_nus0.push_back(rs.SampleFlat());
	}
	m_post.SetDims(m_n, m_k);
	

	m_Elam.clear();
	m_Epi.clear();
	for(int k = 0; k < m_k; k++){
		m_Elam.push_back(0.);
		m_Epi.push_back(0.);
	}

	//init parameters with sub0 values
	for(int k = 0; k < m_k; k++){
		m_alphas.push_back(m_alphas0[k]);
		m_betas.push_back(m_betas0[k]);
		m_means.push_back(m_means0[k]);
		m_Ws.push_back(m_Ws0[k]);
		m_nus.push_back(m_nus0[k]);
	}

	//init responsibility statistics
	for(int k = 0; k < m_k; k++){
		m_norms.push_back(0.);
		m_mus.push_back(Matrix(m_dim,1));
		m_covs.push_back(Matrix(m_dim,m_dim));
	}


};


//for use in ELBO + E-step
//posterior calculation uses these expectation values - need to calculate first
void VarGaussianMixture::CalculateExpectations(){
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




//E-step
//E[z_nk] = r_nk
//(10.46) ln(rho_nk) = E[ln(pi_k)] + 1/2E[ln|lambda_k|] - D/2*ln(2pi) - 1/2E[(x_n - mu_k)Tlambda_k(x_n - mu_k)]
//(10.49) r_nk = rho_nk/sum_k rho_nk
//(10.64) ln(rho_nk) = psi(alpha_k) - psi(alpha_hat) + 1/2(sum^d_i psi( (nu_k + 1 - i) /2) + d*ln2 + ln|W_k| - D/2*ln(2pi) - 1/2( D*beta_k^inv + nu_k*(x_n - m_k)T*W_k*(x_n - m_k) )
void VarGaussianMixture::CalculatePosterior(){
	cout << "CALCULATE POSTERIOR - E STEP" << endl;
	//calculate necessary expectation values for E-step and ELBO
	CalculateExpectations();
cout << "nu:" << endl;
for(int k = 0; k < m_k; k++) cout << m_nus[k] << endl;

cout << "beta:" << endl;
for(int k = 0; k < m_k; k++) cout << m_betas[k] << endl;
	double E_mu_lam, post, norm;
	vector<double> post_norms;
	Matrix x_mat, x_min_m, x_min_mT;
	for(int n = 0; n < m_n; n++){
		norm = 0;
		for(int k = 0; k < m_k; k++){
			//nu_k*(x_n - m_k)T*W_k*(x_n - m_k)
			x_min_m = Matrix(m_dim,1);
			x_min_m.mult(m_means[k],-1.);
			x_mat = Matrix(m_x.at(n).Value());
			x_min_m.add(x_mat);		
			x_min_mT = Matrix(1, m_dim);
			x_min_mT.transpose(x_min_m);
			//full term
			Matrix transp_W = Matrix(1,m_dim);
			transp_W.mult(x_min_mT,m_Ws[k]);
			Matrix full = Matrix(1,1);
			full.mult(transp_W,x_min_m);
			E_mu_lam = m_dim/m_betas[k] + m_nus[k]*full.at(0,0);	

			//gives ln(rho_nk)
			post = m_Epi[k] + 0.5*m_Elam[k] - (m_dim/2.)*log(2*acos(-1)) - 0.5*E_mu_lam;
			post = exp(post);
		
			norm += post;
			//need to normalize
			m_post.SetEntry(post, n, k);
		}
		post_norms.push_back(norm);
	}
cout << "posterior pre-norm" << endl;	
	m_post.Print();
	//normalize
	for(int n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			m_post.SetEntry(m_post.at(n,k)/post_norms[n],n,k);
		}
	}
cout << "posterior normed" << endl;	
	m_post.Print();
cout << "\n" << endl;
};



//M-step
void VarGaussianMixture::UpdateParameters(){
cout << "UPDATE PARAMETERS - M STEP" << endl;
	//responsibility statistics
	//this is for N_k (Bishop eq. 10.51) - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms[k] = 0;
		for(int n = 0; n < m_n; n++){
			m_norms[k] += m_post.at(n,k);
		}
	}

	//this is for x_k (eq. 10.52) - k dx1 matrices
	for(int k = 0; k < m_k; k++){
		//clear and overwrite m_mu[k]
		m_mus[k].clear();
		m_mus[k].InitEmpty();
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
	}
	
	//this is for S_k (eq.10.53) - k dxd matrices
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
	}


	//now update variational distribution parameters
	//updating based on first step (sub-0 params)
//	m_betas.clear();
//	m_alphas.clear();
//	m_nus.clear();
	double prefactor;
	for(int k = 0; k < m_k; k++){
		//betas - eq. 10.60
		m_betas[k] = m_betas0[k] + m_norms[k];
	
		//alphas - eq. 10.58 (all the same in vector)
		m_alphas[k] = m_alpha0[k] + m_norms[k];
	cout << "old nu: " << m_nus[k] << endl;
		cout << "N_k: " << m_norms[k] << endl;
		//nus - eq. 10.63
		m_nus[k] = m_nus0[k] + m_norms[k];
	cout << "new nu: " << m_nus[k] << endl;

		//means - eq. 10.61
		m_means[k].clear();
		m_means[k].InitEmpty();
		//N_k*bar{x}_k
		m_means[k].mult(m_mus[k],m_norms[k]);
		//m_meanBeta0 + N_k*bar{x}_k
		m_means[k].add(m_meanBetas0[k]);
		//normalize by 1/beta
		m_means[k].mult(m_means[k],1./m_betas[k]);

		//Ws - eq. 10.62
		//caluclated for W inverse
		//beta0*N_k/(beta0 + N_k)*(bar{x}_k - m0)(bar{x}_k - m0)T
		m_Ws[k].clear();
		m_Ws[k].InitEmpty();		
		//bar{x}_k - m0
		Matrix x_min_mean = Matrix(m_dim, 1);
		x_min_mean.mult(m_means0[k],-1.);
		x_min_mean.add(x_min_mean,m_mus[k]);
		Matrix x_min_meanT = Matrix(1, m_dim);
		x_min_meanT.transpose(x_min_mean);
		//(bar{x}_k - m0)(bar{x}_k - m0)T
		m_Ws[k].mult(x_min_mean, x_min_meanT);
		prefactor = m_betas0[k]*m_norms[k]/(m_betas0[k] + m_norms[k]);
		m_Ws[k].mult(m_Ws[k], prefactor);
		//N_k*S_k
		Matrix scaledS = Matrix(m_dim, m_dim);
		scaledS.mult(m_covs[k],m_norms[k]);
		//add first two terms to last term
		m_Ws[k].add(scaledS);
		//add W0inv
		m_Ws[k].add(m_Ws0inv[k]);
		//invert (calculated for W_k inverse)
		m_Ws[k].invert(m_Ws[k]);

	}
cout << "AFTER M-STEP" << endl;
cout << "nu:" << endl;
for(int k = 0; k < m_k; k++) cout << m_nus[k] << endl;

cout << "beta:" << endl;
for(int k = 0; k < m_k; k++) cout << m_betas[k] << endl;
cout << "\n" << endl;
};



//calculates ELBO
//(10.70) ELBO = E[ln(p(X|Z,mu,lam))] + E[ln(p(Z|pi)] + E[ln(p(pi))] + E[ln(p(mu,lam))] - E[ln(q(Z))] - E[ln(q(pi))] - E[ln(q(mu,lam))]
double VarGaussianMixture::EvalLogL(){
	double E_p_all, E_p_Z, E_p_pi, E_p_muLam, E_q_Z, E_q_pi, E_q_muLam, E_lam;
	//E[ln p(X|Z,mu,lam)] = 0.5*sum_k( N_k*(ln~lam_k - m_dim/beta_k - nu_k*Tr(S_k*W_k) - nu_k*(mus_k - m_k)T*W_k*(mu_k - m_k) - D*log(2*pi) ))
	E_p_all = 0;
	for(int k = 0; k < m_k; k++){
		//(x_n - m_k)
		Matrix mu_min_m = Matrix(m_dim,1);
		mu_min_m.mult(m_means[k],-1.);
		mu_min_m.add(m_mus[k]);		
		Matrix mu_min_mT = Matrix(1, m_dim);
		mu_min_mT.transpose(mu_min_m);
		//(x_n - m_k)T*W_k*(x_n - m_k)
		Matrix muT_x_W = Matrix(1,m_dim);
		muT_x_W.mult(mu_min_mT,m_Ws[k]);
		Matrix full = Matrix(1,1);
		full.mult(muT_x_W,mu_min_m);
		//tr(S_k*W_k)
		//S_k*W_k = dxd matrix
		Matrix tmp_S_W = Matrix(m_dim, m_dim);
		tmp_S_W.mult(m_covs[k],m_Ws[k]);
			
		E_p_all += m_norms[k]*(m_Elam[k] - m_dim/m_betas[k] - m_nus[k]*tmp_S_W.trace()  - m_nus[k]*full.at(0,0) - m_dim*log(2*acos(-1)));

	}
	E_p_all *= 0.5;

	//E[ln(p(Z|pi))] = sum_n sum_k r_nk*ln(m_Epi[k])
	E_p_Z = 0;
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++)
			E_p_Z += m_post.at(n,k)*m_Epi[k];
	
	//E[ln(p(pi))] = ln(C(alpha)) + (alpha_0 - 1)*sum_k m_Epi[k]
	E_p_pi = 0;
	for(int k = 0; k < m_k; k++)
		E_p_pi += m_Epi[k];
	E_p_pi *= (m_alpha0[0] - 1);
	E_p_pi += log( Dir_C(m_alpha0) );


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
	E_p_muLam += log(Wish_B(m_W0,m_nu0));

	//E[ln(q(Z)]
	E_q_Z = 0;
	for(int n = 0; n < m_n; n++)
		for(int k = 0; k < m_k; k++)
			E_q_Z += m_post.at(n,k)*log(m_post.at(n,k));	

	//E[ln(q(pi))]
	E_q_pi = 0;
	for(int k = 0; k < m_k; k++){
		E_q_pi += (m_alphas[k] - 1)*m_Epi[k];
	}
	E_q_pi += log(Dir_C(m_alphas));
	
	//E[ln(q(mu, lam))]
	E_q_muLam = 0;
	double H;
	for(int k = 0; k < 1; k++){
//	cout << "k: " << k << endl;
//	cout << "	Elam: " << m_Elam[k] << endl;
//	cout << "	beta: " << m_betas[k] << endl;
	H = Wish_H(m_Ws[k], m_nus[k]);
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
	


}; //end ELBO




