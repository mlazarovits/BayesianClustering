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
using std::isnan;
//k = # of clusters (cols)
//n = # of data pts (rows)
GaussianMixture::GaussianMixture(){ 
	m_k = 0;
	m_n = 0;
	m_dim = 0;
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
}

GaussianMixture::GaussianMixture(int k) : BasePDFMixture(k){
	for(int k = 0; k < m_k; k++){
		m_model[k] = new Gaussian();
		m_model[k]->SetPrior(new NormalWishart());
		//init data stat matrices here bc m_dim has been set from data
		_xbar.push_back(Matrix());
		_Sbar.push_back(Matrix());
		m_Elam.push_back(0.);
		m_Epi.push_back(0.);
	}
	//beta > 0
	m_beta0 = 1e-3;
	//m > 0
	m_mean0 = Matrix(m_dim,1);
	//choose m_0 = 0 by symmetry (see Bishop eq. 10.40)
	m_mean0.InitEmpty();
	//nu > d - 1 (degrees of freedom)
	m_nu0 = m_dim;// - 1) + 1e-3;
	
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
}

void GaussianMixture::InitParameters(map<string, Matrix> priors, vector<map<string, Matrix>> prev_posteriors, unsigned long long seed){
//cout << "data" << endl; m_data->Print();
//cout << "m_n " << m_n << " m_k " << m_k << endl;
	//if no previously defined posteriors, use randomly initialized kmeans 
	//to seed data statistics + responsibilities (hard assignments), then calculate posterior parameters
	m_post.SetDims(m_n, m_k);
	if(!priors.empty())
		InitPriorParameters(priors);
	for(int k = 0; k < m_k; k++){
		m_model[k]->SetDim(m_dim);
		m_model[k]->GetPrior()->SetDim(m_dim);
		_xbar[k] = Matrix(m_dim, 1);
		_Sbar[k] = Matrix(m_dim, m_dim);
	}

	//this bypasses kmeans and the first M step (step 0) and directly sets the posteriors from a previous model(s)
	if(!prev_posteriors.empty()){
		bool kmeans = false;
		int n_starting_params = prev_posteriors.size();
		for(auto post : prev_posteriors){
			if(post.empty()) n_starting_params--;
		}
		//# of clusters should be set to be the sum of the given models (ie size of prev_posteriors)
		if(m_k < n_starting_params){
			cout << "Error: # of clusters initialized is " << m_k << " but given " << n_starting_params << " models to initialize posteriors with. The number of starting clusters must be the same as the number of given posterior models. Defaulting to randomly initialized k-means" << endl;
			kmeans = true;
		}
		if(!kmeans){
			//set posterior parameters from previous model(s)
			if(_verb > 7) cout << "Initializing mixture model with parameters from " <<  n_starting_params << " previous posterior(s)" << endl;
			int skip = 0;
			for(int k = 0; k < prev_posteriors.size(); k++){
				map<string, Matrix> params = prev_posteriors[k];
				if(params.empty()){
					//cout << "skipping empty map # " << k << endl;
					skip++; continue;}
				//cout << "k " << k << " skip " << skip << " # xbars " << _xbar.size() << " # Sbars " << _Sbar.size() << " # norms " << m_norms.size() << endl;
				m_model[k-skip]->GetPrior()->SetParameter("dof", params["dof"]);
				m_model[k-skip]->GetPrior()->SetParameter("scale", params["scale"]);
				///Matrix scalemat(3,3);
				///scalemat.InitIdentity();
				///scalemat.mult(scalemat,0.1);
				//could set scalemat to inverse (properly scaled with dof and weighted with r_nk's) of empirical covariance?
				m_model[k-skip]->GetPrior()->SetParameter("scalemat", m_W0);
				//cout << "scalemat " << endl; m_W0.Print();
				m_model[k-skip]->GetPrior()->SetParameter("mean", params["m"]);
				m_model[k-skip]->SetParameter("mean", params["m"]);
				//cout << "mean" << endl; params["m"].Print();
				m_alphas[k-skip] = params["alpha"].at(0,0);	
				//cout << "alpha " << params["alpha"].at(0,0) << endl;
	
				//_xbar[k-skip] = params["m"];
				//Matrix SbarInv = params["scalemat"];
				//SbarInv.mult(SbarInv, m_nu0);
				//_Sbar[k-skip].invert(SbarInv);
				//m_norms[k-skip] = params["alpha"].at(0,0) - m_alpha0;
				
				//cout << " cov " << endl; 
				//Matrix cov = m_model[k]->GetPrior()->GetParameter("scalemat");
				//cov.mult(cov,m_model[k]->GetPrior()->GetParameter("dof").at(0,0));
				//cov.Print();

			}
//cout << "m_k " << m_k << " # starting params " << n_starting_params << endl;
			//if m_k > prev_posteriors, add extra clusters as initialized to center with var of 1, set posteriors to priors
			//if(m_k > n_starting_params)
			//	cout << "adding " << m_k - n_starting_params << " starting pts - mk " << m_k << " # starting params " << n_starting_params <<  endl;
			PointCollection initpts = m_data->SelectPoints(m_k - n_starting_params,seed);
			for(int k = n_starting_params; k < m_k; k++){
				//seed posterior + data statistics means randomly
				//also need to set data stats for initial logLH eval
				m_model[k]->GetPrior()->SetParameter("dof", m_nu0);
				m_model[k]->GetPrior()->SetParameter("scale", m_beta0);
				m_model[k]->GetPrior()->SetParameter("scalemat", m_W0);
				m_model[k]->GetPrior()->SetParameter("mean", Matrix(initpts.at(k - prev_posteriors.size())));	
				m_model[k]->SetParameter("mean", Matrix(initpts.at(k - prev_posteriors.size())));
				m_alphas[k] = 0.01*m_data->Sumw();	
			
				//if a model isn't seeded by a previous posterior in this scheme, it is considered a "ghost" subcluser to study the IR safety of the subcluster scale
				m_model[k]->SetGhost(true);
				
				//_xbar[k] = Matrix(initpts.at(k - prev_posteriors.size())); //center of system
				//_Sbar[k].InitIdentity();
				//m_norms[k] = m_alpha0;
				
				//cout << "k " << k << " alpha " << m_alphas[k] << " mean "  << endl; m_model[k]->GetPrior()->GetParameter("mean").Print(); 
				//cout << " cov " << endl; 
				//Matrix cov = m_model[k]->GetPrior()->GetParameter("scalemat");
				//cov.mult(cov,m_model[k]->GetPrior()->GetParameter("dof").at(0,0));
				//cov.Print();

			}
			int nghosts = 0;
			int nreal = 0;
			for(int k = 0; k < m_model.size(); k++){
				if(m_model[k]->IsGhost()) nghosts++;
				else nreal++;
			}
			//cout << "GaussianMixture - starting with " << nreal << " real models and " << nghosts << " ghost models and " << m_data->GetNPoints() << " points" << endl;
			return;
		}
	}
	//cout << "Initializing mixture model with randomly seeded K-means" << endl;

	//randomly initialize mean, covariance + mixing coeff.
	RandomSample randy(seed);
	randy.SetRange(0.,1.);
	double coeff_norm = 0;
	for(int k = 0; k < m_k; k++){
		m_coeffs[k] = randy.SampleFlat();
		//make sure sum_k m_coeffs[k] = 1
		coeff_norm += m_coeffs[k];
	}
	//make sure sum_k m_coeffs[k] = 1
	for(int k = 0; k < m_k; k++) m_coeffs[k] /= coeff_norm;
	

	//for(int k = 0; k < m_k; k++) cout << "k: " << k << " Nk: " << m_norms[k] << endl;
	//init means
	KMeansCluster kmc = KMeansCluster(m_data, m_k);
	//if have locations from previous models, use those to initialize kmeans
	if(!prev_posteriors.empty()){
		PointCollection start_means;
		for(auto params : prev_posteriors) start_means += params["m"].MatToPoints();
		//add extra random starting positions if more clusters than previous models are given
		if(m_k > prev_posteriors.size()){
			PointCollection initpts = m_data->SelectPoints(m_k - prev_posteriors.size(),seed);	
			start_means += initpts;
		}
		kmc.Initialize(start_means);
	}
	//else initial randomly
	else{
		kmc.Initialize(seed);
	}
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
	kmc.GetMeans(_xbar);

//cout << "initial kmeans params" << endl;
	vector<int> assigns;
	kmc.GetAssignments(assigns);
	double ptnorm = 0;
	for(int k = 0; k < m_k; k++){
		//cout << "k " << k << " xbar " << endl; _xbar[k].Print();
		m_model[k]->SetParameter("mean",_xbar[k]);
		double wtot = 0;
		m_norms[k] = 0;
		for(int n = 0; n < m_n; n++){
			if(assigns[n] == k){ 
				for(int i = 0; i < m_dim; i++){
					for(int j = 0; j < m_dim; j++){
						double disti = m_data->at(n).at(i) - _xbar[k].at(i,0);
						double distj = m_data->at(n).at(j) - _xbar[k].at(j,0);
						_Sbar[k].SetEntry(_Sbar[k].at(i,j) + m_data->at(n).w()*disti*distj,i,j);
					}
				}
				wtot += m_data->at(n).w();
				m_post.SetEntry(m_data->at(n).w(),n,k);
			}
			m_norms[k] += m_post.at(n,k);
		}
		if(wtot != 0) _Sbar[k].mult(_Sbar[k],1/wtot);
		if(m_k == m_n) _Sbar[k].InitIdentity();
		m_model[k]->SetParameter("cov",_Sbar[k]);
		//set "ghosts" from some thresh with k-means clusters
		if(m_norms[k] < 1.){
			m_model[k]->SetGhost(true);
		}
		//cout << "cluster #" << k << " norm " << m_norms[k] << endl; 
		//cout << "k " << k << " mean" << endl; _xbar[k].Print(); cout << " cov" << endl; _Sbar[k].Print(); 

	}
	UpdatePosteriorParameters(); 
	//remove any kmeans initial clusters without any assigned points
	UpdateMixture(0);
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
		for(int n = 0; n < m_n; n++){
			m_norms[k] += m_post.at(n,k);
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
		//cout << "UpdateParameters" << endl;
		//mu.Print();
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



//returns expected value of LH parameters from posterior distributions
map<string, Matrix> GaussianMixture::GetLikelihoodParameters(int k){ 
	map<string, Matrix> p;
	if(k >= m_k) return p;
	if(k < 0) return p;
	//expected value of N(mu_k | mu_k, beta_k*lam_k) = m_k
	p["mean"] = m_model[k]->GetParameter("mean");
	//expected value of Wishart(lam_k | nu_k, W_k) = nu_k*W_k
	p["cov"] = m_model[k]->GetParameter("cov");
	p["pi"] = Matrix((m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_data->Sumw()));
	return p;
};

//returns posterior values of parameters, with Gaussian model parameters set by expected values
map<string, Matrix> GaussianMixture::GetLHPosteriorParameters(int k) const{ 
//cout << "GaussianMixture::GetLHPosteriorParameters - getting params for cluster " << k << " of " << m_k << " " << m_model.size() << " total" << " dim " << m_dim << endl;
	map<string, Matrix> p;
	if(k >= m_k) return p;
	if(k < 0) return p;
	//expected value of N(mu_k | mu_k, beta_k*lam_k) = m_k
	p["mean"] = m_model[k]->GetParameter("mean");
//cout << "model #" << k << " mean has dim " << m_model[k]->GetParameter("mean").GetDims()[0] << " weight " << m_norms[k] << endl;
//m_model[k]->GetParameter("mean").Print();
	p["cov"] = m_model[k]->GetParameter("cov");
	p["pi"] = Matrix((m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_data->Sumw()));
	p["scalemat"] = m_model[k]->GetPrior()->GetParameter("scalemat");
	p["m"] = m_model[k]->GetPrior()->GetParameter("mean");
	p["scale"] = m_model[k]->GetPrior()->GetParameter("scale");
	p["dof"] = m_model[k]->GetPrior()->GetParameter("dof");
	p["alpha"] = Matrix(m_alphas[k]);
//cout << "# xbars " << _xbar.size() << " # Sbar " << _Sbar.size() << " m_Elam " << m_Elam.size() << endl;
	//include data stats for initialization
	p["xbar"] = _xbar[k];
	p["Sbar"] = _Sbar[k];
	//returns double (1./0. bool) for is/isn't a ghost subcluster
	p["ghost"] = Matrix((double)m_model[k]->IsGhost());

	return p;
};

//returns data statistics
map<string, Matrix> GaussianMixture::GetDataStatistics(int k) const{ 
	map<string, Matrix> p;
	if(k >= m_k) return p;
	if(k < 0) return p;
	p["mean"] = _xbar[k];
	p["cov"] = _Sbar[k];
	p["pi"] = Matrix(m_coeffs[k]);
	return p;
};

map<string, Matrix> GaussianMixture::GetOnlyPosteriorParameters(int k){ 
	map<string, Matrix> p;
	if(k >= m_k) return p;
	if(k < 0) return p;
	p["pi"] = Matrix((m_alpha0 + m_norms[k])/(m_k*m_alpha0 + m_data->Sumw()));
	p["scalemat"] = m_model[k]->GetPrior()->GetParameter("scalemat");
	p["mean"] = m_model[k]->GetPrior()->GetParameter("mean");
	p["scale"] = m_model[k]->GetPrior()->GetParameter("scale");
	p["dof"] = m_model[k]->GetPrior()->GetParameter("dof");
	p["alpha"] = Matrix(m_alphas[k]);
	
	return p;
};

//initializes prior parameters to default values 
void GaussianMixture::InitPriorParameters(map<string, Matrix> params){
//cout << "INIT PRIOR PARAMS - start" << endl;
	if(m_dim == 0){
		cout << "VarGaussianMixture Initialize - Error: data has not been set." << endl;
		return;
	}

	m_Elam.clear();
	m_Epi.clear();
	for(int k = 0; k < m_k; k++){
		m_Elam.push_back(0.);
		m_Epi.push_back(0.);
	}


	
	if(params.count("dof") != 0){
		m_nu0 = params["dof"].at(0,0);
	}
	if(params.count("scale") != 0){
		m_beta0 = params["scale"].at(0,0);
	}
	if(params.count("mean") != 0){
		if(params["mean"].GetDims()[0] == m_dim && params["mean"].GetDims()[1] == 1){
			m_mean0 = params["mean"];
		}
	}
	m_meanBeta0.mult(m_mean0, m_beta0);
	if(params.count("scalemat") != 0){
		if(params["scalemat"].GetDims()[0] == m_dim && params["scalemat"].GetDims()[1] == m_dim){
			m_W0 = params["scalemat"];
			m_W0inv.invert(m_W0);
		}
	}
	//init parameters from alpha0
	//assuming a Dirichlet prior on the multinomial (categorical) assignment distribution (over latent variable z - sets pis)
	for(int k = 0; k < m_k; k++){
		m_alphas[k] = m_alpha0;
	}
	//to init prior parameters without calculating Rstats from posterior
	//set to default prior parameters
	for(int k = 0; k < m_k; k++){
		m_model[k]->GetPrior()->SetParameter("dof",Matrix(m_nu0));
		m_model[k]->GetPrior()->SetParameter("scalemat",Matrix(m_W0));
		m_model[k]->GetPrior()->SetParameter("mean",Matrix(m_mean0));
		m_model[k]->GetPrior()->SetParameter("scale",Matrix(m_beta0));
	}

//cout << "INIT PRIOR PARAMS - end" << endl;
 
}






void GaussianMixture::CalculateExpectations(){
	//cout << "CALC EXPECTATIONS - start" << endl;
	//cout << "m_k " << m_k << endl;
	//calculate alpha_hat
	double alpha_hat = 0.;
	//alpha_hat = sum_k alpha_k
	for(int k = 0; k < m_k; k++)
		alpha_hat += m_alphas[k];
//cout << "did alpha hat " << alpha_hat << endl;	
	//calculate Elam (10.65) and Epi (10.66)
	double digam, dof, det;
	Matrix scalemat = Matrix(m_dim, m_dim);
	for(int k = 0; k < m_k; k++){
		scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
 		dof = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		digam = 0;
		for(int d = 1; d < m_dim+1; d++){
			digam += digamma((dof + 1 - d)/2.);
		}
		//cout << "k: " << k << " digam: " << digam << " d*ln2: " << m_dim*log(2) << " lndet: " << log(scalemat.det()) << " det: " << scalemat.det() << endl; scalemat.Print();
		//if model contains points projected to infinity, the determinant calculation will be difficult
		//and may result in a very small, yet negative number
		//if this is the case, just flip the sign for numerical calculations
		det = scalemat.det();
		if(fabs(det) < 1e-10 && det < 0){
			//cout << "flipping sign for det " << det << " from mat" << endl; scalemat.Print();
			det = -det;
		}
		m_Elam[k] = digam + m_dim*log(2) + log(det);
		m_Epi[k] = digamma(m_alphas[k]) - digamma(alpha_hat);
		//cout << "calc expectations - k: " << k << " alpha: " << m_alphas[k] << " dof: " << dof << " Elam: " << m_Elam[k] << " Epi: " << m_Epi[k] << " detW[k]: " << scalemat.det() << " W[k]: " << endl;
		//scalemat.Print();
		if(isnan(m_Elam[k])){ cout << "NAN!!!!! k: " << k << " alpha: " << m_alphas[k] << " dof: " << dof << " Elam: " << m_Elam[k] << " Epi: " << m_Epi[k] << " detW[k]: " << det << " W[k]: " << endl;
		scalemat.Print(); 
			cout << "points" << endl; m_data->Print();
		}
		if(std::isinf(m_Elam[k])){ cout << "INF!!!!! k: " << k << " alpha: " << m_alphas[k] << " dof: " << dof << " Elam: " << m_Elam[k] << " Epi: " << m_Epi[k] << " detW[k]: " << det << " W[k]: " << endl;
		scalemat.Print(); cout << "W0" << endl; m_W0.Print();}
	}	
	//cout << "CALC EXPECTATIONS - end" << endl;
	
}





//E-step - calculate posterior
//E[z_nk] = r_nk
//(10.46) ln(rho_nk) = E[ln(pi_k)] + 1/2E[ln|lambda_k|] - D/2*ln(2pi) - 1/2E[(x_n - mu_k)Tlambda_k(x_n - mu_k)]
//(10.49) r_nk = rho_nk/sum_k rho_nk
//(10.64) ln(rho_nk) = psi(alpha_k) - psi(alpha_hat) + 1/2(sum^d_i psi( (nu_k + 1 - i) /2) + d*ln2 + ln|W_k| - D/2*ln(2pi) - 1/2( D*beta_k^inv + nu_k*(x_n - m_k)T*W_k*(x_n - m_k) )
void GaussianMixture::CalculateVariationalPosterior(){
///if(m_n == 108) cout << "CALCULATE POSTERIOR - E STEP - start" << endl;
	//calculate necessary expectation values for E-step and ELBO
	//for(int k = 0; k < m_k; k++){
	//	Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
 	//	double dof = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
	//	//cout << " e-step - k: " << k << " alpha: " << m_alphas[k] << " dof: " << dof << " Elam: " << m_Elam[k] << " Epi: " << m_Epi[k] << " detW[k]: " << scalemat.det() << " W[k]: " << endl;
	//	//scalemat.Print();
	//}
	//cout << "calcpost - start to calc expectations" << endl;
	CalculateExpectations();
	//cout << "calcpost - done calc expectations" << endl;
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
			//if(n == 96 && m_n == 108){cout << "n: " << n << " k " << k << " x" << endl; x_mat.Print(); cout << "m[k]" << endl; mean.Print(); cout << "x_min_m" << endl; x_min_m.Print();}
			x_min_m.minus(x_mat,mean);		
			//if(n == 96 && m_n == 108){cout << "x - m[k]" << endl;x_min_m.Print();}
			x_min_mT = Matrix(1, m_dim);
			x_min_mT.transpose(x_min_m);
			//full term
			//if(n == 146){cout << "W[k]:" << endl; scalemat.Print();}
			Matrix transp_W = Matrix(1,m_dim);
			transp_W.mult(x_min_mT,scalemat);
			//if(n == 146){cout << "(x - m[k])T*W[k]" << endl; transp_W.Print();}
			Matrix full = Matrix(1,1);
			full.mult(transp_W,x_min_m);
			//if(n == 146){cout << "(x - m[k])T*W[k]*(x - m[k])" << endl; full.Print();}
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
			//need to normalize - done after
			//if(n == 146) cout << "n " << n << " k " << k << " pre norm post " << post << endl;
			m_post.SetEntry(post, n, k);
			//if(k == 1){ cout << std::setprecision(10) << "n: " << n << " k: " << k << " scale: " << scale << " dof: " << dof << " Elam: " << m_Elam[k] << " E_pi: " << m_Epi[k] << " E_mu_lam: " << E_mu_lam << " post: " << post << " mat post: " << m_post.at(n,k) << " full: " << full.at(0,0) << endl;}
			//if(m_n == 108 && n == 96){
			//	cout << "n " << n << " k " << k << " post " << post << " exp[pi] " << m_Epi[k] << " exp[lam] " << m_Elam[k] << " exp[(x - mu)lam(x - mu)] " << E_mu_lam << " dof " << dof << " (x_n - m_k)W_k(x_n - m_k) " << full.at(0,0) << " (x_n - m_k) " << endl; x_min_m.Print(); cout << " m_k " << endl; mean.Print(); cout << " W_k " << endl; scalemat.Print(); cout << "point #" << n << endl; x_mat.Print(); 
			//}
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
	int idx = 0;
	int idx1 = 1;
	//cout << "w->at(" << idx << ") " << m_data->at(idx).w() << endl;
	//if(m_n > 1)
	//	cout << "w->at(" << idx1 << ") " << m_data->at(idx1).w() << endl;
	double testnorm = 0;
	double testnorm1 = 0;
	for(int n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			//will lead to nan
//cout << "post_norms_adj: " << post_norms_adj[n] << " post_n_max: " << post_n_max[n] << " mu_k: " << endl; m_model[k]->GetPrior()->GetParameter("mean").Print();  
			//if(post_norms[n] == 0){ cout << "Entry at n: " << n << " k: " << k << " is " << m_post.at(n,k) << " weight - " << m_data->at(n).w() << " point  " << endl; m_data->at(n).Print(); cout << "m_k: " << endl; m_model[k]->GetPrior()->GetParameter("mean").Print(); } 
			//weight by data weight and adjusted by max ln(p_nk)
			m_post.SetEntry(m_data->at(n).w()*exp(m_post.at(n,k) - post_n_max[n])/post_norms_adj[n],n,k);
			//if(m_n == 108){
			//	cout << "Entry at n: " << n << " k: " << k << " is " << m_post.at(n,k) << " weight - " << m_data->at(n).w() << " unweighted post " << m_post.at(n,k)/m_data->at(n).w() << " point  " << endl; m_data->at(n).Print();
			//}
			
	
			//put in safeguard for computer precision for doubles (~1e\pm308)/rounding
			if(m_post.at(n,k) < 1e-308) m_post.SetEntry(0.,n,k);

			//if(n == idx){
			//	cout << "n " << n << " k " << k << " post " << m_post.at(n,k) << endl;
			//	testnorm += m_post.at(n,k);		
			//}
			//if(n == idx1){
			//	cout << "n " << n << " k " << k << " post " << m_post.at(n,k) << endl;
			//	testnorm1 += m_post.at(n,k);		
			//}
			//if(m_post.at(n,k) > 0 && isinf(1/m_post.at(n,k))){ cout << "Entry at n: " << n << " k: " << k << " is " << m_post.at(n,k) << " weight - " << m_data->at(n).w() << " point  " << endl; m_data->at(n).Print(); cout << "post_norms_adj: " << post_norms_adj[n] << " post_n_max: " << post_n_max[n] << " mu_k: " << endl; m_model[k]->GetPrior()->GetParameter("mean").Print(); } 


			//m_post.SetEntry(m_data->at(n).w()*m_post.at(n,k)/post_norms[n],n,k);
			//uncomment here to check posterior values
		//	cout << "m post = " << m_post.at(n,k) << endl;
			//if(k == 1) cout << "k: " << k << " n: " << n << " post: " << m_post.at(n,k) << " norm: " << post_norms[n] << endl;
		}
	}
	//cout << "sum k r_" << idx << ",k = " << testnorm << endl;
	//if(m_n > 1)
	//	cout << "sum k r_" << idx1 << ",k = " << testnorm1 << endl;

	//cout << "posterior normed" << endl;	
	//m_post.Print();
//cout << "\n" << endl;
//if(m_n == 108) cout << "CALCULATE POSTERIOR - E STEP - end" << endl;
};




void GaussianMixture::CalculateRStatistics(){
//cout << "Calculate RStats - start" << endl;
//cout << "points" << endl; m_data->Print();
	//responsibility statistics
	//this is for N_k (Bishop eq. 10.51) - k entries in this vector
	for(int k = 0; k < m_k; k++){
		m_norms[k] = 0;
		for(int n = 0; n < m_n; n++){
			//if(k == 0) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
			//if(n == 3) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << endl;
			
			//weighted in posterior calculation - need to unweight
			m_norms[k] += m_post.at(n,k);
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
			//cout << "x" << endl; x.Print();
			//weighted by posterior value gamma(z_nk),
			//if(k == 1) cout << "n: " << n << " k: " << k << " post: " << m_post.at(n,k) << " data weight: " << m_data->at(n).w() << endl;
			x.mult(x,m_post.at(n,k));
			//cout << "post*x" << endl; x.Print();
			//cout << "mu" << endl; mu.Print();
			//add to new mu for cluster k
			mu.add(x);
		}
		//cout << "k: " << k << " sum_n post*x" << endl; mu.Print();
		//normalized by sum of posteriors for cluster k
		mu.mult(mu,1./m_norms[k]);
		//cout << "k: " << k << " m_norm: " << m_norms[k] << endl;
//cout << "# xbars " << _xbar.size() << endl;
		//cout << "1/N[k]*(sum_n post*x)" << endl; _xbar[k].Print();
		_xbar[k] = mu;
//cout << "set new xbar" << endl;
		//m_model[k]->SetParameter("mean",mu);		
	}
	


//		cout << std::setprecision(10) << endl;	
	//this is for S_k (eq.10.53) - k dxd matrices
	for(int k = 0; k < m_k; k++){
		//create (x_n - mu)*(x_n - mu)T matrices for each data pt
		Matrix S = Matrix(m_dim,m_dim);
		Matrix mu = _xbar[k];//m_model[k]->GetParameter("mean");
		//to avoid nans, if there are no effective points in a cluster, the update equations dictate that its parameters are just the priors
		//however in calculating the r statistics, N_k = 0 can lead to nans
		if(m_norms[k] == 0){
			m_model[k]->SetParameter("cov",S);		
			continue;	
		}

		//cout << "k: " << k << " mu" << endl; mu.Print();
		double t_sig2;
		for(int n = 0; n < m_n; n++){
			//construct x - mu
			Matrix x_mat = Matrix(m_data->at(n));
		//	cout << "n: " << n << " k: " << k << endl;
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
		
			Matrix S_n = Matrix(m_dim, m_dim);	
			//(x_n - mu_k)*(x_n - mu_k)T
			S_n.mult(x_min_mu,x_min_muT);
			//cout << "(x - mu)*(x - mu)T" << endl;	
			//S_n.Print();
			//if resolution based smearing is set for any dimension of the covariance,
			//alter smearing matrix accordingly
			if(_smear){
				S_n.add(_data_cov);
			}
			//weighting by posterior r_nk
//if(m_data->at(n).at(0) > 1e70 && m_data->at(n).at(1) > 1e70) cout << "k " << k << " r_nk " << m_post.at(n,k) << endl;
			S_n.mult(S_n,m_post.at(n,k));
		//	if(n == 3) cout << "cov calc - post: " << m_post.at(n,k) << endl;
			//cout << "post*(x - mu)*(x - mu)T" << endl; S_n.Print();
			//sum over n
			//cout << "S" << endl; S.Print(); cout << "S_n" << endl; S_n.Print();
			S.add(S_n);
		}	
		//cout << "sum_n post*(x - mu)*(x - mu)T" << endl; S.Print();
		//normalize by N_k
		S.mult(S,1./m_norms[k]);
		//cout << "m_norm: " << m_norms[k] << endl;
		//cout << "k: " << k << " norm: " << m_norms[k] << " alpha: " << m_alphas[k] << " (1/N[k])*sum_n post*(x - mu)*(x - mu)T" << endl; S.Print();
		//if data smear is specified - provides lower bound on covariance -> regularization and provides nonzero covariance in single point case
		//N_k = 0 clusters do not get a smear because there are no points to smear
		//m_model[k]->SetParameter("cov",S);
		_Sbar[k] = S;
//cout << "k " << k << " norm "<< m_norms[k] << " xbar" << endl; _xbar[k].Print(); cout << " Sk " << endl; _Sbar[k].Print();
		//if(k == 1){ cout << "CalculateRStats - cov" << endl; m_model[k]->GetParameter("cov").Print();}
	}

//cout << "Calculate RStats - end" << endl;
}


//M-step
void GaussianMixture::UpdateVariationalParameters(){
	CalculateRStatistics();
	UpdatePosteriorParameters();

}

//M-step
void GaussianMixture::UpdatePosteriorParameters(){
//cout << "UPDATE PARAMETERS - M STEP - start" << endl;

	//now update variational distribution parameters
	//updating based on first step (sub-0 params)
	for(int k = 0; k < m_k; k++){
//cout << "k " << k << " norm " << m_norms[k] << endl;
		//to avoid nans, if there are no effective points in a cluster, the update equations dictate that its parameters are just the priors
		//however in calculating the r statistics, N_k = 0 can lead to nans
		//if N_k == 0, set the parameters of this cluster to initial values
		if(m_norms[k] == 0){
			m_model[k]->GetPrior()->SetParameter("scale", Matrix(m_beta0));
			m_model[k]->GetPrior()->SetParameter("dof", Matrix(m_nu0));
			m_model[k]->GetPrior()->SetParameter("mean", m_mean0);
			m_model[k]->GetPrior()->SetParameter("scalemat", m_W0);
			continue;	
		}

		//already calculated - bar{x}k and Sk
		Matrix mu = _xbar[k];//m_model[k]->GetParameter("mean");
		Matrix cov = _Sbar[k];//m_model[k]->GetParameter("cov");
		//alphas - eq. 10.58 (all the same in vector)
		m_alphas[k] = m_alpha0 + m_norms[k];
		//cout << "k: " << k << " norm: " << m_norms[k] << " alpha: " << m_alphas[k] << endl;	
	
		//betas - eq. 10.60
		double new_scale = m_beta0 + m_norms[k];
		//cout << "beta0 " << m_beta0 << " Nk " << m_norms[k] << " new_scale " << new_scale << endl;	
		m_model[k]->GetPrior()->SetParameter("scale", Matrix(new_scale));
		//means - eq. 10.61 - m_k = (Nk*bar{x}_k + beta0*m0)/beta_k
		//N_k*bar{x}_k
		Matrix new_mean = Matrix(m_dim, 1);
		new_mean.mult(mu,m_norms[k]);
//cout << "Nk " << m_norms[k] << " xbar " << endl; mu.Print(); 
//cout << "newmean = Nk*xbar" << endl; new_mean.Print();
      	//m_meanBeta0 + N_k*bar{x}_k
      	new_mean.add(m_meanBeta0);
//cout << "mean0*beta0 " << endl; m_meanBeta0.Print();		
//cout << "newmean = m_meanBeta0 + N_k*bar{x}_k" << endl; new_mean.Print();
      	//normalize by 1/beta
      	//m_k = N_k*bar{x}_k + beta0*m0
      	new_mean.mult(new_mean,1./new_scale);
//cout << "newscale " << new_scale << endl;
//cout << "newmean = (N_k*bar{x}_k + beta0*m0)/betak" << endl; new_mean.Print();
		if(!_smear){
			//add measurement error covariance - \sum_n r_nk lambda*_n
			Matrix rLamStar(m_dim, m_dim);
			Matrix rLamStar_x(m_dim, 1);
			for(int n = 0; n < m_n; n++){
				Matrix lamStar(m_dim, m_dim);
				lamStar.mult(_lamStar[n],m_post.at(n,k));
				rLamStar.add(lamStar);	

				if(_verb > 6){
					cout << "cluster #" << k << " data pt #" << n << " w " << m_data->at(n).w() << " r " << m_post.at(n,k) << " unweighted r " << m_post.at(n,k)/m_data->at(n).w() << " lamStar " << endl; _lamStar[n].Print();
				}

				Matrix lamStar_x(m_dim, 1);
				lamStar_x.mult(lamStar,Matrix(m_data->at(n)));
				//cout << "for data pt " << endl; m_data->at(n).Print(); cout << " with r_nk " << m_post.at(n,k) << " r_nk * lam*_n * x_n is " << endl; lamStar_x.Print();
				rLamStar_x.add(lamStar_x);	
			}
			//setting lam_k = nu_k*W_k (expected values of Wishart) with old parameters
			Matrix lamExp = m_model[k]->GetPrior()->GetParameter("scalemat");
//cout << "# models " << m_model.size() << endl;
//cout << "model prior mean?" << endl; m_model[k]->GetPrior()->GetParameter("mean").Print();
			lamExp.mult(lamExp,m_model[k]->GetPrior()->GetParameter("dof").at(0,0));
//cout << "lamExp (scalemat*nu)" << endl; lamExp.Print();
//cout << "sum r_nk*lam_n" << endl; rLamStar.Print();
//cout << "sum r_nk*lam_n*x_n" << endl; rLamStar_x.Print();
			//need inverse for m*_k = lam*_k^-1(lam_k*(N_k*bar{x}_k + beta0*m0) + sum_n r_nk*lam*_n*x_n)
			//			= lam*_k^-1(lam_k*m_k + sum_n r_nk*lam*_n*x_n)
			//lam*_k = (lam_k*(N_k + beta0) + sum_n r_nk*lam_n)
			//	 = (lam_k*beta_k + sum_n r_nk*lam_n)
			//with lam_k = nu_k*W_k replaced by expected value of lam_k ~ Wishart(nu_k, W_k) from *previous posterior parameters (before update)*
			//lam*_k and lam*_k^-1
			Matrix lamStar(m_dim, m_dim);
			//cout << "1 - lamStar dims " << lamStar.GetDims()[0] << " " << lamStar.GetDims()[1] << endl;
			lamStar.mult(lamExp,new_scale);
//cout << "lam_k*beta_k " << endl; lamStar.Print();
			//cout << "2 - lamStar dims " << lamStar.GetDims()[0] << " " << lamStar.GetDims()[1] << endl;
			//cout << "rLamStar dims " << rLamStar.GetDims()[0] << " " <<  rLamStar.GetDims()[1] << endl;
		
//cout << "sum_n r_nk*lam_n " << endl; rLamStar.Print();
			lamStar.add(rLamStar);

//cout << "lam_*k = lam_k*beta_k + sum_n r_nk*lam_n" << endl; lamStar.Print();
			if(_verb > 6){
				cout << "scalemat " << endl; m_model[k]->GetPrior()->GetParameter("scalemat").Print();
				cout << "newscale " << new_scale << endl;
				cout << "lamExp" << endl; lamExp.Print();
				cout << "rLamStar" << endl; rLamStar.Print();
				cout << "lamStar" << endl; lamStar.Print();
				cout << "rlamStar_x" << endl; rLamStar_x.Print();
			}
//cout << "total lamStar " << endl; lamStar.Print();
			Matrix lamStarInv(m_dim, m_dim);
			lamStarInv.invert(lamStar);
//cout << "lam*_k^-1 x lam_k " << endl;
//Matrix test(m_dim, m_dim);
//test.mult(lamStarInv,lamExp);
//test.Print();
//test.InitEmpty();
//cout << "lam*_k^-1 x sum r_nk*lam*_n*x_n " << endl;
//Matrix test = Matrix(m_dim,1);
//test.mult(lamStarInv,rLamStar_x);
//test.Print();
			if(_verb > 6){
				cout << "lamStar - 2" << endl; lamStar.Print();
				cout << "lamStarInv" << endl; lamStarInv.Print();
			}

			//lam_k*m_k
			Matrix lam_mean(m_dim, 1);
			lam_mean.mult(lamExp,new_mean);
//cout << "lam_k " << endl; lamExp.Print();
//cout << "(lam_k*(N_k*bar{x}_k + beta0*m0) " << endl; lam_mean.Print();
			//lam_k*m_k + sum_n r_nk*lam_n*x_n
//cout << "sum_n r_nk*lam*_n*x_n " << endl; rLamStar_x.Print();
			lam_mean.add(rLamStar_x);
//cout << "(lam_k*(N_k*bar{x}_k + beta0*m0) + sum_n r_nk*lam*_n*x_n " << endl; lam_mean.Print();
//cout << "lam*_k^-1 " << endl; lamStarInv.Print();

//cout << "initial new mean" << endl; new_mean.Print();
			//redefine m_k (aka new_mean)
			new_mean = Matrix(m_dim, 1);
			new_mean.mult(lamStarInv, lam_mean);	
		}
//cout << "setting posterior mean to " << endl; new_mean.Print();
		m_model[k]->GetPrior()->SetParameter("mean", new_mean);

		//cout << "k: " << k << " scale: " << m_model[k]->GetPrior()->GetParameter("scale").at(0,0) << endl;	
		//nus - eq. 10.63
		double new_dof = m_nu0 + m_norms[k];
		m_model[k]->GetPrior()->SetParameter("dof", Matrix(new_dof));
		//cout << "k: " << k << " dof: " << m_model[k]->GetPrior()->GetParameter("dof").at(0,0) << endl;	
		//cout << "nu0 " << m_nu0 << " Nk " << m_norms[k] << " new_dof " << new_dof << endl;	
		

		//Ws - eq. 10.62
		//caluclated for W inverse
		//beta0*N_k/(beta0 + N_k)*(bar{x}_k - m0)(bar{x}_k - m0)T
		Matrix new_scalemat = Matrix(m_dim, m_dim);
		//bar{x}_k - m0
		Matrix x_min_mean = Matrix(m_dim, 1);
		x_min_mean.minus(mu, m_mean0);
//cout << "xbar" << endl; mu.Print(); cout << "m_mean0" << endl; m_mean0.Print();
//cout << "xbar - m0" << endl; x_min_mean.Print();
//cout << "xbar" << endl;
//mu.Print();
		Matrix x_min_meanT = Matrix(1, m_dim);
		x_min_meanT.transpose(x_min_mean);
//cout << "(xbar - m0)T" << endl; x_min_meanT.Print();
		//(bar{x}_k - m0)(bar{x}_k - m0)T
		new_scalemat.mult(x_min_mean, x_min_meanT);
//cout << "(bar{x}_k - m0)(bar{x}_k - m0)T" << endl; new_scalemat.Print();
		double prefactor = m_beta0*m_norms[k]/(m_beta0 + m_norms[k]);
//cout << "prefactor " << prefactor << " m_beta0 " << m_beta0 << " m_norms[k] " << m_norms[k] << endl;
		new_scalemat.mult(new_scalemat, prefactor);
//cout << "(b*N)/(b + N)*xxT" << endl; new_scalemat.Print();
		//N_k*S_k
		Matrix scaledS = Matrix(m_dim, m_dim);
//cout << "Sbar " << endl; cov.Print();
		scaledS.mult(cov,m_norms[k]);
		//add first two terms to last term
		new_scalemat.add(scaledS);
//cout << "Sk" << endl; cov.Print();
//cout << "Nk*Sk" << endl; scaledS.Print();
//cout << "N*S + (b*N)/(b + N)*xxT" << endl;new_scalemat.Print();
		//add W0inv
		new_scalemat.add(m_W0inv);
//cout << "W0^-1" << endl; m_W0inv.Print();
//cout << std::setprecision(25) << "W-1 + N*S + (b*N)/(b + N)*xxT" << endl;new_scalemat.Print();
//cout << "det " << new_scalemat.det() << endl;
		//invert (calculated for W_k inverse)
		new_scalemat.invert(new_scalemat);
//cout << "k " << k << " inverted new_scalemat det " << new_scalemat.det() << " new mat" << endl; new_scalemat.Print();
		if(isnan(new_scalemat.at(0,0)) || std::isinf(new_scalemat.at(0,0))){
			cout << "W IS NAN/INF!!!!! for cluster " << k << " m_norms: " << m_norms[k] << endl;
			cout << "new scalemat " << endl; new_scalemat.Print();
			cout << "xbar" << endl; mu.Print(); cout << "S" << endl; cov.Print();
			cout << "data" << endl; m_data->Print();
		}	
		m_model[k]->GetPrior()->SetParameter("scalemat", new_scalemat);
		//cout << "k: " << k << " cov: " << endl;
		//cov.Print();
		//cout << "k: " << k << " scalemat: " << endl; 
		//m_model[k]->GetPrior()->GetParameter("scalemat").Print();	

	
		//set individual model parameters as expectation values (replacing data stat values)
		m_model[k]->SetParameter("mean",new_mean);
		Matrix new_cov = new_scalemat;
		new_cov.mult(new_cov, new_dof);
		new_cov.invert(new_cov); //new_cov^-1 = cov (sigma)
		m_model[k]->SetParameter("cov",new_cov);
		//if(_verb > 6){
		//	cout << "k: " << k << " mean: " << endl; 
		//	m_model[k]->GetParameter("mean").Print();	
		//	cout << "cov " << endl; m_model[k]->GetParameter("cov").Print();
		//}


	}
//cout << "UPDATE PARAMETERS - M STEP - end" << endl;
};






//calculates ELBO
//(10.70) ELBO = E[ln(p(X|Z,mu,lam))] + E[ln(p(Z|pi)] + E[ln(p(pi))] + E[ln(p(mu,lam))] - E[ln(q(Z))] - E[ln(q(pi))] - E[ln(q(mu,lam))]
double GaussianMixture::EvalVariationalLogL(){
//cout << "GaussianMixture::EvalVarLogL - start # clusters " << m_k << endl;
	double E_p_all, E_p_Z, E_p_pi, E_p_muLam, E_q_Z, E_q_pi, E_q_muLam, E_lam;
	//E[ln p(X|Z,mu,lam)] = 0.5*sum_k( N_k*(ln~lam_k - m_dim/beta_k - nu_k*Tr(S_k*W_k) - nu_k*(mus_k - m_k)T*W_k*(mu_k - m_k) - D*log(2*pi) ))
	E_p_all = 0;
	//recalculate expectation values E_lam + E_pi with updated parameters
	CalculateExpectations();
	//cout << "E_p_all" << endl;
	for(int k = 0; k < m_k; k++){
		double scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		double nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		Matrix mean = m_model[k]->GetPrior()->GetParameter("mean");
		Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		Matrix mu = _xbar[k];//m_model[k]->GetParameter("mean");
		Matrix cov = _Sbar[k]; //m_model[k]->GetParameter("cov");
	//cout << "k: " << k << " scale: " << scale << " dof: " << nu << " norm: " << m_norms[k] << " alpha: " << m_alphas[k] << endl;
	//cout << "m" << endl;
	//mean.Print();
	//cout << "W" << endl;
	//scalemat.Print();
	
	//cout << "mean" << endl;
	//m_model[k]->GetParameter("mean").Print();
	//cout << "cov" << endl;
	//cov.Print();
		
		//(x_n - m_k)
		//m_xbars[k].Print();
		//cout << "Sbar for cluster " << k << endl; _Sbar[k].Print();
		//cout << "cov" << endl; cov.Print();
		Matrix xbar_min_m = Matrix(m_dim,1);
		xbar_min_m.minus(mu, mean);
	//cout << "x_n - m_k" << endl;
	//xbar_min_m.Print();	
		Matrix xbar_min_mT = Matrix(1, m_dim);
		xbar_min_mT.transpose(xbar_min_m);
		//(x_n - m_k)T*W_k*(x_n - m_k)
		Matrix xbarT_x_W = Matrix(1,m_dim);
		xbarT_x_W.mult(xbar_min_mT,scalemat);
	//cout << "(x_n - m_k)T*W_k" << endl;
	//xbarT_x_W.Print();
		Matrix full = Matrix(1,1);
		full.mult(xbarT_x_W,xbar_min_m);
	//cout << "(x_n - m_k)T*W_k*(x_n - m_k) " << full.at(0,0) << endl;
		//tr(s_k*w_k)
		//S_k*W_k = dxd matrix
		Matrix tmp_S_W = Matrix(m_dim, m_dim);
		tmp_S_W.mult(cov,scalemat);
	//cout << "k " << k << " Sk*Wk" << endl;
	//tmp_S_W.Print();
	//cout << "trace " << tmp_S_W.trace() << endl;
	//cout << "k " << k << " Nk " << m_norms[k] << " Elam " << m_Elam[k] << " dim " << m_dim << " scale " << scale << " nu " << nu << " trace " << tmp_S_W.trace() << " full " << full.at(0,0) << " cov " << endl; cov.Print(); cout << "scalemat " << endl; scalemat.Print();
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
	//cout << "E[ln(p(mu,lambda))]" << endl;
	for(int k = 0; k < m_k; k++){
		//(m_k - m_0)
		double scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		double nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		Matrix mean = m_model[k]->GetPrior()->GetParameter("mean");
		Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
	//cout << "k " << k << " beta_k " << scale << " nu_k " << nu << endl;
	//cout << "m0" << endl; m_mean0.Print();
	//cout << "mk" << endl; mean.Print();
	//cout << "W" << endl; scalemat.Print();
	
		Matrix mk_min_m0 = Matrix(m_dim,1);
		mk_min_m0.minus(mean,m_mean0);		
	//cout << "(mk - m0)" << endl; mk_min_m0.Print();
		Matrix mk_min_m0T = Matrix(1, m_dim);
		mk_min_m0T.transpose(mk_min_m0);
		//(m_k - m_0)T*W_k*(m_k - m_0)
		Matrix mT_x_W = Matrix(1,m_dim);
		mT_x_W.mult(mk_min_m0T,scalemat);
	//cout << "(m_k - m0)T*Wk" << endl;
	//mT_x_W.Print();
		Matrix full = Matrix(1,1);
		full.mult(mT_x_W,mk_min_m0);
	//cout << "(m_k - m0)T*Wk*(m_k - m0) " << full.at(0,0) << endl;
		half_sum += (m_dim*log(m_beta0/(2*acos(-1))) + m_Elam[k] - m_dim*m_beta0/scale - m_beta0*nu*full.at(0,0));
		lam_sum += m_Elam[k];
		
		Matrix tr = Matrix(m_dim,m_dim);
		tr.mult(m_W0inv,scalemat);	
//cout << "k " << k << " half sum contribution " << (m_dim*log(m_beta0/(2*acos(-1))) + m_Elam[k] - m_dim*m_beta0/scale - m_beta0*nu*full.at(0,0)) << " Elam " << m_Elam[k] << " (m_dim*log(m_beta0/(2*acos(-1))) " << m_dim*log(m_beta0/(2*acos(-1))) << " m_dim*m_beta0/scale " << m_dim*m_beta0/scale << " m_beta0*nu*full.at(0,0) " << m_beta0*nu*full.at(0,0)  <<  " nu " << nu << " full " << full.at(0,0) <<  endl;	
	//cout << "m_W0inv" << endl; m_W0inv.Print();
	//cout << "m_W0inv*Wk" << endl; tr.Print();
	//cout << "trace " << tr.trace() << endl;	
		tr_sum += nu*tr.trace();
	}
//cout << "half sum " << half_sum << " lam sum "  << lam_sum << " trace sum " << tr_sum << endl;
	E_p_muLam = 0.5*half_sum + ((m_nu0 - m_dim - 1)/2.)*lam_sum - 0.5*tr_sum;
	Wishart* wish = new Wishart(m_W0, m_nu0);
	E_p_muLam += m_k*wish->lnB();
//cout << "nu0 " << m_nu0 << " W0" << endl; m_W0.Print();
//cout << "wish lnB " << wish->lnB() << endl;
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
	//cout << "E[ln(q(mu, lam))]" << endl;
	for(int k = 0; k < m_k; k++){
		Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
		double nu = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
		double scale = m_model[k]->GetPrior()->GetParameter("scale").at(0,0);
		Wishart* wish_k = new Wishart(scalemat, nu);
		H = wish_k->H();
	//cout << "k " << k << " Elam " << m_Elam[k] << " betak " << scale << " H " << H << " nu " << nu << " scalemat " << endl; scalemat.Print();
		E_q_muLam += 0.5*m_Elam[k] + m_dim/2.*log(scale/(2*acos(-1))) - m_dim/2. - H;

	}
	if(std::isnan(E_p_all)) cout << "E_p_all: " << E_p_all << endl;
	if(std::isnan(E_p_Z)) cout << "E_p_Z: " << E_p_Z << endl;
	if(std::isnan(E_p_pi)) cout << "E_p_pi: " << E_p_pi << endl;
	if(std::isnan(E_p_muLam)) cout << "E_p_muLam: " << E_p_muLam << endl;
	if(std::isnan(E_q_Z)) cout << "E_q_Z: " <<  E_q_Z << endl;
	if(std::isnan(E_q_pi)) cout << "E_q_pi: " << E_q_pi << endl;
	if(std::isnan(E_q_muLam)) cout << "E_q_muLam: " <<  E_q_muLam << endl;
	if(std::isinf(E_p_all)) cout << "E_p_all: " << E_p_all << endl;
	if(std::isinf(E_p_Z)) cout << "E_p_Z: " << E_p_Z << endl;
	if(std::isinf(E_p_pi)) cout << "E_p_pi: " << E_p_pi << endl;
	if(std::isinf(E_p_muLam)) cout << "E_p_muLam: " << E_p_muLam << endl;
	if(std::isinf(E_q_Z)) cout << "E_q_Z: " <<  E_q_Z << endl;
	if(std::isinf(E_q_pi)) cout << "E_q_pi: " << E_q_pi << endl;
	if(std::isinf(E_q_muLam)) cout << "E_q_muLam: " <<  E_q_muLam << endl;
	//for(int k = 0; k < m_k; k++){
	//	Matrix scalemat = m_model[k]->GetPrior()->GetParameter("scalemat");
 	//	double dof = m_model[k]->GetPrior()->GetParameter("dof").at(0,0);
	//	cout << " evalVarLogL - k: " << k << " alpha: " << m_alphas[k] << " dof: " << dof << " Elam: " << m_Elam[k] << " Epi: " << m_Epi[k] << " detW[k]: " << scalemat.det() << " W[k]: " << endl;
	//	scalemat.Print();
	//}
	//cout << "GaussianMixture::EvalVarLogL - end" << endl;
	//if(E_p_all+ E_p_Z + E_p_pi+ E_p_muLam - E_q_Z - E_q_pi - E_q_muLam > 0){
		//cout << "E_p_all: " << E_p_all << endl;
		//cout << "E_p_Z: " << E_p_Z << endl;
		//cout << "E_p_pi: " << E_p_pi << endl;
		//cout << "E_p_muLam: " << E_p_muLam << endl;
		//cout << "E_q_Z: " <<  E_q_Z << endl;
		//cout << "E_q_pi: " << E_q_pi << endl;
		//cout << "E_q_muLam: " <<  E_q_muLam << endl;

	//}
//cout << "GaussianMixture::EvalVarLogL - end " << endl;
	return E_p_all+ E_p_Z + E_p_pi+ E_p_muLam - E_q_Z - E_q_pi - E_q_muLam;
	


}; //end ELBO


double GaussianMixture::EvalNLLTerms(){
	double E_p_all, E_p_Z, E_p_pi, E_p_muLam;
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

//	cout << "E_p_all: " << E_p_all << endl;
//	cout << "E_p_Z: " << E_p_Z << endl;
//	cout << "E_p_pi: " << E_p_pi << endl;
//	cout << "E_p_muLam: " << E_p_muLam << endl;
//	cout << "m_post" << endl;
	return E_p_all+ E_p_Z + E_p_pi+ E_p_muLam;

}; //end NLL

double GaussianMixture::EvalEntropyTerms(){
	double E_q_Z, E_q_pi, E_q_muLam;
	//recalculate expectation values E_lam + E_pi with updated parameters
	CalculateExpectations();
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
//	cout << "E_q_Z: " <<  E_q_Z << endl;
//	cout << "E_q_pi: " << E_q_pi << endl;
//	cout << "E_q_muLam: " <<  E_q_muLam << endl;
//	cout << "m_post" << endl;
	return -E_q_Z - E_q_pi - E_q_muLam;

}; //end Entropy


