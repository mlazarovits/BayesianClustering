#include "GaussianMixture.hh"
#include "RandomSample.hh"

#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;


GaussianMixture::GaussianMixture(){ 
	m_k = 0;
	m_n = 0;
}

GaussianMixture::GaussianMixture(int k){
	m_k = k;
	m_n = 0;
	m_post.SetDims(m_k, m_n);
}

//don't forget to include weights (eventually)
void GaussianMixture::Initialize(){
	//randomly initialize mean, covariance + mixing coeff.
	RandomSample randy;
	//TODO: should use data to set the range on the possible parameter values
	randy.SetRange(0.,1.);
	for(int k = 0; k < m_k; k++){
		m_mus.push_back(randy.SampleFlat());	
		m_covs.push_back(Matrix(m_k));
		m_covs[k].InitRandom(0.,100.);
		m_coeffs.push_back(randy.SampleFlat());
		}
}


//E-step
void GaussianMixture::CalculatePosterior(){
	//each kth vector is a cluster of n data points
	double val;
	vector<vector<double>> gaus;
	vector<double> norms;
	//calculate norms for each data pt (summed over k)
	for(int	n = 0; n < m_n; n++){
		double norm = 0.;
		for(int k = 0; k < m_k; k++){
		//fill posterior with values according to Bishop eq. (9.23)	
			gaus.push_back(Gaus(m_x[n],m_mu[k],m_cov[k]));
			norm += m_coeff[k]*gaus[n][k];
		}
		norms.push_back(norm);
	}
	//calculate posterior for each data pt n in each cluster k
	for(int	n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			val = m_coeff[k]*gaus[n][k]/norms[n];
			m_post.SetEntry(val,n,k);		
		}
	}
}


//M-step
//equations were derived from maximizing the posterior calculated in the E-step
void GaussianMixture::UpdateParameters(){
	//calculate normalization
	vector<double> norms;
	for(int k = 0; k < m_k; k++){
		norms.push_back(0.);
		for(int n = 0; n < m_n; n++){
			norms[k] += m_post[n][k];

		}
	}
	//set new means 	
	for(int k = 0; k < m_k; k++){
		m_mu[k] = 0;
		for(int d = 0; d < m_dim; d++){
			for(int n = 0; n < m_n; n++){
				m_mu[k][d] += m_post[n][k]*x[n][d]/norms[k];
			}
			m_coeff[k] = norms[k]/m_n;
		}
	}

	//set new cov matrix
	for(int k = 0; k < m_k; k++){
		for(int d0 = 0; d0 < m_dim; d0++){
			for(int d1 = 0; d1 < m_dim; d1++){
				for(int n = 0; n < m_n; n++){
					m_cov.SetEntry(m_post[n][k]*(x[n][d0] - m_mu[k])*(x[n][d1] - m_mu[k][d1])/norms[k],d0,d1);

				}
			}

		}
	}



}




void GaussianMixture::EvalLogLH(){
	double logLH = 0;
	double sum_k;
	for(int n = 0; n < m_n; n++){
		sum_k = 0;
		for(int k = 0; k < m_k; k++){
			sum_k += m_coeff[k]*Gaus(x[n],m_mu[k],m_cov[k]);
		} 
		
		logLH += log(sum_k);
	}
	return logLH;
}



}


//consider doing this for whole class or moving Gaus function to more general class
//gaussian for one data point
double GaussianMixture:Gaus(vector<double> x, vector<double> mu, Matrix cov){
	double det = cov.det();
	int dim = x.size();
	if(x.size() != mu.size(){
		cout << "Error: x and mu length do not match." << x.size() << " " << mu.size() << endl;
	}
	Matrix mu = Matrix();
	mu.SetDims(1,m_Xdim);
	for(int d = 0; d < dim; d++){
		mu.SetEntry(x[d] - mu,0,d);
	}
	Matrix muT = mu.transpose();
	double coeff = 1/(pow(det,0.5)*pow(2*acos(-1),0.5*dim));
	//should only be 1 element matrix
	double expon = muT.mult(cov.mult(mu)).GetEntry(0,0);

	return coeff*exp(-0.5*exp);
}

