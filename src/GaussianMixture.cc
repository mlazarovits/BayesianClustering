#include "GaussianMixture.hh"
#include "RandomSample.hh"

#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

//k = # of clusters
//n = # of data pts
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
		m_mus.push_back(Matrix(m_dim,1));
		m_mus[k].InitRandom(-5,5);
		m_covs.push_back(Matrix(m_k));
		m_covs[k].InitRandom(0.,100.);
		m_coeffs.push_back(randy.SampleFlat());
	}
}


//E-step
void GaussianMixture::CalculatePosterior(){
	//each kth vector is a cluster of n data points
	double val;
	Matrix gaus;
	vector<double> norms;
	//calculate norms for each data pt (summed over k)
	for(int	n = 0; n < m_n; n++){
		double norm = 0.;
		for(int k = 0; k < m_k; k++){
		//fill posterior with values according to Bishop eq. (9.23)	
			gaus.SetEntry(Gaus(m_x.at(n),m_mus[k],m_covs[k]),n,k);
			norm += m_coeffs[k]*gaus.at(n,k);
		}
		norms.push_back(norm);
	}
	//calculate posterior for each data pt n in each cluster k
	for(int	n = 0; n < m_n; n++){
		for(int k = 0; k < m_k; k++){
			val = m_coeffs[k]*gaus.at(n,k)/norms[n];
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
			norms[k] += m_post.at(n,k);

		}
	}
	//set new means 	
	double mu;
	for(int k = 0; k < m_k; k++){
		for(int d = 0; d < m_dim; d++){
			mu = 0;
			for(int n = 0; n < m_n; n++){
				mu += m_post.at(n,k)*m_x.at(n).Value(d)/norms[k];
			}
			m_mus[k].SetEntry(mu, d, 0);
			m_coeffs[k] = norms[k]/m_n;
		}
	}

	//set new cov matrix
	for(int k = 0; k < m_k; k++){
		for(int d0 = 0; d0 < m_dim; d0++){
			for(int d1 = 0; d1 < m_dim; d1++){
				for(int n = 0; n < m_n; n++){
			/* redo this - needs to be set as a matrix that is weighted by the post
					m_covs[k].SetEntry(m_post.at(n,k)*(m_x.at(n).Value(d0) - m_mus[k].at(d1,0))*(m_x.at(n).Value(d1) - m_mus[k].at(d1,0)/norms[k],d0,d1),);
			*/

				}
			}

		}
	}



}




void GaussianMixture::EvalLogL(){
//debug
/*
	double logLH = 0;
	double sum_k;
	for(int n = 0; n < m_n; n++){
		sum_k = 0;
		for(int k = 0; k < m_k; k++){
			sum_k += m_coeffs[k]*Gaus(m_x[n],m_mus[k],m_covs[k]);
		} 
		
		logLH += log(sum_k);
	}
	return logLH;
*/
}






//consider doing this for whole class or moving Gaus function to more general class
//gaussian for one data point
double GaussianMixture::Gaus(const Point& x, const Matrix& mu, const Matrix& cov){
	double det = cov.det();
	int dim = x.Dim();
	if(dim != mu.GetDims()[0]){
		cout << "Error: x and mu length do not match." << dim << " " << mu.GetDims()[0] << endl;
		return 0;
	}
	if(mu.GetDims()[0] != cov.GetDims()[0]){
		cout << "Error: covariance and mu dimensions do not match." << cov.GetDims()[0] << " " << mu.GetDims()[0] << endl;
		return 0;
	}
	if(!cov.square()){
		cout << "Error: non-square covariance matrix." << endl;
		return 0;
	}
	Matrix muT;
	muT.transpose(mu);
	double coeff = 1/(pow(det,0.5)*pow(2*acos(-1),0.5*dim));
	//should only be 1 element matrix
	double expon = muT.mult(cov.mult(mu)).at(0,0);

	return coeff*exp(-0.5*expon);
}

