#include "GaussianMixture.hh"
#include "RandomSample.hh"

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
	randy.SetRange(-1,1);
	for(int k = 0; k < m_k; k++){
		m_mus.push_back(randy.SampleFlat());	
		m_covs.push_back(Matrix(m_k));
		m_covs[k].InitRandom();
	}
}


//E-step
void GaussianMixture::CalculatePosterior(){
	//each kth vector is a cluster of n data points
	for(int k = 0; k < m_k; k++){
		for(int	n = 0; n < m_n; n++){
		//fill posterior with values according to Bishop eq. (9.23)	
		}
	}


}

//M-step
//equations were derived from maximizing the posterior calculated in the E-step
void GaussianMixture::UpdateParameters(){


}
