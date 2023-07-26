#include "EMCluster.hh"
#include "RandomSample.hh"
#include "Gaussian.hh"
#include "BasePDFMixture.hh"

#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

//k = # of clusters (cols)
//n = # of data pts (rows)



//E-step - calculate posterior
void EMCluster::Estimate(){
	if(m_dim == 0){
		cout << "EMCluster Initialize - Error: data has not been set." << endl;
		return;
	}
	m_pdfmix->CalculatePosterior();
}



//M-step - update parameters
//equations were derived from maximizing the posterior calculated in the E-step
void EMCluster::Update(){
	m_pdfmix->UpdateParameters();


}


//want to maximize
//this is used to test for algorithm convergence
//ln[p(X | mu, sigma, pi) = sum_n( ln[sum_k(pi_k*Gaus(x_n | mu_k, sigma_k))] )
double EMCluster::EvalLogL(){
	return m_pdfmix->EvalLogL();

}

double EMCluster::EvalVariationalLogL(const Matrix& post){
	return m_pdfmix->EvalVariationalLogL();
}
