#include "VarEMCluster.hh"




//E-step - calculate posterior
void VarEMCluster::Estimate(){
	m_pdfmix->CalculateVariationalPosterior();

}



//M-step - update parameters
//equations were derived from maximizing the posterior calculated in the E-step
void VarEMCluster::Update(){
	m_pdfmix->UpdateVariationalParameters();


}


//want to maximize
//this is used to test for algorithm convergence
//ln[p(X | mu, sigma, pi) = sum_n( ln[sum_k(pi_k*Gaus(x_n | mu_k, sigma_k))] )
double VarEMCluster::EvalLogL(){
	double l = m_pdfmix->EvalVariationalLogL();
	//cout << "CHECK for unnecessary clusters" << endl;
	m_pdfmix->UpdateMixture(_thresh);
	m_k = m_pdfmix->GetNClusters();
	return l;
}

