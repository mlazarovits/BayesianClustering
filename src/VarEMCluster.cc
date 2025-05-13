#include "VarEMCluster.hh"




//E-step - calculate posterior
void VarEMCluster::Estimate(){
//cout << "VarEMCluster - CalculateVarPost start" << endl;
	m_pdfmix->CalculateVariationalPosterior();
//cout << "VarEMCluster - CalculateVarPost end" << endl;

}



//M-step - update parameters
//equations were derived from maximizing the posterior calculated in the E-step
void VarEMCluster::Update(){
//cout << "VarEMCluster - UpdateVarPost start" << endl;
	m_pdfmix->UpdateVariationalParameters();
//cout << "VarEMCluster - UpdateVarPost end" << endl;

}


//want to maximize
//this is used to test for algorithm convergence
//ln[p(X | mu, sigma, pi) = sum_n( ln[sum_k(pi_k*Gaus(x_n | mu_k, sigma_k))] )
double VarEMCluster::EvalLogL(){
//cout << "VarEMCluster - EvalVarLogL start" << endl;
	double l = m_pdfmix->EvalVariationalLogL();
	//cout << "CHECK for unnecessary clusters" << endl;
	//only do this in algorithm, not initial logL eval
//cout << "thresh " << _thresh << " clustering_start " << clustering_start << endl;
	if(_thresh > 0 && clustering_start){
		m_pdfmix->UpdateMixture(_thresh);
		m_k = m_pdfmix->GetNClusters();
	}
//cout << "VarEMCluster - m_k " << m_k << " " << m_pdfmix->GetNClusters() << endl;
//cout << "VarEMCluster - EvalVarLogL end" << endl;
	//cout << "varEM - # clusters " << m_pdfmix->GetNClusters() << endl;
	return l;
}

