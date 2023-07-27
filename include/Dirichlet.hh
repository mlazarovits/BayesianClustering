#ifndef Dirichlet_HH
#define Dirichlet_HH

#include "BasePDF.hh"
#include "Matrix.hh"

class Dirichlet : public BasePDF{
	public:
		Dirichlet(){ m_dim = 0;}
		Dirichlet(vector<double> alphas){ m_alphas = alphas; m_dim = (int)m_alphas.size(); }
		virtual ~Dirichlet(){ };

		void InitParameters();
		void UpdateParameters(){ for(int i = 0; i < m_dim; i++) m_alphas.push_back(m_params["alpha"].at(i,0)); }
		
		void SetAlphas(vector<double> alphas){ m_alphas = alphas; }
		
		double Prob(const Point& x);

		double lnC();
		double ConjugateEvidence(const Point& x){ };

	private:
		//parameters
		vector<double> m_alphas;

	

};
#endif
