#ifndef Dirichlet_HH
#define Dirichlet_HH

#include "BasePDF.hh"
#include "Matrix.hh"

class Dirichlet : public BasePDF{
	public:
		Dirichlet(){ m_dim = 0;}
		Dirichlet(vector<double> alphas){ m_alphas = alphas; m_dim = (int)m_alphas.size(); }
		virtual ~Dirichlet(){ };


		void SetParameters(vector<double> alphas){ m_alphas = alphas; m_dim = (int)m_alphas.size();}
		void GetParameters(vector<double>& alphas){ alphas.clear(); alphas = m_alphas; }	
		
		double Prob(const Point& x);

		double lnC();

	private:
		//parameters
		vector<double> m_alphas;

	

};
#endif
