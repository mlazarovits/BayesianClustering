#ifndef Dirichlet_HH
#define Dirichlet_HH

#include "BasePDF.hh"
#include "Matrix.hh"

class Dirichlet : public BasePDF{
	public:
		Dirichlet(){ m_dim = 0;}
		Dirichlet(vector<double> alphas){ m_alphas = alphas; m_dim = (int)m_alphas.size(); }
		virtual ~Dirichlet(){ };
		std::unique_ptr<BasePDF> clone() const override{
			return std::unique_ptr<BasePDF>(new Dirichlet(*this));
		}

		void InitParameters(unsigned long long seed = 123) override;
		void UpdateParameters() override{ for(int i = 0; i < m_dim; i++) m_alphas.push_back(m_params["alpha"].at(i,0)); }
		
		void SetAlphas(vector<double> alphas){ m_alphas = alphas; }
		
		BasePDF* mult(BasePDF* p1){ return nullptr; }
		
		double Prob(const BayesPoint& x) override;
		double Prob(const PointCollection& x) override{ return -1.; }
		BasePDF* Posterior(){ return nullptr; }

		double lnC();

		void SetDim(int d) override{
			if(m_dim != 0) return;
			BaseSetDim(d);
			for(int i = 0; i < d; i++) m_alphas.push_back(0.);
		} 

	private:
		//parameters
		vector<double> m_alphas;

	

};
#endif
