#ifndef GAUSSIAN_HH
#define GAUSSIAN_HH
#include "BasePDF.hh"
#include "Matrix.hh"
#include "NormalWishart.hh"

class Gaussian : public BasePDF{
	public:
		Gaussian();
		Gaussian(int d);
		Gaussian(BayesPoint mu, Matrix cov);
		Gaussian(Matrix mu, Matrix cov);
		virtual ~Gaussian(){ };
		
		std::unique_ptr<BasePDF> clone() const override{
			return std::unique_ptr<BasePDF>(new Gaussian(*this));
		}

		void InitParameters(unsigned long long seed = 123) override;
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		void UpdateParameters() override{ m_mu = m_params["mean"]; m_cov = m_params["cov"]; }	
		double Prob(const BayesPoint& x) override;
		double Prob(const PointCollection& x) override;
		NormalWishart* Posterior();

		std::unique_ptr<Gaussian> mult(Gaussian* p1); 

		void SetDim(int d) override{
			if(m_dim != 0) return;
			BaseSetDim(d);
			m_mu = Matrix(d, 1);
			m_cov = Matrix(d,d);
		}

	private:
		Matrix m_mu;
		Matrix m_cov;

	

};
#endif
