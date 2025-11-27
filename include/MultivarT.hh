#ifndef MultivarT_HH
#define MultivarT_HH

#include "BasePDF.hh"
#include "Matrix.hh"

class MultivarT : BasePDF{
	public:
		MultivarT(){ };
		MultivarT(int d);
		MultivarT(Matrix mean, Matrix scale, double dof);
		virtual ~MultivarT(){ };
		std::unique_ptr<BasePDF> clone() const override{
			return std::unique_ptr<BasePDF>(new MultivarT(*this));
		}

		void InitParameters(unsigned long long seed = 123) override;
		double Prob(const BayesPoint& x) override;
		double Prob(const PointCollection& x) override{ return -1.; }
		void UpdateParameters() override{ m_mean = m_params["mean"]; m_scale = m_params["scale"]; m_dof = m_params["m_dof"].at(0,0); }	

		BasePDF* Posterior(){ return nullptr; }

		BasePDF* mult(BasePDF* p1){ return nullptr; }

		void SetDim(int d) override{ }

	private:
		Matrix m_mean;
		Matrix m_scale;
		double m_dof;

};
#endif
