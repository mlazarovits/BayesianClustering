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

		void InitParameters(unsigned long long seed = 123);
		double Prob(const BayesPoint& x);
		double Prob(const PointCollection& x){ return -1.; }
		void UpdateParameters(){ m_mean = m_params["mean"]; m_scale = m_params["scale"]; m_dof = m_params["m_dof"].at(0,0); }	

		BasePDF* Posterior(){ return nullptr; }

		BasePDF* mult(BasePDF* p1){ return nullptr; }

	private:
		Matrix m_mean;
		Matrix m_scale;
		double m_dof;

};
#endif
