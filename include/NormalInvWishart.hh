#ifndef NormalInvWishart_HH
#define NormalInvWishart_HH

#include "BasePDF.hh"


class NormalInvWishart : public BasePDF{
	public:
		NormalInvWishart();
		NormalInvWishart(int d);
		NormalInvWishart(Matrix mean, Matrix scalemat, double dof, double scale);
		virtual ~NormalInvWishart(){ };

		double Prob(const Point& x);
		double Prob(const Matrix& mu, const Matrix& cov);		

		void InitParameters();

		void UpdateParameters();

		double ConjugateEvidence(){ };
	private:
		//params
		Matrix m_mean;
		Matrix m_scalemat;
		double m_dof;
		double m_scale;

		double multidim_gam(double x){
			double pi = acos(-1);
			double prod = 1;
			for(int i = 0; i < m_dim; i++) prod *= tgamma(x + (1 - i)/2.);
			return pow(pi,m_dim*(m_dim - 1)/4)*prod;
		}	
};
#endif
