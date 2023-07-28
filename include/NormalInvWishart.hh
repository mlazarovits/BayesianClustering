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

		BasePDF* mult(BasePDF* p1){ };
		BasePDF* Posterior(){ };
		void InitParameters();

		void UpdateParameters();

		double ConjugateEvidence(){ };
	private:
		//params
		Matrix m_mean;
		Matrix m_scalemat;
		double m_dof;
		double m_scale;

};
#endif
