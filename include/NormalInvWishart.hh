#ifndef NormalInvWishart_HH
#define NormalInvWishart_HH

#include "BasePDF.hh"


class NormalInvWishart : public BasePDF{
	public:
		NormalInvWishart();
		NormalInvWishart(int d);
		NormalInvWishart(Matrix scalemat, Matrix mean, double dof, double scale);
		virtual ~NormalInvWishart(){ };

		double Prob(const Point& x);
		double Prob(const PointCollection& x){ return -1; }
		double Prob(const Matrix& mu, const Matrix& cov);		

		BasePDF* mult(BasePDF* p1){ return nullptr; }
		BasePDF* Posterior(){ return nullptr; }
		void InitParameters();

		void UpdateParameters();

		double ConjugateEvidence(const Point& x){ return -1;}
		double ConjugateEvidence(const PointCollection& x){ return -1;}
	private:
		//params
		Matrix m_mean;
		Matrix m_scalemat;
		double m_dof;
		double m_scale;

};
#endif
