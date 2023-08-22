#ifndef NormalWishart_HH
#define NormalWishart_HH

#include "BasePDF.hh"


class NormalWishart : public BasePDF{
	public:
		NormalWishart();
		NormalWishart(int d);
		NormalWishart(Matrix scalemat, Matrix mean, double dof, double scale);
		virtual ~NormalWishart(){ };

		double Prob(const Point& x);
		double Prob(const PointCollection& x){ return -1; }
		double Prob(const Matrix& mu, const Matrix& precision);		

		BasePDF* mult(BasePDF* p1){ return nullptr; }
		BasePDF* Posterior(){ return nullptr; }
		void InitParameters(unsigned long long seed = 123);

		void UpdateParameters();

	private:
		//params
		Matrix m_mean;
		Matrix m_scalemat;
		double m_dof;
		double m_scale;

};
#endif
