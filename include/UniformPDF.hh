#ifndef UNIFORMPDF_HH
#define UNIFORMPDF_HH

#include "BasePDF.hh"

class UniformPDF : public BasePDF{
	public:
		UniformPDF();
		UniformPDF(int d);
		UniformPDF(double a, double b){ _a = a; _b = b; 
			m_params["a"] = Matrix(BayesPoint(_a));
			m_params["b"] = Matrix(BayesPoint(_b));
		}
		virtual ~UniformPDF(){ };
		std::unique_ptr<BasePDF> clone() const override{
			return std::unique_ptr<BasePDF>(new UniformPDF(*this));	
		}

		double Prob(double x);
		double Prob(const BayesPoint& x) override;
		double Prob(const PointCollection& x) override;

		void InitParameters(unsigned long long seed = 123) override;

		void UpdateParameters() override{ _a = m_params["a"].at(0,0); _b = m_params["b"].at(0,0);
			m_params["a"] = Matrix(BayesPoint(_a));
			m_params["b"] = Matrix(BayesPoint(_b));
		 }	

		void SetDim(int d) override{
			if(m_dim != 0) return;
			BaseSetDim(d);
		}

	private:
		double _a, _b;


};

#endif
