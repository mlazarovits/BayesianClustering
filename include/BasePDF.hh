#ifndef BASEPDF_HH
#define BASEPDF_HH

#include "Point.hh"
#include "Matrix.hh"
#include <string>
#include <memory>
using std::string;

class BasePDF{
	public:
		BasePDF(){ m_prefactor = 1; _isGhost = false;}
		BasePDF(int d){ m_dim = d; m_prefactor = 1;  _isGhost = false;}
		BasePDF(const BasePDF& pdf) :
			m_prefactor(pdf.m_prefactor),
			m_dim(pdf.m_dim),
			_isGhost(pdf._isGhost)
		{
			m_params = pdf.m_params;
			if(pdf.m_prior)
				m_prior = pdf.m_prior->clone();
		}
		virtual ~BasePDF(){ }	
		virtual std::unique_ptr<BasePDF> clone() const = 0;

		virtual double Prob(const BayesPoint& x) = 0;
		virtual double Prob(const PointCollection& x) = 0;

		virtual void InitParameters(unsigned long long seed = 123) = 0;
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		virtual void GetParameters(map<string, Matrix>& p){ p.clear(); p = m_params; }
		void SetParameter(string name, const Matrix& param){
			map<string, Matrix>::iterator it = m_params.find(name);
			if(it != m_params.end()){
				it->second = param;
			}
			else{
				cout << "Parameter " << name << " not valid for this PDF." << endl;
				return;
			}
			//cout << "updated map need to update parameters" << endl;
			UpdateParameters();
		}
		//to call once the SetParameters function has been called
		virtual void UpdateParameters() = 0;
		void GetParameter(string name, Matrix& param){
			if(m_params.empty()) param = Matrix(); 
			if(m_params.find(name) != m_params.end()) param = m_params[name]; 
			else{
				cout << "Parameter " << name << " not valid for this PDF." << endl;
				param = Matrix();
			}
		}
		void BaseSetDim(int d){ m_dim = d; }
		virtual void SetDim(int d) = 0;
	
		int Dim(){ return m_dim; }
		int m_dim;
		
		void SetPrefactor(double p){ m_prefactor = p; }	

		void SetPrior(std::unique_ptr<BasePDF> p){ m_prior = std::move(p); }
		BasePDF* GetPrior(){ return m_prior.get(); }

		map<string, Matrix> m_params;

		std::unique_ptr<BasePDF> m_prior;

		double m_prefactor;

		double multidim_gam(double x){
			double pi = acos(-1);
			double prod = 1;
			for(int i = 1; i < m_dim+1; i++) prod *= tgamma(x + (1 - i)/2.); 
			return pow(pi,m_dim*(m_dim - 1)/4.)*prod;
		}	
		
		void SetGhost(bool g){ _isGhost = g; }
		bool IsGhost(){return _isGhost; }
		bool _isGhost; //use for studying IR safety of subcluster scale

};

#endif
