#ifndef BASEPDF_HH
#define BASEPDF_HH

#include "Point.hh"
#include "Matrix.hh"
#include <string>
using std::string;


class BasePDF{
	public:
		BasePDF(){ m_prefactor = 1; }
		BasePDF(int d){ m_dim = d; m_prefactor = 1; }
		virtual ~BasePDF(){ }	

		virtual double Prob(const Point& x) = 0;
		virtual double Prob(const PointCollection& x) = 0;

		virtual void InitParameters() = 0;
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		virtual map<string, Matrix> GetParameters(){ return m_params; }
//		virtual void SetParameters(map<string, Matrix> params) = 0;
		void SetParameter(string name, const Matrix& param){
			map<string, Matrix>::iterator it = m_params.find(name);
			if(it != m_params.end()){
				it->second = param;
			}
			else{
				cout << "Parameter " << name << " not valid for this PDF." << endl;
			}
			UpdateParameters();
		}
		//to call once the SetParameters function has been called
		virtual void UpdateParameters() = 0;
		Matrix GetParameter(string name){ if(m_params.find(name) != m_params.end()) return m_params[name]; else return Matrix(); }
		void SetDim(int d){ m_dim = d; }
	
		int Dim(){ return m_dim; }
		int m_dim;

		virtual BasePDF* mult(BasePDF* p1) = 0;
	
		void SetPrefactor(double p){ m_prefactor = p; }	

		void SetPrior(BasePDF* p){ m_prior = p; }
		BasePDF* GetPrior(){ return m_prior; }
		virtual double ConjugateEvidence(const Point& x) = 0;

		map<string, Matrix> m_params;

		BasePDF* m_prior;

		double m_prefactor;



};

#endif
