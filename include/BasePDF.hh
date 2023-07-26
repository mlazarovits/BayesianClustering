#ifndef BASEPDF_HH
#define BASEPDF_HH

#include "Point.hh"
#include "Matrix.hh"
#include <string>
using std::string;


class BasePDF{
	public:
		BasePDF(){ };
		BasePDF(int d){ m_dim = d; }
		virtual ~BasePDF(){ }	

		virtual double Prob(const Point& x) = 0;

		virtual void InitParameters() = 0;
		//returns a map from string name of parameter to vector (1 per cluster) of parameter value
		virtual map<string, vector<Matrix>> GetParameters() = 0;
//		virtual void SetParameters(map<string, Matrix> params) = 0;
		void SetParameter(string name, Matrix param){
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
		Matrix GetParameter(string name){ return m_params[name]; }
		void SetDim(int d){ m_dim = d; }
	
		int Dim(){ return m_dim; }
		int m_dim;

		map<string, Matrix> m_params;







};

#endif
