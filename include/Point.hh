#ifndef BayesPoint_HH
#define BayesPoint_HH

#include "Tools.hh"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
using std::cout;
using std::endl;


class BayesPoint{
	public:
		//BayesPoint() = default;
		
		BayesPoint(){
			_nDim = 0;
			_weight = 1.;
			_skipDim = -1;
			_userIdx = -1;
		}

		BayesPoint(const int d){
			_nDim = d;
			for(int i = 0; i < _nDim; i++) _value.push_back(-1);	
			_weight = 1.;
			_skipDim = -1;
			_userIdx = -1;
		}
		
		//copy constructor
		BayesPoint(const BayesPoint &p){
			_nDim = p._nDim;
			_value.clear();
			_value = p._value;
			_weight = p._weight;
			_skipDim = p._skipDim;
			_userIdx = p._userIdx;
		}

		BayesPoint(const vector<double>& x){
			_nDim = (int)x.size();
			for(int i = 0; i < _nDim; i++) _value.push_back(x[i]);	
			_weight = 1.;
			_skipDim = -1;
			_userIdx = -1;
		}
		
		
		BayesPoint& operator =(const BayesPoint& p){
			_nDim = p._nDim;
			_value.clear();
			_value = p._value;
			_weight = p._weight;
			_skipDim = p._skipDim;
			_userIdx = p._userIdx;
			return *this;
		}
		bool operator == (const BayesPoint& pt2) const{
			return !(*this != pt2);
		}

		
		bool operator != (const BayesPoint& pt2) const{
			if(_nDim != pt2.Dim()) return true;
			if(_weight != pt2.w()) return true;
			if(_skipDim != pt2.GetSkipDim()) return true;
			if(_userIdx != pt2.GetUserIdx()) return true;
			for(int i = 0; i < _nDim; i++){
				//round within 10 places
				double val = _value[i];
				double val2 = pt2.Value(i);
			        if(i == 1){ //account for small changes in phi due to 2pi modulo
					if( fabs(val - val2) > 1e-10) return true;
				}
				else{
					if( val != val2 ) return true;
				}	
				//if(val != val2) return true;
			}
			return false;
		}
		
		~BayesPoint(){
			_value.clear();
		}
		

		void Value(vector<double>& val) const{
			val.clear();
			val = _value;
		}
		//return value at dimension d
		double Value(int d) const{
			if(d < _value.size()){
				return _value[d];
			}
			else{
				cout << "Error: attempting to retrieve value at dimension " << d << " for point of dimension " << _nDim << endl;
				return -999;
			}
		}
		double at(int d) const{
			if(d < _value.size()){
				return _value[d];
			}
			else{
				cout << "Error: attempting to retrieve value at dimension " << d << " for point of dimension " << _nDim << endl;
				return -999;
			}
		}
		
		void SetValue(double v, int d){
			if(d >= _value.size()){
				cout << "Error: initial values not set." << endl;
				return;
			}
			_value[d] = v;
			return;
		}
		
		void SetValue(vector<double>& v){
			if(v.size() != _nDim){
				cout << "Error: length of vector " << v.size() << " does not match BayesPoint dimension " << _nDim << endl;
				return;
			}
			_value.clear();
			_value = v;
			return;
		}
		
		void SetWeight(double w){ _weight = w; }
		double Weight() const{ return _weight; }	
		double w() const{ return _weight; }

		int Dim() const{ return _nDim; }

		//compare along axis
		bool eq(BayesPoint pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return(_value[d] == pt2.Value(d));
		}
		
		//compare along axis
		bool neq(BayesPoint pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return(_value[d] != pt2.Value(d));
		}
		
		//compare along axis
		bool ge(BayesPoint pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return (_value[d] > pt2.Value(d));
		}
		//compare along axis
		bool le(BayesPoint pt2, int d = 0){
			return(_value[d] < pt2.Value(d));
		}

		double size(){
			return _value.size();
		}

		void Print() const{
			std::string out = "(";
			std::stringstream stream;
			stream.str("");
			stream << std::setprecision(5) << "(";
			for(int d = 0; d < _nDim; d++){ 
					stream << _value[d];
				if( d < _nDim-1){
					stream << ",";
				}
				else
					stream << ")";
			}	
			cout << std::setprecision(5) << stream.str() << " w = " << _weight << endl;

		}


	
	void Scale(double s){
		for(int d = 0; d < _nDim; d++){
				_value[d] *= s;
			}
		}


	void Invert(){
		for(int d = 0; d < _nDim; d++){
				_value[d] = 1./_value[d];
			}
		}
		

	void Put02pi(int d){
		double pi = 4*atan(1);
		int nit = 0;
		double ogval = _value[d];
		while(!(_value[d] >= 0 && _value[d] < 2*pi)){
			//if pt is negative
			if(_value[d] < -1e-16) _value[d] = _value[d] + 2*pi;
			//if pt is geq 2*pi
			else _value[d] = _value[d] - 2*pi;
			nit++;
			if(nit > 100){
				cout << "Error: BayesPoint::Put02pi for dim " << d << " not converging -  setting to value " << _value[d] << " with original value " << ogval << endl;
				break;
			}
		}
		
	}
	void CircularTranslate(double t, int d){
		double pi = acos(-1);
		double theta = atan2( sin(_value[d] - t), cos(_value[d] - t));
		_value[d] = theta;
		//round to 0
		if(fabs(_value[d]) < 1e-16){ _value[d] = 0.;}
		if(fabs(fabs(_value[d]) - 2*acos(-1)) < 1e-10){ _value[d] = 0.;}
	}
	void Translate(double t, int d){
		_value[d] = _value[d] - t;
		//round to 0
		if(fabs(_value[d]) < 1e-16) _value[d] = 0.;
	}


	void SetSkipDim(int s){ _skipDim = s; }
	int GetSkipDim() const{ return _skipDim; }
	bool Skip() const{ return _skipDim != -1; }

	void SetUserIdx(int i){ _userIdx = i; }
	int GetUserIdx() const{return _userIdx; }

	private:
		int _nDim;
		vector<double> _value;
		double _weight;
		int _skipDim;
		int _userIdx;
	







};
#endif
