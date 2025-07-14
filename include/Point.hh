#ifndef BayesPoint_HH
#define BayesPoint_HH

#include "Tools.hh"
#include <iostream>
#include <iomanip>
#include <sstream>
using std::cout;
using std::endl;


class BayesPoint{
	public:
		//BayesPoint() = default;
		
		BayesPoint(){
			_nDim = 0;
			_weight = 1.;
		}

		BayesPoint(const int d){
			_nDim = d;
			for(int i = 0; i < _nDim; i++) _value.push_back(-1);	
			for(int i = 0; i < _nDim; i++) _rank.push_back(-1.);	
			_weight = 1.;
		}
		
		//copy constructor
		BayesPoint(const BayesPoint &p){
			_nDim = p.Dim();
			_value.clear();
			_rank.clear();
			_value = p.Value();
			_rank = p.Rank();	
			_weight = p._weight;
		}

		BayesPoint(const vector<double>& x){
			_nDim = (int)x.size();
			for(int i = 0; i < _nDim; i++) _value.push_back(x[i]);	
			for(int i = 0; i < _nDim; i++) _rank.push_back(-1.);	
			_weight = 1.;
		}
		
		
		BayesPoint& operator =(const BayesPoint& p){
			_nDim = p.Dim();
			_value.clear();
			_rank.clear();
			_value = p.Value();
			_rank = p.Rank();	
			_weight = p.Weight();
			return *this;
		}
		bool operator == (const BayesPoint& pt2) const{
			return !(*this != pt2);
		}

		
		bool operator != (const BayesPoint& pt2) const{
			if(_nDim != pt2.Dim()) return true;
			if(_weight != pt2.w()) return true;
			for(int i = 0; i < _nDim; i++){
				if(_value[i] != pt2.Value(i)) return true;
			}
			return false;
		}
		
		~BayesPoint(){
			_value.clear();
			_rank.clear();
		}
		

		vector<double> Value() const{return _value;}
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

		vector<double> Rank() const{return _rank;}
		//return rank at dimension d
		double Rank(int d) const{return _rank[d];}
		
		int Dim() const{ return _nDim; }

		void SetRank(vector<double>& r){
			if(r.size() != _nDim){
				cout << "Error: length of vector " << r.size() << " does not match BayesPoint dimension " << _nDim << endl;
				return;
			}
			_rank.clear();
			_rank = r;
			return;
		}

		void SetRank(double r, int d){
			if(d > _rank.size()){
				cout << "Error: initial values not set." << endl;
				return;
			}
			if(d > _nDim){
				cout << "Error: dimension: " << d << " inaccessible for BayesPoint of dimension: " << _nDim << endl;
				return;
			}
			_rank[d] = r;
			return;
		}


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
		double pi = acos(-1);
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


	private:
		int _nDim;
		vector<double> _value;
		vector<double> _rank;
		double _weight;










};
#endif
