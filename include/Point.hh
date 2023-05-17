#ifndef Point_HH
#define Point_HH

#include "Tools.hh"
#include <iostream>

using std::cout;
using std::endl;


class Point{
	public:
		Point() = default;
		
		//copy constructor
		Point(const Point &p){
			_nDim = p.Dim();
			_value.clear();
			_rank.clear();
			_value = p.Value();
			_rank = p.Rank();	
		}

		Point(const vector<double>& x){
		_nDim = (int)x.size();
			for(int i = 0; i < _nDim; i++) _value.push_back(x[i]);	
			for(int i = 0; i < _nDim; i++) _rank.push_back(-1.);	
		}
		
		
		Point& operator =(const Point& p){
			_nDim = p.Dim();
			_value.clear();
			_rank.clear();
			_value = p.Value();
			_rank = p.Rank();	
			return *this;
		}
		bool operator == (const Point& pt2) const{
			return !(*this != pt2);
		}

		
		bool operator != (const Point& pt2) const{
			if(_nDim != pt2.Dim()) return true;
			for(int i = 0; i < _nDim; i++){
				if(_value[i] != pt2.Value(i)) return true;
			}
			return false;
		}
		
		~Point(){
			_value.clear();
			_rank.clear();
		}
		

		vector<double> Value() const{return _value;}
		//return value at dimension d
		double Value(int d) const{return _value[d];}
		
		void SetValue(double v, int d){
			if(d > _value.size()){
				cout << "Error: initial values not set." << endl;
				return;
			}
			_value[d] = v;
			return;
		}
		
		void SetValue(vector<double>& v){
			if(v.size() != _nDim){
				cout << "Error: length of vector " << v.size() << " does not match Point dimension " << _nDim << endl;
				return;
			}
			_value.clear();
			_value = v;
			return;
		}
		
		vector<double> Rank() const{return _rank;}
		//return rank at dimension d
		double Rank(int d) const{return _rank[d];}
		
		int Dim() const{ return _nDim; }

		void SetRank(vector<double>& r){
			if(r.size() != _nDim){
				cout << "Error: length of vector " << r.size() << " does not match Point dimension " << _nDim << endl;
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
				cout << "Error: dimension: " << d << " inaccessible for Point of dimension: " << _nDim << endl;
				return;
			}
			_rank[d] = r;
			return;
		}


		//compare along axis
		bool eq(Point pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return(_value[d] == pt2.Value(d));
		}
		
		//compare along axis
		bool neq(Point pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return(_value[d] != pt2.Value(d));
		}
		
		//compare along axis
		bool ge(Point pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return (_value[d] > pt2.Value(d));
		}
		//compare along axis
		bool le(Point pt2, int d = 0){
			return(_value[d] < pt2.Value(d));
		}

		double size(){
			return _value.size();
		}

		void Print() const{
			std::string out = "(";
			for(int d = 0; d < _nDim; d++){ 
				if( d < _nDim-1)
					out += std::to_string(_value[d])+",";

				else
					out += std::to_string(_value[d])+")";
			}	
			cout << out << endl;

		}

	private:
		int _nDim;
		vector<double> _value;
		vector<double> _rank;











};
#endif
