#ifndef Pt_HH
#define Pt_HH

#include "Tools.hh"
#include <iostream>

using std::cout;
using std::endl;


class Pt{
	public:
		//Pt() = default;
		
		Pt(){
			_nDim = 0;
			_weight = 1.;
		}

		Pt(const int d){
			_nDim = d;
			for(int i = 0; i < _nDim; i++) _value.push_back(-1);	
			for(int i = 0; i < _nDim; i++) _rank.push_back(-1.);	
			_weight = 1.;
		}
		
		//copy constructor
		Pt(const Pt &p){
			_nDim = p.Dim();
			_value.clear();
			_rank.clear();
			_value = p.Value();
			_rank = p.Rank();	
			_weight = p._weight;
		}

		Pt(const vector<double>& x){
			_nDim = (int)x.size();
			for(int i = 0; i < _nDim; i++) _value.push_back(x[i]);	
			for(int i = 0; i < _nDim; i++) _rank.push_back(-1.);	
			_weight = 1.;
		}
		
		
		Pt& operator =(const Pt& p){
			_nDim = p.Dim();
			_value.clear();
			_rank.clear();
			_value = p.Value();
			_rank = p.Rank();	
			_weight = p.Weight();
			return *this;
		}
		bool operator == (const Pt& pt2) const{
			return !(*this != pt2);
		}

		
		bool operator != (const Pt& pt2) const{
			if(_nDim != pt2.Dim()) return true;
			for(int i = 0; i < _nDim; i++){
				if(_value[i] != pt2.Value(i)) return true;
			}
			return false;
		}
		
		~Pt(){
			_value.clear();
			_rank.clear();
		}
		

		vector<double> Value() const{return _value;}
		//return value at dimension d
		double Value(int d) const{return _value[d];}
		double at(int d) const{return _value[d];}
		
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
				cout << "Error: length of vector " << v.size() << " does not match Pt dimension " << _nDim << endl;
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
				cout << "Error: length of vector " << r.size() << " does not match Pt dimension " << _nDim << endl;
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
				cout << "Error: dimension: " << d << " inaccessible for Pt of dimension: " << _nDim << endl;
				return;
			}
			_rank[d] = r;
			return;
		}


		//compare along axis
		bool eq(Pt pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return(_value[d] == pt2.Value(d));
		}
		
		//compare along axis
		bool neq(Pt pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return(_value[d] != pt2.Value(d));
		}
		
		//compare along axis
		bool ge(Pt pt2, int d = 0){
		if(_value.size() < 1) cout << "Error: size of point is less than 1: " << _value.size() << endl;
			return (_value[d] > pt2.Value(d));
		}
		//compare along axis
		bool le(Pt pt2, int d = 0){
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
			cout << out << " w = " << _weight << endl;

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
		


	private:
		int _nDim;
		vector<double> _value;
		vector<double> _rank;
		double _weight;










};
#endif
