#ifndef PointCollection_HH
#define PointCollection_HH

#include "Point.hh"
#include "RandomSample.hh"

#include <vector>
#include <iostream>
using std::vector;
using std::cout;
using std::endl;

class PointCollection{
	public:
		PointCollection() = default;
		PointCollection(Point pt){
			_pts.push_back(pt);
			_nDim = pt.Dim();
		}
		PointCollection(const vector<Point>& pts){
			_nDim = pts[0].Dim();
			for(auto p : pts){
				if(p.Dim() != _nDim){
					cout << "Error: uneven dimensions in points: " << _nDim << " vs " << p.Dim() << endl;
					break;
				}	
				_pts.push_back(p);
			}
		}

		//copy constructor
		PointCollection(const PointCollection& pts){
			_nDim = pts.Dim();
			_pts.clear();
			_pts = pts._pts;
			//for(int i = 0; i < (int)pts.GetNPoints(); i++) _pts.push_back(pts.at(i));
		}
		
		virtual ~PointCollection(){
			_pts.clear();
		}


//		PointCollection& operator =(const PointCollection& pts){
//			_pts.clear();
//			_nDim = pts.Dim();
//			for(int i = 0; i < pts.GetNPoints(); i++) _pts.push_back(pts[i]);
//			return *this;
//		}

		//add only unique points
		PointCollection& operator +=(PointCollection pts){
			//cout << "PointCollection+=PointCollection" << endl;
			for(int i = 0; i < pts.GetNPoints(); i++){
				if(!this->Has(pts.at(i))){
					_pts.push_back(pts.at(i));
				}
			}
			if(_nDim == 0) _nDim = pts.Dim();
			return *this;
		}
		
		//add only unique points
		PointCollection& operator +=(Point& pt){
			//cout << "PointCollection+=Point" << endl;
				if(!this->Has(pt)){
					_pts.push_back(pt);
				}
			if(_nDim == 0) _nDim = pt.Dim();
			return *this;
		}
		
/*
		PointCollection& operator +=(Point& pt){
//			cout << "PointCollection+=Point&" << endl;
				if(!this->Has(pt)){
//					cout << "unique point:" << endl;
//					cout << pt.Value(0) << endl;
					_pts.push_back(&pt);
				}
				else
				return *this;
			if(_nDim == 0) _nDim = pt.Dim();
			return *this;
		}
	
		Point operator[](int i) const{
			return _pts[i];
		}	

*/
	
		Point& at(int i){
			return _pts.at(i);
		}
		

		void Clear(){
			_pts.clear();
		}

		void Remove(int i){
			_pts.erase(_pts.begin()+i);	
		}	
		
		int GetNPoints() const{
			return (int)_pts.size();
		}

		int Dim() const{ return _nDim; }
		
		//check if point is in this collection
		bool Has(Point& pt){		
		//cout << "PointCollection::Has" << endl;
		//cout << "size for for loop:" << (int)_pts.size() << endl;
			for(int p = 0; p < (int)_pts.size(); p++){
				if(_pts[p] == pt){
					return true;
				}
			}
			return false;
		}
/*
		bool Has(Point* pt) const{		
		cout << "PointCollection::Has" << endl;
		cout << "size for for loop:" << (int)_pts.size() << endl;
			for(int p = 0; p < (int)_pts.size(); p++){
				bool eq = (_pts[p] == pt);
				cout << "looking at point #" << p << " match? " << eq << endl;
				if(_pts[p] == pt){
				cout << "same point: " << pt->Value(0) << " " << _pts[p]->Value(0) << endl;
					return true;
				}
			}
			return false;
		}
*/
		//add only unique points
		void AddPoint(Point& pt){
		//cout << "addpoint: " << &pt << endl;
		//pt.Print();
		if(_nDim == 0) _nDim = pt.Dim();
			if(_nDim != 0 && pt.Dim() != _nDim){
				cout << "Error: point dimension " << pt.Dim() << " does not match collection dim: " << _nDim << endl; 
				return;
			}
			if(!this->Has(pt)){
				_pts.emplace_back(pt);
			}
		//cout << "addpoint pts: " << &_pts[0] << endl;
		//_pts[0].Print();
		}

		//add only unique points
		void AddPoints(PointCollection& pts){
			if(_pts.size() != 0){
				if(pts.Dim() == _nDim){
					for(int i = 0; i < pts.GetNPoints(); i++){
						if(!this->Has(pts.at(i))){ 
						_pts.push_back(pts.at(i));}
					}
				}
				else{
					cout << "Error: PointCollection dimension given " << pts.Dim() << " does not match PointCollection dimension " << _nDim << "." << endl;
				}
			}
			else{
				for(int i = 0; i < pts.GetNPoints(); i++) _pts.push_back(pts.at(i));
				_nDim = pts.Dim();
			}
		}

		void AddPoints(vector<Point> pts){
			_nDim = pts[0].Dim();
			for(auto p : pts){
				if(p.Dim() != _nDim){
					cout << "Error: uneven dimensions in points: " << _nDim << " vs " << p.Dim() << endl;
					break;
				}	
				_pts.push_back(p);
			}
		}

		void Print(){
			for(int i = 0; i < (int)_pts.size(); i++) _pts[i].Print();

		}		




	//Randomly selects points from in array, fills out array with selected points
	PointCollection SelectPoints(int nIn, int nOut, unsigned long long seed = 123){
		RandomSample rs(seed);
		rs.SetRange(0,nIn);
		PointCollection out;
		for(int i = 0; i < nOut; i++){
			int idx = rs.SampleFlat();
			out += _pts[idx];
		}
		return out;
	
	}


	void Sort(int d, unsigned long long seed = 123){
		int N = _pts.size();
		if(N < 2) return;
		RandomSample rs(seed);	
		
		PointCollection low;
		PointCollection same;
		PointCollection high;
	
		rs.SetRange(0,N);
		int idx = rs.SampleFlat();		
		Point pivot = _pts[idx];
		if(pivot.Value().size() < 1){
			cout << "Error: no value for pivot point #" << idx << " found." << endl;
			return;
		}
		for(int i = 0; i < N; i++){
			if(_pts[i].size() < 1){
				cout << "Error: point #" << i << " has no value set." << endl;
				return;
			}
			if( pivot.ge(_pts[i], d) )
				low += _pts[i];
			else if( pivot.eq(_pts[i], d) )
				same += _pts[i];
			else
				high += _pts[i];
		}
		low.Sort(d);
		high.Sort(d);
	
		_pts.clear();
		for(int i = 0; i < low.GetNPoints(); i++) _pts.push_back(low.at(i));
		for(int i = 0; i < same.GetNPoints(); i++) _pts.push_back(same.at(i));
		for(int i = 0; i < high.GetNPoints(); i++) _pts.push_back(high.at(i));
	
	}
/*

	void FillGaus1D(vector<double> mean, vector<double> sigma, int Nsample, int dim = 1, unsigned long long seed = 123){
		if((int)mean.size() != dim || (int)sigma.size() != dim || mean.size() != sigma.size()){
			cout << "Please provide consistent parameters and dimensionality. Mean: " << mean.size() << " sigma: " << sigma.size() << " dimensions: " << dim << endl;
			return;
		}
		_pts.clear();
		
		if(_nDim != 0 && dim != _nDim){
			cout << "Error: dimensions provided " << dim << " are incompatible with existing dimensions " << _nDim << "." << endl;
		return;
		}
		if(_nDim == 0) _nDim = dim;

		RandomSample rs(seed);
		rs.SetRange(0.,1.);
		vector<vector<double>> dims(dim);
		//sample each dimension independently
		for(int d = 0; d < dim; d++){
			if(sigma[d] < 0){
				cout << "Please input valid sigma > 0" << endl;
				return;
			}
			dims[d] = rs.SampleGaussian(mean[d], sigma[d], Nsample);
		}
		for(int i = 0; i < Nsample; i++){
			VD x;
			for(int d = 0; d < dim; d++){
				x += dims[d][i];
			}
			_pts.push_back(Point(x));
			x.clear();
		
		}	
	}	

*/
	private:
		int _nDim = 0;
		vector<Point> _pts;



};
#endif
