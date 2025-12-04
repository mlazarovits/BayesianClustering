#ifndef PointCollection_HH
#define PointCollection_HH

#include "Point.hh"
#include "RandomSample.hh"
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
using std::vector;
using std::cout;
using std::endl;


class PointCollection{
	public:
		PointCollection() = default;
		PointCollection(const BayesPoint& pt){
			_pts.push_back(pt);
			_nDim = pt.Dim();
			for(int d = 0; d < _nDim; d++)
				_infs.push_back((1+(d)/10.)*1e70);
			if(_skipDim == -1) _skipDim = pt.GetSkipDim();
			AllSkip();
		}
		PointCollection(const vector<BayesPoint>& pts){
			_nDim = pts[0].Dim();
			for(auto p : pts){
				if(p.Dim() != _nDim){
					cout << "Error: uneven dimensions in points: " << _nDim << " vs " << p.Dim() << endl;
					break;
				}	
				_pts.push_back(p);
			for(int d = 0; d < _nDim; d++)
				_infs.push_back((1+(d)/10.)*1e70);
			}
			if(_skipDim == -1) _skipDim = pts[0].GetSkipDim();
			AllSkip();
		}

		//copy constructor
		PointCollection(const PointCollection& pts){
			_nDim = pts.Dim();
			for(int d = 0; d < _nDim; d++)
				_infs.push_back((1+(d)/10.)*1e70);
			_pts.clear();
			//_pts = pts._pts;
			//TODO: make sure copy is decoupled from original pts (vector::assign?)
			for(int i = 0; i < (int)pts.GetNPoints(); i++) _pts.push_back(pts.at(i));
			_skipDim = pts._skipDim;
			AllSkip();
		}
		
		virtual ~PointCollection(){
			_pts.clear();
			_infs.clear();
		}


		//add only unique points
		PointCollection& operator +=(const PointCollection& pts){
			//cout << "PointCollection+=PointCollection" << endl;
			for(int i = 0; i < pts.GetNPoints(); i++){
				if(!this->Has(pts.at(i))){
					_pts.push_back(pts.at(i));
				}
				if(!pts.at(i).Skip()) _allskip = false; 
			}
			if(_nDim == 0) _nDim = pts.Dim();
			//add infs if not already there
			if(_infs.size() == 0){
				for(int d = 0; d < _nDim; d++)
					_infs.push_back((1+(d)/10.)*1e70);
			}
			if(_skipDim == -1) _skipDim = pts.at(0).GetSkipDim();
			return *this;
		}
		
		//add only unique points
		PointCollection& operator +=(const BayesPoint& pt){
			//cout << "PointCollection+=Point" << endl;
			if(!this->Has(pt)){
				_pts.push_back(pt);
			}
			if(!pt.Skip()) _allskip = false; 
			if(_nDim == 0) _nDim = pt.Dim();
			if(_infs.size() == 0){
				for(int d = 0; d < _nDim; d++)
					_infs.push_back((1+(d)/10.)*1e70);
			}
			if(_skipDim == -1) _skipDim = pt.GetSkipDim();
			return *this;
		}
		
		const BayesPoint& at(int i) const{
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
		bool Has(const BayesPoint& pt) const{		
			for(int p = 0; p < (int)_pts.size(); p++){
				if(_pts[p] == pt){
					return true;
				}
			}
			return false;
		}
		
		//add only unique points
		void AddPoint(const BayesPoint& pt){
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
			if(_infs.size() == 0){
				for(int d = 0; d < _nDim; d++)
					_infs.push_back((1+(d)/10.)*1e70);
			}
			if(_skipDim == -1) _skipDim = pt.GetSkipDim();
			if(!pt.Skip()) _allskip = false; 
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
						if(!pts.at(i).Skip()) _allskip = false; 
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
			if(_infs.size() == 0){
				for(int d = 0; d < _nDim; d++)
					_infs.push_back((1+(d)/10.)*1e70);
			}
			if(_skipDim == -1) _skipDim = pts.at(0).GetSkipDim();
		}

		void AddPoints(const vector<BayesPoint>& pts){
			_nDim = pts[0].Dim();
			for(auto p : pts){
				if(p.Dim() != _nDim){
					cout << "Error: uneven dimensions in points: " << _nDim << " vs " << p.Dim() << endl;
					break;
				}	
				_pts.push_back(p);
				if(p.Skip()) _allskip = false; 
			}
			if(_infs.size() == 0){
				for(int d = 0; d < _nDim; d++)
					_infs.push_back((1+(d)/10.)*1e70);
			}
			if(_skipDim == -1) _skipDim = pts[0].GetSkipDim();
		}

		void Print() const{
			for(int i = 0; i < (int)_pts.size(); i++) _pts[i].Print();

		}		




	//Randomly selects points 
	PointCollection SelectPoints(int nOut, unsigned long long seed = 123){
		RandomSample rs(seed);
		int nIn = (int)_pts.size(); 
		if(nOut > nIn) return *this;
		rs.SetRange(0,nIn);
		PointCollection out;
		while(out.GetNPoints() < nOut){
			int idx = rs.SampleFlat();
			out += _pts[idx];
		}
		return out;
	
	}


	void Sort(int d, unsigned long long seed = 123){
		if(d >= _nDim){
			cout << "Error: dimension " << d << " not valid for PointCollection with max dimension " << _nDim << " " << this->Dim() << endl;
			return;
		}
		int N = _pts.size();
		if(N < 2) return;
		RandomSample rs(seed);	
		
		PointCollection low;
		PointCollection same;
		PointCollection high;
		rs.SetRange(0,N);
		int idx = rs.SampleFlat();		
		BayesPoint pivot = _pts[idx];
		if(pivot.Dim() < 1){
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
		if(low.GetNPoints() > 1) low.Sort(d, seed);
		if(high.GetNPoints() > 1) high.Sort(d, seed);
	
		_pts.clear();
		for(int i = 0; i < low.GetNPoints(); i++) _pts.push_back(low.at(i));
		for(int i = 0; i < same.GetNPoints(); i++) _pts.push_back(same.at(i));
		for(int i = 0; i < high.GetNPoints(); i++) _pts.push_back(high.at(i));
	
	}
	
	//sort by weight
	void Sort(unsigned long long seed = 123){
		int N = _pts.size();
		if(N < 2) return;
		RandomSample rs(seed);	
		
		PointCollection low;
		PointCollection same;
		PointCollection high;
	
		rs.SetRange(0,N);
		int idx = rs.SampleFlat();		
		BayesPoint pivot = _pts[idx];
		if(pivot.Dim() < 1){
			cout << "Error: no value for pivot point #" << idx << " found." << endl;
			return;
		}
		for(int i = 0; i < N; i++){
			if(_pts[i].size() < 1){
				cout << "Error: point #" << i << " has no value set." << endl;
				return;
			}
			if( pivot.w() > _pts[i].w() )
				low += _pts[i];
			else if( pivot.w() == _pts[i].w() )
				same += _pts[i];
			else
				high += _pts[i];
		}
		low.Sort(seed);
		high.Sort(seed);
	
		_pts.clear();
		for(int i = 0; i < low.GetNPoints(); i++) _pts.push_back(low.at(i));
		for(int i = 0; i < same.GetNPoints(); i++) _pts.push_back(same.at(i));
		for(int i = 0; i < high.GetNPoints(); i++) _pts.push_back(high.at(i));
	
	}
	
	void AllSkip(){
		for(int i = 0; i < GetNPoints(); i++){
			if(!_pts[i].Skip()) _allskip = false;
		}
		_allskip = true;
	}

	double max(int d = 0){
		vector<double> pts;
		for(int i = 0; i < GetNPoints(); i++){
			//only include skipped pts if the pointcollection is made up of entirely skipped pt
			if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
			pts.push_back(_pts[i].Value(d));
		}
		return *std::max_element(pts.begin(), pts.end());
	}

	double min(int d = 0){
		vector<double> pts; 
		for(int i = 0; i < GetNPoints(); i++){
			//only include skipped pts if the pointcollection is made up of entirely skipped pt
			if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
			pts.push_back(_pts[i].Value(d));
		}
		return *std::min_element(pts.begin(), pts.end());
	}

	

	//center all dimensions independently
	//shifts average to zero
	BayesPoint Center(){
		//includes weights
		BayesPoint avg = mean();//BayesPoint(_nDim);//mean();
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				//only include skipped pts if the pointcollection is made up of entirely skipped pt
				if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
				_pts[i].SetValue(_pts[i].at(d) - avg.at(d),d);
				//round to 0
				if(fabs(_pts[i].at(d)) < 1e-16) _pts[i].SetValue(0,d);
			}
		}
		return avg;

	}


	//shifts min to zero
	BayesPoint MinCenter(){
		//vector<double> min;
		BayesPoint min = BayesPoint(_nDim);
		for(int d = 0; d < _nDim; d++){
			//min.push_back(this->min(d));
			min.SetValue(this->min(d),d);
			for(int i = 0; i < (int)_pts.size(); i++){
				//only include skipped pts if the pointcollection is made up of entirely skipped pt
				if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
				_pts[i].SetValue(_pts[i].at(d) - min.at(d),d);
			}
		}
		return min;
	}
	void Translate(const BayesPoint& t){
		if(t.Dim() != _nDim) return;
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				_pts[i].SetValue(_pts[i].at(d) - t.at(d),d);
				//round to 0
				if(fabs(_pts[i].at(d)) < 1e-16) _pts[i].SetValue(0,d);
			}
		}

	}
	void Translate(const vector<double>& t){
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				_pts[i].SetValue(_pts[i].at(d) - t[d],d);
				//round to 0
				if(fabs(_pts[i].at(d)) < 1e-16) _pts[i].SetValue(0,d);
			}
		}

	}
	void Translate(double t){
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				_pts[i].SetValue(_pts[i].at(d) - t,d);
				//round to 0
				if(fabs(_pts[i].at(d)) < 1e-16) _pts[i].SetValue(0,d);
			}	
		}        	

	}
	
	void Translate(double t, int d){
		for(int i = 0; i < (int)_pts.size(); i++){
			_pts[i].SetValue(_pts[i].at(d) - t,d);
			//round to 0
			if(fabs(_pts[i].at(d)) < 1e-16) _pts[i].SetValue(0,d);
		}

	}

	//commented out lines don't always hold - regular translate is being called in BasePDFMixture for shifting data	
	void CircularTranslate(double t, int d){
		double pi = acos(-1);
		for(int i = 0; i < (int)_pts.size(); i++){
		//turned off acos(cos()) from subcluster shifts	
			//if(fabs(_pts[i].at(d) - t) > acos(-1)){
			//	//if(cos(_pts[i].at(d) - t) > 0) _pts[i].SetValue(acos(cos(_pts[i].at(d) - t)),d);
			//	//else _pts[i].SetValue(-acos(cos(_pts[i].at(d) - t)),d);
			//	_pts[i].SetValue(-acos(cos(_pts[i].at(d) - t)), d);
			//}
			//else{
			//	_pts[i].SetValue(_pts[i].at(d) - t,d);
			//}
				//double theta = _pts[i].at(d) - t;
				////do wraparound - DOESNT WORK WHEN ORIGINAL POINTS ARE OUTSIDE [0, 2pi] (event with put02pi)
				//but don't force points on 02pi bc the translation factor may have been derived in a different frame
				//theta = fmod((theta + pi), 2*pi) - pi;
				//think of this as a rotation by -t to get t at 0
				double theta = atan2( sin(_pts[i].at(d) - t), cos(_pts[i].at(d) - t));
				_pts[i].SetValue(theta,d);
			//round to 0
			if(fabs(_pts[i].at(d)) < 1e-16){ _pts[i].SetValue(0,d);}
			if(fabs(fabs(_pts[i].at(d)) - 2*acos(-1)) < 1e-10){ _pts[i].SetValue(0,d);}
		}	
	}


	void EtaToTheta(int d){
		for(int i = 0; i < (int)_pts.size(); i++){
			double eta = _pts[i].at(d);
			_pts[i].SetValue(2*atan(exp(-eta)),d);
		}

	}
	
	void ThetaToEta(int d){
		for(int i = 0; i < (int)_pts.size(); i++){
			double theta = _pts[i].at(d);
			_pts[i].SetValue(-log(tan(theta/2)),d);
		}


	}


	void AngleToPlaneProject(int d){
		//for small theta, theta ~ sin(theta) ~ tan(theta), cos(theta) ~ 1
		//tan(theta/2) = sin(theta)/(1 + cos(theta))
		//for small theta => theta/2
		//use tan(theta/2) to get range to be (-pi,pi) as max deviation then multiply by 2 to get back original (small) theta
		for(int i = 0; i < (int)_pts.size(); i++){
			//_pts[i].SetValue(tan(_pts[i].at(d)),d);
			//if within range [-pi/2, pi/2], map to plane
			if(fabs(_pts[i].at(d)) < acos(-1)/2.){
				_pts[i].SetValue(tan(_pts[i].at(d)),d);
			}
			else{ //map to infinity - add original value to preserve unique coordinate
				if(_pts[i].at(d) > 0){
					_pts[i].SetValue(_infs[d],d);
				}
				else{	
					_pts[i].SetValue(-_infs[d],d);
				}
			}
		}
	}

	int GetInfVal(int d){
		return _infs[d];
	}


	bool HasInf(int d){
		for(int i = 0; i < (int)_pts.size(); i++){
			//if(fabs(_pts[i].at(d)) == _infs[d]){ cout << "inf here " << endl; _pts[i].Print(); return true; }
			if(fabs(_pts[i].at(d)) == _infs[d]){ return true; }
		}
		return false;
	}
	

	//domain only defined for [-pi/2,pi/2] (max deviation allowed)
	void PlaneToAngleProject(int d){
		//for small theta, theta ~ sin(theta) ~ tan(theta), cos(theta) ~ 1
		//tan(theta/2) = sin(theta)/(1 + cos(theta))
		//for small theta => theta/2
		//use tan(theta/2) to get range to be (-pi,pi) as max deviation then multiply by 2 to get back original (small) theta
		for(int i = 0; i < (int)_pts.size(); i++){
			double val = atan2(_pts[i].at(d),1);
			//put on [0,2pi] range - which may not be OG value if mirrored point...
			if(val < 0) val += 2*acos(-1);
			if(val >= 2*acos(-1)) val -= 2*acos(-1);
			_pts[i].SetValue(val,d);
		}

	}


	void Put02pi(int d){
		for(int i = 0; i < (int)_pts.size(); i++){
			_pts[i].Put02pi(d);
		}
		//double pi = acos(-1);
		//for(int i = 0; i < (int)_pts.size(); i++){
		//	int nit = 0;
		//	double ogval = _pts[i].at(d);
		//	while(!(_pts[i].at(d) >= 0 && _pts[i].at(d) < 2*pi)){
		//		//if pt is negative
		//		if(_pts[i].at(d) < -1e-16) _pts[i].SetValue(_pts[i].at(d) + 2*pi,d);
		//		//if pt is geq 2*pi
		//		else _pts[i].SetValue(_pts[i].at(d) - 2*pi, d);
		//		nit++;
		//		if(nit > 100){
		//			cout << "Error: PointCollection::Put02pi for dim " << d << " not converging for pt #" << i << " setting to value " << _pts[i].at(d) << " with original value " << ogval << endl;
		//			break;
		//		}
		//	}
		//}
	}

	//normalize all dimensions indepedently
	//vector<double> Normalize(){
	BayesPoint Normalize(){
		//translate first
		MinCenter();
		double min, max;
		//vector<double> scale;
		BayesPoint scale = BayesPoint(_nDim);
		for(int d = 0; d < _nDim; d++){
			min = this->min(d);
			max = this->max(d);
			for(int i = 0; i < (int)_pts.size(); i++){
				if(_pts[i].Skip() && _pts[i].GetSkipDim() == d) continue;
				_pts[i].SetValue(_pts[i].at(d)/(max - min),d);
			}
			//scale.push_back(1./(max - min));
			scale.SetValue((max - min),d);
		}
		return scale;
	}
	
	void Scale(const BayesPoint& s){
		if(s.Dim() != _nDim) return;
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				_pts[i].SetValue(_pts[i].at(d)*(s.at(d)),d);
			}
		}

	}
	void Scale(const vector<double>& s){
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				_pts[i].SetValue(_pts[i].at(d)*(s[d]),d);
			}
		}

	}
	void Scale(double s){
		for(int d = 0; d < _nDim; d++){
			for(int i = 0; i < (int)_pts.size(); i++){
				_pts[i].SetValue(_pts[i].at(d)*(s),d);
			}
		}

	}
	

	BayesPoint mean() const{
		BayesPoint m = BayesPoint(_nDim);
		for(int d = 0; d < _nDim; d++) m.SetValue(0.,d);
		for(int i = 0; i < (int)_pts.size(); i++){
			for(int d = 0; d < _nDim; d++){
				//only include skipped pts if the pointcollection is made up of entirely skipped pt
				if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
				m.SetValue(m.at(d)+_pts[i].at(d),d);
			}
		}
		for(int d = 0; d < _nDim; d++) m.SetValue(m.at(d)/(double)_pts.size(),d);
		return m;
	}

	BayesPoint stddev() const{
		BayesPoint sd = BayesPoint(_nDim);
		BayesPoint m = mean();
		double val;
		for(int d = 0; d < _nDim; d++) m.SetValue(0.,d);
		for(int i = 0; i < (int)_pts.size(); i++){
			for(int d = 0; d < _nDim; d++){
				//only include skipped pts if the pointcollection is made up of entirely skipped pt
				if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
				val = pow(_pts[i].at(d) - m.at(d),2);
				sd.SetValue(sd.at(d)+val,d);
			}
		}
		return sd;
		for(int d = 0; d < _nDim; d++) sd.SetValue(sqrt(sd.at(d)/(double)_pts.size()),d);
	}

	void GetWeights(vector<double>& w){
		w.clear();
		for(int i = 0; i < (int)_pts.size(); i++)
			w.push_back(_pts[i].w());
	}

	void SetWeights(vector<double>& w){
		for(int i = 0; i < (int)_pts.size(); i++)
			_pts[i].SetWeight(w[i]);
	}

	//unweighted circular mean	
	double CircularMean(int d) const{
		//make new PointCollection that will store angles on [0, 2pi]
		PointCollection wrapped = PointCollection(*this);
		wrapped.Put02pi(d);
		double avg_s = 0;
		double avg_c = 0;
		double max = -1e308;
		double min = 1e308;
		for(int i = 0; i < (int)_pts.size(); i++){
			//only include skipped pts if the pointcollection is made up of entirely skipped pt
			if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
			avg_s += sin(wrapped.at(i).at(d));
			avg_c += cos(wrapped.at(i).at(d));
			if(wrapped.at(i).at(d) > max) max = wrapped.at(i).at(d);
			if(wrapped.at(i).at(d) < min) min = wrapped.at(i).at(d);
		}
		avg_s /= (double)_pts.size();
		avg_c /= (double)_pts.size();
		double mean = atan2(avg_s,avg_c);
		//put atan2 in [0,twopi)
		double pi = 4*atan(1);
		//2pi is only defined up to 10 decimal places in Dnn2piCylinder (which is reasonable) so as not to double count 0 + 2pi, round here
		if(fabs(mean) < 1e-10) mean = 0;
		//if mean is truly negative (ie on [-pi, 0)) then reflect across to match wrapped [0, 2pi) interval
		if(mean < 0) mean += 2*pi;
		//get raw mean of original data
		double euclid_mean = this->mean().at(d);	
		//find best unwrapped version of mean
		double best = mean;
		//if points cross 0--2pi seam, just use atan2 wrapped mean
		//explicity check for seam crossings
		//if((max > 3*pi/2 && max <= 2*pi) && (min < pi/2 && min >= 0.)){	
		if(max - min > pi){	
			if(euclid_mean < 0 && mean > 0){ return mean - 2*pi; }
			//if(euclid_mean > 2*pi && mean < 2*pi){ cout << " euclid_mean " << euclid_mean << " mean " << mean << " would return " << mean + 2*pi << " points " <<  endl; this->Print();}
			if(euclid_mean > 2*pi && mean < 2*pi){ return mean + 2*pi; }
			return mean;
		}
		double min_diff = fabs(mean - euclid_mean);
		for(int k = -4; k <=4; k++){
			double cand = mean + k*2*pi;
			double diff = fabs(cand - euclid_mean);
			if(diff < min_diff){
				best = cand;
				min_diff = diff;
			}
		}
		//cout << "CircularMean() - wrapped mean " << mean << " raw mean " << euclid_mean << " unwrapped mean " << best << endl;
		return best;
	};

	//weighted circular mean	
	double CircularCentroid(int d) const{
		//make new PointCollection that will store angles on [0, 2pi]
		PointCollection wrapped;
		for(int i = 0; i < (int)_pts.size(); i++){
			wrapped += _pts[i];
		}		
		wrapped.Put02pi(d);
		double avg_s = 0;
		double avg_c = 0;
		double max = -1e308;
		double min = 1e308;
		for(int i = 0; i < (int)_pts.size(); i++){
			//only include skipped pts if the pointcollection is made up of entirely skipped pt
			if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
			avg_s += sin(wrapped.at(i).at(d))*wrapped.at(i).w();
			avg_c += cos(wrapped.at(i).at(d))*wrapped.at(i).w();
			if(wrapped.at(i).at(d) > max) max = wrapped.at(i).at(d);
			if(wrapped.at(i).at(d) < min) min = wrapped.at(i).at(d);
		}
		avg_s /= (double)_pts.size();
		avg_c /= (double)_pts.size();
		double mean = atan2(avg_s,avg_c);
		//put atan2 in [0,twopi)
		double pi = 4*atan(1);
		if(mean < 0) mean += 2*pi;
		//get raw mean of original data
		double euclid_mean = this->Centroid(d);
		//find best unwrapped version of mean
		double best = mean;
		//if points cross 0--2pi seam, just use atan2 wrapped mean
		//explicity check for seam crossings
		if(max - min > pi){	
			if(euclid_mean < 0 && mean > 0){ return mean - 2*pi; }
			//if(euclid_mean > 2*pi && mean < 2*pi){ cout << " euclid_mean " << euclid_mean << " mean " << mean << " would return " << mean + 2*pi << " points " <<  endl; this->Print();}
			if(euclid_mean > 2*pi && mean < 2*pi){ return mean + 2*pi; }
			return mean;
		}
		double min_diff = fabs(mean - euclid_mean);
		for(int k = -4; k <=4; k++){
			double cand = mean + k*2*pi;
			double diff = fabs(cand - euclid_mean);
			if(diff < min_diff){
				best = cand;
				min_diff = diff;
			}
		}
		//cout << "CircularCentroid() - wrapped mean " << mean << " raw mean " << euclid_mean << " unwrapped mean " << best << endl;
		return best;
	};

	//weighted mean	
	double Centroid(int d) const{
		double cent = 0;
		double sum = 0;
		for(int i = 0; i < (int)_pts.size(); i++){
			//only include skipped pts if the pointcollection is made up of entirely skipped pt
			if(_pts[i].Skip() && _pts[i].GetSkipDim() == d && !_allskip) continue;
			cent += _pts[i].at(d)*_pts[i].w();
			sum += _pts[i].w();
		}
		cent /= sum;
		return cent;
	};

	double Sumw() const{
		double ret = 0;
		for(int i = 0; i < (int)_pts.size(); i++) ret += _pts[i].w();
		return ret;
	}

	void SetPointSkip(int n, int d){
		if(_skipDim != -1)
			_pts[n].SetSkipDim(d);
	}
	void GetSkipDims(vector<int>& ds) const{
		ds.clear();
		for(int n = 0; n < _pts.size(); n++){
			ds.push_back(_pts[n].GetSkipDim());
		}
	}
	void GetUserIdxs(vector<int>& ds) const{
		ds.clear();
		for(int n = 0; n < _pts.size(); n++){
			ds.push_back(_pts[n].GetUserIdx());
		}
	}

	//want to check points with wraparound but save points in original phi
	void GetUniquePoints(const vector<PointCollection>& collections){
		Clear();
		PointCollection uniquepts;
		for(int i = 0; i < collections.size(); i++){
			//put phi on 02pi
			if(false) cout << "collection #" << i << " has " << collections[i].GetNPoints() << " points" << endl;
			for(int j = 0; j < collections[i].GetNPoints(); j++){
				if(false) cout << "checking pt # " << j << " for collection # " << i << endl;
				BayesPoint pt = collections[i]._pts[j];
				pt.Put02pi(1);
				bool has = uniquepts.Has(pt);
				if(false){ cout << std::setprecision(20) << "Checking pt" << endl; pt.Print(); cout << "original pt" << endl; collections[i]._pts[j].Print(); }
				if(false) cout << "uniquepts has pt? " << has << endl;
				if(false) cout << "checking pt # " << j << " for collection # " << i << " has? " << has << endl;
				if(uniquepts.Has(pt)) continue;
				this->AddPoint(collections[i]._pts[j]);
				uniquepts.AddPoint(pt);
			}
		}
	}

	private:
		int _nDim = 0;
		vector<BayesPoint> _pts = {};
		vector<double> _infs = {};
		int _skipDim = -1;
		bool _allskip = false;	

};
#endif
