#ifndef Tools_HH
#define Tools_HH


#include<vector>
#include<iostream>

using std::cout;
using std::endl;
using std::vector;
using std::pair;

class VD : public vector<double> {
public:
  VD() {}

  VD(const vector<double>& vd){
    for(auto d : vd)
      *this += d;
  }

  VD(int d){
	  for(int i = 0; i < d; i++) (*this).push_back(0.);
  }

  VD& a(double d){
    (*this).push_back(d);
    return *this;
  }
  
  virtual ~VD() {}

  VD& operator += (double d){
    this->push_back(d);
    return *this;
  }

  VD& operator += (const VD& list){
    for(int i = 0; i < int(list.size()); i++)
      this->push_back(list[i]);
  
    return *this;
  }

  bool operator ==(const VD& l2){
	  if((*this).size() != l2.size()){
		  return false;
	  }
	  for(int i = 0; i < (int)(*this).size(); i++){
		  if( (*this)[i] != l2[i])
			  return false;
	  }
	  return true;
  }

  bool operator !=(const VD& l2){
	  if((*this).size() != l2.size()){
		  return true;
	  }
	  for(int i = 0; i < (int)(*this).size(); i++){
		  if( (*this)[i] != l2[i])
			  return true;
	  }
	  return false;
  }
 };


/*
class VpairD : public vector<pair<double,double>> {
public:
  VpairD() {}
  
  VpairD(int N) {
  	for(int i = 0; i < N; i++)
		*this += 0.;
  
  }

  VpairD(const vector<pair<double,double>>& vd){
    for(auto d : vd)
      *this += d;
  }

  VpairD& a(pair<double,double> d){
    (*this).push_back(d);
    return *this;
  }

  VD GetFirst(){
	VD ret;
	for(int i = 0; i < this->size(); i++) ret += (*this)[i].first;
  	return ret;
  }
  
  VD GetSecond(){
	VD ret;
	for(int i = 0; i < this->size(); i++) ret += (*this)[i].second;
  	return ret;
  }

  //quick sort by first element in pair
  void Sort(unsigned long long seed = 111){
	  int N = (int)this->size();
	if(N < 2) return;

	VpairD low;
	VpairD same;
	VpairD high;

	RandomSample rs(seed);
	rs.SetRange(0,N);

	double pivot = (*this)[(int)rs.SampleFlat()].first;

	for(int i = 0; i < N; i++){
		// Elements that are smaller than the `pivot` go to
    		// the `low` list. Elements that are larger than
    		// `pivot` go to the `high` list. Elements that are
    		// equal to `pivot` go to the `same` list.
		if((*this)[i].first < pivot)
			low += (*this)[i];
		else if((*this)[i].first == pivot)
			same += (*this)[i];
		else
			high += (*this)[i];

	}

	//continue recursively
	low.Sort();
	high.Sort();
	
	//recreate original vector
	(*this).clear();
	for(auto val : low)
		(*this) += val;
	for(auto val : same)
		(*this) += val;
	for(auto val : high)
		(*this) += val;

  }


  
  virtual ~VpairD() {}

  VpairD& operator += (pair<double,double> d){
    this->push_back(d);
    return *this;
  }

  VpairD& operator += (const VpairD& list){
    for(int i = 0; i < int(list.size()); i++)
      this->push_back(list[i]);
  
    return *this;
  }
};
*/
#endif
