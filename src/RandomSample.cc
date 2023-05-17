#include "RandomSample.hh"
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
using std::pair;

RandomSample::RandomSample(unsigned long long seed){
	_b = 4101842887655102017LL;
	_c = 1;
	_a = seed ^ _b;

	int64();
	_b = _a;
	int64();
	_c = _b;
	int64();


}


RandomSample::~RandomSample() { }

//returns pseudo-random 64bit integer
unsigned long long RandomSample::int64(){
	_a = _a * 2862933555777941757LL + 7046029254386353087LL;
	_b ^= _b >> 17;
	_b ^= _b << 31;
	_b ^= _b >> 8;
	_c = 4294957665U*(_c & 0xffffffff) + (_c >> 32);

	unsigned long long x = _a ^ (_a << 21);
	x ^= x >> 35;
	x ^= x << 4;
	return (x + _b) ^ _c;
}


//return random 32bit integer
unsigned int RandomSample::int32(){
	return (unsigned int)int64();
}

//return random double between (0,1) (uniform distribution)
double RandomSample::rand(){
	return 5.42101086242752217E-20 * int64();
}

//flat distribution scaled to Gaussian max
double RandomSample::FlatGausScaled(){
	return 1./sqrt(2.*acos(-1));
}

//normal distribution with specified mean and sigma
double RandomSample::Gaussian(double x, double mu, double sigma){
	double x_scaled = x - mu;
	return 1./(sigma*sqrt(2.*acos(-1)))*exp(-(x_scaled*x_scaled)/(2.*sigma*sigma));
}

//get random x value according to flat distribuion
double RandomSample::SampleFlat(){
	return _xmin + (_xmax - _xmin)*rand();
}

void RandomSample::SetRange(double xmin, double xmax){
	if(xmin >= xmax) return;
	_xmin = xmin;
	_xmax = xmax;
}



//return random double (xmin, xmax) according to Gaussian distribuion
vector<double> RandomSample::SampleGaussian(double mean, double sigma, int Nsample){
	vector<double> samples;
	if(sigma < 0){
		cout << "Please input valid sigma > 0" << endl;
		return samples;
	}
	double X, ran, R;
	//double samples[Nsample];
	int Ntrial = 0;
	for(int i = 0; i < Nsample; i++){
		Ntrial += 1;
		X = SampleFlat();
		R = Gaussian(X,mean,sigma)/FlatGausScaled();
		ran = rand();
		if(ran > R){ // reject if "die thrown" (random number) is above distribution
			i -= 1; // decrease i and try again
			continue;
		} else // accept if "die" is below (within) distribution
			samples.push_back(X);
	}
	return samples;
}
/*
//return random double (xmin, xmax) according to Gaussian distribuion
PointCollection RandomSample::SampleGaussian(vector<double> mean, vector<double> sigma, int Nsample, int dim){
	PointCollection samples;
	if((int)mean.size() != dim || (int)sigma.size() != dim || mean.size() != sigma.size()){
		cout << "Please provide consistent parameters and dimensionality. Mean: " << mean.size() << " sigma: " << sigma.size() << " dimensions: " << dim << endl;
		return samples;
	}

	vector<vector<double>> dims(dim);
	//sample each dimension independently
	for(int d = 0; d < dim; d++){
		if(sigma[d] < 0){
			cout << "Please input valid sigma > 0" << endl;
			return samples;
		}
		dims[d] = SampleGaussian(mean[d], sigma[d], Nsample);
	}
	vector<Point> pts;
	for(int i = 0; i < Nsample; i++){
		VD pt;
		for(int d = 0; d < dim; d++){
			pt += dims[d][i];
		}
		pts.push_back(Point(pt));
	}	
	samples = PointCollection(pts);
	return samples;
}

//Randomly selects points from in array, fills out array with selected points
vector<double> RandomSample::SelectPoints(vector<double> in, int nIn, int nOut){
	SetRange(0,nIn);
	vector<double> out;
	for(int i = 0; i < nOut; i++){
		int idx = SampleFlat();
		out.push_back(in[idx]);
	}
	return out;

}
//Randomly selects points from in array, fills out array with selected points

vector<pair<double,double>> RandomSample::SelectPairs(vector<pair<double,double>> in, int nIn, int nOut){
	SetRange(0,nIn);
	vector<pair<double,double>> out;
	for(int i = 0; i < nOut; i++){
		//choose a random index
		int idx = SampleFlat();
		out.push_back(in[idx]);
	}
	return out;

}

vector<map<double,double>> RandomSample::SelectMultiDimPairs(vector<map<double,double>> in, int nIn, int nOut){
	cout << "SelectMultiDimPoints" << endl;
	SetRange(0,nIn);
	vector<map<double,double>> out;
	int nDim = in.size();
	map<double,double> map1D;
	for(int d = 0; d < nDim; d++){
		for(int i = 0; i < nOut; i++){
			//choose a random rank (key)
			double key = SampleFlat()/nIn;
			map1D[key] = in[d][key];
			cout << "chose point: (" << key << ", " << map1D[key] << ")" << endl;
		}
		out.push_back(map1D);
		map1D.clear();
	}
	return out;

}

*/



