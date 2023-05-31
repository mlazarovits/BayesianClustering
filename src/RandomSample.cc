#include "RandomSample.hh"
#include "PointCollection.hh"
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


RandomSample::~RandomSample() { 

}

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
	unsigned long long x = int64();
	double ret = 5.42101086242752217E-20 * x;
	return ret;
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
	double x = _xmin + (_xmax - _xmin)*rand();
	return x;
}

//get random n-dim x value according to flat distribuion
Point RandomSample::SampleNDimFlat(int dim){
	Point x = Point(dim);
	for(int i = 0; i < dim; i++)
		x.SetValue(_xmin + (_xmax - _xmin)*rand(), i);
	return x;
}

void RandomSample::SetRange(double xmin, double xmax){
	if(xmin >= xmax) return;
	_xmin = xmin;
	_xmax = xmax;
}



//return random double (xmin, xmax) according to Gaussian distribuion
vector<double> RandomSample::SampleGaussian(double mean, double sigma, int Nsample){
	cout << "RandomSample::SampleGaussian with " << Nsample << " samples" << endl;
	vector<double> samples;
	//vector<double> test = {0., 1.};
	//test.emplace_back(0);
	//cout << "emplaced back 1" << endl;
	//test.emplace_back(1);
	//cout << "emplaced back 2" << endl;
	//test.emplace_back(2);
	//cout << "pushed back 3" << endl;
/*
	if(sigma < 0){
		cout << "Please input valid sigma > 0" << endl;
		return samples;
	}
	double X, ran, R;
	int Ntrial = 0;

	for(int i = 0; i < Nsample; i++){
		Ntrial += 1;
cout << "i: " << i << endl;
		X = i;//SampleFlat();
cout << "X: " << &X << " " << &samples << endl;
		R = Gaussian(X,mean,sigma)/FlatGausScaled();
		ran = rand();
		if(ran > R){ // reject if "die thrown" (random number) is above distribution
			i -= 1; // decrease i and try again
			continue;
		} else{ // accept if "die" is below (within) distribution
		cout << "accept1: " << samples.size() << " X: " << &ran << endl;
			samples.push_back(X);
		cout << "accept2: " << samples.size()  << " i: " << i << " " << &samples[i] << endl;
		}
cout << "\n" << endl;
	}
*/
	cout << "RandomSample::SampleGaussian - end" << endl;

	return samples;
}
