#include "RandomSample.hh"
#include "Matrix.hh"
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
//n-dim normal distribution with specified mean (n x 1 matrix) and sigma (n x n matrix) for one n-dim point x
double RandomSample::NDimGaussian(Point x, Matrix mu, Matrix sigma){
	/*
	int dim = x.Dim();
	Matrix mu_scaled = Matrix(x.Dim(),1);
	for(int d = 0; d < dim; d++)
		mu_scaled.SetEntry(x.Value(d) - mu.at(d,0), d, 0);
	double det = sigma.det();
	//invert sigma then multiply by x_scaled*sigma*x_scaledT
	Matrix sigInv;
	sigInv.invert(sigma);
	Matrix muT;
	muT.transpose(mu_scaled);
	//should only be 1x1 matrix - one entry
	Matrix expon = muT.mult(sigInv.mult(mu_scaled));
	
	return 1./(det*sqrt(2.*acos(-1)))*exp(-0.5*expon.at(0,0));
	*/
	return 0.;
}


//get random x value according to flat distribuion
double RandomSample::SampleFlat(){
	return _xmin + (_xmax - _xmin)*rand();
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
//return Nsample random d-dim points (xmin, xmax) according to n-dim Gaussian distribuion
Matrix RandomSample::SampleNDimGaussian(Matrix mean, Matrix sigma, int Nsample){
	//need cholesky decomposition: takes sigma = LL^T
	//returns Nsample normally distributed n-dim points (one point = x)
	//x = L*y + mu for vectors x, y, mu and matrix L
	Matrix L = sigma.cholesky();	
	int d = mean.GetDims()[0];


	//aggregate points
	PointCollection pc;
	for(int n = 0; n < Nsample; n++){
		//y_i ~ N(0,1)
		vector<double> y = SampleGaussian(0,1,d);
		Matrix mat_y(y);
		//L*y	
		Matrix Ly = Matrix(d,1);
		Ly.mult(L,mat_y);	
		//L*y + mu
		//for(int i = 0; i < d; i++)
		//pc += pt;
	}

/*
	//samples should be an n x d (n rows of d-dim data pts) matrix
	int dim = mean.GetDims()[0];
	Matrix samples = Matrix(Nsample,dim);
//	if(sigma < 0){
//		cout << "Please input valid sigma > 0" << endl;
//		return samples;
//	}
	Point X;
	double ran, R;
	for(int i = 0; i < Nsample; i++){
		int Ntrial = 0;
		for(int d = 0; d < dim; d++){ 
			Ntrial += 1;
			X = SampleFlat();
			R = NDimGaussian(X,mean,sigma)/FlatGausScaled();
			ran = rand();
			if(ran > R){ // reject if "die thrown" (random number) is above distribution
				i -= 1; // decrease i and try again
				continue;
			} else // accept if "die" is below (within) distribution
				samples.SetEntry(X,n,d);
		}
	}
	return samples;
*/
	return L;
}
