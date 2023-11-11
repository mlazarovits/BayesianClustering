#ifndef RandomSample_HH
#define RandomSample_HH
#include <cmath>
#include <vector>
#include <map>
#include "Point.hh"

using std::map;
using std::vector;
using std::pair;

class RandomSample{
	public:
		RandomSample(unsigned long long seed = 111);
		virtual ~RandomSample();
		unsigned long long int64();
		unsigned int int32();
		double rand();
		double FlatGausScaled();
		double Gaussian(double x, double mu, double sigma);
		double Poisson(int x, double rate);
		double SampleFlat();
		Point SampleNDimFlat(int dim);
		void SetRange(double xmin, double xmax);
		vector<double> SampleGaussian(double mean, double sigma, int Nsample);
		vector<double> SamplePoisson(double rate, int Nsample);
		int SampleCategorical(const vector<double>& weights);

		double _xmax = 5;
		double _xmin = -5;

	private:
		unsigned long long _a;
		unsigned long long _b;
		unsigned long long _c;





};







#endif
