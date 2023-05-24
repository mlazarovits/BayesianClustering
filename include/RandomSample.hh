#ifndef RandomSample_HH
#define RandomSample_HH
#include <cmath>
#include <vector>
#include <map>
//#include "PointCollection.hh"
#include "Point.hh"
#include "Matrix.hh"

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
		double SampleFlat();
		Point RandomSample::SampleNDimFlat(int dim)
		void SetRange(double xmin, double xmax);
		std::vector<double> SampleGaussian(double mean, double sigma, int Nsample);
		double NDimGaussian(Point x, Matrix mu, Matrix sigma);
		Matrix SampleNDimGaussian(Matrix mean, Matrix sigma, int Nsample);
		//PointCollection SampleGaussian(vector<double> mean, vector<double> sigma, int Nsample, int dim = 1);
		//std::vector<double> SelectPoints(vector<double> in, int nIn, int nOut);	
		//vector<pair<double,double>> SelectPairs(vector<pair<double,double>> in, int nIn, int nOut);
		//vector<map<double,double>> SelectMultiDimPairs(vector<map<double,double>> in, int nIn, int nOut);

		double _xmax = 5;
		double _xmin = -5;

	private:
		unsigned long long _a;
		unsigned long long _b;
		unsigned long long _c;





};







#endif
