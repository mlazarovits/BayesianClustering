#ifndef EMCluster_HH
#define EMCluster_HH
#include "BaseCluster.hh"
#include "BasePDF.hh"
#include "Matrix.hh"
#include "PointCollection.hh"
#include <vector>
using std::vector;


class EMCluster : public BaseCluster{
	//import all inherited ctors
	using BaseCluster::BaseCluster;
	public:
		virtual ~EMCluster(){ };
	
		void Initialize(unsigned long long seed = 111);
		//E-step - calculates posterior
		void Estimate();
		//M-step - updates parameters
		void Update();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();
		//eval - returns log-likelihood value at given iteration
		double EvalVariationalLogL(const Matrix& post);
		


};


#endif
