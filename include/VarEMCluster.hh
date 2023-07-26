#ifndef VarEMCluster_HH
#define VarEMCluster_HH

#include "BaseCluster.hh"

class VarEMCluster : public BaseCluster{
	//import all inherited ctors
	using BaseCluster::BaseCluster;
	public:
		virtual ~VarEMCluster(){ };	

		//E-step - calculates posterior
		void Estimate();
		//M-step - updates parameters
		void Update();
		//eval - returns log-likelihood value at given iteration
		double EvalLogL();
		





};
#endif
