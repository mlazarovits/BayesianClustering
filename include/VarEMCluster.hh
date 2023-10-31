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
		
		//set threshold for E[pi] cutoff
		void SetThresh(double t){ _thresh = t; }
		//this is done automatically in full algo but not in the GMM/EM subroutine (ie for photons)
		//so need to do this by hand
		void SetClusterStart(){ clustering_start = true; }

		
	private:
		double _thresh = 0;



};
#endif
