#ifndef JETCLUSTERIZER_HH
#define JETCLUSTERIZER_HH


#include "Jet.hh"
#include "VarGaussianMixture.hh"
#include "GaussianMixture.hh"

using std::vector;

//this class is a wrapper for the clustering algorithms (see: fastjet::ClusterSequence)
//should b a recursive implementation
class JetClusterizer{
	public:
		JetClusterizer();
		virtual ~JetClusterizer();











};
#endif
