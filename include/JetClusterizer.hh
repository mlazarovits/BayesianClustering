#ifndef JETCLUSTERIZER_HH
#define JETCLUSTERIZER_HH

#include "PointCollection.hh"
#include "Jet.hh"
#include "JetTree.hh"
#include "GaussianMixture.hh"

using std::vector;

//this class is a wrapper for the clustering algorithms (see: fastjet::ClusterSequence)
//should b a recursive implementation
class JetClusterizer{
	public:
		JetClusterizer();
		JetClusterizer(vector<Jet> jets);
		JetClusterizer(vector<JetPoint> rhs);
		JetClusterizer(Jet jet);
		virtual ~JetClusterizer();

		//call recursively - runs everything (varGMM + BHC)
		void Cluster();

		GaussianMixture* FindSubjets(PointCollection* points, double LogLthresh, int maxNit, int maxK, bool viz);
		//just runs varGMM over given jets
		vector<Jet> FindSubjets_XYZ(Jet jet, double LogLthresh = 0.0001, int maxNit = 1, int maxK = 10, bool viz = false);
		//just runs varGMM over given jets
		vector<Jet> FindSubjets_etaPhi(Jet jet, double LogLthresh = 0.0001, int maxNit = 1, int maxK = 10, bool viz = false);

//		void SetMaxNClusters(int k){ m_maxK = k; }	

		JetTree GetTree() const{ return m_tree; }

	private:
		//BHC object
		JetTree m_tree;
		//jets to points
		PointCollection m_points;
		//running number of jets
		int m_nJets;
		vector<Jet> m_newJets;
		vector<Jet> m_oldJets;




};
#endif
