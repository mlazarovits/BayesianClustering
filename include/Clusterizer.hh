#ifndef CLUSTERIZER_HH
#define CLUSTERIZER_HH

#include "PointCollection.hh"
#include "Jet.hh"
#include "JetTree.hh"
#include "GaussianMixture.hh"
#include <string>
using std::string;
using std::vector;

//this class is a wrapper for the clustering algorithms (see: fastjet::ClusterSequence)
//should b a recursive implementation
class Clusterizer{
	public:
		Clusterizer();
		Clusterizer(vector<Jet> jets);
		Clusterizer(vector<JetPoint> rhs);
		Clusterizer(Jet jet);
		virtual ~Clusterizer();

		//set point smear
		void SetDataSmear(const Matrix& cov){ _data_smear = cov; }
		void SetPriorParameters(map<string, Matrix> params){ _params = params; }

		//runs everything (varGMM + BHC)
		//change to run over generic points (or vector of rhs)
		void Cluster(Jet jet, double alpha = 0.1, double thresh = 1., bool viz = false, int verb = 0, string fname = "");

		//just runs varGMM over given jets
		void FindSubjets(Jet jet, double alpha = 0.1, double thresh = 1., bool viz = false, int verb = 0, int maxK = 2, string fname = "");


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
		Matrix _data_smear;

		map<string, Matrix> _params;

};
#endif