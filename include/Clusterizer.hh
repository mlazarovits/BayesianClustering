#ifndef CLUSTERIZER_HH
#define CLUSTERIZER_HH

#include "PointCollection.hh"
#include "Jet.hh"
#include "JetTree.hh"
#include "GaussianMixture.hh"
#include "BaseTree.hh"
#include <string>
using std::string;
using std::vector;
using node = BaseTree::node;

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
		void SetDataSmear(const Matrix& cov){ _data_smear = cov; _smeared = true; }
		void SetPriorParameters(map<string, Matrix> params){ _params = params; }
		void SetThresh(double t){_thresh = t; }
		void SetVerbosity(int v){ _verb = v; }
		void SetMaxNClusters(int k){ _maxK = k; }
		void SetWeighted(bool w){_weighted = w; }
		void SetDistanceConstraint(bool d){ _distconst = d;}	

		//runs everything (varGMM + BHC)
		//change to run over generic points (or vector of rhs)
		vector<node*> Cluster(Jet jet, string fname = "");

		//just runs varGMM over given jets
		GaussianMixture* FindSubjets(Jet jet, string fname = "");


//		void SetMaxNClusters(int k){ m_maxK = k; }	

		JetTree GetTree() const{ return m_tree; }

		//set alpha for BHC
		void SetClusterAlpha(double a){_bhcAlpha = a;}
		//set alpha for EM
		void SetSubclusterAlpha(double a){_emAlpha = a;}

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

		double _bhcAlpha;
		double _emAlpha;
		double _thresh;
		int _verb;
		int _maxK;
		bool _weighted;
		bool _smeared;
		bool _distconst;

};
#endif
