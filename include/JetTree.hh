#ifndef JETTREE_HH
#define JETTREE_HH

#include "Jet.hh"

using std::vector;

//this class keeps track of the cluster (merge) history
class JetTree{
	public:
		JetTree();
		virtual ~JetTree();

		void AddLayer(vector<Jet> jets);	
	
		//jets at depth i
		vector<Jet> at(int i);

		//jet j at depth i
		Jet at(int i, int j);

		int GetNLayers();
		int GetNJets();

		//removal functions?

		//check clustertree is valid
		//parents are in layer before, less or equal jets in successive layer, etc.
		bool check_jet_in_layer(Jet& jet, int i);
		void check_validity();

	private:
		//m_tree[i] = vector<jet> jets at depth i -> m_tree[i][j] = jet j at depth i
		vector<vector<Jet>> m_tree;
















};
#endif
