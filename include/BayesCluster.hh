#ifndef BAYESCLUSTER_HH
#define BAYESCLUSTER_HH

//The structure of this class is respectfully repurposed from ClusterSequence in FastJet (Cacciari, Salam, Soyez).

// This work was modified from its original form by Margaret Lazarovits on October 2, 2023. 
// The original version of this work was released
// under version 2 of the GNU General Public License. As of v3 of GNU GPL,
// any conditions added in Section 7 also apply. 
//----------------------------------------------------------------------
// Copyright (c) 2005-2021, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//----------------------------------------------------------------------


#include "Dnn2piCylinder.hh"
#include "Jet.hh"
#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;


class BayesCluster{
	public:
		BayesCluster(const vector<Jet>& pseudojets){ 
			// this will ensure that we can point to jets without difficulties
			// arising.
			_jets.reserve(pseudojets.size()*2);
			
			// insert initial jets this way so that any type L that can be
			// converted to a pseudojet will work fine (basically PseudoJet
			// and any type that has [] subscript access to the momentum
			// components, such as CLHEP HepLorentzVector).
			for (unsigned int i = 0; i < pseudojets.size(); i++) {
			  _jets.push_back(pseudojets[i]);}
	
			_initialize_and_run();
		};
		virtual ~BayesCluster(){ _jets.clear(); };


	protected:
		//need to typedef some stuff to build probability map used for determining cluster pairs
		typedef std::pair<int,int> verts;
		typedef std::pair<double,verts> RkEntry;
		//use a multimap so multiple keys can have the same value
		//also it's automatically sorted
		typedef std::multimap<double,verts> ProbMap;

		struct history_element{
			/// index in _history where first parent of this jet
			/// was created (InexistentParent if this jet is an
			/// original particle)
			int parent1; 
 
			/// index in _history where second parent of this jet
			/// was created (InexistentParent if this jet is an
			/// original particle); BeamJet if this history entry
			/// just labels the fact that the jet has recombined
			/// with the beam)
			int parent2; 
 
			/// index in _history where the current jet is
			            /// recombined with another jet to form its child. It
			            /// is Invalid if this jet does not further
			            /// recombine.
			int child;   
 
			/// index in the _jets vector where we will find the
			/// PseudoJet object corresponding to this jet
			/// (i.e. the jet created at this entry of the
			/// history). NB: if this element of the history
			/// corresponds to a beam recombination, then
			/// jetp_index=Invalid.
			int jetp_index; 
 
			/// the distance corresponding to the recombination
			/// at this stage of the clustering.
			double dij;  
 
			/// the largest recombination distance seen
			/// so far in the clustering history.
			double max_dij_so_far;


			/// maximum value of merge posterior rk  
			/// at this stage of the clustering.
			/// if max_rk < 0.5, clustering can stop
			double max_rk;
		};


		void _initialize_and_run(){ 
			_fill_initial_history();
			// don't run anything if the event is empty
			if (n_particles() == 0) return;	

			this->_cluster();
		};


		void _fill_initial_history () {
			//if (_jets.size() == 0) {throw Error("Cannot run jet-finder on empty event");}
			// reserve sufficient space for everything
			_jets.reserve(_jets.size()*2);
			_history.reserve(_jets.size()*2);
		 
			_Qtot = 0;
		 
			for (int i = 0; i < static_cast<int>(_jets.size()) ; i++) {
				history_element element;
				element.parent1 = InexistentParent;
				element.parent2 = InexistentParent;
				element.child   = Invalid;
				element.jetp_index = i;
				element.dij     = 0.0;
				element.max_dij_so_far = 0.0;
				element.max_rk = -1.0;	
		 
				_history.push_back(element);
			
				// not sure if i need the stuff below
				/*	
				// do any momentum preprocessing needed by the recombination scheme
				_jet_def.recombiner()->preprocess(_jets[i]);
		 
				// get cross-referencing right from PseudoJets
				_jets[i].set_cluster_hist_index(i);
				_set_structure_shared_ptr(_jets[i]);
		 		*/
				// determine the total energy in the event
				_Qtot += _jets[i].E();
			}
			_initial_n = _jets.size();
		}

	private:
		vector<Jet> _jets;
		vector<history_element> _history;
		double _Qtot;
		enum JetType {Invalid=-3, InexistentParent = -2, BeamJet = -1};
		int _initial_n;

		void _cluster();
		int n_particles() const;
	 
};
#endif

