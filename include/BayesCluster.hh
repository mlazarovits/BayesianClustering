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
#include "Matrix.hh"
#include "GaussianMixture.hh"
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
			
			_thresh = -999;
			_alpha = 0.1;
			_subalpha = 0.5;
			
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);

			_smear = Matrix();
			_verb = 0;
			_trees = {};

			_check_merges = false;
			
			//add initial jets (rec hits) to merge tree as leaves
			int n = (int)_jets.size();
			vector<double> weight;
       			for (int i = 0; i < n; i++) {
				//if jet has no rec hits - skip
				if(_jets[i].GetNRecHits() < 1) continue;
				PointCollection pc;
				BayesPoint pt(3);
				pt.SetValue(_jets[i].rap(), 0);
				pt.SetValue(_jets[i].phi_02pi(), 1);
				pt.SetValue(_jets[i].time(), 2);
				_jets[i].GetWeights(weight);
				pt.SetWeight(weight[0]); 
			//	//make sure phi is in the right range - [0,2pi]
				_sanitize(pt);
				if(!(pt.at(1) >= 0.0 && pt.at(1) < 2*acos(-1))) cout << "i: " << i << " bad phi: " << pt.at(1) << endl;
				assert(pt.at(1) >= 0.0 && pt.at(1) < 2*acos(-1));
				pc.AddPoint(pt);
				_points.push_back(pc);
			}
			_cell = acos(-1)/180;
			_tresCte = 0.133913;
			_tresStoch = 1.60666; 
			_tresNoise = 0.00691415;

		};


		virtual ~BayesCluster(){
			////for(auto n : _trees) delete n;
			for(auto n : _trees) n = nullptr;
			//_trees.clear();
			_jets.clear();
			_points.clear();
			_history.clear();
			_prior_params.clear();
		};


		//for jets - BHC for clusters + GMM EM for subclusters
		const vector<node*>& NlnNCluster(){
			return this->_delauney_cluster();
		}

		const vector<node*>& N2Cluster(){
			return this->_naive_cluster();
		}

		//for photons - subclusters only
		GaussianMixture* SubCluster(string oname = ""){
			return this->_subcluster(oname);
		}

		void SetThresh(double t){ _thresh = t; }
		void SetDataSmear(const Matrix& cov){ _smear = cov; }		
		void SetAlpha(double a){ _alpha = a; }
		void SetSubclusterAlpha(double a){ _subalpha = a; }
		void SetPriorParameters(map<string, Matrix> params){ _prior_params = params; }
		void SetVerbosity(int v){ _verb = v; }
		bool _check_merges;
		void CheckMerges(bool t){ _check_merges = t; }

		double _cell, _tresCte, _tresStoch, _tresNoise;
		void SetMeasErrParams(double spatial, double tresCte, double tresStoch, double tresNoise){ _cell = spatial; _tresCte = tresCte; _tresStoch = tresStoch; _tresNoise = tresNoise; 
		//cout << "BayesCluster set - Using _tresCte = " << _tresCte << " ns, _tresStoch = " << _tresStoch << " ns and _tresNoise = " << _tresNoise << " ns" << endl;
		}

	protected:
		//need to typedef some stuff to build probability map used for determining cluster pairs
		typedef std::pair<int,int> verts;
		typedef std::pair<double,verts> RkEntry;
		typedef std::pair<double,verts> CompEntry;
		typedef std::pair<int,pair<int,double>> InvCompEntry;
		//use a multimap so multiple keys can have the same value
		//also it's automatically sorted
		typedef std::multimap<double,verts> CompareMap;
		typedef std::map<int,pair<int,double>> InvCompareMap;
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

			/// the merge probability corresponding to the recombination
			/// at this stage of the clustering.
			double rk;  

			/// maximum value of merge posterior rk  
			/// at this stage of the clustering.
			/// if max_rk < 0.5, clustering can stop
			double max_rk_so_far;
		};
		void _sanitize(BayesPoint pt){
			if(pt.Dim() < 2) return;
  			double twopi = 6.28318530717;
			double second = pt.at(1);
			if (second <  0)     second += twopi; 
			if (second >= twopi) second -= twopi;
			pt.SetValue(second,1);
		}



		//void _initialize_and_run(){ 
		//	_fill_initial_history();
		//	// don't run anything if the event is empty
		//	if (n_particles() == 0) return;	

		//	this->_cluster();
		//};
		void _initialize(){ 
			_fill_initial_history();
			// don't run anything if the event is empty
			if (n_particles() == 0) return;	

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
				element.rk = -1.0;	
				element.max_rk_so_far = -1.0;	
		 
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


		void _do_ij_recombination_step(
				const int jet_i, const int jet_j,
				const double rk,
				int& newjet_k){
			Jet newjet;
			if(_verb > 0){
				cout << "recombining jets " << jet_i << ": ";
				_jets[jet_i].Print(); cout << "eta " << _jets[jet_i].eta() << " phi " << _jets[jet_i].phi() << " time " << _jets[jet_i].time() << endl;
				cout << "and " << jet_j << ": ";
				_jets[jet_j].Print(); cout << "eta " << _jets[jet_j].eta() << " phi " << _jets[jet_j].phi() << " time " << _jets[jet_j].time() << endl;
			}

			recombine(_jets[jet_i], _jets[jet_j], newjet);
			_jets.push_back(newjet);
			newjet_k = _jets.size() - 1;

			//do history stuff
			int newstep_k = _history.size();
			//provide jet with history info
			_jets[newjet_k].SetHistoryIndex(newstep_k);

			//sort history
			int hist_i = _jets[jet_i].GetHistoryIndex();
			int hist_j = _jets[jet_j].GetHistoryIndex();

			//add_step_to_history(min(hist_i, hist_j), max(hist_i, hist_j), newjet_k, rk);
			
		}

		//the straight up addition of 3-mom and E is the "E_scheme" recombination
		//option in FastJet
		void recombine(const Jet& pa, const Jet& pb, Jet& pab) const{ 
			pab = Jet(pa);
			pab.add(pb);
		}


		// initialise the history in a standard way
		void _add_step_to_history (
			             const int parent1, 
			             const int parent2, const int jetp_index,
			             const double rk) {
 
			history_element element;
			element.parent1 = parent1;
			element.parent2 = parent2;
			element.jetp_index = jetp_index;
			element.child = Invalid;
			element.rk   = rk;
			element.max_rk_so_far = std::max(rk,_history[_history.size()-1].max_rk_so_far);
			_history.push_back(element);
 
			int local_step = _history.size()-1;
			// sanity check: make sure the particles have not already been recombined
			//
			// Note that good practice would make this an assert (since this is
			// a serious internal issue). However, we decided to throw an
			// InternalError so that the end user can decide to catch it and
			// retry the clustering with a different strategy.
 
			assert(parent1 >= 0);
			if (_history[parent1].child != Invalid){
			  cout << "trying to recomine an object that has previsously been recombined" << endl;
			}
			_history[parent1].child = local_step;
			if (parent2 >= 0) {
			  if (_history[parent2].child != Invalid){
			  cout << "trying to recomine an object that has previsously been recombined" << endl;
			  }
			  _history[parent2].child = local_step;
			}
 
			// get cross-referencing right from PseudoJets
			if (jetp_index != Invalid) {
			  assert(jetp_index >= 0);
			  _jets[jetp_index].SetHistoryIndex(local_step);
			}
 
		}
		void _add_entry_to_maps(const int i, CompareMap& inmap, const Dnn2piCylinder* DNN, bool prob = true);
		void _add_entry_to_maps(const int i, InvCompareMap& inmap, const Dnn2piCylinder* DNN);
		void _phi_wraparound(PointCollection& pc);

	private:
		vector<Jet> _jets;
		vector<PointCollection> _points; //to pass to merge tree in cluster function
		vector<node*> _trees;
		vector<history_element> _history;
		double _Qtot;
		enum JetType {Invalid=-3, InexistentParent = -2, BeamJet = -1};
		int _initial_n;
		double _thresh, _alpha, _subalpha;
		map<string, Matrix> _prior_params;
		Matrix _smear;
		int _verb;

		const vector<node*>& _delauney_cluster();
		const vector<node*>& _naive_cluster();
		GaussianMixture* _subcluster(string oname = "");
		int n_particles() const{ return _initial_n; }


	 
};
#endif

