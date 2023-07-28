#ifndef MergeTree_HH
#define MergeTree_HH

#include "BaseTree.hh"
#include "BasePDFMixture.hh"
#include "MultivarT.hh"

using node = BaseTree::node;

class MergeTree : public BaseTree{
	public:
		MergeTree(){ _alpha = 0; }

		MergeTree(vector<double> vals, PointCollection* pc){
			//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < (int)vals.size(); i++){
				AddLeaf(vals[i], &pc->at(i));	
			}
		}

		virtual ~MergeTree(){ _nodes.clear(); }

		node* Get(int i){ return _nodes[i]; }

		void Merge(node *l, node *r){
			//double v = CalculateMerge(l, r): calculate p_lr (posterior from these two nodes)
			//combine points from l + r into one pc
			//construct new node x, assign v (posterior) as val and pc as points, calculate and set d_k, need to set prob_tk
				//pass by reference from CalculateMerge? return and recalculate posterior here?
			//remove nodes l and r
			//insert x into tree
		
		}

		//assuming Dirichlet Process Model (sets priors)
		double CalculateMerge(node *l, node* r);
		

		void RemoveMerge(node *l, node *r){
			//remove nodes l and r (that combine merge x)
		}

		void InsertMerge(node *x){
			//insert node into tree (vector)
			//call insert sort
		}

		void SetAlpha(double alpha){ _alpha = alpha; }	

		void SetModel(BasePDF* m){ _model = m; _prior = m->GetPrior(); }


	protected:
		void AddLeaf(double val, const Point* pt = nullptr){
			node* x = (node*)malloc(sizeof *x);
			x->l = _z; x->r = _z;
			x->val = val;
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			x->d = _alpha; 
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)
			x->prob_tk = _model->ConjugateEvidence(*pt);
			if(pt != nullptr) x->points = new PointCollection(*pt);
			_nodes.push_back(x);
		
		}

	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _nodes;
		//Dirichlet prior parameter
		double _alpha;
		//mixture model
		BasePDFMixture* _mix_model;
		BasePDF* _model;
		BasePDF* _prior;	


};
#endif
