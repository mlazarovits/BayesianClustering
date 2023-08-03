#ifndef MergeTree_HH
#define MergeTree_HH

#include "BaseTree.hh"
#include "BasePDFMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "NodeList.hh"

using node = BaseTree::node;
class MergeTree : public BaseTree{
	public:
		MergeTree() : BaseTree(){ _alpha = 0; }

		MergeTree(BasePDF* model) : BaseTree(){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			_model = model; _prior = model->GetPrior(); 
		}

		//copy constructor
		MergeTree(const MergeTree& tree){
			_head = tree._head;
			_z = tree._z;
			_t = tree._t;
			_alpha = tree._alpha;
			_clusters = tree._clusters;
			//mixture model
			_model = tree._model;
			_prior = tree._prior;	
		}

		virtual ~MergeTree(){  }

		void AddData(PointCollection* pc){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));
		//		_clusters[i]->name = std::to_string(i);
			}
		//	NodeSort(_clusters);
		}
		
		node* Get(int i){ return _clusters[i]; }

		vector<node*> GetClusters(){ return _clusters; }

		node* Merge(node *l, node *r){
			//calculate p_lr (posterior from these two nodes)
			struct node* x = CalculateMerge(l, r); 	
			//remove nodes l and r
			Remove(l);
			Remove(r);
			//insert x into tree
			_clusters.push_back(x);
			return x;
		}

		void Insert(node* x){
			_clusters.push_back(x);
		}

		//assuming Dirichlet Process Model (sets priors)
		node* CalculateMerge(node *l, node* r);
		

		void Remove(node *l){
			//remove nodes l and r (that combine merge x)
			auto it = find(_clusters.begin(), _clusters.end(), l);
			int idx;
			// If element was found
			if (it != _clusters.end()){ 
				idx = it - _clusters.begin();
			}
			else
				return;
			_clusters.erase(_clusters.begin()+idx);
		}

		void SetAlpha(double alpha){ _alpha = alpha; }	

		int GetNClusters(){ return (int)_clusters.size(); }	



	protected:
		void AddLeaf(const Point* pt = nullptr){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			node* x = (node*)malloc(sizeof *x);
			x->l = _z; x->r = _z;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			x->val = 1.;	
			x->d = _alpha; 
			//////initialize probability of subtree to be null hypothesis for leaf
			//////p(D_k | T_k) = p(D_k | H_1^k)
			if(_model == nullptr) cout << "null" << endl;
			x->prob_tk = _model->ConjugateEvidence(*pt);
			if(pt != nullptr) x->points = new PointCollection(*pt);
			_clusters.push_back(x);
		}


	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _clusters;
		//Dirichlet prior parameter
		double _alpha;
		BasePDF* _model;
		BasePDF* _prior;	


};
#endif
