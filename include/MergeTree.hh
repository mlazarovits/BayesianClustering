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
			_roots = tree._roots;
			//mixture model
			_model = tree._model;
			_prior = tree._prior;	
		}

		virtual ~MergeTree(){  }

		void AddData(PointCollection* pc){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));
			}
		//	NodeSort(_roots);
		}
		
		node* Get(int i){ return _roots[i]; }

		node* Merge(node *l, node *r){
			//calculate p_lr (posterior from these two nodes)
			struct node* x = CalculateMerge(l, r); 		
			//set index to lower branch
			x->idx = x->l->idx;
			//remove nodes l and r
			Remove(l);
			Remove(r);
			//insert x into tree
			_roots.push_back(x);
			return x;
		}

		void Insert(node* x){
			_roots.push_back(x);
		}

		//assuming Dirichlet Process Model (sets priors)
		node* CalculateMerge(node *l, node* r);
		

		void Remove(node *l){
			//remove nodes l and r (that combine merge x)
			auto it = find(_roots.begin(), _roots.end(), l);
			int idx;
			// If element was found
			if (it != _roots.end()){ 
				idx = it - _roots.begin();
			}
			else
				return;
			_roots.erase(_roots.begin()+idx);
		}

		void SetAlpha(double alpha){ _alpha = alpha; }	

		int GetNPoints(){ return (int)_roots.size(); }	



	protected:
		void AddLeaf(const Point* pt = nullptr){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			node* x = (node*)malloc(sizeof *x);
			x->l = _z; x->r = _z;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			x->val = 1.;	
			x->d = _alpha; 
			x->idx = -1;	
			//////initialize probability of subtree to be null hypothesis for leaf
			//////p(D_k | T_k) = p(D_k | H_1^k)
			if(_model == nullptr) cout << "null" << endl;
			x->prob_tk = _model->ConjugateEvidence(*pt);
			if(pt != nullptr) x->points = new PointCollection(*pt);
			_roots.push_back(x);
		}


/*
		//quicksort for nodes
		void NodeSort(vector<node*> nodes, unsigned long long seed = 123){
			int N = _roots.size();
			if(N < 2) return;
			RandomSample rs(seed);	
			
			vector<node*> low;
			vector<node*> same;
			vector<node*> high;
			rs.SetRange(0,N);
			int idx = rs.SampleFlat();		
			node* pivot = nodes[idx];
			double v;
			for(int i = 0; i < N; i++){
				v = pivot->val;
				if( v > _roots[i]->val )
					low.push_back(nodes[i]);
				if( v == _roots[i]->val )
					same.push_back(nodes[i]);
				else
					high.push_back(nodes[i]);
			}
			NodeSort(low);
			NodeSort(same);
			NodeSort(high);
	
			_roots.clear();
			for(int i = 0; i < (int)low.size(); i++) _roots.push_back(low.at(i));
			for(int i = 0; i < (int)same.size(); i++) _roots.push_back(same.at(i));
			for(int i = 0; i < (int)high.size(); i++) _roots.push_back(high.at(i));

		}
	*/

	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _roots;
		//Dirichlet prior parameter
		double _alpha;
		BasePDF* _model;
		BasePDF* _prior;	


};
#endif
