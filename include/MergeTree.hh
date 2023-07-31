#ifndef MergeTree_HH
#define MergeTree_HH

#include "BaseTree.hh"
#include "BasePDFMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"


using node = BaseTree::node;

class MergeTree : public BaseTree{
	public:
		MergeTree(){ _alpha = 0; }

		MergeTree(PointCollection* pc){
			//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));	
			}
			NodeSort(_roots);
		}

		//copy constructor
		MergeTree(const MergeTree& tree){
			_alpha = tree._alpha;
			_roots = tree._roots;
			//mixture model
			_mix_model = tree._mix_model;
			_model = tree._model;
			_prior = tree._prior;	
		}

		virtual ~MergeTree(){ _roots.clear(); }

		node* Get(int i){ return _roots[i]; }

		void Merge(node *l, node *r){
			//calculate p_lr (posterior from these two nodes)
			double p = CalculateMerge(l, r); 		
			//combine points from l + r into one pc
			PointCollection* newpts;
			newpts->AddPoints(*l->points);
			newpts->AddPoints(*r->points);
			//construct new node x
			struct node *x = (struct node*)malloc(sizeof *x);
			//assign v (posterior) as val and pc as points
			x->val = p;
			x->points = newpts;
			//calculate and set d_k
			x->d = _alpha*tgamma(newpts->GetNPoints()) + l->d*r->d;
			//need to set prob_tk
				//pass by reference from CalculateMerge? return and recalculate posterior here?
			x->l = l;
			x->r = r;
			//remove nodes l and r
			Remove(l, r);
			//insert x into tree
			Insert(x);
		}

		//assuming Dirichlet Process Model (sets priors)
		double CalculateMerge(node *l, node* r);
		

		void Remove(node *l, node *r){
			//remove nodes l and r (that combine merge x)
			auto it = find(_roots.begin(), _roots.end(), l);
			int idx;
			// If element was found
			if (it != _roots.end()) 
				idx = it - _roots.begin();
			else
				idx = -1;
		
			_roots.erase(_roots.begin(),_roots.begin()+idx);
			it = find(_roots.begin(), _roots.end(), r);
			// If element was found
			if (it != _roots.end()) 
				idx = it - _roots.begin();
			else
				idx = -1;
			_roots.erase(_roots.begin(),_roots.begin()+idx);
		}

		void Insert(node *x){
			//insert node into tree (vector)
			//call insert sort
			_roots.push_back(x);
			InsertSort(x, _roots);
		}

		void SetAlpha(double alpha){ _alpha = alpha; }	

		void SetModel(BasePDF* m){ _model = m; _prior = m->GetPrior(); }


	protected:
		void AddLeaf(const Point* pt = nullptr){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			node* x = (node*)malloc(sizeof *x);
			x->l = _z; x->r = _z;
			//x->val = val;
			//if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			x->val = 1.;	
			x->d = _alpha; 
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)
			x->prob_tk = _model->ConjugateEvidence(*pt);
			if(pt != nullptr) x->points = new PointCollection(*pt);
			_roots.push_back(x);
		
		}

		//TEST
		//insert sort for adding one node
		void InsertSort(node* x, vector<node*> nodes){
			int n = (int)nodes.size();
			double v = x->val;
			for(int i = 1; i < n; i++){
				int j = i;
				while(nodes[j-1]->val > v){
					nodes[j] = nodes[j-1]; j--;
				}
				nodes[j] = x;
				break;
			}
		}

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

	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _roots;
		//Dirichlet prior parameter
		double _alpha;
		//mixture model
		BasePDFMixture* _mix_model;
		BasePDF* _model;
		BasePDF* _prior;	


};
#endif
