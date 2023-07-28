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

		MergeTree(vector<double> vals, PointCollection* pc){
			//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < (int)vals.size(); i++){
				AddLeaf(vals[i], &pc->at(i));	
			}
		}

		virtual ~MergeTree(){ _roots.clear(); }

		node* Get(int i){ return _roots[i]; }

		void Merge(node *l, node *r){
			//calculate p_lr (posterior from these two nodes)
			double p = CalculateMerge(l, r); 		
			//combine points from l + r into one pc
			PointCollection* newpts;
			newpts->add(*l->points);
			newpts->add(*r->points);
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
			_roots.push_back(x);
		
		}

		void InsertSort(node* x, vector<node*> nodes){
			int n = (int)nodes.size();
			double v = x->val;
			for(int i = 0; i < n; i++){
				int j = i;
				while(nodes[j-1]->val > v){
					nodes[j] = nodes[j-1]; j--;
				}
				nodes[j] = x;
			}
		}

		//quicksort for nodes
		void NodeSort(unsigned long long seed = 123){
			int N = _roots.size();
			if(N < 2) return;
			RandomSample rs(seed);	
			
			vector<node*> low;
			vector<node*> same;
			vector<node*> high;
/*	
			rs.SetRange(0,N);
			int idx = rs.SampleFlat();		
			node* pivot = _roots[idx];
			for(int i = 0; i < N; i++){
				if( pivot.ge(_pts[i], d) )
					low += _pts[i];
				else if( pivot.eq(_pts[i], d) )
					same += _pts[i];
				else
					high += _pts[i];
			}
			low.Sort(d);
			high.Sort(d);
	
			_pts.clear();
			for(int i = 0; i < low.GetNPoints(); i++) _pts.push_back(low.at(i));
			for(int i = 0; i < same.GetNPoints(); i++) _pts.push_back(same.at(i));
			for(int i = 0; i < high.GetNPoints(); i++) _pts.push_back(high.at(i));
*/

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
