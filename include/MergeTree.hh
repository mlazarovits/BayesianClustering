#ifndef MergeTree_HH
#define MergeTree_HH

#include "BaseTree.hh"
#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "NodeList.hh"

using node = BaseTree::node;
class MergeTree : public BaseTree{
	public:
		MergeTree() : BaseTree(){ _alpha = 0; }

		MergeTree(double alpha) : BaseTree(){
			_alpha = alpha;
		}

		//copy constructor
		MergeTree(const MergeTree& tree){
			_head = tree._head;
			_z = tree._z;
			_t = tree._t;
			_alpha = tree._alpha;
			_clusters = tree._clusters;
		}

		virtual ~MergeTree(){  }

		void AddData(PointCollection* pc){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));
			}
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
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)
			x->prob_tk = Evidence(x);//_model->ConjugateEvidence(*pt);
			if(pt != nullptr) x->points = new PointCollection(*pt);
			_clusters.push_back(x);
		}

		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
			int k;
			//if leaf node (ie r == _z && l == _z) -> set k = 1
			if(x->l == _z && x->r == _z) k = 1;
			//number of clusters in node x = k_l + k_r for left and right nodes
			else k = x->l->model->GetNClusters() + x->r->model->GetNClusters();

			x->model = new GaussianMixture(k); //p(x | theta)
			x->model->SetData(x->points);
			x->model->InitParameters();
			x->model->InitPriorParameters();
			
			VarEMCluster* algo = new VarEMCluster(x->model, k);	
			algo->SetThresh(1.);
			//cluster
			double oldLogL = algo->EvalLogL();
			double LogLThresh = 0.01;
			double newLogL, dLogL; 
			while(dLogL < LogLThresh){
				newLogL = algo->Cluster();
				dLogL = newLogL - oldLogL;
				oldLogL = newLogL;
			}
			return newLogL;
		}


	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _clusters;
		//Dirichlet prior parameter
		double _alpha;


};
#endif
