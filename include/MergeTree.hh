#ifndef MergeTree_HH
#define MergeTree_HH

#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "NodeStack.hh"

//using node = BaseTree::node;
class MergeTree : BDATA
(0.504315,0.882420,1.097047)
(0.249523,-0.069064,-0.501077)
(0.799247,1.714346,2.391901)
(0.314927,0.822628,1.321417)
(0.085644,0.038067,0.008893)
(0.580332,1.479676,2.369276)aseTree{
	public:
		MergeTree(){ 
		cout << "empty ctor - z val: " << _z->val << ": " <<  _z << endl;
		_alpha = 0; _c = -1; 
		cout << "empty ctor end - z val: " << _z->val << ": " <<  _z << endl;
		}

		MergeTree(double alpha){
		cout << "alpha ctor - z val: " << _z->val << ": " <<  _z << endl;
			_alpha = alpha;
			_c = -1;
		cout << "alpha ctor end - z val: " << _z->val << ": " << _z << endl;
		}

		//copy constructor
		MergeTree(const MergeTree& tree){
			_head = tree._head;
			_z = tree._z;
			_t = tree._t;
			_alpha = tree._alpha;
			_clusters = tree._clusters;
			_c = tree._c;
		}

		virtual ~MergeTree(){ _c = 0; }


		double zval(){ cout << _z << " "; return _z->val; }

		void AddData(PointCollection* pc){
		cout << "AddData - z val: " << _z->val << endl;
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));
			}
		}
		
		node* Get(int i){ return _clusters[i]; }

		vector<node*> GetClusters(){ return _clusters; }
/*
//not actually called
		node* Merge(node *l, node *r){
			//calculate p_lr (posterior from these two nodes)
			struct node* x = CalculateMerge(l, r);
		//increment color counter
			_c++; 	
			//remove nodes l and r
			Remove(l);
			Remove(r);
			//insert x into tree
			_clusters.push_back(x);
			return x;
		}
	*/
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
			cout << "z val: " << _z->val << endl;
			x->l = _z; x->r = _z;
			cout << "add leaf" << endl;
			pt->Print();
			cout << "add leaf with left val: " << x->l->val << endl;
			cout << "add leaf with right val: " << x->r->val << endl;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			x->val = 1.;	
			x->d = _alpha; 
			if(pt != nullptr) x->points = new PointCollection(*pt);
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)
			x->prob_tk = Evidence(x);//_model->ConjugateEvidence(*pt);
			x->color = _c;
			_clusters.push_back(x);
		}

		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
		//	cout << "MergeTree::Evidence" << endl;
			int k;
			//if leaf node (ie r == _z && l == _z) -> set k = 1
			if(x->l == _z && x->r == _z) k = 1;
			//number of clusters in node x = k_l + k_r for left and right nodes
			else k = x->l->model->GetNClusters() + x->r->model->GetNClusters();
			//cout << "max k: " << k << endl;
			x->model = new GaussianMixture(k); //p(x | theta)
			x->model->SetData(x->points);
			x->model->InitParameters();
			x->model->InitPriorParameters();
			
		//	cout << "data" << endl;
		//	x->model->GetData()->Print();

			VarEMCluster* algo = new VarEMCluster(x->model, k);	
			algo->SetThresh(1.);
			
			//cluster
			double oldLogL = algo->EvalLogL();
			double LogLThresh = 0.01;
			double newLogL;
			double dLogL = 999; 
			int it = 0;
			while(dLogL > LogLThresh){
				newLogL = algo->Cluster();
		//cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << " old ELBO: " << oldLogL << " new ELBO: " << newLogL << endl;
				dLogL = newLogL - oldLogL;
				oldLogL = newLogL;
				it++;
			}
			return newLogL;
		}


	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _clusters;
		//Dirichlet prior parameter
		double _alpha;
		int _c;

};
#endif
