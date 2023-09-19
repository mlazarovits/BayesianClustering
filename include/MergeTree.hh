#ifndef MergeTree_HH
#define MergeTree_HH

#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "NodeStack.hh"

//using node = BaseTree::node;
class MergeTree : BaseTree{
	public:
		MergeTree(){ 
		_alpha = 0; _thresh = 0.; _verb = 0; 
		}

		MergeTree(double alpha){
			_alpha = alpha;
			_thresh = 1.; _verb = 0;
		}

		//copy constructor
		MergeTree(const MergeTree& tree){
			_head = tree._head;
			_z = tree._z;
			_t = tree._t;
			_alpha = tree._alpha;
			_clusters = tree._clusters;
			_thresh = tree._thresh;
			_verb = tree._verb;
		}

		virtual ~MergeTree(){ }

		void SetThresh(double t){ _thresh = t; }

		void AddData(PointCollection* pc){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));
			}
		}
		
		node* Get(int i){ return _clusters[i]; }

		vector<node*> GetClusters(){ return _clusters; }
		
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


		void SetDataSmear(const Matrix& cov){ _data_smear = cov; }

		void SetVerbosity(int v){ _verb = v; }

		void SetPriorParameters(map<string, Matrix> params){ _params = params; }

	protected:
		void AddLeaf(const Point* pt = nullptr){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			node* x = (node*)malloc(sizeof *x);
			//cout << "z val: " << _z->val << endl;
			x->l = _z; x->r = _z;
			//cout << "add leaf" << endl;
			//pt->Print();
			//cout << "add leaf with left val: " << x->l->val << endl;
			//cout << "add leaf with right val: " << x->r->val << endl;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			x->val = 1.;	
			x->d = _alpha; 
			if(pt != nullptr) x->points = new PointCollection(*pt);
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)		
			x->prob_tk = exp(Evidence(x)); //Evidence = ELBO \approx log(LH)
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
			if(_verb != 0) x->model->SetVerbosity(_verb-1);
			x->model->SetAlpha(_alpha);
			x->model->SetData(x->points);
			if(!_data_smear.empty()) x->model->SetDataSmear(_data_smear);
			x->model->InitParameters();
			x->model->InitPriorParameters();

			if(!_params.empty()) x->model->SetPriorParameters(_params);


			/*
			cout << "Initial prior parameters" << endl;
			map<string, Matrix> params;
			for(int i = 0; i < x->model->GetNClusters(); i++){
				params = x->model->GetPriorParameters(i);	
				cout << "mean " << i << endl;
				params["m"].Print();
				cout << "Gaus scale " << i << ": " << params["scale"].at(0,0) << endl;
				cout << "dof " << i << ": " << params["dof"].at(0,0) << endl;
				cout << "scalemat " << i << endl;
				params["scalemat"].Print();
				params.clear();
		}
			*/
		//	cout << "data" << endl;
		//	x->model->GetData()->Print();


			VarEMCluster* algo = new VarEMCluster(x->model, k);	
			if(x->points->GetNPoints() > 2) algo->SetThresh(_thresh);
			
			//cluster
			double oldLogL = algo->EvalLogL();
			double LogLThresh = 0.01;
			double newLogL;
			double dLogL = 999; 
			int it = 0;
			while(dLogL > LogLThresh){
				newLogL = algo->Cluster();
		if(isnan(newLogL)) cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << " old ELBO: " << oldLogL << " new ELBO: " << newLogL << endl;
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
		//threshold on variational EM
		double _thresh;

		Matrix _data_smear;
		int _verb;
		map<string, Matrix> _params;
};
#endif
