#ifndef MergeTree_HH
#define MergeTree_HH

#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "NodeStack.hh"

class MergeTree : BaseTree{
	public:
		MergeTree(){ 
		_alpha = 0; _thresh = 0.; _verb = 0;
		_constraint_a = 0; _constraint_b = acos(-1)/2.; 
		_constraint = false;
		}

		MergeTree(double alpha){
			_alpha = alpha;
			_thresh = 1.; _verb = 0;
			_constraint_a = 0; _constraint_b = acos(-1)/2.;
			_constraint = false;
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
			_constraint_a = tree._constraint_a; _constraint_b = tree._constraint_b;
			_constraint = tree._constraint;
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

		const vector<node*>& GetClusters() const{ 
			return _clusters;
		}
		
		void Insert(node* x){
			_clusters.push_back(nullptr);
			x->idx = (int)_clusters.size() - 1;
			_clusters[(int)_clusters.size() - 1] = x;
		}

		//assuming Dirichlet Process Model (sets priors)
		node* CalculateMerge(node *l, node* r);
		double CalculateMerge(int i, int j);
		node* Merge(node* l, node* r);
		void Merge(int i, int j);

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
			//_clusters.erase(_clusters.begin()+idx);
			//setting the node to null matches the implementation in the FastJet code (see DnnPlane)
			_clusters[idx] = NULL;
		}

		void Remove(int i){ _clusters[i] = NULL; }

		void SetAlpha(double alpha){ _alpha = alpha; }	

		//counts clusters that have been removed
		//use this function when looping over clusters s.t. the element indices line up
		int GetNAllClusters(){ 			
			return (int)_clusters.size();
		}	

		//this does not count removed clusters
		//use this function to get number of actual clusters
		//when no merges have been made (ie clusters are only leaves)
		//this should match GetNAllClusters 
		int GetNClusters(){ 			
			int n = 0;
			for(int i = 0; i < _clusters.size(); i++){
				if(_clusters[i] == nullptr) continue;
				n++;
			}
			return n;
		}	
			

		void SetSubclusterAlpha(double a){ _emAlpha = a; }		

		void SetDataSmear(const Matrix& cov){ _data_smear = cov; }

		void SetVerbosity(int v){ _verb = v; }

		void SetPriorParameters(map<string, Matrix> params){ _params = params; }

		void SetDistanceConstraint(double a, double b){ _constraint = true; _constraint_a = a; _constraint_b = b; }
		//uses a triangular distribution to constraint certain clusterings
		double DistanceConstraint(node* i, node* j){
			_constraint = true;
			double c = (_constraint_a + _constraint_b)/2.;
			double pi = acos(-1);
			//phi
			double cent1 = i->points->Centroid(1);
			double cent2 = j->points->Centroid(1);
			
			//transform deltaPhi to be on [0,pi], wrapped s.t. 0 is close to 2pi (-3 close to 3)
			double d = fabs(cent1 - cent2);
			//update 3D nearest neighbors for mirror point calculation
			double dist3d = _euclidean_3d(i, j);
			if(dist3d < i->nn3dist) i->nn3dist = dist3d;
			if(dist3d < j->nn3dist) j->nn3dist = dist3d;	

	
			double phi = 0;
			double theta = 0;
			if(d >= _constraint_a && d <= _constraint_b) phi = cos(d);
		
	
			//eta	
			cent1 = i->points->Centroid(0); 
			cent2 = j->points->Centroid(0); 
		
			//eta to theta
			cent1 = 2*atan(exp(-cent1));
			cent2 = 2*atan(exp(-cent2));
			//don't need to wrap eta -> only goes from 0 to pi in theta
			d = fabs(cent1 - cent2);	
			if(d >= _constraint_a && d <= _constraint_b) theta = cos(d);
			return theta*phi;		
		}

		void AddLeaf(const Point* pt = nullptr){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			_clusters.push_back(nullptr);
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
			int n = _clusters.size();
			x->idx = n-1;
			_clusters[n-1] = x;
		}


		//note: this ONLY gets called in the N^2 clustering
		//for the delauney (NlnN) clustering this is taken care of
		//by Dnn2piCylinder
		void CreateMirrorNode(node* x);

	protected:

		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
			int k;
  			//if leaf node (ie r == _z && l == _z) -> set k = 1
			if(x->l == _z && x->r == _z) k = 1;
			//number of clusters in node x = k_l + k_r for left and right nodes
			else k = x->l->model->GetNClusters() + x->r->model->GetNClusters();
			
			x->model = new GaussianMixture(k); //p(x | theta)
			if(_verb != 0) x->model->SetVerbosity(_verb-1);
			x->model->SetAlpha(_emAlpha);
			x->model->SetData(x->points);
			if(!_data_smear.empty()) x->model->SetDataSmear(_data_smear);
			x->model->InitParameters();
			x->model->InitPriorParameters();

			if(!_params.empty()) x->model->SetPriorParameters(_params);

			VarEMCluster* algo = new VarEMCluster(x->model, k);	
			if(x->points->Sumw() >= _thresh) algo->SetThresh(_thresh);
			//cluster
			double oldLogL = algo->EvalLogL();
			//if(algo->GetNClusters() < 1){ cout << "Error: threshold too high for successful clustering with " << x->points->GetNPoints() << " points (" << x->points->Sumw() << " weighted). Please adjust accordingly to continue hierarchically clustering" << endl; return -999; }
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


		double _euclidean_3d(node* i, node* j){
			double deta = i->points->mean().at(0) - j->points->mean().at(0);
			double dphi = i->points->mean().at(1) - j->points->mean().at(1);
			double dtime = i->points->mean().at(2) - j->points->mean().at(2);

			return sqrt(deta*deta + dphi*dphi + dtime*dtime);

		}


	private:
		//keep list of nodes since tree is built bottom up
		vector<node*> _clusters;
		//Dirichlet prior parameter
		double _alpha, _emAlpha;
		//threshold on variational EM
		double _thresh;
		//endpoints for distance constraining
		double _constraint_a, _constraint_b;

		Matrix _data_smear;
		int _verb;
		map<string, Matrix> _params;
		bool _constraint;
		void _remap_phi(PointCollection& points);
};
#endif
