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
		_alpha = 0; _thresh = 0.; _verb = 0; _wraparound = false; 
		}

		MergeTree(double alpha){
			_alpha = alpha;
			_thresh = 1.; _verb = 0;_wraparound = false;
		}

		//copy constructor
		MergeTree(const MergeTree& tree){
			_head = tree._head;
			_z = tree._z;
			_t = tree._t;
			_alpha = tree._alpha;
			_clusters = tree._clusters;
			_thresh = tree._thresh;
			_verb = tree._verb;_wraparound = tree._wraparound;
		}

		virtual ~MergeTree(){ }

		void SetThresh(double t){ _thresh = t; }

		void SetPhiWraparound(bool p){_wraparound = p; }

		void AddData(PointCollection* pc){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(&pc->at(i));
			}
		}
		
		node* Get(int i){ return _clusters[i]; }

		vector<node*> GetClusters(){ 
			vector<node*> clusters;
			for(int i = 0; i < _clusters.size(); i++){
				if(_clusters[i] == nullptr) continue;
				clusters.push_back(_clusters[i]);
			} 
			return clusters;
		}
		
		void Insert(node* x){
			_clusters.push_back(x);
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

		int GetNClusters(){ 			
			vector<node*> clusters;
			for(int i = 0; i < _clusters.size(); i++){
				if(_clusters[i] == nullptr) continue;
				clusters.push_back(_clusters[i]);
			}
		return (int)clusters.size(); }	
			

		void SetSubclusterAlpha(double a){ _emAlpha = a; }		

		void SetDataSmear(const Matrix& cov){ _data_smear = cov; }

		void SetVerbosity(int v){ _verb = v; }

		void SetPriorParameters(map<string, Matrix> params){ _params = params; }

		void SetDistanceConstraint(double a, double b){ _constraint = true; _constraint_a = a; _constraint_b = b; }
		//uses a triangular distribution to constraint certain clusterings
		double DistanceConstraint(node* i, node* j){
			double c = (_constraint_a + _constraint_b)/2.;
			double pi = acos(-1);
			//phi
			double cent1 = i->points->Centroid(1);
			double cent2 = j->points->Centroid(1);
			
			//transform deltaPhi to be on [0,pi], wrapped s.t. 0 is close to 2pi (-3 close to 3)
			double d = fabs(cent1 - cent2);
			if(_wraparound){
				if(d > pi) d = 2*pi - d; 
			}	
		
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
		
			return sqrt(theta*phi);//tri->Prob(d)/tri->Prob(c);
		
		}

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



	protected:
		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
			int k;
  			//if leaf node (ie r == _z && l == _z) -> set k = 1
			if(x->l == _z && x->r == _z) k = 1;
			//number of clusters in node x = k_l + k_r for left and right nodes
			else k = x->l->model->GetNClusters() + x->r->model->GetNClusters();
			//center points in this cluster (bucket)
			//this accounts for phi wraparound
			Point transf = Point();
			if(_wraparound) transf = x->points->Center();
			
			//cout << "transformed points" << endl;
			//x->points->Print();
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
			//x - avg -> x - avg + avg = x
			//transform back relevant parameters in each subcluster (just centers - matrices unaffected)
			if(_wraparound){
				Matrix mu, mean;
				map<string, Matrix> params;
				//for(int d = 0; d < transf.Dim(); d++) transf.SetValue(-transf.at(d),d);
				for(int k = 0; k < x->model->GetNClusters(); k++){
					params = x->model->GetPriorParameters(k);
					mu = params["mean"];
					mean = params["m"];


					mu.add(Matrix(transf));
					x->model->GetModel(k)->SetParameter("mean",mu);
	
					mean.add(Matrix(transf));
					x->model->GetModel(k)->GetPrior()->SetParameter("mean",mean);
				}
				for(int d = 0; d < transf.Dim(); d++) transf.SetValue(-transf.at(d),d);
				//transform data back
				x->points->Translate(transf); 
			}

			return newLogL;
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
		bool _wraparound, _constraint;
};
#endif
