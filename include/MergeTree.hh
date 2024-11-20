#ifndef MergeTree_HH
#define MergeTree_HH

#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "HierGaussianMixture.hh"
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

		virtual ~MergeTree(){ 
			for(auto n : _clusters) delete n;
			_clusters.clear();
		}

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

		void Remove(node *l){
			//remove given node and associated mirror node if exists
			//also remove mirror node
			if(l->mirror != nullptr) _clusters[l->mirror->idx] = NULL;
			//setting the node to null matches the implementation in the FastJet code (see DnnPlane)
			_clusters[l->idx] = NULL;
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
		//set a weight-dependent resolution smearing for dimension d
		//ie energy dependent time resolution smearing
		//res = c^2 + n^2/w^2
		void SetWeightBasedResSmear(double c, double n, int d){
			_res_smear_c.push_back(c); _res_smear_n.push_back(n); _res_smear_d.push_back(d);
		} 

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
			double dist2d = _euclidean_2d(i, j);
			if(dist2d < i->nndist) i->nndist = dist2d;
			if(dist2d < j->nndist) j->nndist = dist2d;	

	
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

		void AddLeaf(const BayesPoint* pt = nullptr){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			_clusters.push_back(nullptr);
			node* x = (node*)malloc(sizeof *x);
			x->l = _z; x->r = _z;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			x->val = 1.;	
			x->d = _alpha; 
			x->model = nullptr;
			if(pt != nullptr) x->points = new PointCollection(*pt);
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)		
			int n = _clusters.size();
			x->idx = n-1;
			x->nndist = 1e300;
			x->mirror = nullptr;
			x->prob_tk = exp(Evidence(x)); //Evidence = ELBO \approx log(LH)
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
		

			//transform points into local coordinates
			//for GMM parameter estimation
			//use weighted mean as center to be set to 0 point
			//Point center({x->points->Centroid(0), x->points->Centroid(1), x->points->Centroid(2)});
			//copy points for parameter estimation
			//so original points don't get overwritten
			PointCollection newpts = PointCollection(*x->points);
			BayesPoint center = newpts.Center();//Translate(center);

			x->model = new HierGaussianMixture(k); //p(x | theta)
			if(_verb != 0) x->model->SetVerbosity(_verb-1);
			x->model->SetAlpha(_emAlpha);
			x->model->SetData(&newpts);
			if(!_data_smear.empty()){
				x->model->SetDataSmear(_data_smear);
				if(_res_smear_d.size() > 0){
					for(int i = 0; i < _res_smear_d.size(); i++){
						x->model->SetWeightBasedResSmear(_res_smear_c[i], _res_smear_n[i], _res_smear_d[i]);
					}
				}
			}
			x->model->InitParameters();
			x->model->InitPriorParameters();

			if(!_params.empty()) x->model->SetPriorParameters(_params);

			cout << "EVIDENCE FOR NODE " << x->idx << " WITH " << x->model->GetData()->GetNPoints() << " POINTS AND " << k << " clusters";
			if(x->mirror != nullptr) cout << " and mirror pt " << x->mirror->idx << endl;
			else cout << endl;

			VarEMCluster* algo = new VarEMCluster(x->model, k);
			//make sure sum of weights for points is above threshold
			//otherwise the algorithm will exclude all components	
			if(x->points->Sumw() >= _thresh) algo->SetThresh(_thresh);
			//cluster
			double oldLogL = algo->EvalLogL();
			double LogLThresh = 1e-10;
			double newLogL, entropy, nll;
			double dLogL = 999; 
			int it = 0;
			while(dLogL > LogLThresh){
				newLogL = algo->Cluster();
				entropy = x->model->EvalEntropyTerms();
				nll = x->model->EvalNLLTerms();
			cout << "it " << it << " newLogL " << newLogL << " entropy " << entropy << " NLL " << nll << endl;
		if(std::isnan(newLogL)) cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << " old ELBO: " << oldLogL << " new ELBO: " << newLogL << endl;
				dLogL = newLogL - oldLogL;
				oldLogL = newLogL;
				it++;
			}

			//translate the parameters back into global coordinates
			Matrix matshift;
			matshift.PointToMat(center);
			Matrix mean, priormean;
			map<string, Matrix> params;
			for(int k = 0; k < x->model->GetNClusters(); k++){
				//only need to do this once for prior
				if(k == 0){
					params = x->model->GetOnlyPriorParameters(k);
					priormean = params["mean"];
					priormean.add(matshift);
					params["m"] = priormean;
					x->model->SetPriorParameters(params);
				}
				params = x->model->GetPriorParameters(k);
				//only need to translate mean + mean on prior - stddev doesn't change
				mean = params["mean"];
				mean.add(matshift);
				params["mean"] = mean;
				x->model->GetModel(k)->SetParameter("mean",mean);
			}
			//reset data to original points
			x->model->SetData(x->points);

			return newLogL;
		}



		double _euclidean_2d(node* i, node* j){
			double deta = i->points->mean().at(0) - j->points->mean().at(0);
			double dphi = i->points->mean().at(1) - j->points->mean().at(1);

			return sqrt(deta*deta + dphi*dphi);

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
		vector<double> _res_smear_c, _res_smear_n;
		vector<int> _res_smear_d;
		int _verb;
		map<string, Matrix> _params;
		bool _constraint;
		void _remap_phi(PointCollection& points);
};
#endif
