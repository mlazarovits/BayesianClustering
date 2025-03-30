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
			_cell = acos(-1)/180; //default is CMS ECAL cell size
			_tresCte = 0.133913;
			_tresStoch = 1.60666; 
			_tresNoise = 0.00691415;
		}

		MergeTree(double alpha){
			_alpha = alpha;
			_thresh = 1.; _verb = 0;
			_constraint_a = 0; _constraint_b = acos(-1)/2.;
			_constraint = false;
			_cell = acos(-1)/180; //default is CMS ECAL cell size
			_tresCte = 0.133913;
			_tresStoch = 1.60666; 
			_tresNoise = 0.00691415;
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
			_cell = tree._cell; //default is CMS ECAL cell size
			_tresCte = tree._tresCte;
			_tresStoch = tree._tresStoch; 
			_tresNoise = tree._tresNoise; 
			//above tres params are for gev = 1
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

		double _cell, _tresCte, _tresStoch, _tresNoise;
		void SetMeasErrParams(double spatial, double tresCte, double tresStoch, double tresNoise){
			_cell = spatial;
			_tresCte = tresCte;
			_tresStoch = tresStoch;
			_tresNoise = tresNoise;
		}

	protected:

		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
			int k;
  			//if leaf node (ie r == _z && l == _z) -> set k = 1
			if(x->l == _z && x->r == _z) k = 1; 
			//number of clusters in node x = k_l + k_r for left and right nodes
			else//{
			k = x->l->model->GetNClusters() + x->r->model->GetNClusters();
			//k = x->l->model->GetData()->GetNPoints() + x->r->model->GetData()->GetNPoints();
			//cout << "k L " <<  x->l->model->GetNClusters() << " k R " << x->r->model->GetNClusters() << " k " << k << endl;}
 //cout << "og points" << endl; x->points->Print();
			//cout << "original points" << endl; x->points->Print();
			//do phi wraparound? may already be taken care of in local coords + mirror pts
			PointCollection* newpts = new PointCollection(*x->points);
			//cout << "pts " << endl; newpts->Print();
			BayesPoint center({newpts->Centroid(0), newpts->CircularCentroid(1), newpts->Centroid(2)});
			//cout << "center" << endl; center.Print();

			//scale points s.t. 1 cell ~ 0.0174 = 1 unit in eta-phi
			//x'' = x'/b = (x-a)/b
			//sets relative importance of dimensions
			//decreasing cell -> eta/phi distance more important
			//increasing entry (2,2) -> time distance more important
			//TODO: make configurable externally - or not?? kind of integral to algorithm
			double cell = 4*atan(1)/180;
			Matrix Rscale(3,3);
			Rscale.SetEntry(1/cell,0,0);
			Rscale.SetEntry(1/cell,1,1);
			Rscale.SetEntry(1e-1,2,2); //what should be considered "out of time" for a jet? 1 ns? 2 ns? 10 ns? 
			//for right now, 1e-1 is good, nominal time is in ns, ns*1e-3 = us, ns*1e3 = ps



			x->model = new GaussianMixture(k); //p(x | theta)
			if(_verb != 0) x->model->SetVerbosity(_verb-1);
			x->model->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				x->model->SetDataSmear(_data_smear);
			}

			x->model->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
			
			//in local space, circular coordinates (like phi) can go negative
			x->model->SetData(newpts); //may need to make copy of de-referenced object so as not to change the original points	
			x->model->ShiftData(center);
			//cout << "centroid " << endl; 
			
			//cout << "translated pts" << endl;
			//x->model->GetData()->Print();
			//cout << "scale data + lam*s" << endl;	
			x->model->ScaleData(Rscale);
			
			//cout << "transformed points" << endl;
			//x->model->GetData()->Print();
			//cout << endl;
	
			x->model->InitParameters();
			x->model->InitPriorParameters();

			if(!_params.empty()) x->model->SetPriorParameters(_params);
			//inverse transformations
			Matrix RscaleInv;
			RscaleInv.invert(Rscale);
			center.Scale(-1);


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
			//cout << "it " << it << " newLogL " << newLogL << " entropy " << entropy << " NLL " << nll << endl;
		if(std::isnan(newLogL)) cout << "iteration #" << it << " log-likelihood: " << newLogL << " dLogL: " << dLogL << " old ELBO: " << oldLogL << " new ELBO: " << newLogL << endl;
				dLogL = newLogL - oldLogL;
				oldLogL = newLogL;
				it++;
			}
			//cout << "EVIDENCE FOR NODE " << x->idx << " WITH " << x->model->GetData()->GetNPoints() << " POINTS AND " << k << " max clusters and " << x->model->GetNClusters() << " found clusters - evidence " << exp(newLogL) << " ELBO " << newLogL << " with points in node " << endl; x->points->Print();  

			//transform the parameters back into global coordinates
			//need to unscale first then uncenter since x'' = (x-a)/b (see above)
			//need to unscale data - also unscales lamStar measurement errors 
			//cout << "inverse scale data + lam*s" << endl;	
			x->model->ScaleData(RscaleInv);	
			//cout << "unscaled points" << endl;
			//x->model->GetData()->Print();
			//need to unscale mean + covariances
			x->model->ScaleParameters(RscaleInv);	
			
			//need to put GMM parameters AND points back in detector coordinates (not local)
			//cout << "ismirror " << x->ismirror << " left " << x->l->ismirror << " right " << x->r->ismirror << endl;
			//x->model->ShiftData(center);
			//cout << "untransformed points" << endl;
			//resets data to original points
			x->model->SetData(x->points);
			//to consider: keeping the data in the model as the transformed points that the algorithm actually runs on and the points in the node the original ones in the detector system
			//x->model->GetData()->Print();
			//cout << "end evidence" << endl;
			x->model->ShiftParameters(center);

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
		int _verb;
		map<string, Matrix> _params;
		bool _constraint;
		void _remap_phi(PointCollection& points);
};
#endif
