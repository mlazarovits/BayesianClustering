#ifndef MergeTree_HH
#define MergeTree_HH
#include "BaseTree.hh"
#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "Gaussian.hh"

using std::map;

class MergeTree : BaseTree{
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW //to ensure proper alignment for heap allocation and STL containers
		MergeTree() :
			BaseTree(), 
			_alpha(0),
			_thresh(0.), 
			_verb(0),
			_cell(acos(-1)/180), //default is CMS ECAL cell size
			_tresCte(0.133913),
			_tresStoch(1.60666), 
			_tresNoise(0.00691415),
			_check_merges(false),
			_nGhosts(0),
			_Rscale(3,3),
			_RscaleInv(3,3)
		{ 
			//x'' = x'/b = (x-a)/b
			//sets relative importance of dimensions
			//decreasing cell -> eta/phi distance more important
			//increasing entry (2,2) -> time distance more important
			_Rscale.SetEntry(1/_cell,0,0);
			_Rscale.SetEntry(1/_cell,1,1);
			_Rscale.SetEntry(1,2,2); 
			_RscaleInv.invert(_Rscale);
		}

		MergeTree(double alpha) : 
			BaseTree(), 
			_alpha(alpha),
			_thresh(0.1), 
			_verb(0),
			_cell(acos(-1)/180), //default is CMS ECAL cell size
			_tresCte(0.133913),
			_tresStoch(1.60666), 
			_tresNoise(0.00691415),
			_check_merges(false),
			_nGhosts(0),
			_Rscale(3,3),
			_RscaleInv(3,3)
		{
			//x'' = x'/b = (x-a)/b
			//sets relative importance of dimensions
			//decreasing cell -> eta/phi distance more important
			//increasing entry (2,2) -> time distance more important
			_Rscale.SetEntry(1/_cell,0,0);
			_Rscale.SetEntry(1/_cell,1,1);
			_Rscale.SetEntry(1,2,2); 
			_RscaleInv.invert(_Rscale);
		}

		//copy constructor
		MergeTree(const MergeTree& tree) :
			BaseTree(tree),
			_alpha(tree._alpha),
			_thresh(tree._thresh),
			_verb(tree._verb),
			_cell(tree._cell), //default is CMS ECAL cell size
			_tresCte(tree._tresCte),
			_tresStoch(tree._tresStoch), 
			_tresNoise(tree._tresNoise), 
			_check_merges(tree._check_merges),
			_nGhosts(tree._nGhosts),
			//matrices
			_Rscale(tree._Rscale),
			_RscaleInv(tree._RscaleInv)
			//pointers
			//_clusters()	
		{
			//for(auto& c : tree._clusters) _clusters.push_back(make_unique<node>(*c));
			//above tres params are for gev = 1
		}

		virtual ~MergeTree(){ 
		}
		
		void avg_time(){
			cout << "Average CalcMerge time " << _total_calcmerge_time / _n_calcmerge_calls << " seconds with " << _n_calcmerge_calls << " calls and " << _total_calcmerge_time << " seconds total time" << endl;
			cout << "Average Evidence time " << _total_evidence_time / _n_evidence_calls << " seconds with " << _n_evidence_calls << " calls and " << _total_evidence_time << " seconds total time" << endl;
		}

		void SetThresh(double t){ _thresh = t; }

		void AddData(PointCollection* pc){
		//sort nodes of merge tree once here then add nodes to search tree and merge tree (as leaves)	
			for(int i = 0; i < pc->GetNPoints(); i++){
				AddLeaf(pc->at(i));
			}
		}


		void CheckMerges(bool t){ _check_merges = t; if(_check_merges) cout << "Checking merged models" << endl;}	
	
		//node* Get(int i){ return _clusters[i].get(); }

		//relinquishes ownership of nodes to caller
		//vector<std::shared_ptr<node>> GetClusters(){ 
		//	return std::move(_clusters);
		//}
	
		//make sure children of x are set before Insert	
		//void Insert(std::shared_ptr<node> x){
		////cout << "inserting node with ismirror " << x->ismirror << " and pts " << endl; x->points->Print();
		//	_clusters.push_back(nullptr);
		//	x->idx = (int)_clusters.size() - 1;
		//	_clusters[(int)_clusters.size() - 1] = std::move(x);
		//}

		//assuming Dirichlet Process Model (sets priors)
		std::shared_ptr<node> CalculateMerge(node *l, node* r);
		double CalculateMerge(int i, int j);

		//void Remove(node* l){
		//	//remove given node and associated mirror node if exists
		//	//also remove mirror node
		//	if(l->mirror != nullptr) _clusters[l->mirror->idx] = NULL;
		//	//setting the node to null matches the implementation in the FastJet code (see DnnPlane)
		//	_clusters[l->idx] = NULL;
		//}

		//void Remove(int i){ _clusters[i] = NULL; }

		void SetAlpha(double alpha){ _alpha = alpha; }	

		//counts clusters that have been removed
		//use this function when looping over clusters s.t. the element indices line up
		//int GetNAllClusters(){ 			
		//	return (int)_clusters.size();
		//}	

		//this does not count removed clusters
		//use this function to get number of actual clusters
		//when no merges have been made (ie clusters are only leaves)
		//this should match GetNAllClusters 
		//int GetNClusters(){ 			
		//	int n = 0;
		//	for(int i = 0; i < _clusters.size(); i++){
		//		if(_clusters[i] == nullptr) continue;
		//		n++;
		//	}
		//	return n;
		//}	
			

		void SetSubclusterAlpha(double a){ _emAlpha = a; }		

		void SetDataSmear(const Matrix& cov){ _data_smear = cov; }

		void SetVerbosity(int v){ _verb = v; }

		void SetPriorParameters(const map<string, Matrix>& params){ _params = params; }


		std::shared_ptr<node> MakeLeaf(const BayesPoint& pt){
			//cout << "MergeTree::MakeLeaf - start" << endl;
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			//_clusters.push_back(nullptr);
			std::shared_ptr<node> x = std::make_shared<node>();
			x->l = nullptr; x->r = nullptr;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			//x->val = 1.;	
			//x->d = _alpha;
			x->log_d = log(_alpha);
			x->model = nullptr;
			x->points = std::make_unique<PointCollection>(pt);//std::make_unique<PointCollection>(*pt);
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)		
			//int n = _clusters.size();
			x->idx = -1;
			x->nndist = 1e300;
			x->mirror = nullptr;
			x->ismirror = false;
			double elbo = Evidence(x.get());
			//cpp_bin_float_100 elbo_100 = cpp_bin_float_100(elbo);
			//x->prob_tk = exp(elbo);//p_dk_tk = p_dk_h1 since cannot be divided further
			x->log_prob_tk = elbo;
			//cout << "adding node with evidence " << Evidence(x) << " and weight " << pt->w() << endl;
			//Evidence = ELBO \approx log(LH)
			//_clusters[n-1] = std::move(x);
			//cout << "adding leaf with ismirror " << x->ismirror << endl;
			//cout << "MergeTree::MakeLeaf - end" << endl;
			return x;
		}


		void AddLeaf(const BayesPoint& pt){
			if(_alpha == 0) cout << "MergeTree - need to set alpha" << endl;
			//_clusters.push_back(nullptr);
			unique_ptr<node> x = make_unique<node>();
			x->l = nullptr; x->r = nullptr;
			//////if leaf -> p(Dk | Tk) = p(Dk | H1k) => rk = 1
			//x->val = 1.;	
			//x->d = _alpha;
			x->log_d = log(_alpha);
			x->model = nullptr;
			x->points = std::make_unique<PointCollection>(pt);//std::make_unique<PointCollection>(*pt);
			//initialize probability of subtree to be null hypothesis for leaf
			//p(D_k | T_k) = p(D_k | H_1^k)		
			//int n = _clusters.size();
			x->idx = -1;
			x->nndist = 1e300;
			x->mirror = nullptr;
			x->ismirror = false;
			double elbo = Evidence(x.get());
			//cpp_bin_float_100 elbo_100 = cpp_bin_float_100(elbo);
			//x->prob_tk = exp(elbo);//p_dk_tk = p_dk_h1 since cannot be divided further
			x->log_prob_tk = elbo;
			//cout << "adding node with evidence " << Evidence(x) << " and weight " << pt->w() << endl;
			//Evidence = ELBO \approx log(LH)
			//_clusters[n-1] = std::move(x);
			//cout << "adding leaf with ismirror " << x->ismirror << endl;
		}


		//note: this ONLY gets called in the N^2 clustering
		//for the delauney (NlnN) clustering this is taken care of
		//by Dnn2piCylinder
		void CreateMirrorNode(std::shared_ptr<node> x);

		double _cell, _tresCte, _tresStoch, _tresNoise;
		void SetMeasErrParams(double spatial, double tresCte, double tresStoch, double tresNoise){
			_cell = spatial;
			_tresCte = tresCte;
			_tresStoch = tresStoch;
			_tresNoise = tresNoise;
		}
		void SetNGhosts(int n){ _nGhosts = n; }



		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
			if(_verb > 5) cout << "MergeTree::Evidence - start" << endl;
			clock_t t = clock();
			int k;
			vector<map<string,Matrix>> prev_posts = {};
			map<int,int> left_post_indices, right_post_indices;
  			//if leaf node (ie r == nullptr && l == nullptr) -> set k = 1
			if(x->points->GetNPoints() == 1 || x->points->Sumw() < _thresh){ 
				//cout << "leaf nodes - setting max # clusters == 1" << endl;
		 		k = 1;
			}
			//if the sum of all weights is not enough for 1 subcluster to be "above threshold" they should not be able to be broken into multiple subclusters
			//number of clusters in node x = k_l + k_r for left and right nodes
			else{
				//cout << "point a1" << endl;
				//cout << "left null? " << (x->l == nullptr) << " right null? " << (x->r == nullptr) << endl;
				//cout << "left model null? " << (x->l->model == nullptr) << endl;
				//cout << "left points" << endl; x->l->points->Print();
				//cout << " right model null? " << (x->r->model == nullptr) << endl;
				//cout << "right points" << endl; x->r->points->Print();
				//cout << "total # pts " << x->points->GetNPoints() << endl; x->points->Print();	
				int mincls = x->l->model->GetNClusters() + x->r->model->GetNClusters();
				int npts = x->points->GetNPoints();//x->l->model->GetData()->GetNPoints() + x->r->model->GetData()->GetNPoints();
				k = npts < mincls ? npts : mincls;
				//cout << "not leaf nodes - mincls " << mincls << " npts " << npts << " k " << k << endl;
				//k = npts;
				//k = x->l->model->GetNClusters() + x->r->model->GetNClusters();
				//setting initial posterior parameters from models of previous steps
				vector<double> nnorms; x->l->model->GetNorms(nnorms);
				for(int kk = 0; kk < x->l->model->GetNClusters(); kk++){
					//cout << "left cluster #" << kk << endl;
					prev_posts.emplace_back(map<string, Matrix>());
					x->l->model->GetLHPosteriorParameters(kk,prev_posts[prev_posts.size()-1]);
					if(_verb > 5){
						cout << "left node subcluster # " <<kk << endl;
				       		cout << "left points " << x->l->points->GetNPoints() << " sumw " << x->l->points->Sumw() << endl;
						cout << "with left mean " << endl;prev_posts[prev_posts.size()-1].at("mean").Print();
						cout << "and left norm " << nnorms[kk] << endl;

					}
					left_post_indices[kk] = prev_posts.size()-1;
				} 
				x->r->model->GetNorms(nnorms);
				for(int kk = 0; kk < x->r->model->GetNClusters(); kk++){
					//cout << "right cluster #" << kk << endl;
					prev_posts.emplace_back(map<string, Matrix>());
					x->r->model->GetLHPosteriorParameters(kk, prev_posts[prev_posts.size()-1]);
					if(_verb > 5){
						cout << "right node subcluster # " <<kk << endl; 
				       		cout << "right points" << endl; x->l->points->Print();	
						cout << "with right mean " << endl;prev_posts[prev_posts.size()-1].at("mean").Print();
						cout << "and right norm " << nnorms[kk] << endl;
					}
					right_post_indices[kk] = prev_posts.size()-1;
				}

			}
			x->points->Sort();
			auto newpts = make_unique<PointCollection>(*x->points);
			if(_verb > 5){ cout << newpts->GetNPoints() << " original pts " << endl; newpts->Print();}
			//sort points so random initialization is consistent based on seed
		
			x->model = std::make_unique<GaussianMixture>(k); //p(x | theta)
			if(_verb != 0) x->model->SetVerbosity(_verb-1);
			x->model->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				x->model->SetDataSmear(_data_smear);
			}

			x->model->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
		
			//in local space, circular coordinates (like phi) can go negative
			x->model->SetData(std::move(newpts)); //need to make copy of de-referenced object so as not to change the original point and maintain correct ownership	
			x->model->PutPhi02pi();		


			//change eta to theta BEFORE calculating centroid
			x->model->EtaToTheta();
			if(_verb > 5){cout << "eta to theta" << endl; x->model->GetData()->Print();}
			//make sure all data in model is on [0,2pi] to begin with
			//some mirror nodes will have pts < 0 or > 2 pi
			//but the shifting + plane project depends on all pts initially being on unit circle
			//it doesn't matter if this is done before or after centroid calculation
			//bc of the way the centroid is calculated for phi (ie CircularCentroid)
		

			//since the first dimension is now an angle, its centroid needs to be circularly calculated
			//but since eta is only defined for theta on [-pi/2,pi/2], it shouldn't make that much of a difference
			BayesPoint center({x->model->GetData()->CircularCentroid(0), x->model->GetData()->CircularCentroid(1), x->model->GetData()->Centroid(2)});
			if(_verb > 5){ cout << "center" << endl; center.Print();}
			x->model->ShiftData(center);


			if(_verb > 5){ cout << "translated pts" << endl; x->model->GetData()->Print(); }
			//project phi onto plane
			bool phiInf = x->model->ProjectPhi();
			bool thetaInf = x->model->ProjectTheta();
			if(_verb > 5){ cout << "projected pts" << endl; x->model->GetData()->Print();}
			//check for infinities in eta + phi from measurement error
			if(phiInf || thetaInf){
				if(_verb > 5) cout << "found inf in pts, returning elbo as " << -1e308 << endl;
				//reset model data for further merges
				x->model->SetTestIndex(69);
				x->model->SetData(std::move(make_unique<PointCollection>(*x->points)));
				return -1e308;
			}


			//consider this: if inserting numerical 'infinities' causes numerical instabilities in the model, (especially in the matrix inversions)
			//derive from the ELBO equations that large infinities originating from the plane projection
			//should lead to large negative ELBOs (ie not a good lower bound on the LH)
			//and are probabilistically disfavored
			//so, due to aforementioned numerical instabilities, these quantities are not explicitly calculated
			//and are automatically excluded from merging (ie the ELBO is set to a large negative number
			//s.t. the merge is never chosen) 

			x->model->ScaleData(_Rscale);
			if(_verb > 5){ cout << "scaled pts" << endl; x->model->GetData()->Print();}
		
			//put means from previous step in same transformed space as data	
			for(int kk = 0; kk < prev_posts.size(); kk++){
				auto params = prev_posts[kk];
				Matrix mean = params["m"];
				TransformMean(mean, center);
				prev_posts[kk]["m"] = mean;
				//do for xbar too
				mean = params["xbar"];
				TransformMean(mean, center);
				prev_posts[kk]["xbar"] = mean;

				//don't have to do for model->GetParameter("mean") bc this isn't used in the initial logLH eval
			} 
			//needs to be done after data is set + transformed bc the initialization procedure depends on data
			x->model->InitParameters(_params,prev_posts);
	
			//calculate list of inner products of subcluster PDFs to compare for merges
			//only do if for # subclusters in node l = n and # subclusters in node l = m
			//n > 1 && m > 1 => n + m > 2
			map<double, pair<int, int>, std::greater<double>> prodmap;
			if(_check_merges && x->l != nullptr && x->r != nullptr && x->model->GetNClusters() > 1){
				//construct PDF inner product map to pairs of PDFs
				//get vectors of PDFs from left (first) then right (second)
				vector<BasePDF*> l_pdfs, r_pdfs;
				for(int n = 0; n < x->l->model->GetNClusters(); n++){
					l_pdfs.push_back(x->l->model->GetModel(n));
				}
				
				for(int n = 0; n < x->r->model->GetNClusters(); n++){
					r_pdfs.push_back(x->r->model->GetModel(n));
				}
				_match_pdfs(l_pdfs, r_pdfs, prodmap);
			}
			//nominal GMM (no merges)
			auto algo = std::make_unique<VarEMCluster>(x->model.get(), k);
			double thresh = 1e-1;//0.01*x->points->Sumw(); //1% of total weight
			algo->SetThresh(thresh);
			//cluster
			double oldLogL = algo->EvalLogL();
			double LogLThresh = 1e-3;
			double newLogL = 0;
			double entropy = 0;
			double nll = 0;
			double dLogL = 1e308; 
			int it = -1;
			while(dLogL > LogLThresh*fabs(oldLogL) || it == 0){
				it++;
				newLogL = algo->Cluster();
				//ELBO should maximizing LH -> therefore newLogL > oldLogL if both are < 0	
				dLogL = newLogL - oldLogL;
				if(std::isnan(newLogL)) cout << std::setprecision(10) << "it " << it << " new logl " << newLogL << " oldlogl " << oldLogL << " dlogl " << dLogL << " # clusters " << x->model->GetNClusters() << endl;
				oldLogL = newLogL;
			}
			if(_verb > 5) cout << "EVIDENCE FOR NODE " << x->idx << " WITH " << x->model->GetData()->GetNPoints() << " POINTS AND " << k << " max clusters and " << x->model->GetNClusters() << " found clusters - evidence " << exp(newLogL) << " ELBO " << newLogL << endl;
			//compare nominal model to merges IN ORDER and with memory (ie if merge1 > nom -> check merge1 & merge2)
			if(_check_merges && x->l != nullptr && x->r != nullptr && x->model->GetNClusters() > 1){
				double merge_elbo;
				//compare merges in its respective projected data space
				vector<map<string, Matrix>> merge_posteriors;
				vector<pair<int,int>> merge_pairs;
				for(auto prodit = prodmap.begin(); prodit != prodmap.end(); prodit++){
					if(prodit->second.first == -1 || prodit->second.second == -1) continue;
					merge_pairs.push_back(std::make_pair(left_post_indices[prodit->second.first], right_post_indices[prodit->second.second]));
				}
				//this function compares all the "local" (ie projected) subcluster merges to not merging, then sets the starting params for the merged model from these merge decisions
				unique_ptr<GaussianMixture> merged_model = _merge_model(x, prev_posts, merge_pairs, merge_elbo);
				if(merge_elbo > newLogL){
					if(_verb > 1) cout << "replacing nominal model with " << x->model->GetNClusters() << " subclusters and elbo = " << newLogL << " with merged model with " << merged_model->GetNClusters() << " subclusters and elbo = " << merge_elbo << endl;
					x->model = std::move(merged_model);
					newLogL = merge_elbo;	
				}
				if(_verb > 1) cout << endl;
			}	
			

			//inverse transformations - for after fitting
			center.Scale(-1);
			//transform the parameters back into global coordinates
			//need to unscale first then uncenter since x'' = (x-a)/b (see above)
			//need to unscale data - also unscales lamStar measurement errors 
			
			x->model->ScaleData(_RscaleInv);	 
			
			//need to unscale mean + covariances
			x->model->ScaleParameters(_RscaleInv);

			//need to put GMM parameters AND points back in detector coordinates (not local)

			//only does for mean rn - not sure how to do for cov...
			x->model->UnprojectPhi_params();
			x->model->UnprojectTheta_params();
      			x->model->ShiftParameters(center);

			x->model->PutPhi02pi_params(); //does for data and parameters - do after data + parameters shift so the [0,2pi] transformation doesn't get shifted
			x->model->ThetaToEta_params();
			bool isnan = false;
			for(int k = 0; k < x->model->GetNClusters(); k++){
				map<string, Matrix> params;
				x->model->GetLHPosteriorParameters(k, params);
				//if( (0 < phi && phi < 1.57) || (4.7 < phi && phi < 2*acos(-1))){
				//	cout << "final center #" << k << endl; params["m"].Print(); 
				//}
				if(std::isnan(params["mean"].at(0,0))){
					isnan = true;
				}
			}

			if(isnan){
				for(int k = 0; k < x->model->GetNClusters(); k++){
					map<string, Matrix> params;
					x->model->GetLHPosteriorParameters(k, params);
					cout << "theta to eta params" << endl;
					cout << "cluster #" << k << " with weight " << params["alpha"].at(0,0) - _emAlpha << endl;
					cout << "cluster #" << k << endl;
					cout << "mean" << endl;
					params["mean"].Print();
				}
			}
			//resets data to original points
			//x->model->SetData(x->points.get());
			//x->model->SetData(x->points);
			t = clock() - t;
			_total_evidence_time += (double)t/CLOCKS_PER_SEC;
			_n_evidence_calls++;
			if(_verb > 5) cout << "MergeTree::Evidence - end" << endl;
			return newLogL;
		}
		
	protected:
		void TransformMean(Matrix& mat, BayesPoint& center){
				//eta to theta
				double eta = mat.at(0,0);
				mat.SetEntry(2*atan(exp(-eta)),0,0);
				
				//put 02pi
				PointCollection mat_pt; mat.MatToPoints(mat_pt);
				mat_pt.Put02pi(1);	

				//shift
				mat_pt.CircularTranslate(center.at(0),0);
				mat_pt.CircularTranslate(center.at(1),1);
				mat_pt.Translate(center.at(2),2);
				//project phi
				mat_pt.AngleToPlaneProject(1);	

				//project eta
				mat_pt.AngleToPlaneProject(0);	
				mat = Matrix(mat_pt);

				//scale
				mat.mult(_Rscale,mat);	
		}



		double _euclidean_2d(node* i, node* j){
			double deta = i->points->mean().at(0) - j->points->mean().at(0);
			double dphi = i->points->CircularMean(1) - j->points->CircularMean(1);
			dphi = acos(cos(dphi));

			return sqrt(deta*deta + dphi*dphi);

		}
		double _euclidean_3d(node* i, node* j){
			double deta = i->points->mean().at(0) - j->points->mean().at(0);
			double dphi = i->points->CircularMean(1) - j->points->CircularMean(1);
			dphi = acos(cos(dphi));
			double dtime = i->points->mean().at(2) - j->points->mean().at(2);
			//cout << "deta " << deta << " dphi " << dphi << " dtime " << dtime << endl;
			return sqrt(deta*deta + dphi*dphi + dtime*dtime);

		}
		
		double _euclidean_3d_fromCentroid(node* x){
			if(x->l == nullptr && x->r == nullptr) return 0;
			BayesPoint center({x->points->Centroid(0), x->points->CircularCentroid(1), x->points->Centroid(2)});
			double deta = x->l->points->mean().at(0) - center.at(0);
			double dphi = x->l->points->CircularMean(1) - center.at(1);
			dphi = acos(cos(dphi));
			double dtime = x->l->points->mean().at(2) - center.at(2);
			//cout << "l: deta " << deta << " dphi " << dphi << " dtime " << dtime << endl;
			deta = x->r->points->mean().at(0) - center.at(0);
			dphi = x->r->points->CircularMean(1) - center.at(1);
			dphi = acos(cos(dphi));
			dtime = x->r->points->mean().at(2) - center.at(2);
			//cout << "r: deta " << deta << " dphi " << dphi << " dtime " << dtime << endl;
			return 0;//sqrt(deta*deta + dphi*dphi + dtime*dtime);

		}

		void _run_model(GaussianMixture* model, double& elbo){
			auto algo = std::make_unique<VarEMCluster>(model, model->GetNClusters());
			double thresh = 1e-1;//0.01*model->GetData()->Sumw(); //1% of total weight
			algo->SetThresh(thresh);
			//cluster
			double oldLogL = algo->EvalLogL();
			double LogLThresh = 1e-3;
			double newLogL = 0;
			double dLogL = 1e308; 
			int it = -1;
			while(dLogL > LogLThresh*fabs(oldLogL) || it == 0){
				it++;
				newLogL = algo->Cluster();
				//ELBO should maximizing LH -> therefore newLogL > oldLogL if both are < 0	
				dLogL = newLogL - oldLogL;
		if(std::isnan(newLogL)) cout << std::setprecision(10) << "it " << it << " new logl " << newLogL << " oldlogl " << oldLogL << " dlogl " << dLogL << " # clusters " << model->GetNClusters() << endl;
				oldLogL = newLogL;
			}

			elbo = newLogL;

		}

		//project points using nominal model responsibilities to isolate subcluster merges
		//need to project data from each node for node subcluster
		unique_ptr<PointCollection> _project_data(node* x, int cl_left, int cl_right){
			Matrix r_nk_l; x->l->model->GetPosterior(r_nk_l);
			if(cl_left >= r_nk_l.nCols() || cl_left < 0){
				cout << "Error: accessing cluster " << cl_left << " with posterior of " << r_nk_l.nCols() << " # cols # clusters in this model " << x->l->model->GetNClusters() << endl;
				return nullptr;
			}
			Matrix r_nk_r; x->r->model->GetPosterior(r_nk_r);
			if(cl_right >= r_nk_r.nCols() || cl_right < 0){
				cout << "Error: accessing cluster " << cl_right << " with posterior of " << r_nk_r.nCols() << " # cols # clusters in this model " << x->r->model->GetNClusters() << endl;
				return nullptr;
			}
			auto newpts = make_unique<PointCollection>();
			for(int i = 0; i < x->model->GetData()->GetNPoints(); i++){
				double w_i = x->model->GetData()->at(i).w();
				bool matched = false;
				for(int ll = 0; ll < x->l->model->GetData()->GetNPoints(); ll++){
					if(x->l->model->GetData()->at(ll).w() == w_i){
						BayesPoint pt = x->model->GetData()->at(i); 
						pt.SetWeight(r_nk_l.at(ll,cl_left));
						newpts->AddPoint(pt);
						break;
					}
				}
				if(matched) continue;
				for(int rr = 0; rr < x->r->model->GetData()->GetNPoints(); rr++){
					if(x->r->model->GetData()->at(rr).w() == w_i){
						matched = true;
						//get unweighted responsibility
						BayesPoint pt = x->model->GetData()->at(i);
						pt.SetWeight(r_nk_r.at(rr,cl_right));
						newpts->AddPoint(pt);
						break;
					}
				}
	
			}
			newpts->Sort();
			return newpts;
		}



		//if passing all starting params, then merge_pair idxs are in "global" frame, need in local
		std::unique_ptr<GaussianMixture> _compare_projected_models(node* x, const vector<map<string,Matrix>>& starting_params, const pair<int, int>& merge_pair){
			//project data in x->l and x->r according to given merge pair (0: left, 1: right)
			int cl1 = merge_pair.first;
			int cl2 = merge_pair.second;

			//cout << "_compare_projected_models - comparing subcls " << cl1 << " and " << cl2 << " of " << starting_params.size() << " starting params " << endl;
			//cout << "cl1 local idx " << cl1 << " of " << x->l->model->GetNClusters() << " clusters and cl2 local idx " << cl2 - x->l->model->GetNClusters() << " of " << x->r->model->GetNClusters() << " clusters" << endl;


			//get centroid of starting params
			map<string, Matrix> new_cl;
			Matrix new_center(3,1);
			double w1 = starting_params[cl1].at("alpha").at(0,0) - _emAlpha;
			double w2 = starting_params[cl2].at("alpha").at(0,0) - _emAlpha;
			PointCollection cents, cents2; starting_params[cl1].at("m").MatToPoints(cents, {w1});
			starting_params[cl2].at("m").MatToPoints(cents2, {w2});
			cents += cents2;
			//centroid bw clusters 1 and 2 - parameters have been transformed + projected to be on a plane
			BayesPoint center({cents.Centroid(0), cents.Centroid(1), cents.Centroid(2)});
			new_cl["m"] = center;
			//need to also give Nk (norms) and cov to seed posteriors
			new_cl["alpha"] = Matrix(_emAlpha + (w1 + w2));
			//cov - average of two
			Matrix new_cov(3,3);
			Matrix cov1 = starting_params[cl1].at("cov");
			cov1.mult(cov1,w1);
			Matrix cov2 = starting_params[cl2].at("cov");
			cov2.mult(cov2,w2);
			new_cov.add(cov1,cov2);
			new_cov.mult(new_cov,1./(w1+w2));
			//mult by nu0 to undo this in InitParameters
			new_cov.mult(new_cov,_params.at("dof").at(0,0));
			//need to invert bc pulls as scalemat
			new_cov.invert(new_cov);
			new_cl["scalemat"] = new_cov;
			new_cl["scale"] = Matrix(_params.at("scale").at(0,0) + (w1 + w2));
			new_cl["dof"] = Matrix(_params.at("dof").at(0,0) + (w1 + w2));
			

			//project_data needs cl1 and cl2 in their "local" values
			unique_ptr<PointCollection> mergepts = _project_data(x,cl1,cl2 - x->l->model->GetNClusters());
			if(mergepts == nullptr) return nullptr;
			//cout << "MergeTree::_compare_project_models - merge data" << endl; mergepts->Print();
			//cout << "# total mergepts " << mergepts->GetNPoints() << endl;
			auto seppts = make_unique<PointCollection>(*mergepts);		
			//cout << "MergeTree::_compare_project_models - sep data" << endl; seppts->Print();
	
			//do merge model
			//run GMM with projected data with 2 initial starting subclusters (ie starting_params[merge_pair.first], starting_params[merge_pair.second])
			auto mergemodel = std::make_unique<GaussianMixture>(1);
			if(_verb != 0) mergemodel->SetVerbosity(_verb-1);
			mergemodel->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				mergemodel->SetDataSmear(_data_smear);
			}
			mergemodel->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
			mergemodel->SetData(std::move(mergepts));
			mergemodel->InitParameters(_params,{new_cl});
				// save ELBO
			double elbo_pair;
			_run_model(mergemodel.get(), elbo_pair);
	//cout << "elbo for proj merge " << elbo_pair << " # final clusters " << mergemodel->GetNClusters() << endl;	
	
			//do separate (original) model
			auto sepmodel = std::make_unique<GaussianMixture>(2);
			if(_verb != 0) sepmodel->SetVerbosity(_verb-1);
			sepmodel->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				sepmodel->SetDataSmear(_data_smear);
			}
			sepmodel->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
			sepmodel->SetData(std::move(seppts));
			sepmodel->InitParameters(_params,{starting_params[cl1], starting_params[cl2]});
				// save ELBO
			double elbo_sep;
			_run_model(sepmodel.get(), elbo_sep);
	//cout << "elbo for proj sep " << elbo_sep << " # final clusters " << sepmodel->GetNClusters() << endl;	
			//compare 2 subcluster to 1 subcluster model
			//take ending subclusters of better model to end_params
			//then will use end_params as starting point of full "merged" model to compare to full nominal
		
	//cout << "starting place for merge model" << endl; new_cl["m"].Print();
	//cout << "ending place for merge model with " << mergemodel->GetNClusters() << " subclusters" << endl; mergemodel->GetLHPosteriorParameters(0)["m"].Print();

	//cout << "starting places for sep model" << endl;
	//cout << "cluster #0" << endl; starting_params[cl1]["m"].Print();
	//cout << "cluster #1" << endl; starting_params[cl2]["m"].Print();
	//cout << "ending place for sep model with " << sepmodel->GetNClusters() << " subclusters " << endl; 
	//for(int i = 0; i < sepmodel->GetNClusters(); i++){
	//	cout << "cluster #" << i << endl;
	//	sepmodel->GetLHPosteriorParameters(i)["m"].Print();
	//}
			if(elbo_pair > elbo_sep)
				return mergemodel;
			else
				return sepmodel;

		}

		//will run model seeded from projected merges
		//need to pass full node so we can project the data from the potential merge (ie from the x->l and x->r nodes) by their own posterior responsibilities
		std::unique_ptr<GaussianMixture> _merge_model(node* x, const vector<map<string,Matrix>>& starting_params, const vector<pair<int, int>>& merge_pairs, double& merge_elbo){
			vector<map<string, Matrix>> merge_starting_params;
			//for each merge 
			for(int m = 0; m < merge_pairs.size(); m++){
				//do projected merge - return better model (merge or separate)
				//cout << "looking at merge for left subcl " << merge_pairs[m].first << " and right subcl " << merge_pairs[m].second << endl;
				std::unique_ptr<GaussianMixture> proj_model = _compare_projected_models(x, starting_params, merge_pairs[m]);
				if(proj_model == nullptr){
					cout << "Error: projected model null." << endl;
					break;
				}
				//cout << "# clusters in proj model " << proj_model->GetNClusters() << endl;
				//get final clusters of proj_model - these will be starting place for full merge model
				for(int k = 0; k < proj_model->GetNClusters(); k++){
					map<string, Matrix> params;
					proj_model->GetLHPosteriorParameters(k, params);
					merge_starting_params.emplace_back(params);
				}
			}
			//loop through starting params and add those that aren't in a merge model
			for(int i = 0; i < starting_params.size(); i++){
				bool inPair = false;
				for(int m = 0; m < merge_pairs.size(); m++){
					if(i == merge_pairs[m].first || i == merge_pairs[m].second){
						inPair = true;
						break;
					}
				}
				if(inPair) continue;
				merge_starting_params.push_back(starting_params[i]);

			}

			//starting from merge_starting_params, fit full data
			//don't forget to include parameters that DONT have a merge pair
			//cout << "merge model starts with " << merge_starting_params.size() << endl;
			//cout << "starting place for merge model" << endl; for(int m = 0; m < merge_starting_params.size(); m++){ cout << "cluster #" << m << endl; merge_starting_params[m]["m"].Print();}
			int mincls = merge_starting_params.size() + _nGhosts; 
			int npts = x->model->GetData()->GetNPoints();
			int k = npts < mincls ? npts : mincls;

			//cout << "# ghosts " << _nGhosts << " total k " << merge_starting_params.size()+_nGhosts << " # merge starting params " << merge_starting_params.size() << " # pts " << npts << " k " << k << endl; 
			auto mergemodel = std::make_unique<GaussianMixture>(k);
			if(_verb != 0) mergemodel->SetVerbosity(_verb-1);
			mergemodel->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				mergemodel->SetDataSmear(_data_smear);
			}
			mergemodel->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
			auto newpts = make_unique<PointCollection>(*x->model->GetData());
			mergemodel->SetData(std::move(newpts));
			mergemodel->InitParameters(_params,merge_starting_params);
			int nghosts = 0;
			int nreal = 0;
			for(int k = 0; k < mergemodel->GetNClusters(); k++){
				map<string, Matrix> params;
				mergemodel->GetLHPosteriorParameters(k, params);
				if((bool)params["ghost"].at(0,0)) nghosts++;
				else nreal++;
			}
			//cout << "MergeTree::_merge_model - starting with " << nreal << " real models and " << nghosts << " ghost models" << endl;
			_run_model(mergemodel.get(),merge_elbo);
			nghosts = 0;
			nreal = 0;
			for(int k = 0; k < mergemodel->GetNClusters(); k++){
				map<string, Matrix> params;
				mergemodel->GetLHPosteriorParameters(k, params);
				if((bool)params["ghost"].at(0,0)) nghosts++;
				else nreal++;
			}
			//cout << "MergeTree::_merge_model - ended with " << nreal << " real models and " << nghosts << " ghost models" << endl;
			//cout << "ending place for merge model" << endl; for(int m = 0; m < mergemodel->GetNClusters(); m++){ cout << "cluster #" << m << endl; mergemodel->GetLHPosteriorParameters(m)["m"].Print();} 
			return mergemodel; 
		}


		//find pdf that most closely is inner product-matched to another pdf (exclusively) - map automatically sorts in descending order
		//matching bw subcls in pseudojet n and subcls in pseudojet m with
		//i = # subclusters in pseudojet n and j = # subclusters in pseudojet m
		//if i > j then there *will* be "no match" cases
		//returning indices s.t. the pair is (idx_i, idx_j) for subcluster from pseudojet i with index idx_i and subcluster idx_j from pseudojet j
		void _match_pdfs(const vector<BasePDF*>& inpdfs_n, const vector<BasePDF*>& inpdfs_m, map<double,pair<int,int>, std::greater<double>>& bestMatchIdxs){
			//cout << "matching pdfs - " << inpdfs_n.size() << " to " << inpdfs_m.size() << endl;
			//loop through pdfs
			double bestProd, prod;
			bestMatchIdxs.clear();
			bestMatchIdxs = {};
		
			if(inpdfs_n.size() < 1 || inpdfs_m.size() < 1) return;

			//prods[i][j] = prod for jet i and jet j
			vector<vector<double>> prods;
			vector<int> idxs;
			for(int j = 0; j < inpdfs_n.size(); j++){
				prods.push_back({});
				for(int g = 0; g < inpdfs_m.size(); g++){
					prods[j].push_back(-1e308);
					//cout << "product bw pdf " << j << " of pdfs_n and " << g << " of pdfs_m "; 
					prod = _gaussian_L2_inner_product(inpdfs_n[j], inpdfs_m[g]);	
					prods[j][g] = prod;
				}
				//cout << "size of prod " << j << " is " << prods[j].size() << endl;
			}
			vector<int> best_idxs; //one per jet
			vector<double> best_prods; //one per jet
			int otherPdf, thismatchidx, thisPdf;
			//go back through pdfs can check to see if there are overlapping matches
			for(int j = 0; j < inpdfs_n.size(); j++){
				//for(int g = 0; g < nMatch; g++){
				//	cout << "jet " << j << " and match jet " << g << " have prod " << drs[j][g] << endl;
				//}
				//cout << "pdf " << j << " has best prod " << *max_element(prods[j].begin(), prods[j].end()) << " at match pdf " << find(prods[j].begin(), prods[j].end(), *max_element(prods[j].begin(), prods[j].end())) - prods[j].begin() << endl;
				double maxprod = *max_element(prods[j].begin(), prods[j].end());
				int matchidx = find(prods[j].begin(), prods[j].end(), maxprod) - prods[j].begin();
				if(maxprod == -1e308) matchidx = -1; //no match found (ie no available match jet for best match)
				best_idxs.push_back(matchidx);
				best_prods.push_back(maxprod);
				thismatchidx = matchidx;
				thisPdf = j;
				//if other jets have the same genidx matched, go through and disambiguate until there is only 1 instance of genidx
				while(count(best_idxs.begin(), best_idxs.end(), matchidx) > 1 && matchidx != -1){
					otherPdf = find(best_idxs.begin(), best_idxs.end(), matchidx) - best_idxs.begin();
					//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
					//skip otherPdf (ie thisPdf) in this case and look at all other jets
					if(otherPdf == thisPdf){
						otherPdf = find(best_idxs.begin()+otherPdf+1, best_idxs.end(), matchidx) - best_idxs.begin();
					}
					//for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
					//cout << " found another match at pdf " << otherPdf << " at matchidx " << matchidx << " with other prod " << prods[otherPdf][matchidx] << " against this pdf " << thisPdf << endl;
					//if other prod is better than current mindr
					if(prods[otherPdf][matchidx] > maxprod){
						//set this prod to 999 (is invalid), find new min for this jet, reset genidx to this index
						prods[thisPdf][matchidx] = -1e308;
						maxprod = *max_element(prods[thisPdf].begin(), prods[thisPdf].end());
						if(maxprod == -1e308) matchidx = -1;
						else matchidx = find(prods[thisPdf].begin(), prods[thisPdf].end(), maxprod) - prods[thisPdf].begin();
						best_idxs[thisPdf] = matchidx;
						best_prods[thisPdf] = maxprod;
						//cout << " reset match of this pdf " << thisPdf << " to pdf " << matchidx << " with prod " << maxprod << endl;
			
					}
					//if this prod is not as good as current prod
					else{
						//set other prod to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
						prods[otherPdf][matchidx] = -1e308;
						thismatchidx = matchidx;
						//for(int p = 0; p < prods[otherPdf].size(); p++) cout << "prod #" << p << ": " << prods[otherPdf][p] << " for pdf " << otherPdf << endl; 
						maxprod = *max_element(prods[otherPdf].begin(), prods[otherPdf].end());
						//cout << "maxprod " << maxprod << endl;
						if(maxprod == -999) matchidx = -1;
						else matchidx = find(prods[otherPdf].begin(), prods[otherPdf].end(), maxprod) - prods[otherPdf].begin();
						thisPdf = otherPdf;
						best_idxs[thisPdf] = matchidx;
						best_prods[thisPdf] = maxprod;
						//cout << " reset match of other pdf " << otherPdf << " to pdf " << matchidx << " with prod " << maxprod << endl;
					}	
					//cout << "matchidx is now " << matchidx << " with count " << count(best_idxs.begin(), best_idxs.end(), matchidx) << " for pdf " << thisPdf << endl;

				}
				
				//cout << "pdf " << thisPdf << " has best exclusive match with " << matchidx << " " << best_idxs[j] << " " << best_idxs[thisPdf] << " with inner prod " << maxprod << "\n" << endl;
			}
			for(int i = 0; i < best_idxs.size(); i++){
				//cout << "best_idxs for pdf #" << i << " is " << best_idxs[i] << " with prod " << best_prods[i] << endl;
				bestMatchIdxs[best_prods[i]] = std::make_pair(i, best_idxs[i]);	
			}
		}

		double _gaussian_dist(BasePDF* pdf1, BasePDF* pdf2){
			Matrix mu1, mu2;
			pdf1->GetParameter("mean", mu1);
			pdf2->GetParameter("mean", mu2);

			double d1 = mu1.at(0,0) - mu2.at(0,0);
			double d2 = mu1.at(1,0) - mu2.at(1,0);
			double d3 = mu1.at(2,0) - mu2.at(2,0);
		
			return sqrt(d1*d1 + d2*d2 + d3*d3);
		}


		double _gaussian_L2_inner_product(BasePDF* pdf1, BasePDF* pdf2){
			//d = integral N(x | mu1, lam1)*N(x | mu2, lam2) dx
			//analytically this integral becomes
			//exp(-muT*lam*mu + mu_1T*lam_1*mu_1 + mu_2T*lam_2*mu_2)
			//with some coefficients that i'm gonna ignore since we're comparing inner products
			//where lam = (lam_1 + lam_2)
			//mu = lam^-1(lam_1*mu_1 + lam_2*mu_2)
			//returns LOG of inner product
			Matrix mu1; pdf1->GetParameter("mean", mu1);
			Matrix mu1T(1,3);
			mu1T.transpose(mu1);
			Matrix cov1; pdf1->GetParameter("cov", cov1);
			Matrix lam1(3,3);
			lam1.invert(cov1);
			
			Matrix mu2; pdf2->GetParameter("mean", mu2);
			Matrix mu2T(1,3);
			mu2T.transpose(mu2);
			Matrix cov2; pdf2->GetParameter("cov", cov2);
			Matrix lam2(3,3);
			lam2.invert(cov2);
		
			Matrix lam;
			lam.add(lam1, lam2);
			Matrix lamInv(3,3);
			lamInv.invert(lam);			
		
			Matrix mu(3,1);
			mu.mult(lam1,mu1);
			Matrix mu_pt2(3,1);
			mu_pt2.mult(lam2,mu2);
			mu.add(mu,mu_pt2);
			mu.mult(lamInv,mu);
			Matrix muT(1,3);
			muT.transpose(mu);
		
			Matrix muLam_pt1(1,3);
			muLam_pt1.mult(muT,lam);
			Matrix muLam(1,1);
			muLam.mult(muLam_pt1,mu);
		
			Matrix muLam1_pt1(1,3);
			muLam1_pt1.mult(mu1T,lam1);
			Matrix muLam1(1,1);
			muLam1.mult(muLam1_pt1,mu1);
			
			Matrix muLam2_pt1(1,3);
			muLam2_pt1.mult(mu1T,lam1);
			Matrix muLam2(1,1);
			muLam2.mult(muLam2_pt1,mu1);
		
			double p = -0.5*(-muLam.at(0,0) + muLam1.at(0,0) + muLam2.at(0,0));	
			//determinant factors
			double det1 = lam1.det();
			double det2 = lam2.det();
			double det = lam.det();
			//cout << "det1 " << det1 << " det2 " << det2 << " det " << det << " p " << p << " ";
			p = log(sqrt(det1)*sqrt(det2)/sqrt(det)) + p;
		
			return p;
		
		}


	private:
		//keep list of nodes since tree is built bottom up
		//vector<std::shared_ptr<node>> _clusters;
		//Dirichlet prior parameter
		double _alpha, _emAlpha;
		//threshold on variational EM
		double _thresh;

		Matrix _data_smear;
		int _verb;
		map<string, Matrix> _params;
		void _remap_phi(PointCollection& points);
		Matrix _Rscale;
		Matrix _RscaleInv;
		bool _check_merges;
		int _nGhosts;
		double _log_sum_exp(vector<double>& terms);
		double _total_calcmerge_time = 0;
		double _n_calcmerge_calls = 0;
		double _total_evidence_time = 0;
		double _n_evidence_calls = 0;
};
#endif
