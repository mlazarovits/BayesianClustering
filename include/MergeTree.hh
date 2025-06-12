#ifndef MergeTree_HH
#define MergeTree_HH

#include "VarEMCluster.hh"
#include "GaussianMixture.hh"
#include "MultivarT.hh"
#include "RandomSample.hh"
#include "NodeStack.hh"
#include "Gaussian.hh"

class MergeTree : BaseTree{
	public:
		MergeTree(){ 
			_alpha = 0; _thresh = 0.; _verb = 0;
			_cell = acos(-1)/180; //default is CMS ECAL cell size
			_tresCte = 0.133913;
			_tresStoch = 1.60666; 
			_tresNoise = 0.00691415;
			//x'' = x'/b = (x-a)/b
			//sets relative importance of dimensions
			//decreasing cell -> eta/phi distance more important
			//increasing entry (2,2) -> time distance more important
			_Rscale = Matrix(3,3);
			_RscaleInv = Matrix(3,3);
			_Rscale.SetEntry(1/_cell,0,0);
			_Rscale.SetEntry(1/_cell,1,1);
			_Rscale.SetEntry(1,2,2); 
			_RscaleInv.invert(_Rscale);

			_check_merges = false;
		}

		MergeTree(double alpha){
			_alpha = alpha;
			_thresh = 1.; _verb = 0;
			_cell = acos(-1)/180; //default is CMS ECAL cell size
			_tresCte = 0.133913;
			_tresStoch = 1.60666; 
			_tresNoise = 0.00691415;
			//x'' = x'/b = (x-a)/b
			//sets relative importance of dimensions
			//decreasing cell -> eta/phi distance more important
			//increasing entry (2,2) -> time distance more important
			_Rscale = Matrix(3,3);
			_RscaleInv = Matrix(3,3);
			_Rscale.SetEntry(1/_cell,0,0);
			_Rscale.SetEntry(1/_cell,1,1);
			_Rscale.SetEntry(1,2,2); 
			_RscaleInv.invert(_Rscale);
			_check_merges = false;
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
			_cell = tree._cell; //default is CMS ECAL cell size
			_tresCte = tree._tresCte;
			_tresStoch = tree._tresStoch; 
			_tresNoise = tree._tresNoise; 
			//above tres params are for gev = 1
		
			_Rscale = tree._Rscale;
			_RscaleInv = tree._RscaleInv;
			_check_merges = tree._check_merges;
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


		void CheckMerges(bool t){ _check_merges = t; if(_check_merges) cout << "Checking merged models" << endl;}	
	
		node* Get(int i){ return _clusters[i]; }

		const vector<node*>& GetClusters() const{ 
			return _clusters;
		}
		
		void Insert(node* x){
		//cout << "inserting node with ismirror " << x->ismirror << " and pts " << endl; x->points->Print();
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
			x->ismirror = false;
			x->prob_tk = exp(Evidence(x));//p_dk_tk = p_dk_h1 since cannot be divided further
			//cout << "adding node with evidence " << Evidence(x) << " and weight " << pt->w() << endl;
			//Evidence = ELBO \approx log(LH)
			_clusters[n-1] = x;
			//cout << "adding leaf with ismirror " << x->ismirror << endl;
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

		void TransformMean(Matrix& mat, BayesPoint& center){
				//eta to theta
				double eta = mat.at(0,0);
				mat.SetEntry(2*atan(exp(-eta)),0,0);
				
				//put 02pi
				PointCollection mat_pt = mat.MatToPoints();
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


		//runs varEM to get Evidence (ELBO) for given GMM
		double Evidence(node* x){
			int k;
			vector<map<string,Matrix>> prev_posts = {};
			map<int,int> left_post_indices, right_post_indices;	
  			//if leaf node (ie r == _z && l == _z) -> set k = 1
			if(x->l == _z && x->r == _z || x->points->Sumw() < _thresh){ 
				//cout << "leaf nodes - setting max # clusters == 1" << endl;
		 		k = 1;
			//cout << "# clusters allowed " << k << " # pts " << x->points->GetNPoints() << " weight " << x->points->Sumw() << endl;
			}
			//if the sum of all weights is not enough for 1 subcluster to be "above threshold" they should not be able to be broken into multiple subclusters
			//number of clusters in node x = k_l + k_r for left and right nodes
			else{
				//cout << "not leaf nodes - setting according to max # clusters" << endl;
				int mincls = x->l->model->GetNClusters() + x->r->model->GetNClusters();
				int npts = x->l->model->GetData()->GetNPoints() + x->r->model->GetData()->GetNPoints();
				k = npts < mincls ? npts : mincls;
				//k = npts;
				//k = x->l->model->GetNClusters() + x->r->model->GetNClusters();
				//setting initial posterior parameters from models of previous steps
				for(int kk = 0; kk < x->l->model->GetNClusters(); kk++){
					prev_posts.push_back(x->l->model->GetLHPosteriorParameters(kk));
					left_post_indices[kk] = prev_posts.size()-1;
				} 
				for(int kk = 0; kk < x->r->model->GetNClusters(); kk++){
					prev_posts.push_back(x->r->model->GetLHPosteriorParameters(kk));
					right_post_indices[kk] = prev_posts.size()-1;
				}

				//int npts_abovethresh = 0;
				//for(int n = 0; n < x->points->GetNPoints(); n++){
				//	if(x->points->at(n).w() > _thresh) npts_abovethresh++;
				//}
 
			//cout << "# clusters allowed " << k << " # pts " << x->points->GetNPoints() << " # clusters from l " << x->l->model->GetNClusters() << " from r " << x->r->model->GetNClusters() << " weight " << x->points->Sumw() << " with " << npts_abovethresh << " pts above thresh: " << _thresh <<  endl;
			}
			//k = x->l->model->GetData()->GetNPoints() + x->r->model->GetData()->GetNPoints();
 //cout << "og points" << endl; x->points->Print();
			//cout << "original points" << endl; x->points->Print();
			x->points->Sort();
			PointCollection* newpts = new PointCollection(*x->points);
			if(_verb > 5){ cout << newpts->GetNPoints() << " original pts " << endl; newpts->Print();}
			//sort points so random initialization is consistent based on seed




			x->model = new GaussianMixture(k); //p(x | theta)
			if(_verb != 0) x->model->SetVerbosity(_verb-1);
			x->model->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				x->model->SetDataSmear(_data_smear);
			}

			x->model->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
			
			//in local space, circular coordinates (like phi) can go negative
			x->model->SetData(newpts); //may need to make copy of de-referenced object so as not to change the original points	
	
			//change eta to theta BEFORE calculating centroid
			x->model->EtaToTheta();
			if(_verb > 5){cout << "eta to theta" << endl; x->model->GetData()->Print();}
			//make sure all data in model is on [0,2pi] to begin with
			//some mirror nodes will have pts < 0 or > 2 pi
			//but the shifting + plane project depends on all pts initially being on unit circle
			//it doesn't matter if this is done before or after centroid calculation
			//bc of the way the centroid is calculated for phi (ie CircularCentroid)
		
			x->model->PutPhi02pi();		
			
			//since the first dimension is now an angle, its centroid needs to be circularly calculated
			//but since eta is only defined for theta on [-pi/2,pi/2], it shouldn't make that much of a difference
			BayesPoint center({x->model->GetData()->CircularCentroid(0), x->model->GetData()->CircularCentroid(1), x->model->GetData()->Centroid(2)});
			if(_verb > 5){ cout << "center" << endl; center.Print();}
		//cout << "theta circular centroid " << x->model->GetData()->CircularCentroid(0) << " regular centroid " << x->model->GetData()->Centroid(0) << endl;	
//cout << "original data (phi on 02pi)" << endl; x->model->GetData()->Print();
			x->model->ShiftData(center);



			if(_verb > 5){ cout << "translated pts" << endl; x->model->GetData()->Print(); }
			//project phi onto plane
			bool phiInf = x->model->ProjectPhi();
			bool thetaInf = x->model->ProjectTheta();
			if(_verb > 5){ cout << "projected pts" << endl; x->model->GetData()->Print();}
			//check for infinities in eta + phi from measurement error
			if(phiInf || thetaInf){
				if(_verb > 5) cout << "found inf, returning elbo as " << -1e308 << endl;
				//reset model data for further merges
				x->model->SetData(x->points);
				return -1e308;
			}


			//consider this: if inserting numerical 'infinities' causes numerical instabilities in the model, (especially in the matrix inversions)
			//derive from the ELBO equations that large infinities originating from the plane projection
			//should lead to large negative ELBOs (ie not a good lower bound on the LH)
			//and are probabilistically disfavored
			//so, due to aforementioned numerical instabilities, these quantities are not explicitly calculated
			//and are automatically excluded from merging (ie the ELBO is set to a large negative number
			//s.t. the merge is never chosen) 

			//check for infinities in eta + phi here
			if(x->model->GetData()->HasInf(0) || x->model->GetData()->HasInf(1)){
			//cout << "found inf, returning elbo as " << -1e308 << endl;
				//reset model data for further merges
				x->model->SetData(x->points);
				return -1e308;
			}

			//cout << "scale data + lam*s" << endl;	
			x->model->ScaleData(_Rscale);
			if(_verb > 5){ cout << "scaled pts" << endl; x->model->GetData()->Print();}
		
			//put means from previous step in same transformed space as data	
			for(int kk = 0; kk < prev_posts.size(); kk++){
				auto params = prev_posts[kk];
				Matrix mean = params["m"];
				
				//cout << "pre-transform - k " << kk << " alpha " << params["alpha"].at(0,0) << " mean "  << endl; mean.Print();
				TransformMean(mean, center);
				/*		
				if(fabs(mean.at(0,0)) > 1e20 || fabs(mean.at(1,0)) > 1e20){
					cout << "MEAN AT INFINITY" << endl; mean.Print();
					cout << "original mean" << endl; params["m"].Print();
					cout << "left node has " << x->l->model->GetNClusters() << " subclusters" << " and " << x->l->points->GetNPoints() << " pts" << endl;
					//cout << "left node pts" << endl; x->l->points->Print();
					cout << "centroid in phi of left node pts " << x->l->points->CircularCentroid(1) << " in eta " << x->l->points->Centroid(0) << endl;
					cout << "left means" << endl;
					for(int kkk = 0; kkk < x->l->model->GetNClusters(); kkk++){
					auto paramsl = x->l->model->GetLHPosteriorParameters(kkk);
						cout << "left subcl #" << kkk << " weight " << paramsl["alpha"].at(0,0) << " mean" << endl; paramsl["m"].Print();
					}
				
	
					cout << "right node has " << x->r->model->GetNClusters() << " subclusters" << endl;
					//cout << "right node pts" << endl; x->r->points->Print();
					//cout << "original data" << endl; x->points->Print();
					//cout << "transf data" << endl; x->model->GetData()->Print();
				}
				*/
				prev_posts[kk]["m"] = mean;
				//cout << "post-transform - k " << kk << " alpha " << params["alpha"].at(0,0) << " mean "  << endl; mean.Print();

				//do for xbar too
				mean = params["xbar"];
				//cout << " xbar " << endl; mean.Print();
				TransformMean(mean, center);
				prev_posts[kk]["xbar"] = mean;

				//don't have to do for model->GetParameter("mean") bc this isn't used in the initial logLH eval
			} 
		
			//cout << "scaled points" << endl; x->model->GetData()->Print();
			//needs to be done after data is set + transformed bc the initialization procedure depends on data
			x->model->InitParameters(_params,prev_posts);
			//x->model->InitParameters(_params);
	
			//calculate list of inner products of subcluster PDFs to compare for merges
			//only do if for # subclusters in node l = n and # subclusters in node l = m
			//n > 1 && m > 1 => n + m > 2
			map<double, pair<int, int>, std::greater<double>> prodmap;
			if(_check_merges && x->l != _z && x->r != _z && x->model->GetNClusters() > 1){
				if(x->l->model->GetNClusters() + x->r->model->GetNClusters() > 2){
				//cout << "data in x" << endl; x->model->GetData()->Print();
//cout << "data in x->l" << endl; x->l->model->GetData()->Print();
//cout << "data in x->r" << endl; x->r->model->GetData()->Print();
				//cout << "Checking merges bw left node with " << x->l->model->GetNClusters() << " clusters and right node with " << x->r->model->GetNClusters() << " clusters" << endl;
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
			}
			//nominal GMM (no merges)
			VarEMCluster* algo = new VarEMCluster(x->model, k);
			//make sure sum of weights for points is above threshold
			//otherwise the algorithm will exclude all components	
			//if(x->points->Sumw() >= _thresh){
			//	algo->SetThresh(_thresh);
			//}
			//else
			//	algo->SetThresh(x->points->Sumw());
			algo->SetThresh(1e-3);
			//cluster
			double oldLogL = algo->EvalLogL();
		 //cout << std::setprecision(10) << " it -1  firstlogl " << oldLogL <<  " # clusters " << x->model->GetNClusters() << endl;
			double LogLThresh = 1e-3;
			double newLogL = 0;
			double entropy = 0;
			double nll = 0;
			double dLogL = 1e308; 
			int it = -1;
			//cout << "oldLogL (initial) " << oldLogL << endl;
			while(dLogL > LogLThresh*fabs(oldLogL) || it == 0){
				it++;
				newLogL = algo->Cluster();
				//entropy = x->model->EvalEntropyTerms();
				//nll = x->model->EvalNLLTerms();
			//cout << "it " << it << " newLogL " << newLogL << " entropy " << entropy << " NLL " << nll << endl;
			
				//ELBO should maximizing LH -> therefore newLogL > oldLogL if both are < 0	
				dLogL = newLogL - oldLogL;
		if(std::isnan(newLogL)) cout << std::setprecision(10) << "it " << it << " new logl " << newLogL << " oldlogl " << oldLogL << " dlogl " << dLogL << " # clusters " << x->model->GetNClusters() << endl;
		 //cout << std::setprecision(10) << " it " << it << " new logl " << newLogL << " oldlogl " << oldLogL << " dlogl " << dLogL << " # clusters " << x->model->GetNClusters() << endl;
				oldLogL = newLogL;
			}
//cout << "finished in " << it << " iterations with final dLogL " << dLogL << " and final logL " << newLogL << endl;
			if(_verb > 5) cout << "EVIDENCE FOR NODE " << x->idx << " WITH " << x->model->GetData()->GetNPoints() << " POINTS AND " << k << " max clusters and " << x->model->GetNClusters() << " found clusters - evidence " << exp(newLogL) << " ELBO " << newLogL << endl;
//cout << " with points in node model " << endl; x->model->GetData()->Print();  
//cout << "original points" << endl; x->points->Print();
//cout << "x is mirror? " << x->ismirror << " x->l ismirror? " << x->l->ismirror << " x->r ismirror? " << x->r->ismirror << endl;
	//cout << "model has " << x->model->GetData()->GetNPoints() << " points" << endl;
	//if(x->model->GetData()->GetNPoints() < 3){ cout << "model pts" << endl; x->model->GetData()->Print(); cout << "node pts" << endl; x->points->Print(); }
	//cout << "node x means pre scale" << endl;
	if(newLogL > 0){
cout << "from model" << endl;
	for(int k = 0; k < x->model->GetNClusters(); k++){
		auto params = x->model->GetLHPosteriorParameters(k);
		cout << "cluster #" << k << " mean with weight " << params["alpha"].at(0,0) - _emAlpha << endl;
		params["mean"].Print();
		//cout << " cov" << endl;
		//params["cov"].Print();
	}
x->model->GetData()->Print();
	}
	//if(x->model->GetNClusters() <= 2){
		//cout << "everything found to be in 1 cluster" << endl;
		//cout << "from model" << endl;
		//for(int k = 0; k < x->model->GetNClusters(); k++){
		//	auto params = x->model->GetLHPosteriorParameters(k);
		//	cout << "cluster #" << k << " mean with weight " << params["alpha"].at(0,0) - _emAlpha << endl;
		//	params["mean"].Print();
		//	//cout << " cov" << endl;
		//	//params["cov"].Print();
		//}
		//cout << endl;
	//}
	//compare nominal model to merges IN ORDER and with memory (ie if merge1 > nom -> check merge1 & merge2)
	if(_check_merges && x->l != _z && x->r != _z && x->model->GetNClusters() > 2){
		double merge_elbo;
		//vector<map<string, Matrix>> starting_params = prev_posts;
		//cout << "x->model currently has " << x->model->GetNClusters() << " clusters" << endl;
		/*
		///cout << "from model" << endl;
		//for(int k = 0; k < x->model->GetNClusters(); k++){
		//	auto params = x->model->GetLHPosteriorParameters(k);
		//	cout << "cluster #" << k << " mean with weight " << params["alpha"].at(0,0) - _emAlpha << endl;
		//	params["mean"].Print();
		//	//cout << " cov" << endl;
		//	//params["cov"].Print();
		//}
		

		//compare merges in its respective projected data space
		vector<map<string, Matrix>> merge_posteriors;
		for(auto prodit = prodmap.begin(); prodit != prodmap.end(); prodit++){
			if(prodit->second.first == -1 || prodit->second.second == -1) continue;
			
			//compare this merge in projected data to no merge
			//project data ...
			//GaussianMixture* merged_model = _merge_model(x, prev_posts, left_post_indices[prodit->second.first], right_post_indices[prodit->second.second], merge_elbo);
	
			//if merge is better, remove associated posteriors from starting parameters and replace with their centroid


		}
		//if any merges are good, rerun whole model (unprojected) with new set of starting posteriors according to "local" merges
		*/
		/*	
		GaussianMixture* merged_model_tot = new GaussianMixture();
		merged_model_tot->SetData(newpts); 
		for(auto prodit = prodmap.begin(); prodit != prodmap.end(); prodit++){
			if(prodit->second.first == -1 || prodit->second.second == -1) continue;
			//fit GMM with merge prodit->second pair (if valid)
			cout << "best match bw clusters " << prodit->second.first << " in left " << left_post_indices[prodit->second.first] << " and " << prodit->second.second << " in right " << right_post_indices[prodit->second.second] << " with prod " << prodit->first << endl;
			GaussianMixture* merged_model = _merge_model(x, prev_posts, left_post_indices[prodit->second.first], right_post_indices[prodit->second.second], merge_elbo);
			//GaussianMixture* merged_model = _merge_model(x->model, starting_params, left_post_indices[prodit->second.first], right_post_indices[prodit->second.second], merge_elbo);
			cout << "merge ELBO " << merge_elbo << " nominal " << newLogL << endl;
			if(merge_elbo - newLogL > 0){ //ie merge model is better
				if(_verb > 1) cout << "model merging clusters " << prodit->second.first << " and " << prodit->second.second << " with product " << prodit->first << " and # clusters " << merged_model->GetNClusters() << " and ELBO " << merge_elbo << " better than nominal " << newLogL << " updating model" << endl;
				cout << "model merging clusters " << prodit->second.first << " and " << prodit->second.second << " with product " << prodit->first << " and # clusters " << merged_model->GetNClusters() << " and ELBO " << merge_elbo << " better than nominal " << newLogL << " updating model" << endl;
				
				//merged_model_tot->add(merged_model);

				//x->model = merged_model;
				//cout << "x->model now has " << x->model->GetNClusters() << " clusters" << " compared to " << prev_posts.size() << " # prev post subclusters " << endl;
				//newLogL = merge_elbo;
				//set starting params to params of merged model
				//starting_params.clear();
				//for(int n = 0; n < x->model->GetNClusters(); n++) starting_params.push_back(x->model->GetLHPosteriorParameters(n));

			}
		}
		*/
		cout << "final node model has " << x->model->GetNClusters() << " with centers" << endl;
		for(int k = 0; k < x->model->GetNClusters(); k++){
			auto params = x->model->GetLHPosteriorParameters(k);
			cout << "cluster #" << k << " mean with weight " << params["alpha"].at(0,0) - _emAlpha << endl;
			params["mean"].Print();
			cout << " cov" << endl;
			params["cov"].Print();
		}
		cout << endl;
	}	
	

	//inverse transformations - for after fitting
	center.Scale(-1);
//cout << "transf data" << endl;
//x->model->GetData()->Print();
////cout << " # clusters found " << x->model->GetNClusters() << endl;
//for(int k = 0; k < x->model->GetNClusters(); k++){
//		cout << "cluster #" << k << " transf mean " << endl;
//		auto params = x->model->GetLHPosteriorParameters(k);
//		params["mean"].Print();
//}
			//transform the parameters back into global coordinates
			//need to unscale first then uncenter since x'' = (x-a)/b (see above)
			//need to unscale data - also unscales lamStar measurement errors 
			//cout << "inverse scale data + lam*s" << endl;	
			
			x->model->ScaleData(_RscaleInv);	 
			
			//cout << "unscaled points" << endl;
			//x->model->GetData()->Print();
			//need to unscale mean + covariances
			x->model->ScaleParameters(_RscaleInv);

//cout << "unscaled" << endl;	
// x->model->GetData()->Print();
//	for(int k = 0; k < x->model->GetNClusters(); k++){
//		cout << "cluster #" << k << " unscaled mean" << endl;
//		auto params = x->model->GetLHPosteriorParameters(k);
//		params["mean"].Print();
//}
		
	
			//need to put GMM parameters AND points back in detector coordinates (not local)
			//cout << "ismirror " << x->ismirror << " left " << x->l->ismirror << " right " << x->r->ismirror << endl;
			//x->model->ShiftData(center);
			//cout << "untransformed points" << endl;

			//only does for mean rn - not sure how to do for cov...
			x->model->UnprojectPhi_params();
			x->model->UnprojectTheta_params();

	//cout << "unprojected" << endl;
	//for(int k = 0; k < x->model->GetNClusters(); k++){
	//	cout << "cluster #" << k << " unproj mean" << endl;
	//	auto params = x->model->GetLHPosteriorParameters(k);
	//	//auto params = x->model->GetDataStatistics(k);
	//	cout << "mean" << endl;
	//	params["mean"].Print();
	//	//cout << "cov" << endl;
	//	//params["cov"].Print();
	//}
      		x->model->ShiftParameters(center);


	//cout << "unshifted" << endl;
	//for(int k = 0; k < x->model->GetNClusters(); k++){
	//	cout << "cluster #" << k << " unshifted means" << endl;
	//	auto params = x->model->GetLHPosteriorParameters(k);
	//	cout << "mean" << endl;
	//	params["mean"].Print();
	//	//cout << "cov" << endl;
	//	//params["cov"].Print();
	//}
			x->model->PutPhi02pi_params(); //does for data and parameters - do after data + parameters shift so the [0,2pi] transformation doesn't get shifted
//cout << "put params and data on 02pi" << endl;
			x->model->ThetaToEta_params();
//cout << "theta to eta params" << endl;
	for(int k = 0; k < x->model->GetNClusters(); k++){
		auto params = x->model->GetLHPosteriorParameters(k);
		if(isnan(params["mean"].at(0,0){
			cout << "cluster #" << k << endl;
			cout << "mean" << endl;
			params["mean"].Print();
		}
		//cout << "cov" << endl;
		//params["cov"].Print();
	}
			//resets data to original points
			x->model->SetData(x->points);
			//x->model->GetData()->Print();
			//cout << "end evidence" << endl;
			//to consider: keeping the data in the model as the transformed points that the algorithm actually runs on and the points in the node the original ones in the detector system
		//cout << "original allowed # subclusters " << k << " found subclusters " << x->model->GetNClusters() << endl;
	//cout << "node x means post transformation" << endl;
	//vector<double> norms;
	//x->model->GetNorms(norms);
	//for(int k = 0; k < x->model->GetNClusters(); k++){
	//	cout << "cluster #" << k << " with weight " << norms[k] << endl;
	//	auto params = x->model->GetLHPosteriorParameters(k);
	//	//auto params = x->model->GetDataStatistics(k);
	//	cout << "mean" << endl;
	//	params["mean"].Print();
	//	//cout << "cov" << endl;
	//	//params["cov"].Print();
	//}
	//cout << endl;
	//cout << std::setprecision(5) << endl;
			return newLogL;
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
			if(x->l == _z && x->r == _z) return 0;
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

		//project points using nominal model responsibilities to isolate subcluster merges
		//need to project data from each node for node subcluster
		PointCollection* _project_data(node* x, int cl_left, int cl_right){
cout << "projecting data into left subcl " << cl_left << " and right subcl " << cl_right << endl;
cout << "og data " << endl; x->points->Print();
			Matrix r_nk_l = x->l->model->GetPosterior();
			if(cl_left >= r_nk_l.GetDims()[1]){
				cout << "Error: accessing cluster " << cl_left << " with posterior of " << r_nk_l.GetDims()[1] << " # cols" << endl;
				return nullptr;
			}
			PointCollection* newpts = new PointCollection();
			PointCollection* pts_l = new PointCollection(*x->l->model->GetData()); //data should already be transformed in common frame
			for(int n = 0; n < pts_l->GetNPoints(); n++){
				//get unweighted responsibility
				BayesPoint pt = pts_l->at(n);
				pt.SetWeight(r_nk_l.at(n,cl_left));
				newpts->AddPoint(pt);
			}
			Matrix r_nk_r = x->r->model->GetPosterior();
			if(cl_right >= r_nk_r.GetDims()[1]){
				cout << "Error: accessing cluster " << cl_right << " with posterior of " << r_nk_r.GetDims()[1] << " # cols" << endl;
				return nullptr;
			}
			PointCollection* pts_r = new PointCollection(*x->r->model->GetData()); //data should already be transformed in common frame
			for(int n = 0; n < pts_r->GetNPoints(); n++){
				//get unweighted responsibility
				BayesPoint pt = pts_r->at(n);
				pt.SetWeight(r_nk_l.at(n,cl_right));
				newpts->AddPoint(pt);
			}
			newpts->Sort();
cout << "newpts " << endl; newpts->Print();
			return newpts;
		}


		//runs EM algorithm to fit GMM with cl1 and cl2 merged initially in model
		//cl1 is cluster idx in left node
		//cl2 is cluster idx in right node
		//GaussianMixture* _merge_model(BasePDFMixture* model, vector<map<string,Matrix>> starting_params, int cl1, int cl2, double& merge_elbo){
		GaussianMixture* _merge_model(node* x, vector<map<string,Matrix>> starting_params, int cl1, int cl2, double& merge_elbo){
			//cout << "starting means for merge model check with " << starting_params.size() << " starting params" << endl;
			if(starting_params.size() == 0){
				merge_elbo = -1e308;
				return nullptr;
			}
			for(int s = 0; s < starting_params.size(); s++){
				if(starting_params[s].empty()) continue;
				cout << "subcluster #" << s << " has mean with weight " << starting_params[s]["alpha"].at(0,0) - _emAlpha << endl; starting_params[s]["m"].Print();
			}
			cout << "# initial clusters " << x->model->GetNClusters() << " # starting params " << starting_params.size() << " cl1 " << cl1 << " cl2 " << cl2 << endl;
			//get initial parameters and data + number of subclusters from model	
			//create new starting center from centroid of cl1 + cl2
			map<string, Matrix> new_cl;
			Matrix new_center(3,1);
			double w1 = starting_params[cl1]["alpha"].at(0,0) - _emAlpha;
			double w2 = starting_params[cl2]["alpha"].at(0,0) - _emAlpha;
			PointCollection cents = starting_params[cl1]["m"].MatToPoints({w1});
			cents += starting_params[cl2]["m"].MatToPoints({w2});
			//centroid bw clusters 1 and 2 - parameters have been transformed + projected to be on a plane
			BayesPoint center({cents.Centroid(0), cents.Centroid(1), cents.Centroid(2)});
			new_cl["m"] = center;
			//need to also give Nk (norms) and cov to seed posteriors
			new_cl["alpha"] = _emAlpha + (w1 + w2);
			//cov - average of two
			Matrix new_cov(3,3);
			Matrix cov1 = starting_params[cl1]["cov"];
			cov1.mult(cov1,w1);
			Matrix cov2 = starting_params[cl2]["cov"];
			cov2.mult(cov2,w2);
			new_cov.add(cov1,cov2);
			new_cov.mult(new_cov,1./(w1+w2));
			//mult by nu0 to undo this in InitParameters
			new_cov.mult(new_cov,_params["dof"].at(0,0));
			//need to invert bc pulls as scalemat
			new_cov.invert(new_cov);
			new_cl["scalemat"] = new_cov;


			//remove clusters specified in args - in correct order s.t. removing a subcluster doesn't affect a later subcluster
			cout << "removing clusters " << cl1 << " and " << cl2 << endl;
			starting_params[cl1].clear();
			starting_params[cl2].clear();
			//add new centroid
			starting_params.push_back(new_cl);
			int nstart = 0;
			for(int i = 0; i < starting_params.size(); i++){
				if(starting_params[i].empty()) continue;
				cout << "starting subcluster #" << i << " with weight " << starting_params[i]["alpha"].at(0,0) - _emAlpha << endl; starting_params[i]["m"].Print();
				nstart++;
			}
			GaussianMixture* mergemodel = new GaussianMixture(nstart);
			if(_verb != 0) mergemodel->SetVerbosity(_verb-1);
			mergemodel->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				mergemodel->SetDataSmear(_data_smear);
			}

			mergemodel->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 

			//project data
			PointCollection* mergepts = _project_data(x,cl1,cl2);
			//PointCollection* mergepts = new PointCollection(*model->GetData()); //data should already be transformed in common frame
			mergemodel->SetData(mergepts);
			mergemodel->InitParameters(_params,starting_params);
		
			VarEMCluster* algo = new VarEMCluster(mergemodel, nstart);
			//make sure sum of weights for points is above threshold
			//otherwise the algorithm will exclude all components	
			//if(mergepts->Sumw() >= _thresh){
			//	algo->SetThresh(_thresh);
			//}
			//else
			//	algo->SetThresh(mergepts->Sumw());
				algo->SetThresh(1e-3);
			//cluster
			double oldLogL = algo->EvalLogL();
		 //cout << std::setprecision(10) << " it -1  firstlogl " << oldLogL <<  " # clusters " << x->model->GetNClusters() << endl;
			double LogLThresh = 1e-3;
			double newLogL = 0;
			double entropy = 0;
			double nll = 0;
			double dLogL = 1e308; 
			int it = -1;
			//cout << "oldLogL (initial) " << oldLogL << endl;
			while(dLogL > LogLThresh*fabs(oldLogL) || it == 0){
				it++;
				newLogL = algo->Cluster();
				//ELBO should maximizing LH -> therefore newLogL > oldLogL if both are < 0	
				dLogL = newLogL - oldLogL;
		if(std::isnan(newLogL)) cout << std::setprecision(10) << "it " << it << " new logl " << newLogL << " oldlogl " << oldLogL << " dlogl " << dLogL << " # clusters " << mergemodel->GetNClusters() << endl;
				oldLogL = newLogL;
			}

			merge_elbo = newLogL;
			cout << "merge model has " << mergemodel->GetData()->GetNPoints() << " points and " << mergemodel->GetNClusters() << " final clusters with means " << endl;
			for(int k = 0; k < mergemodel->GetNClusters(); k++){
				auto params = mergemodel->GetLHPosteriorParameters(k);
				cout << "cluster #" << k << " mean with weight " << params["alpha"].at(0,0) - _emAlpha << endl;
				params["mean"].Print();
			}
			return mergemodel;
		}


		//find pdf that most closely is inner product-matched to another pdf (exclusively) - map automatically sorts in descending order
		//matching bw subcls in pseudojet n and subcls in pseudojet m with
		//i = # subclusters in pseudojet n and j = # subclusters in pseudojet m
		//if i > j then there *will* be "no match" cases
		//returning indices s.t. the pair is (idx_i, idx_j) for subcluster from pseudojet i with index idx_i and subcluster idx_j from pseudojet j
		void _match_pdfs(vector<BasePDF*>& inpdfs_n, vector<BasePDF*> inpdfs_m, map<double,pair<int,int>, std::greater<double>>& bestMatchIdxs){
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
			Matrix mu1 = pdf1->GetParameter("mean");
			Matrix mu2 = pdf2->GetParameter("mean");

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
			Matrix mu1 = pdf1->GetParameter("mean");
			Matrix mu1T(1,3);
			mu1T.transpose(mu1);
			Matrix cov1 = pdf1->GetParameter("cov");
			Matrix lam1(3,3);
			lam1.invert(cov1);
			
			Matrix mu2 = pdf2->GetParameter("mean");
			Matrix mu2T(1,3);
			mu2T.transpose(mu2);
			Matrix cov2 = pdf2->GetParameter("cov");
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
		vector<node*> _clusters;
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
};
#endif
