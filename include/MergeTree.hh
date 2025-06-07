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


		void CheckMerges(bool t){ _check_merges = t; }	
	
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
				} 
				for(int kk = 0; kk < x->r->model->GetNClusters(); kk++){
					prev_posts.push_back(x->r->model->GetLHPosteriorParameters(kk));
				}

				//int npts_abovethresh = 0;
				//for(int n = 0; n < x->points->GetNPoints(); n++){
				//	if(x->points->at(n).w() > _thresh) npts_abovethresh++;
				//}
 
			//cout << "# clusters allowed " << k << " # pts " << x->points->GetNPoints() << " # clusters from l " << x->l->model->GetNClusters() << " from r " << x->r->model->GetNClusters() << " weight " << x->points->Sumw() << " with " << npts_abovethresh << " pts above thresh: " << _thresh <<  endl;
			}
			//k = x->l->model->GetData()->GetNPoints() + x->r->model->GetData()->GetNPoints();
			//cout << "k L " <<  x->l->model->GetNClusters() << " k R " << x->r->model->GetNClusters() << " k " << k << endl;}
 //cout << "og points" << endl; x->points->Print();
			//cout << "original points" << endl; x->points->Print();
			//do phi wraparound? may already be taken care of in local coords + mirror pts
			PointCollection* newpts = new PointCollection(*x->points);
			newpts->Sort(0);
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
			if(_check_merges){
				if(x->l->model->GetNClusters() + x->r->model->GetNClusters() > 2){
				      //construct PDF inner product map to pairs of PDFs
					map<double, pair<int, int>, std::greater<double>> prodmap;
				      //get vectors of PDFs from left (first) then right (second)
					vector<BasePDF*> l_pdfs, r_pdfs;
					for(int n = 0; n < x->l->model->GetNClusters(); n++)
						l_pdfs.push_back(x->l->model->GetModel(n));	
					for(int n = 0; n < x->r->model->GetNClusters(); n++)
						r_pdfs.push_back(x->r->model->GetModel(n));
					_match_pdfs(l_pdfs, r_pdfs, prodmap);
				}
			}

			//nominal GMM (no merges)
			VarEMCluster* algo = new VarEMCluster(x->model, k);
			//make sure sum of weights for points is above threshold
			//otherwise the algorithm will exclude all components	
			if(x->points->Sumw() >= _thresh){
				algo->SetThresh(_thresh);
			}
			else
				algo->SetThresh(x->points->Sumw());
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
		cout << "cluster #" << k << " mean " << endl;
		auto params = x->model->GetLHPosteriorParameters(k);
		params["mean"].Print();
		cout << " cov" << endl;
		params["cov"].Print();
	}
x->model->GetData()->Print();
	}

	//compare nominal model to merges IN ORDER and with memory (ie if merge1 > nom -> check merge1 & merge2)
	//	for(auto prodit = prodmap.begin(); prodit != prodmap.end(); prodit++){
			//merge_elbo = fit GMM with merge prodit->second pair
	//		if(merge_elbo - newLogL > 0){ (ie merge model is better)
	//			x->model = merge_model
	//			newLogL = merge_elbo
	//		}
	//	}	
	//}


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
//	for(int k = 0; k < x->model->GetNClusters(); k++){
//		cout << "cluster #" << k << endl;
//		auto params = x->model->GetLHPosteriorParameters(k);
//		cout << "mean" << endl;
//		params["mean"].Print();
//		//cout << "cov" << endl;
//		//params["cov"].Print();
//	}
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


		//runs EM algorithm to fit GMM with cl1 and cl2 merged initially in model
		BasePDFMixture* _merge_model(BasePDFMixture* model, vector<map<string,Matrix>>& starting_params, int cl1, int cl2, double& merge_elbo){
			//get initial parameters and data + number of subclusters from model	
			int merge_k = model->GetNClusters() - 1;
			//create new starting center from centroid of cl1 + cl2
			map<string, Matrix> new_cl;
			Matrix new_center(3,1);
			double w1 = starting_params[cl1]["alpha"].at(0,0) - _emAlpha;
			double w2 = starting_params[cl2]["alpha"].at(0,0) - _emAlpha;
			PointCollection cents = starting_params[cl1]["m"].MatToPoints({w1});
			cents += starting_params[cl2]["m"].MatToPoints({w2});
			//centroid bw clusters 1 and 2
			BayesPoint center({cents.CircularCentroid(0), cents.CircularCentroid(1), cents.Centroid(2)});
			new_cl["m"] = center;
			//only using posterior centers of subclusters to seed parameters for now (see GaussianMixture::InitParameters())

			//remove clusters specified in args - in correct order s.t. removing a subcluster doesn't affect a later subcluster
			if(cl1 > cl2){
				starting_params.erase(starting_params.begin() + cl1);
				starting_params.erase(starting_params.begin() + cl2);
			}
			else{
				starting_params.erase(starting_params.begin() + cl2);
				starting_params.erase(starting_params.begin() + cl1);
			}
			//add new centroid
			starting_params.push_back(new_cl);

			BasePDFMixture* mergemodel = new GaussianMixture(merge_k);
			if(_verb != 0) mergemodel->SetVerbosity(_verb-1);
			mergemodel->SetAlpha(_emAlpha);
			if(!_data_smear.empty()){
				mergemodel->SetDataSmear(_data_smear);
			}

			mergemodel->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 

			PointCollection* mergepts = new PointCollection(*mergemodel->GetData()); //data should already be transformed in common frame
			mergemodel->SetData(mergepts);
			mergemodel->InitParameters(_params,starting_params);
		
			//nominal GMM (no merges)
			VarEMCluster* algo = new VarEMCluster(mergemodel, merge_k);
			//make sure sum of weights for points is above threshold
			//otherwise the algorithm will exclude all components	
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
			return mergemodel;
		}


		//find pdf that most closely is inner product-matched to another pdf (exclusively) - map automatically sorts in descending order
		//matching bw subcls in pseudojet n and subcls in pseudojet m with
		//i = # subclusters in pseudojet n and j = # subclusters in pseudojet m
		//if i > j then there *will* be "no match" cases
		//returning indices s.t. the pair is (idx_i, idx_j) for subcluster from pseudojet i with index idx_i and subcluster idx_j from pseudojet j
		void _match_pdfs(vector<BasePDF*>& inpdfs_n, vector<BasePDF*> inpdfs_m, map<double,pair<int,int>, std::greater<double>>& bestMatchIdxs){
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
					prods[j].push_back(-999);
			
					prod = _gaussian_L2_inner_product(inpdfs_n[j], inpdfs_m[g]);	
					prods[j][g] = prod;
				}
			}
			vector<int> best_idxs; //one per jet
			int otherPdf, thismatchidx, thisPdf;
			//go back through pdfs can check to see if there are overlapping matches
			for(int j = 0; j < inpdfs_n.size(); j++){
				//for(int g = 0; g < nMatch; g++){
				//	cout << "jet " << j << " and match jet " << g << " have prod " << drs[j][g] << endl;
				//}
				//cout << "jet " << j << " has best prod " << *min_element(drs[j].begin(), drs[j].end()) << " at match jet " << find(drs[j].begin(), drs[j].end(), *min_element(drs[j].begin(), drs[j].end())) - drs[j].begin() << endl;
				double maxprod = *max_element(prods[j].begin(), prods[j].end());
				int matchidx = find(prods[j].begin(), prods[j].end(), maxprod) - prods[j].begin();
				if(maxprod == -999) matchidx = -1; //no match found (ie no available match jet for best match)
				best_idxs.push_back(matchidx);
				thismatchidx = matchidx;
				thisPdf = j;
				//if other jets have the same genidx matched, go through and disambiguate until there is only 1 instance of genidx
				while(count(best_idxs.begin(), best_idxs.end(), matchidx) > 1 && matchidx != -1){
					otherPdf = find(best_idxs.begin(), best_idxs.end(), matchidx) - best_idxs.begin();
					//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
					//skip otherJet (ie thisJet) in this case and look at all other jets
					if(otherPdf == thisPdf){
						otherPdf = find(best_idxs.begin()+otherPdf+1, best_idxs.end(), matchidx) - best_idxs.begin();
					}
					//for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
					//cout << " found another match at jet " << otherJet << " with other prod " << drs[otherJet][genidx] << " against this jet " << thisJet << endl;
					//if other prod is less than current mindr
					if(prods[otherPdf][matchidx] < maxprod){
						//set this prod to 999 (is invalid), find new min for this jet, reset genidx to this index
						prods[thisPdf][matchidx] = -999;
						maxprod = *max_element(prods[thisPdf].begin(), prods[thisPdf].end());
						if(maxprod == -999) matchidx = -1;
						else matchidx = find(prods[thisPdf].begin(), prods[thisPdf].end(), maxprod) - prods[thisPdf].begin();
						best_idxs[thisPdf] = matchidx;
						//cout << " reset gen match of this jet " << thisJet << " to gen jet " << genidx << " with prod " << mindr << endl;
			
					}
					//if this prod is less than (or equal to) current mindr
					else{
						//set other prod to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
						prods[otherPdf][matchidx] = -999;
						thismatchidx = matchidx;
						maxprod = *max_element(prods[otherPdf].begin(), prods[otherPdf].end());
						if(maxprod == -999) matchidx = -1;
						else matchidx = find(prods[otherPdf].begin(), prods[otherPdf].end(), maxprod) - prods[otherPdf].begin();
						thisPdf = otherPdf;
						best_idxs[thisPdf] = matchidx;
						//cout << " reset match of other jet " << otherJet << " to  jet " << idx << " with prod " << mindr << endl;
					}	
					//cout << "matchidx is now " << matchidx << " with count " << count(best_idxs.begin(), best_idxs.end(), matchidx) << " for jet " << thisJet << endl;

				}
				//cout << "jet " << j << " has best exclusive match with " << best_idxs[j] << "\n" << endl;
				bestMatchIdxs[maxprod] = std::make_pair(thisPdf,matchidx);
			}
		}

		double _gaussian_L2_inner_product(BasePDF* pdf1, BasePDF* pdf2){
			//d = integral N(x | mu1, lam1)*N(x | mu2, lam2) dx
			//analytically this integral becomes
			//exp(-muT*lam*mu + mu_1T*lam_1*mu_1 + mu_2T*lam_2*mu_2)
			//with some coefficients that i'm gonna ignore since we're comparing inner products
			//where lam = (lam_1 + lam_2)
			//mu = lam^-1(lam_1*mu_1 + lam_2*mu_2)
				
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
		
			double p = exp(-0.5*(-muLam.at(0,0) + muLam1.at(0,0) + muLam2.at(0,0)));	
		
			//determinant factors
			double det1 = lam1.det();
			double det2 = lam2.det();
			double det = lam.det();
		
			p = (sqrt(det1)*sqrt(det2)/sqrt(det))*p;
		
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
