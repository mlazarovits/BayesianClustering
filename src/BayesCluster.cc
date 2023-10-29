#include "BayesCluster.hh"
#include "DynamicNearestNeighbours.hh"
#include "FullViz3D.hh"
#include "VarClusterViz3D.hh"

// The structure of this method is respectfully repurposed from ClusterSequence_Delaunay in FastJet (Cacciari, Salam, Soyez).
// This work was modified from its original form by Margaret Lazarovits on October 2, 2023. 
// The original version of this work was released
// under version 2 of the GNU General Public License. As of v3 of GNU GPL,
// any conditions added in Section 7 also apply. 

//----------------------------------------------------------------------
// Copyright (c) 2005-2021, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//----------------------------------------------------------------------
/// Run the clustering using a Hierarchical Delaunay triangulation and
/// STL maps to achieve O(N*ln N) behaviour.
///
/// There may be internally asserted assumptions about absence of
/// points with coincident eta-phi coordinates.
const vector<node*>& BayesCluster::_delauney_cluster(){
	//the 2D Delauney triangulation used in FastJet will be used to seed the 3D clustering
	MergeTree* mt = new MergeTree(_alpha);
	mt->SetSubclusterAlpha(_subalpha);
	mt->SetVerbosity(_verb);
	//set data smear
	if(!_smear.empty()) mt->SetDataSmear(_smear);
	if(_thresh != -999) mt->SetThresh(_thresh);
	//set distance constraint
	mt->SetDistanceConstraint(0,acos(-1)/2.);
	int n = _points.size();	
	for (int i = 0; i < n; i++) {
		//should only be one point per entry in points
		if(_points[i].GetNPoints() != 1){
			cerr << "BayesCluster - Error: multiple points in one collection of starting vector." << endl;
			return _trees;
		} 
		mt->AddLeaf(&_points[i].at(0));
	}
	const bool verbose = false;
	const bool ignore_nearest_is_mirror = true; //based on _Rparam < twopi, should always be true for this 
	if(_verb > 0) cout << "BayesCluster - # clusters: " << mt->GetNClusters() << endl;
	Dnn2piCylinder* DNN = new Dnn2piCylinder(_points, ignore_nearest_is_mirror, mt, verbose);

	//need to make a distance map like in FastJet, but instead of clustering
	//based on geometric distance, we are using merge probability (posterior) from BHC
	//all three dimensions will go into calculating the probabilityu
	//but the map will be built only in 2D space
	//structure is typdef'ed in header
	CompareMap ProbMap, DistMap;
	InvCompareMap InvDistMap;
	//fill map with initial potential clusterings
	for(int i = 0; i < n; i++){
		_add_entry_to_maps(i, ProbMap, DNN);
		_add_entry_to_maps(i, DistMap, DNN, false);	
		_add_entry_to_maps(i, InvDistMap, DNN);	
	}
	bool done = false;
	std::pair<std::multimap<double,verts>::iterator, std::multimap<double,verts>::iterator> ret;	
	//run the clustering
	for(int i = 0; i < n; i++){
		// find largest rk value in map (last entry)
		double BestRk;
		verts BestRkPair;
		std::multimap<double,verts>::iterator map_it;
		int jet_i, jet_j;
		bool Valid2;
		double mindist = 999;
		double dist;
		do{
			map_it = ProbMap.end();
			map_it--;
			BestRk = map_it->first;
			BestRkPair = map_it->second;
			//check for equal rks (more than three - may be two for i,j and j,i), break tie with 3d distance
			//right now - only equal rks are for equivalent merges
			if(ProbMap.count(BestRk) > 2){
				ret = ProbMap.equal_range(BestRk);
				for(std::multimap<double,verts>::iterator it = ret.first; it != ret.second; ++it){
					jet_i = it->second.first;
					jet_j = it->second.second;
					if(_verb > 1) cout << "rk: " << BestRk << " with points: " << jet_i << ", " << jet_j << endl;
					dist = InvDistMap[jet_i].second;
					if(dist < mindist) mindist = dist;

				}
			if(_verb > 1) cout << "mindist: " << mindist << " jet_i: " << DistMap.find(mindist)->second.first << " jet_j: " << DistMap.find(mindist)->second.second  << endl;
			jet_i = DistMap.find(mindist)->second.first; jet_j = DistMap.find(mindist)->second.second;
		
			}else{
				jet_i = BestRkPair.first;
				jet_j = BestRkPair.second;
			}

			if (_verb > 0) cout << "BayesCluster found recombination candidate: " << jet_i << " " << jet_j << " " << BestRk << " " << ProbMap.size() << endl; // GPS debugging
 			//also need to erase any impossible merges from map too
			ProbMap.erase(map_it);
			Valid2 = DNN->Valid(jet_j);
			if (_verb > 1) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << endl;
		} while((!DNN->Valid(jet_i) || !Valid2) && ProbMap.size() >= 1); //this is what checks to see if merges are still allowed or if they include points that have already been merged
		//if point matches to itself (mirror point), find best geo match
		if((jet_i == jet_j) && (ProbMap.size() > 1)){
			// find largest rk value in map (last entry)
			double BestDist;
			verts BestDistPair;
			std::multimap<double,verts>::iterator map_it_dist;
			for(map_it_dist = DistMap.begin(); map_it_dist != DistMap.end(); ++map_it_dist){
				BestDist = map_it_dist->first;
				BestDistPair = map_it_dist->second;
				jet_i = BestDistPair.first;
				jet_j = BestDistPair.second;
				if (_verb > 0) cout << "BayesCluster found distance recombination candidate: " << jet_i << " " << jet_j << " " << BestRk << endl; // GPS debugging
			if (_verb > 1) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << endl;
				//check validity (to see if it has been clustered before)
				//if best match is between mirrored points AND there are no points that haven't been geometrically clustered yet,
				//we're done
				if(!DNN->Valid(jet_i) || !DNN->Valid(jet_j)){ done = true; break; }
				//also need to erase any impossible merges from map too
				DistMap.erase(map_it_dist);
			break;
			} 
		}
		if(done){ cout << "done clustering" << endl; break;}
if(_verb > 1) cout << "get best rk for jet " << jet_i << " and " << jet_j << endl;
		//if max rk < 0.5, can stop clustering
		if(BestRk < 0.5){ done = true; break; }	


		int nn;
		if(_verb > 0) cout << "BayesCluster call _do_ij_recomb: " << jet_i << " " << jet_j << " " << BestRk << endl; // GPS debug
		//do_ij_recomb - this should be the same as in the OG code (except rk instead of dij)
      		_do_ij_recombination_step(jet_i, jet_j, BestRk, nn);
		// exit the loop because we do not want to look for nearest neighbours
		// etc. of zero partons
		if (i == n-1) {break;}//get eta phi (and time!) of new point - centroid of points in combined vertices
		int pt3;
		vector<int> updated_neighbors;
		//update DNN with RemoveCombinedAddCombination
		//this should also update the merge tree - RemoveAndAddPoints in DnnPlane does
		vector<JetPoint> jps = _jets[_jets.size() - 1].GetJetPoints();
		PointCollection newpts = PointCollection();
		for(int i = 0; i < (int)jps.size(); i++){
			Point pt = Point({jps[i].eta(), jps[i].phi_02pi(), jps[i].t()});
			newpts += pt;
		}
		DNN->RemoveCombinedAddCombination(jet_i, jet_j,
							newpts, pt3, updated_neighbors);
	if(_verb > 1) cout << "updating map" << endl;
		//update map
		vector<int>::iterator it = updated_neighbors.begin();
		for(; it != updated_neighbors.end(); ++it){
			int ii = *it;
			_add_entry_to_maps(ii, ProbMap, DNN);
			_add_entry_to_maps(ii, DistMap, DNN, false);	
			_add_entry_to_maps(ii, InvDistMap, DNN);
		}
		if(_verb > 0) cout << "\n" << endl;
	}

	//MergeTree* endmerge = DNN->GetMergeTree();
	_trees = mt->GetClusters();
	cout << mt->GetNClusters() << " final clusters" << endl;
	double nmirror = 0;
	for(int i = 0; i < _trees.size(); i++){
		//ignore removed (clustered) points and don't plot mirror points
		if(_trees[i] == nullptr) continue;
		if(_trees[i]->points->mean().at(1) > 2*acos(-1) || _trees[i]->points->mean().at(1) < 0){
			nmirror++; 
			//trees.erase(trees.begin()+i);
			//i--;
			continue; }
		cout << "getting " << _trees[i]->points->GetNPoints() << " points in cluster #" << i << endl;
		_trees[i]->points->Print();
		//cout << trees[i]->l->points->GetNPoints() << " in left branch " << trees[i]->r->points->GetNPoints() << " in right branch" << endl;
	}
	if(_verb > 0) cout << nmirror << " mirror points." << endl;
	cout << mt->GetNClusters() - nmirror << " jets found." << endl;
	return _trees;
	
	
}



const vector<node*>& BayesCluster::_naive_cluster(){
        node* di; node* dj;
        double rk;
        double rk_max;
	NodeStack list; //for sorting nodes by merge probability and tracking potential merges
	
	//set up merge tree that will calculate merges and track nodes
	MergeTree* mt = new MergeTree(_alpha); //keeps track of merges made + merge history
	mt->SetSubclusterAlpha(_subalpha);
	mt->SetVerbosity(_verb);
	//set data smear
	if(!_smear.empty()) mt->SetDataSmear(_smear);
	if(_thresh != -999) mt->SetThresh(_thresh);
	//set distance constraint
	mt->SetDistanceConstraint(0,acos(-1)/2.);
	int n = _points.size();
	if(_verb > 0) cout << "Creating rk list for all pairings (this might take a minute...)" << endl;	
	for (int i = 0; i < n; i++) {
		//should only be one point per entry in points
		if(_points[i].GetNPoints() != 1){
			cerr << "BayesCluster - Error: multiple points in one collection of starting vector. " << _points[i].GetNPoints() << " points in collection " << i << "." << endl;
			return _trees;
		} 
		mt->AddLeaf(&_points[i].at(0));
		mt->CreateMirrorNode(mt->Get(i));
	}

	int it = 0;
	//loop through all possible pairs to create list of merge probabilites r_{i+j} = r_k
	//O(N^2)
	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n; j++){
			di = mt->Get(i);
			dj = mt->Get(j);
			node* x = mt->CalculateMerge(di, dj);
		if(_verb > 2){
		cout << "checking nodes " << i << ": ";
		di->points->Print();
		cout << " and " << j << ": ";
		dj->points->Print();
		cout << "this rk: " << x->val <<  "\n" << endl;
		}	
			list.insert(x);
		}
	}
	//do clustering
	double maxrk = 0.5;
	while(mt->GetNClusters() > 1){
                if(_verb > 0) cout << "---------- iteration: " << it << " ----------" << endl;
                //if(_verb > 0) list.Print();
		
		list.sort();
	
		//get node with maximum merge probability + remove impossible merges because of max node
		node* max = list.fullpop();
		//check that popped node is not null
		if(max == nullptr) break;
		if(_verb > 2){ cout << "max merge is " << max->l->idx << " and " << max->r->idx << " rk: " << max->val << " " << endl;
		max->points->Print();
		cout << "\n" << endl;
		}
		//if maximum merge probability is less than 0.5, stop clustering (not probabilistically favored merges left)
                if(max->val < maxrk){
                        cout << "reached min rk = " << max->val << " <  " << maxrk << " - final iteration: " << it <<  " - " << mt->GetNClusters() << " clusters" << endl;
                        break;
                }
		
		//update merge tree with selected merge
		mt->Insert(max);
		//make new mirror node for max if necessary
		mt->CreateMirrorNode(max);
		mt->Remove(max->l);
		mt->Remove(max->r);

		//if this merge clusters all avaible points, stop clustering
		if(max->points->GetNPoints() == n){
                        cout << "All points clustered. Exiting clustering..." << endl;
                        break;
		}	
		//print remaining possible merges
		if(_verb > 2){
		cout << "remaining possible merges" << endl;
		list.Print(1);
		cout << "\n" << endl; 
		}

		//create new merges with the remaining nodes and the newly formed cluster
		//this operation is O(N)
		if(_verb > 2) cout << "new merges - " << mt->GetNClusters() << " remaining clusters" << endl;
		node* dl;
		for(int i = 0; i < mt->GetNAllClusters(); i++){
			dl = mt->Get(i);
			//make sure this node exists
			if(dl == nullptr) continue;
			if(_verb > 2) cout << "checking pair " << i << " or " << dl->idx << " and max " << max->idx << endl;
			//make sure this is not the max merge cluster - doesn't need to be paired with itself
			if(dl == max) continue;
			//calculate merge for node i and max node
			node* x = mt->CalculateMerge(dl, max);
		if(_verb > 2){
		cout << "checking nodes " << i << ": ";
		dl->points->Print();
		cout << " and max " << max->idx << ": ";
		max->points->Print();
		cout << "this rk: " << x->val << "\n" << endl;
		}
			//make sure rk is not nan
			if(isnan(x->val)) break;
			//add new potential merge to list
			list.insert(x); 
		}


		it++;

	}

	if(_verb > 0) cout << "---------- final iteration: " << it <<  " - " << mt->GetNClusters() << " clusters ----------" << endl;
        _trees = mt->GetClusters();
	cout << mt->GetNClusters() << " final clusters" << endl;
	double nmirror = 0;
	for(int i = 0; i < _trees.size(); i++){
		//ignore removed (clustered) points and don't plot mirror points
		if(_trees[i] == nullptr) continue;
		if(_trees[i]->points->mean().at(1) > 2*acos(-1) || _trees[i]->points->mean().at(1) < 0){
			nmirror++; 
			//trees.erase(trees.begin()+i);
			//i--;
			continue; }
		cout << "getting " << _trees[i]->points->GetNPoints() << " points in cluster #" << i << endl;
		_trees[i]->points->Print();
		//cout << trees[i]->l->points->GetNPoints() << " in left branch " << trees[i]->r->points->GetNPoints() << " in right branch" << endl;
	}
	if(_verb > 0) cout << nmirror << " mirror points." << endl;
	cout << mt->GetNClusters() - nmirror << " jets found." << endl;
	return _trees;

}



void BayesCluster::_subcluster(string oname){
	//create GMM model
	PointCollection* points = new PointCollection();
	for(auto point : _points){
		points->AddPoints(point);
	}
	int maxK = points->GetNPoints();
	GaussianMixture* gmm = new GaussianMixture(maxK);
	
	gmm->SetData(points);
	gmm->SetAlpha(_subalpha);
	gmm->SetVerbosity(_verb);
	gmm->InitParameters();
	gmm->InitPriorParameters();

	gmm->SetDataSmear(_smear);

	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,maxK);
	algo->SetThresh(_thresh);
	

	map<string, vector<Matrix>> params;
	bool viz = false;
	if(!oname.empty()){
		viz = true;
	}

	VarClusterViz3D cv3D(algo);
	if(viz){
		cv3D.SetVerbosity(_verb);
		cv3D.UpdatePosterior();
		cv3D.WriteJson(oname+"/it0");
	}
	//loop
	double dLogL, newLogL;
	double LogLthresh = 0.01;
	double oldLogL = algo->EvalLogL();
	////////run EM algo////////
	//maximum of 50 iterations
	for(int it = 0; it < 50; it++){
	
		//E step
		algo->Estimate();
		//M step
		algo->Update();
		
		//Plot
		if(viz){
			cv3D.UpdatePosterior();
			cv3D.WriteJson(oname+"/it"+std::to_string(it+1));
		}
		//Check for convergence
		newLogL = algo->EvalLogL();
		if(isnan(newLogL)){
			cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << endl;
			return;
		}
		dLogL = oldLogL - newLogL;
		if(_verb > 0) cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(fabs(dLogL) < LogLthresh){// || dLogL > 0){
			if(_verb > 0){
				cout << "Reached convergence at iteration " << it+1 << endl;
			}
			break;
		}
		oldLogL = newLogL;
	}
	if(_verb > 1){
		cout << "Estimated parameters" << endl;
		map<string, Matrix> params;
		for(int i = 0; i < gmm->GetNClusters(); i++){
			params = gmm->GetPriorParameters(i);	
			cout << "weight " << i << ": " << params["pi"].at(0,0) << " alpha: " << params["alpha"].at(0,0) << endl;
			cout << "mean " << i << endl;
			params["mean"].Print();
			cout << "cov " << i << endl;
			params["cov"].Print();
			params.clear();
		}

	}

	return;








}



void BayesCluster::_add_entry_to_maps(const int i, InvCompareMap& inmap, const Dnn2piCylinder* DNN){
		double dist;
		int j;
		dist = DNN->NearestNeighbourDistance(i);
		j = DNN->NearestNeighbourIndex(i);
//cout << "adding entry " << i << " " << j << " with dist " << dist << " to inv map" << endl;
		inmap.insert(InvCompEntry(i,std::make_pair(j,dist)));
}

//need to add idx corresponding to plane index because those are the nodes stored in merge tree
void BayesCluster::_add_entry_to_maps(const int i, CompareMap& inmap, const Dnn2piCylinder* DNN, bool prob){
		bool verbose = false;
		double comp;
		int j, cyl_j;
		if(prob){
			comp = DNN->NearestNeighbourProb(i);
			j = DNN->NearestNeighbourProbIndex(i,cyl_j);
			if(i == j){ j = cyl_j;} //don't want to match to own (mirrored) point
			if(_verb > 1) cout << "adding entry " << i << " with best probability " << comp << " pair " << j << " cyl index: " << cyl_j << endl;
		}
		else{
			comp = DNN->NearestNeighbourDistance(i);
			j = DNN->NearestNeighbourIndex(i);
	if(verbose) cout << "adding entry " << i << " with best distance " << comp << " pair " << j << endl;
		}
		inmap.insert(CompEntry(comp,verts(i,j)));
}

