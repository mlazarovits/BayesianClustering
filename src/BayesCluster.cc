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
	vector<double> ws;
	_jets[0].GetWeights(ws);
	double gev = ws[0]/_jets[0].E();
	//the 2D Delauney triangulation used in FastJet will be used to seed the 3D clustering
	MergeTree* mt = new MergeTree(_alpha);
	mt->SetSubclusterAlpha(_subalpha);
	mt->SetVerbosity(_verb);
	//set data smear
	if(!_smear.empty()) mt->SetDataSmear(_smear);
	//set time resolution smear: c^2 + n^2/e^2
	//remember time is already in ns
	//e = w/gev
	if(_thresh != -999) mt->SetThresh(_thresh);
	cout << "BayesCluster thresh " << _thresh << " alpha " << _alpha << " EM alpha " << _subalpha << endl;
	cout << "beta0" << endl;
	_prior_params["scale"].Print();
	cout << "mean0" << endl;
	_prior_params["mean"].Print();
	cout << "nu0" << endl;
	_prior_params["dof"].Print();
	cout << "W0" << endl;
	_prior_params["scalemat"].Print();

	cout << "BayesCluster delauney cluster - Using tresCte " << _tresCte << " _tresStoch " << _tresStoch << " _tresNoise " << _tresNoise << " gev " << gev << " _cell " << _cell << endl;
	mt->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
 
	mt->SetPriorParameters(_prior_params);
	//set distance constraint
	//mt->SetDistanceConstraint(0,acos(-1)/2.);
	int n = _points.size();	
cout << "n starting pts " << n << endl;
if(_verb > 3)cout <<  "original pts " << endl;
	for (int i = 0; i < n; i++) {
		//should only be one point per entry in points
		if(_points[i].GetNPoints() != 1){
			cerr << "BayesCluster - Error: multiple points in one collection of starting vector." << endl;
			cout << "return" << endl;
			return _trees;
		}
		if(_verb > 3){cout << i <<" "; _points[i].Print();}
		mt->AddLeaf(&_points[i].at(0));
	}
	if(_verb > 0) cout << "--------------------------------\nBayesCluster - # clusters: " << mt->GetNClusters() << endl;
	const bool verbose = false;
	const bool ignore_nearest_is_mirror = true; //based on _Rparam < twopi, should always be true for this 
	Dnn2piCylinder* DNN = new Dnn2piCylinder(_points, ignore_nearest_is_mirror, mt, verbose);
	//cout << "post mirror # clusters " << mt->GetNAllClusters() << endl;
	
	//need to make a distance map like in FastJet, but instead of clustering
	//based on geometric distance, we are using merge probability (posterior) from BHC
	//all three dimensions will go into calculating the probabilityu
	//but the map will be built only in 2D space
	//structure is typdef'ed in header
	CompareMap ProbMap, DistMap;
	InvCompareMap InvDistMap;
	//fill map with initial potential clusterings - including mirror pts added to MergeTree
	for(int i = 0; i < n; i++){
		_add_entry_to_maps(i, ProbMap, DNN);
		_add_entry_to_maps(i, DistMap, DNN, false);	
		_add_entry_to_maps(i, InvDistMap, DNN);	
	}
	bool done = false;
	std::pair<std::multimap<double,verts>::iterator, std::multimap<double,verts>::iterator> ret;	
	//run the clustering
	int faken = n;
	for(int i = 0; i < faken; i++){
		// find largest rk value in map (last entry)
		double BestRk;
		verts BestRkPair;
		std::multimap<double,verts>::iterator map_it;
		int jet_i, jet_j;
		bool Valid2;
		double dist;
		//cout << "Prob map size " << ProbMap.size() << endl;
		//cout << " dnn validity bool " << (!DNN->Valid(jet_i) || !Valid2) << " total bool " << ((!DNN->Valid(jet_i) || !Valid2) && ProbMap.size() > 0) << endl;
		if(ProbMap.size() == 0){done = true; break;}
		do{
			if(_verb > 1) cout << "probmap size " << ProbMap.size() << endl;
			map_it = ProbMap.end();
			map_it--;
			if(_verb > 1)cout << "updated to prob map pair " << map_it->second.first << " " << map_it->second.second << endl;
			BestRk = map_it->first;
			BestRkPair = map_it->second;
			//check for equal rks (more than three - may be two for i,j and j,i), break tie with 3d distance
			//right now - only equal rks are for equivalent merges
			//cout << "number of pairs with same bestRk = " << BestRk << ": " << ProbMap.count(BestRk) << endl;
			double mindist = 999;
			if(ProbMap.count(BestRk) > 2){
				ret = ProbMap.equal_range(BestRk);
				for(std::multimap<double,verts>::iterator it = ret.first; it != ret.second; ++it){
					int jjet_i = it->second.first;
					int jjet_j = it->second.second;
					if(InvDistMap.find(jjet_i) == InvDistMap.end()) continue; //skip points already combined
					dist = InvDistMap[jjet_i].second;
					if(_verb > 1) cout << "rk: " << BestRk << " with points: " << jjet_i << ", " << jjet_j << " and dist " << dist << endl;
					if(dist < mindist){
						mindist = dist;
						map_it = it;
						jet_i = jjet_i;
						jet_j = jjet_j;
					}
				}
			if(_verb > 1) cout << "mindist: " << mindist << " jet_i: " << jet_i << " jet_j: " << jet_j  << endl;
			//jet_i = DistMap.find(mindist)->second.first; jet_j = DistMap.find(mindist)->second.second;
		
			}else{
				jet_i = BestRkPair.first;
				jet_j = BestRkPair.second;
			}

			if (_verb > 1){ cout << "BayesCluster found recombination candidate: " << jet_i << " " << jet_j << " " << BestRk << " " << ProbMap.size() << endl;
			} // GPS debugging
 			//also need to erase any impossible merges from map too
			if(_verb > 1)cout << "erasing from prob map pair " << map_it->second.first << " " << map_it->second.second << endl;
			ProbMap.erase(map_it); //erase from InvDistMap too
			if(InvDistMap.find(jet_i) != InvDistMap.end()) InvDistMap.erase(InvDistMap.find(jet_i));
			if(InvDistMap.find(jet_j) != InvDistMap.end()) InvDistMap.erase(InvDistMap.find(jet_j));
			Valid2 = DNN->Valid(jet_j);
			if (_verb > 1) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << " prob map size " << ProbMap.size() << endl;
		} while((!DNN->Valid(jet_i) || !Valid2) && ProbMap.size() > 0); //this is what checks to see if merges are still allowed or if they include points that have already been merged
		//if point matches to itself (mirror point), find best geo match - this shouldn't happen...there is a safety in SetNearest and SetAndUpdateNearest in DnnPlane to skip calculating probabilties for points + their mirrors...
		if((jet_i == jet_j) && (ProbMap.size() > 1)){
			cout << "Uh oh best probability merger is jet to its mirror point. Returning for debugging..." << endl;
			vector<JetPoint> jps_i = _jets[jet_i].GetJetPoints();
			vector<JetPoint> jps_j = _jets[jet_j].GetJetPoints();
			cout << "jet_i pts" << endl;
			for(int i = 0; i < (int)jps_i.size(); i++){
				BayesPoint pt = BayesPoint({jps_i[i].eta(), jps_i[i].phi_02pi(), jps_i[i].t()});
				pt.SetWeight(jps_i[i].GetWeight());
				pt.Print();
			}
			cout << "jet_j pts" << endl;
			for(int i = 0; i < (int)jps_j.size(); i++){
				BayesPoint pt = BayesPoint({jps_j[i].eta(), jps_j[i].phi_02pi(), jps_j[i].t()});
				pt.SetWeight(jps_j[i].GetWeight());
				pt.Print();
			}
			return _trees;
		/*
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
			if (_verb > 1) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << " validity constraint " << (!DNN->Valid(jet_i) || !DNN->Valid(jet_j)) << " " << (DNN->Valid(jet_i) && DNN->Valid(jet_j)) << endl;
				//check validity (to see if it has been clustered before)
				//if best match is between mirrored points AND there are no points that haven't been geometrically clustered yet,
				//we're done
				//also need to erase any impossible merges from map too
				if(DNN->Valid(jet_i) && DNN->Valid(jet_j)){ done = true; break; 
				DistMap.erase(map_it_dist);}
				else continue;
			} 
		*/
		}
		//cout << "BestRk " << BestRk << endl;
		//if max rk < 0.5, can stop clustering
		if(BestRk < 0.5){ done = true; if(_verb > 0) cout << "stop with BestRk " << BestRk << " for combo " << jet_i << " + " << jet_j << endl; break; }	
		//if either sides of best recombination candidate found is not valid - break
                if((!DNN->Valid(jet_i) || !Valid2)){done = true; if(_verb > 0) cout << "best recomb candidate not valid + prob map exhausted - stop" << endl; break;}

		int nn;
		if(_verb > 1){
			cout << "BayesCluster call _do_ij_recomb: " << jet_i << " " << jet_j << " " << BestRk << endl << " with points " << endl;
			vector<JetPoint> jps_i = _jets[jet_i].GetJetPoints();
			vector<JetPoint> jps_j = _jets[jet_j].GetJetPoints();
			cout << "jet_i pts" << endl;
			for(int i = 0; i < (int)jps_i.size(); i++){
				BayesPoint pt = BayesPoint({jps_i[i].eta(), jps_i[i].phi_02pi(), jps_i[i].t()});
				pt.SetWeight(jps_i[i].GetWeight());
				pt.Print();
			}
			cout << "jet_j pts" << endl;
			for(int i = 0; i < (int)jps_j.size(); i++){
				BayesPoint pt = BayesPoint({jps_j[i].eta(), jps_j[i].phi_02pi(), jps_j[i].t()});
				pt.SetWeight(jps_j[i].GetWeight());
				pt.Print();
			}
		}

		//do_ij_recomb - this should be the same as in the OG code (except rk instead of dij)
      		_do_ij_recombination_step(jet_i, jet_j, BestRk, nn);
		// exit the loop because we do not want to look for nearest neighbours
		// etc. of zero partons
		if(_verb > 1) cout << "BayesCluster - i " << i << " n " << n << endl;
		if (i == n-1) {cout << "i " << i << " n " << n << " i == n-1 - breaking" << endl;break;}//get eta phi (and time!) of new point - centroid of points in combined vertices
		int pt3;
		vector<int> updated_neighbors;
		//update DNN with RemoveCombinedAddCombination
		//this should also update the merge tree - RemoveAndAddPoints in DnnPlane does
		vector<JetPoint> jps = _jets[_jets.size() - 1].GetJetPoints();
		PointCollection newpts = PointCollection();
		for(int i = 0; i < (int)jps.size(); i++){
			//cout << "adding pt to new cluster - eta " << jps[i].eta() << " raw phi " << jps[i].phi() << " time " << jps[i].t() << endl;
			BayesPoint pt = BayesPoint({jps[i].eta(), jps[i].phi_02pi(), jps[i].t()});
			pt.SetWeight(jps[i].GetWeight());
			newpts += pt;
		}
		if(_verb > 1)cout <<"remove combined add combination start" << endl;
		DNN->RemoveCombinedAddCombination(jet_i, jet_j,
							newpts, pt3, updated_neighbors);
		if(_verb > 1)cout <<"remove combined add combination done\n" << endl;
		//cout << "newpts" << endl;
		//newpts.Print(); 
		if(_verb > 1)cout << "\n\n\n" << endl;
		if(_verb > 1) cout << "updating map: adding new cluster " << pt3 << " = " << jet_i << " + " << jet_j << endl;
		//update map
		vector<int>::iterator it = updated_neighbors.begin();
		for(; it != updated_neighbors.end(); ++it){
			int ii = *it;
			_add_entry_to_maps(ii, ProbMap, DNN);
			_add_entry_to_maps(ii, DistMap, DNN, false);	
			_add_entry_to_maps(ii, InvDistMap, DNN);
		}
		if(_verb > 0) cout << mt->GetNClusters() << " current clusters at iteration " << i << endl;
		if(_verb > 0) cout << "\n" << endl;
	}

	vector<node*> trees = mt->GetClusters();
	double nmirror = 0;
	double nnull = 0; 	
	map<string,Matrix> params;
	for(int i = 0; i < trees.size(); i++){
		//ignore removed (clustered) points and don't plot mirror points
		if(trees[i] == nullptr){ nnull++; continue; }
		//if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0){
		if(trees[i]->ismirror){
			nmirror++;
				//cout << "\nMIRROR tree model params " << trees[i]->model->GetNClusters() << " clusters and "<< trees[i]->model->GetData()->GetNPoints() << " points and " << trees[i]->model->GetData()->Sumw() << " weight" << endl;
			//cout << " points" << endl; trees[i]->model->GetData()->Print();
			continue; }
		if(trees[i]->points->GetNPoints() < 2 || trees[i]->points->Sumw() < _thresh) continue;

		//cout << "getting " << _trees[i]->points->GetNPoints() << " " << _trees[i]->model->GetData()->GetNPoints() << " points in cluster #" << i << endl;
		//cout << trees[i]->l->points->GetNPoints() << " in left branch " << trees[i]->r->points->GetNPoints() << " in right branch" << endl;
		_trees.push_back(trees[i]);
		//cout << "\nREAL tree model params " << trees[i]->model->GetNClusters() << " clusters and "<< trees[i]->model->GetData()->GetNPoints() << " points and " << trees[i]->model->GetData()->Sumw() << " weight" << endl;
		//cout << " points" << endl; trees[i]->model->GetData()->Print();
		cout << " tree has subclusters " << endl;
		for(int k = 0; k < trees[i]->model->GetNClusters(); k++){
			params = trees[i]->model->GetLHPosteriorParameters(k);
			cout << " k " << k << " center " << endl; params["mean"].Print();
			//cout << "cov " << endl; params["cov"].Print();
		}
		cout << trees[i]->points->GetNPoints() << " points for jet " << i << " with " << trees[i]->model->GetNClusters() << " subclusters" << endl; trees[i]->model->GetData()->Print();
		if(trees[i]->mirror != nullptr){
			cout << " tree has mirror node with subclusters " << endl;
			for(int k = 0; k < trees[i]->mirror->model->GetNClusters(); k++){
				params = trees[i]->mirror->model->GetLHPosteriorParameters(k);
				cout << " k " << k << " center " << endl; params["mean"].Print();
			}
		}
	}
	//cout << " all points" << endl;
	//for (int i = 0; i < n; i++) {	_points[i].Print(); }
	cout << _trees.size() << " clustered trees " << trees.size() << " found trees " << nnull << " null trees " << nmirror << " mirror trees" << endl;
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
	//set time resolution smear: c^2 + n^2/e^2
	//remember time is already in ns
	//e = w/gev
	if(_thresh != -999) mt->SetThresh(_thresh);
	_prior_params["scale"] = Matrix(1e-3);
	_prior_params["dof"] = Matrix(3);
	Matrix W(3,3);
	double r = 0.4;
	W.SetEntry(1./(r*r/3),0,0);
	W.SetEntry(1./(r*r/3),1,1);
	W.SetEntry(1./(0.1*0.1/3),2,2);
	_prior_params["scalemat"] = W;
	_prior_params["mean"] = Matrix(3,1);
	mt->SetPriorParameters(_prior_params);
	//set distance constraint
	//mt->SetDistanceConstraint(0,acos(-1)/2.);
	int n = _points.size();
	if(_verb > 0) cout << "Creating rk list for all pairings (this might take a minute...)" << endl;	
	for (int i = 0; i < n; i++) {
		//should only be one point per entry in points
		if(_points[i].GetNPoints() != 1){
			cerr << "BayesCluster - Error: multiple points in one collection of starting vector. " << _points[i].GetNPoints() << " points in collection " << i << "." << endl;
			return _trees;
		} 
		mt->AddLeaf(&_points[i].at(0));
	//experiment with not adding mirror points for leaves
	//mirror points are calculated here and CalculateMerge because nndist is set with DistanceConstraint (called in calculate merge)
	//so here nndist for all leaves is set to default
	//in order to set nn3dist for all leaves, that would be another n^2 operation
	//see if cluster phi's are enough so you only do n^2 once
	//	mt->CreateMirrorNode(mt->Get(i));
	}

	int it = 0;
	//loop through all possible pairs to create list of merge probabilites r_{i+j} = r_k
	//O(N^2)
	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n; j++){
			di = mt->Get(i);
			dj = mt->Get(j);
			node* x = mt->CalculateMerge(di, dj);
		if(_verb > 1){
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
		if(_verb > 1){ cout << "max merge is " << max->l->idx << " and " << max->r->idx << " rk: " << max->val << " " << endl;
		max->l->points->Print();
		cout << "and" << endl;
		max->r->points->Print();
	//	cout << "\n" << endl;
		}
		//if maximum merge probability is less than 0.5, stop clustering (not probabilistically favored merges left)
                if(max->val < maxrk){
                        cout << "reached min rk = " << max->val << " <  " << maxrk << " - final iteration: " << it <<  " - " << mt->GetNClusters() << " clusters" << endl;
                        break;
                }
	
	
	//	mt->Merge(max->l, max->r);		
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
		if(_verb > 1){
		cout << "remaining possible merges" << endl;
		//list.Print(1);
		}

		//create new merges with the remaining nodes and the newly formed cluster
		//this operation is O(N)
		if(_verb > 1) cout << "new merges - " << mt->GetNClusters() << " remaining clusters" << endl;
		for(int i = 0; i < mt->GetNAllClusters(); i++){
			node* dl = mt->Get(i);
			//make sure this node exists
			if(dl == nullptr) continue;
			if(_verb > 2) cout << "checking pair " << i << " or " << dl->idx << " and max " << max->idx << endl;
			//make sure this is not the max merge cluster - doesn't need to be paired with itself
			if(dl == max) continue;
			//calculate merge for node i and max node
			node* x = mt->CalculateMerge(dl, max);
		if(_verb > 1){
		cout << std::setprecision(10) << endl;
		cout << "checking nodes " << dl->idx << ": ";
		//dl->points->Print();
		cout << " and max " << max->idx << ": ";
		//max->points->Print();
		cout << "this rk: " << x->val << " with " << dl->points->GetNPoints() + max->points->GetNPoints() << " total points\n" << endl;
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
			continue; }
		if(_verb > 0){
			cout << "getting " << _trees[i]->points->GetNPoints() << " points in cluster #" << i << endl;
			_trees[i]->points->Print();
		}
	}
	if(_verb > 0) cout << nmirror << " mirror points." << endl;
	cout << mt->GetNClusters() - nmirror << " jets found." << endl;
	return _trees;

}



GaussianMixture* BayesCluster::_subcluster(string oname){
	vector<double> ws;
	_jets[0].GetWeights(ws);
	double gev = ws[0]/_jets[0].E();
	//create GMM model
	PointCollection* points = new PointCollection();
	for(auto point : _points){
		points->AddPoints(point);
	}
	int maxK = points->GetNPoints();
	

	GaussianMixture* gmm = new GaussianMixture(maxK);
	//double tresCte = 0.2 * 1e-9; //time resolution for CMS ECAL (s) (200 ps)
	//double tresStoch = 0.34641 * 1e-9; //rate of time res that gives 400 ps at E = 1 GeV (in [GeV*s])
	//tresStoch = 2.4999200e-05; //rate of time res that gives 5 ns at E = 5 GeV (in [GeV*s])
	//tresStoch *= gev;
	//needs to be before SetData bc thats when the measurement errors are set
	if(_verb > 6) cout << "BayesCluster subcluster - Using tresCte " << _tresCte << " _tresStoch " << _tresStoch << " _tresNoise " << _tresNoise << " gev " << gev << " _cell " << _cell << endl;
	gmm->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
	gmm->SetData(points);
	//cout << "1 - set gmm data as" << endl; gmm->GetData()->Print();
	//needs to be before ScaleData() bc this method also scales the smear
	if(!_smear.empty()){ gmm->SetDataSmear(_smear); }
	gmm->SetAlpha(_subalpha);
	gmm->SetVerbosity(_verb);
	gmm->InitParameters();
	gmm->InitPriorParameters();
	//set prior parameters
	gmm->SetPriorParameters(_prior_params);	


	//cout << "old points w/ wraparound" << endl;
	//points->Print();

	//transform points into local coordinates
	//for GMM parameter estimation
	//use weighted mean as center to be set to 0 point
	//zero points by energy centroid
	//x' = x - a
	BayesPoint center({points->Centroid(0), points->CircularCentroid(1), points->Centroid(2)});
	gmm->ShiftData(center);
	//cout << "centroid " << endl; center.Print();
	
	//cout << "translated pts" << endl;
	//cout << "shifted - gmm data is" << endl; gmm->GetData()->Print();
	//gmm->GetData()->at(0).Print();

	//scale points s.t. 1 cell ~ 0.0174 = 1 unit in eta-phi
	//x'' = x'/b = (x-a)/b
	//sets relative importance of dimensions
	//decreasing cell -> eta/phi distance more important
	//increasing entry (2,2) -> time distance more important
	double cell = 4*atan(1)/180;
	Matrix Rscale(3,3);
	Rscale.SetEntry(1/cell,0,0);
	Rscale.SetEntry(1/cell,1,1);
	Rscale.SetEntry(1,2,2);

	//assuming prior parameters are given in shifted + scaled frame
	gmm->ScaleData(Rscale);
	//cout << "scaled - gmm data is" << endl; gmm->GetData()->Print();
	//cout << "scaled points" << endl;
	//gmm->GetData()->at(0).Print();
	
	
	//create EM algo
	VarEMCluster* algo = new VarEMCluster(gmm,maxK);
	algo->SetThresh(_thresh);



	map<string, vector<Matrix>> params;
	bool viz = false;
	if(!oname.empty()){
		viz = true;
	}

	//inverse transformations
	Matrix RscaleInv;
	RscaleInv.invert(Rscale);
	center.Scale(-1);
	
	VarClusterViz3D cv3D(algo);
	cv3D.SetScale(RscaleInv);
	cv3D.SetShift(center);
	if(isnan(gev) || isinf(gev)){ cout << "Energy-weighting scale factor " << gev << " is not valid." << endl; return gmm; }
	if(viz){
		cv3D.SetVerbosity(_verb);
		cv3D.SetTransfFactor(gev);
		cv3D.UpdatePosterior();
		cv3D.WriteJson(oname+"/it0");
	}
	//loop
	double dLogL, newLogL;
	double LogLthresh = 1e-10;
	double oldLogL = algo->EvalLogL();
	////////run EM algo////////
	//maximum of 50 iterations
	for(int it = 0; it < 50; it++){
	
		algo->SetClusterStart();
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
			return gmm;
		}
		dLogL = oldLogL - newLogL;
		if(_verb > 3) cout << "iteration #" << it+1 << " log-likelihood: " << newLogL << " dLogL: " << dLogL << endl;
		if(fabs(dLogL) < LogLthresh){// || dLogL > 0){
			if(_verb > 2){
				cout << "Reached convergence at iteration " << it+1 << endl;
			}
			break;
		}
		oldLogL = newLogL;
	}

	//need to unscale first then uncenter since x'' = (x-a)/b (see above)
	//need to unscale data 
	gmm->ScaleData(RscaleInv);	
	//cout << "unscaled - gmm data is" << endl; gmm->GetData()->Print();
	//need to unscale mean + covariances
	gmm->ScaleParameters(RscaleInv);	
	
	//need to put GMM parameters AND points back in detector coordinates (not local)
	gmm->ShiftData(center);
	//cout << "unshifted - gmm data is" << endl; gmm->GetData()->Print();
	gmm->ShiftParameters(center);
	gmm->PutPhi02pi(); //does for data and parameters - do after data + parameters shift so the [0,2pi] transformation doesn't get shifted
	//cout << "phi02pi - gmm data is" << endl; gmm->GetData()->Print();
	//cout << "center " << endl; center.Print();
	//cout << "predicted center - lead only - nclusters " << gmm->GetNClusters() << endl;
	if(_verb > 2){
		cout << std::setprecision(10) << endl;
		cout << "Estimated parameters" << endl;
		vector<double> norms;
		gmm->GetNorms(norms);
		for(int k = 0; k < gmm->GetNClusters(); k++){
			auto params1 = gmm->GetLHPosteriorParameters(k);
			auto params2 = gmm->GetDataStatistics(k);
			cout << "weight " << k << ": " << params1["pi"].at(0,0) << " alpha: " << params1["alpha"].at(0,0) << " eff evts: " << norms[k] << endl;
			cout << "subcl #" << k << endl;	
			cout << "mean " << endl;
			params1["mean"].Print();
			cout << "r stat mean" << endl; params2["mean"].Print();
			cout << "cov " << endl;
			params1["cov"].Print();
			cout << "r stat cov" << endl; params2["cov"].Print();
			cout << endl;
		}
	}




	return gmm;

}



void BayesCluster::_phi_wraparound(PointCollection& pc){
	BayesPoint mean = pc.mean();
	double pi = acos(-1);

	//calculate max distance from mean for all points (O(n))
	double maxdist = 0;
	double dphi;
	for(int i = 0; i < pc.GetNPoints(); i++){
		dphi = fabs(pc.at(i).at(1) - mean.at(1));
		if(dphi > maxdist) maxdist = dphi;
	} 

	//if maxdist is less than pi/2. return
	if(maxdist < pi/2.) return;

	//else, mirror all points above pi to below 0
	PointCollection mirrorpts;
	for(int i = 0; i < pc.GetNPoints(); i++){
		BayesPoint pt = pc.at(i);
		if(pt.at(1) > pi){ 
			pt.SetValue(2*pi - pt.at(1),1);
		}
		mirrorpts.AddPoint(pt);
	}
	//recalculate max distance from mean
	mean = mirrorpts.mean();
	double newmaxdist = 0;
	for(int i = 0; i < mirrorpts.GetNPoints(); i++){
		dphi = fabs(mirrorpts.at(i).at(1) - mean.at(1));
		if(dphi > newmaxdist) newmaxdist = dphi;
	} 

	//if new maxdist < maxdist, keep mirrored points
	if(newmaxdist < maxdist) pc = mirrorpts;
	//else keep original points
	return;

	//else unmirror points
}

void BayesCluster::_add_entry_to_maps(const int i, InvCompareMap& inmap, const Dnn2piCylinder* DNN){
		double dist;
		int j;
		dist = DNN->NearestNeighbourDistance(i);
		j = DNN->NearestNeighbourIndex(i);
if(_verb > 1)cout << "adding entry " << i << " " << j << " with dist " << dist << " to inv map" << endl;
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
			if(_verb > 1) cout << std::setprecision(20) << "adding entry " << i << " with best probability " << comp << " pair " << j << " cyl index: " << cyl_j << std::setprecision(5) << endl;
		}
		else{
			comp = DNN->NearestNeighbourDistance(i);
			j = DNN->NearestNeighbourIndex(i);
			if(_verb > 1) cout << std::setprecision(20) << "adding entry " << i << " with best distance " << comp << " pair " << j << std::setprecision(5) << endl;
		}
		inmap.insert(CompEntry(comp,verts(i,j)));
}

