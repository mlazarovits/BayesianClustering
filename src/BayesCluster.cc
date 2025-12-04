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
void BayesCluster::_delauney_cluster(vector<std::shared_ptr<BaseTree::node>>& trees){
	//cout << "BayesCluster::_delauney_cluster - start" << endl;
	trees.clear();
	vector<double> ws;
	_jets[0].GetWeights(ws);
	double gev = ws[0]/_jets[0].E();
	//the 2D Delauney triangulation used in FastJet will be used to seed the 3D clustering
	std::unique_ptr<MergeTree> mt = std::make_unique<MergeTree>(_alpha);
	mt->SetSubclusterAlpha(_subalpha);
	mt->SetVerbosity(_verb);
	mt->CheckMerges(_check_merges);
	mt->SetNGhosts(_nGhosts);
	//set data smear
	if(!_smear.empty()) mt->SetDataSmear(_smear);
	//set time resolution smear: c^2 + n^2/e^2
	//remember time is already in ns
	//e = w/gev
	if(_thresh != -999) mt->SetThresh(_thresh);
	if(_verb > -1){
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
	}
	mt->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
	mt->SetPriorParameters(_prior_params);
	int n = _jets.size();
	if(n < 2){
		cout << "ERROR: only have 1 pt for clustering - returning vector of nulls" << endl;
		return;
	}
	if(_verb > -1) cout << "n starting pts " << n << endl;
	for(int j = 0; j < _jets.size(); j++){
		std::shared_ptr<BaseTree::node> x = _jet_to_leaf_node(_jets[j], mt.get());
		_nodes.push_back(std::move(x));
	}
	const bool verbose = false;
	const bool ignore_nearest_is_mirror = true; //based on _Rparam < twopi, should always be true for this 
	std::unique_ptr<Dnn2piCylinder> DNN = std::make_unique<Dnn2piCylinder>(_nodes, ignore_nearest_is_mirror, std::move(mt), verbose);
	//need to make a distance map like in FastJet, but instead of clustering
	//based on geometric distance, we are using merge probability (posterior) from BHC
	//all three dimensions will go into calculating the probabilityu
	//but the map will be built only in 2D space
	//structure is typdef'ed in header
	CompareMap ProbMap, DistMap;
	//fill map with initial potential clusterings - including mirror pts added to MergeTree
	for(int i = 0; i < n; i++){
		_add_entry_to_maps(i, ProbMap, DNN.get());
		_add_entry_to_maps(i, DistMap, DNN.get(), false);	
	}
	bool done = false;
	std::pair<std::multimap<double,verts>::iterator, std::multimap<double,verts>::iterator> ret;	
	//run the clustering
	int faken = n;
	for(int i = 0; i < faken; i++){
		// find largest rk value in map (last entry)
		double BestRk;
		verts BestRkPair;
		std::multimap<double,verts>::iterator map_it, map_it_i, map_it_j;
		int jet_i, jet_j;
		bool Valid2;
		double dist_i, dist_j;
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
			map_it_i = ProbMap.end();
			map_it_j = ProbMap.end();
			if(ProbMap.count(BestRk) > 2){
				ret = ProbMap.equal_range(BestRk);
				if(_verb > 1) cout << ProbMap.count(BestRk) << " number of pairs with the same BestRk = " << BestRk << endl;
				double mindist = 999;
				for(std::multimap<double,verts>::iterator it = ret.first; it != ret.second; ++it){
					int jjet_i = it->second.first;
					int jjet_j = it->second.second;
					if(_verb > 1) cout << std::setprecision(5) << "rk: " << BestRk << " with points: " << jjet_i << " and " << jjet_j << endl;
					double deta = _jets[jjet_i].eta() - _jets[jjet_j].eta();
					double dphi = _jets[jjet_i].phi() - _jets[jjet_j].phi();
					dphi = acos(cos(dphi));
					double dist = sqrt(deta*deta + dphi*dphi);
					if(dist < mindist){
						mindist = dist;
						jet_i = jjet_i;
						jet_j = jjet_j;
						map_it = it;	
					}
				}
			if(_verb > 1) cout << "mindist: " << mindist << " jet_i: " << jet_i << " jet_j: " << jet_j  << endl;
		
			}else{
				jet_i = BestRkPair.first;
				jet_j = BestRkPair.second;
			}

			if(_verb > 1){ cout << "BayesCluster found recombination candidate: " << jet_i << " " << jet_j << " " << BestRk << " " << ProbMap.size() << endl;} // GPS debugging
 			//also need to erase any impossible merges from map too
			if(_verb > 1)cout << "erasing from prob map pair " << jet_i << " " << jet_j << " from map it " << map_it->second.first << " " << map_it->second.second << endl;
			if(map_it != ProbMap.end()){ if(_verb > 1) cout << "erasing it for jets " << map_it->second.first << " " << map_it->second.second << endl; ProbMap.erase(map_it);}
			//not guaranteed that pair that is merged is corresponding pair in prob map
			if(map_it_i != ProbMap.end()){ if(_verb > 1)cout << "erasing jet_i " << jet_i << " from prob map " << endl; ProbMap.erase(map_it_i);}
			if(map_it_j != ProbMap.end()){ if(_verb > 1)cout << "erasing jet_j " << jet_j << " from prob map " << endl; ProbMap.erase(map_it_j);}
			Valid2 = DNN->Valid(jet_j);
			if(_verb > 1) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << " prob map size " << ProbMap.size() << endl;
		} while((!DNN->Valid(jet_i) || !Valid2) && ProbMap.size() > 0); //this is what checks to see if merges are still allowed or if they include points that have already been merged
		//if point matches to itself (mirror point), find best geo match - this shouldn't happen...there is a safety in SetNearest and SetAndUpdateNearest in DnnPlane to skip calculating probabilties for points + their mirrors...
		if((jet_i == jet_j) && (ProbMap.size() > 1)){
			cout << "Uh oh best probability merger is jet to its mirror point. Returning for debugging..." << endl;
			vector<JetPoint> jps_i, jps_j;
			_jets[jet_i].GetJetPoints(jps_i);
			_jets[jet_j].GetJetPoints(jps_j);
			PointCollection jeti_pts, jetj_pts;
			for(int i = 0; i < (int)jps_i.size(); i++){
				BayesPoint pt = BayesPoint({jps_i[i].eta(), jps_i[i].phi_02pi(), jps_i[i].t()});
				pt.SetWeight(jps_i[i].GetWeight());
				pt.SetUserIdx(jps_i[i].rhId());
				if(jps_i[i].InvalidTime())
					pt.SetSkipDim(2);
				_sanitize(pt);
				jeti_pts += pt;
			}
			if(_verb > 3) cout << "# jet_i pts " << jeti_pts.GetNPoints() << endl;
			//jeti_pts.Print();
			BayesPoint jeti_mean = BayesPoint({jeti_pts.mean().at(0), jeti_pts.CircularMean(1), jeti_pts.mean().at(2)});
			if(_verb > 3){ cout << "with mean " << endl; jeti_mean.Print(); }
			for(int i = 0; i < (int)jps_j.size(); i++){
				BayesPoint pt = BayesPoint({jps_j[i].eta(), jps_j[i].phi_02pi(), jps_j[i].t()});
				pt.SetWeight(jps_j[i].GetWeight());
				pt.SetUserIdx(jps_j[i].rhId());
				if(jps_j[i].InvalidTime())
					pt.SetSkipDim(2);
				_sanitize(pt);
				jetj_pts += pt;
			}
			if(_verb > 3) cout << "# jet_j pts " << jetj_pts.GetNPoints() << endl;
			//jetj_pts.Print();
			BayesPoint jetj_mean = BayesPoint({jetj_pts.mean().at(0), jetj_pts.CircularMean(1), jetj_pts.mean().at(2)});
			if(_verb > 3){ cout << "with mean " << endl; jetj_mean.Print();}
	
			return;
		}
		if(_verb > 1) cout << "BayesCluster - current best rk " << BestRk << endl;
		//if max rk < 0.0 (0.5 without computational/numerical trickery), can stop clustering
		if(BestRk <= 0.){ 
			done = true; 
			if(_verb > -1) cout << "stop with BestRk " << BestRk << " for combo " << jet_i << " + " << jet_j << " with phis " << _jets[jet_i].phi() << " and " << _jets[jet_j].phi() << " probmap size " << ProbMap.size() << endl; 
			if(_verb > 3){
				for(auto it = ProbMap.begin(); it != ProbMap.end(); it++){
					int i = it->second.first;
					int j = it->second.second;
					cout << "remaining best rk " << it->first << " for jets " << i << " and " << j << " with phis " << _jets[i].phi_02pi() << " and " << _jets[j].phi_02pi() << " and eta " << _jets[i].eta() << " and " << _jets[j].eta() << " and # rhs " << _jets[i].GetNRecHits() << " and " << _jets[j].GetNRecHits() << endl;

				}
			} 

		//_bestRk = BestRk; 
		break; }	
		//if either sides of best recombination candidate found is not valid - break
                if((!DNN->Valid(jet_i) || !Valid2)){done = true; if(_verb > 0) cout << "best recomb candidate not valid + prob map exhausted - stop" << endl; break;}

		int nn;
		if(_verb > 1){
			cout << "BayesCluster call _do_ij_recomb: " << jet_i << " " << jet_j << " " << BestRk << " with points " << endl;
			/*
			vector<JetPoint> jps_i, jps_j;
			_jets[jet_i].GetJetPoints(jps_i);
			_jets[jet_j].GetJetPoints(jps_j);
			PointCollection jeti_pts, jetj_pts;
			for(int i = 0; i < (int)jps_i.size(); i++){
				BayesPoint pt = BayesPoint({jps_i[i].eta(), jps_i[i].phi_02pi(), jps_i[i].t()});
				pt.SetWeight(jps_i[i].GetWeight());
				pt.SetUserIdx(jps_j[i].rhId());
				if(jps_j[i].InvalidTime())
					pt.SetSkipDim(2);
				_sanitize(pt);
				jeti_pts += pt;
			}
			if(_verb > 1){ cout << "# jet_i pts " << jeti_pts.GetNPoints() << endl; jeti_pts.Print();}
			BayesPoint jeti_mean = BayesPoint({jeti_pts.mean().at(0), jeti_pts.CircularMean(1), jeti_pts.mean().at(2)});
			if(_verb > 3){ cout << "with mean " << endl; jeti_mean.Print();}
			for(int i = 0; i < (int)jps_j.size(); i++){
				BayesPoint pt = BayesPoint({jps_j[i].eta(), jps_j[i].phi_02pi(), jps_j[i].t()});
				pt.SetWeight(jps_j[i].GetWeight());
				pt.SetUserIdx(jps_j[i].rhId());
				if(jps_j[i].InvalidTime())
					pt.SetSkipDim(2);
				_sanitize(pt);
				jetj_pts += pt;
			}
			if(_verb > 1){ cout << "# jet_j pts " << jetj_pts.GetNPoints() << endl; jetj_pts.Print();}
			BayesPoint jetj_mean = BayesPoint({jetj_pts.mean().at(0), jetj_pts.CircularMean(1), jetj_pts.mean().at(2)});
			if(_verb > 3){cout << "with mean " << endl; jetj_mean.Print();}
			*/
		}
		//do_ij_recomb - this should be the same as in the OG code (except rk instead of dij)
		_do_ij_recombination_step(jet_i, jet_j, BestRk, nn, DNN.get());
		// exit the loop because we do not want to look for nearest neighbours
		// etc. of zero partons
		if(_verb > 1) cout << "BayesCluster - i " << i << " n " << n << endl;
		if (i == n-1) {
			if(_verb > -1) cout << "i " << i << " n " << n << " i == n-1 - breaking" << endl;
			break;}//get eta phi (and time!) of new point - centroid of points in combined vertices
		int pt3;
		vector<int> updated_neighbors;
		//update DNN with RemoveCombinedAddCombination
		//this should also update the merge tree - RemoveAndAddPoints in DnnPlane does
		/*
		vector<JetPoint> jps;
		_jets[_jets.size() - 1].GetJetPoints(jps);
		PointCollection newpts = PointCollection();
		for(int i = 0; i < (int)jps.size(); i++){
			BayesPoint pt = BayesPoint({jps[i].eta(), jps[i].phi_02pi(), jps[i].t()});
			_sanitize(pt);
			if(jps[i].InvalidTime())
				pt.SetSkipDim(2);
			pt.SetUserIdx(jps[i].rhId());
			pt.SetWeight(jps[i].GetWeight());
			newpts += pt;
		}
		if(newnode->points->HasInf(0) || newnode->points->HasInf(1)){
			cout << "BayesCluster - ERROR adding jet with inf" << endl;
		}
		*/
		if(_verb > 1)cout <<"remove combined add combination start" << endl;
		DNN->RemoveCombinedAddCombination(jet_i, jet_j,
							_nodes[nn], pt3, updated_neighbors);
		if(_verb > 1)cout <<"remove combined add combination done\n" << endl;
		if(_verb > 1)cout << "\n\n\n" << endl;
		if(_verb > 1) cout << "updating map: adding new cluster " << pt3 << " = " << jet_i << " + " << jet_j << endl;
		//update map
		vector<int>::iterator it = updated_neighbors.begin();
		for(; it != updated_neighbors.end(); ++it){
			int ii = *it;
			_add_entry_to_maps(ii, ProbMap, DNN.get());
			_add_entry_to_maps(ii, DistMap, DNN.get(), false);	
		}
		if(_verb > 0) cout << "\n" << endl;
	}

	vector<std::shared_ptr<BaseTree::node>> mt_trees;
	DNN->GetValidNodes(mt_trees);
	double nmirror = 0;
	double nnull = 0; 	
	map<string,Matrix> params;
	for(int i = 0; i < mt_trees.size(); i++){
		//ignore removed (clustered) points and don't plot mirror points
		if(mt_trees[i] == nullptr){ nnull++; continue; }
		//if(mt_trees[i]->points->mean().at(1) > 2*acos(-1) || mt_trees[i]->points->mean().at(1) < 0){
		if(mt_trees[i]->ismirror){
			nmirror++;
			continue; }
		if(mt_trees[i]->points->GetNPoints() < 2 || mt_trees[i]->points->Sumw() < _thresh) continue;
		if(_verb > -1){
			cout << " tree has subclusters " << endl;
			vector<double> norms;
			mt_trees[i]->model->GetNorms(norms);
			for(int k = 0; k < mt_trees[i]->model->GetNClusters(); k++){
				std::map<string, Matrix> params;
				mt_trees[i]->model->GetLHPosteriorParameters(k, params);
				cout << " k " << k << " weight " << norms[k] << " center " << endl; params["mean"].Print();
				cout << "cov " << endl; params["cov"].Print();
			}
			cout << mt_trees[i]->points->GetNPoints() << " " << mt_trees[i]->model->GetData()->GetNPoints() << " points for jet " << i << " with " << mt_trees[i]->model->GetNClusters() << " subclusters" << endl; mt_trees[i]->points->Print(); 
			BayesPoint center({mt_trees[i]->points->Centroid(0), mt_trees[i]->points->CircularCentroid(1), mt_trees[i]->points->Centroid(2)});
			cout << "with centroid" << endl; center.Print();

			if(mt_trees[i]->mirror != nullptr){
				cout << " tree has mirror node with subclusters " << endl;
				for(int k = 0; k < mt_trees[i]->mirror->model->GetNClusters(); k++){
					std::map<string, Matrix> params;
					mt_trees[i]->mirror->model->GetLHPosteriorParameters(k, params);
					cout << " k " << k << " center " << endl; params["mean"].Print();
				}
			}
			cout << endl;
		}
		trees.push_back(std::move(mt_trees[i]));
	}
	//mt->avg_time();
	//cout << " all points" << endl;
	//for (int i = 0; i < n; i++) {	_points[i].Print(); }
	if(_verb > -1) cout << trees.size() << " clustered trees " << trees.size() << " found trees " << nnull << " null trees " << nmirror << " mirror trees" << endl;
	
	
	//cout << "BayesCluster::_delauney_cluster - end" << endl;
}

/*

void BayesCluster::_naive_cluster(vector<std::shared_ptr<BaseTree::node>>& trees){
	trees.clear();
        node* di; node* dj;
        double rk;
        double rk_max;
	NodeStack list; //for sorting nodes by merge probability and tracking potential merges
	
	//set up merge tree that will calculate merges and track nodes
	std::unique_ptr<MergeTree> mt = std::make_unique<MergeTree>(_alpha); //keeps track of merges made + merge history
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
			return;
		} 
		mt->AddLeaf(_points[i].at(0));
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
			std::shared_ptr<node> x = mt->CalculateMerge(di, dj);
		if(_verb > 1){
		cout << "checking nodes " << i << ": ";
		di->points->Print();
		cout << " and " << j << ": ";
		dj->points->Print();
		cout << "this rk: " << x->val <<  "\n" << endl;
		}	
			list.insert(std::move(x));
		}
	}
	//do clustering
	double maxrk = 0.5;
	while(mt->GetNClusters() > 1){
                if(_verb > 0) cout << "---------- iteration: " << it << " ----------" << endl;
                //if(_verb > 0) list.Print();
		
		list.sort();
	
		//get node with maximum merge probability + remove impossible merges because of max node
		shared_ptr<node> max = list.fullpop();
		//check that popped node is not null
		if(max == nullptr) break;
		//if maximum merge probability is less than 0.5, stop clustering (not probabilistically favored merges left)
                if(max->val < maxrk){
                        cout << "reached min rk = " << max->val << " <  " << maxrk << " - final iteration: " << it <<  " - " << mt->GetNClusters() << " clusters" << endl;
                        break;
                }
	
		//if this merge clusters all avaible points, stop clustering
		if(max->points->GetNPoints() == n){
                        cout << "All points clustered. Exiting clustering..." << endl;
                        break;
		}	
	
		//update merge tree with selected merge
		mt->Remove(max->l.get());
		mt->Remove(max->r.get());
		//make new mirror node for max if necessary
		mt->CreateMirrorNode(max.get());
		mt->Insert(std::move(max));
			

		//print remaining possible merges
		if(_verb > 1){
		cout << "remaining possible merges" << endl;
		//list.Print(1);
		}

		//create new merges with the remaining nodes and the newly formed cluster
		//this operation is O(N)
		if(_verb > 1) cout << "new merges - " << mt->GetNClusters() << " remaining clusters" << endl;
		for(int i = 0; i < mt->GetNAllClusters()-1; i++){ //skip last inserted node (max)
			node* dl = mt->Get(i);
			//make sure this node exists
			if(dl == nullptr) continue;
			if(_verb > 2) cout << "checking pair " << i << " or " << dl->idx << " and max " << max->idx << endl;
			//make sure this is not the max merge cluster - doesn't need to be paired with itself
			//calculate merge for node i and max node
			unique_ptr<node> x = mt->CalculateMerge(dl, mt->Get(mt->GetNAllClusters()-1));
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
			list.insert(std::move(x)); 
		}


		it++;

	}

	if(_verb > 0) cout << "---------- final iteration: " << it <<  " - " << mt->GetNClusters() << " clusters ----------" << endl;
        vector<unique_ptr<BaseTree::node>> mt_trees = mt->GetClusters();
	cout << mt->GetNClusters() << " final clusters" << endl;
	double nmirror = 0;
	for(int i = 0; i < mt_trees.size(); i++){
		//ignore removed (clustered) points and don't plot mirror points
		if(mt_trees[i] == nullptr) continue;
		if(mt_trees[i]->points->mean().at(1) > 2*acos(-1) || mt_trees[i]->points->mean().at(1) < 0){
			nmirror++; 
			continue; }
		if(_verb > 0){
			cout << "getting " << mt_trees[i]->points->GetNPoints() << " points in cluster #" << i << endl;
			mt_trees[i]->points->Print();
		}
		trees.push_back(std::move(mt_trees[i]));
	}
	if(_verb > 0) cout << nmirror << " mirror points." << endl;
	cout << trees.size() << " jets found." << endl;
}
*/



std::unique_ptr<GaussianMixture> BayesCluster::_subcluster(string oname){
	vector<double> ws;
	_jets[0].GetWeights(ws);
	double gev = ws[0]/_jets[0].E();
	//create GMM model
	std::unique_ptr<PointCollection> points = std::make_unique<PointCollection>();
	for(auto point : _points){
		points->AddPoints(point);
	}
	BayesPoint center({points->Centroid(0), points->CircularCentroid(1), points->Centroid(2)});
	//int npts = points->GetNPoints();
	//int mcls = 10;
	//int maxK = npts < mcls ? npts : mcls;
	int maxK = points->GetNPoints();


	std::unique_ptr<GaussianMixture> gmm = std::make_unique<GaussianMixture>(maxK);
	//double tresCte = 0.2 * 1e-9; //time resolution for CMS ECAL (s) (200 ps)
	//double tresStoch = 0.34641 * 1e-9; //rate of time res that gives 400 ps at E = 1 GeV (in [GeV*s])
	//tresStoch = 2.4999200e-05; //rate of time res that gives 5 ns at E = 5 GeV (in [GeV*s])
	//tresStoch *= gev;
	//needs to be before SetData bc thats when the measurement errors are set
	if(_verb > 6) cout << "BayesCluster subcluster - Using tresCte " << _tresCte << " _tresStoch " << _tresStoch << " _tresNoise " << _tresNoise << " gev " << gev << " _cell " << _cell << endl;
	gmm->SetMeasErrParams(_cell, _tresCte, _tresStoch, _tresNoise); 
	gmm->SetData(std::move(points));

	//set # ghosts by points with weight below a threshold - done automatically in InitParameters()
	
	//transform points into local coordinates
	//for GMM parameter estimation
	//use weighted mean as center to be set to 0 point
	//zero points by energy centroid
	//x' = x - a
	center.Put02pi(1);	
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
	//cout << "scaled and shifted - gmm data is" << endl; gmm->GetData()->Print();
	//cout << "scaled points" << endl;
	//gmm->GetData()->Print();
	//gmm->GetData()->at(0).Print();
	
	//cout << "1 - set gmm data as" << endl; gmm->GetData()->Print();
	//needs to be before ScaleData() bc this method also scales the smear
	if(!_smear.empty()){ gmm->SetDataSmear(_smear); }
	gmm->SetAlpha(_subalpha);
	gmm->SetVerbosity(_verb);

	gmm->InitParameters(_prior_params);

	
	//create EM algo
	std::unique_ptr<VarEMCluster> algo = std::make_unique<VarEMCluster>(gmm.get(),maxK);
	double totw = gmm->GetData()->Sumw();
	algo->SetThresh(0.01*totw);



	map<string, vector<Matrix>> params;
	bool viz = false;
	if(!oname.empty()){
		viz = true;
	}

	//inverse transformations
	Matrix RscaleInv;
	RscaleInv.invert(Rscale);
	center.Scale(-1);
	
	VarClusterViz3D cv3D(algo.get());
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
	double dLogL = 999;
	double newLogL = 0;
	double LogLthresh = 1e-3; //use relative error to old logLH
	//cout << "initial # clusters " << gmm->GetNClusters() << endl;
	double oldLogL = algo->EvalLogL();
	//cout << "after first logLH eval " << gmm->GetNClusters() << endl;
	algo->SetClusterStart();
	////////run EM algo////////
	int it = -1;
	while(dLogL > fabs(oldLogL)*LogLthresh || it == 0){
		it++;
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
		//newLogL > oldLogL
		dLogL = newLogL - oldLogL;
		if(_verb > 2)cout << std::setprecision(10) << "it " << it << " new logl " << newLogL << " oldlogl " << oldLogL << " dlogl " << dLogL << " # clusters " << gmm->GetNClusters() << endl;
		oldLogL = newLogL;
	}
		if(_verb > 2){
			cout << "Reached convergence at iteration " << it << endl;
		}
//cout << "finished in " << it << " iterations with final dLogL " << dLogL << " and final logL " << newLogL << endl;

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
	gmm->PutPhi02pi(); //does for data and parameters (not parameters you stupid slut) 
	gmm->PutPhi02pi_params(); //does for data and parameters (not parameters you stupid slut) 
	//- do after data + parameters shift so the [0,2pi] transformation doesn't get shifted
	//cout << "phi02pi - gmm data is" << endl; gmm->GetData()->Print();
	//cout << "center " << endl; center.Print();
	//cout << "predicted center - lead only - nclusters " << gmm->GetNClusters() << endl;
	if(_verb > 2){
		cout << std::setprecision(10) << endl;
		cout << "Estimated parameters" << endl;
		vector<double> norms;
		gmm->GetNorms(norms);
		for(int k = 0; k < gmm->GetNClusters(); k++){
			std::map<string, Matrix> params1, params2;
			gmm->GetLHPosteriorParameters(k, params1);
			gmm->GetDataStatistics(k, params2);
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
//need to add idx corresponding to plane index because those are the nodes stored in merge tree
void BayesCluster::_add_entry_to_maps(const int i, CompareMap& inmap, const Dnn2piCylinder* DNN, bool prob){
		bool verbose = false;
		double comp;
		int j, cyl_j;
		if(prob){
			comp = DNN->NearestNeighbourProb(i);
			j = DNN->NearestNeighbourProbIndex(i);
			if(_verb > 1) cout << std::setprecision(20) << "adding entry " << i << " with best probability " << comp << " pair " << j << std::setprecision(5) << " eta for this jet (i) " << _nodes[i]->points->mean().at(0) << " phi " << _nodes[i]->points->CircularMean(1) << " # rhs " << _nodes[i]->points->GetNPoints() << " eta for pair jet (j) " << _nodes[j]->points->mean().at(0) << " phi " << _nodes[j]->points->CircularMean(1) << " # rhs " << _nodes[j]->points->GetNPoints() << " and total sumw " << _nodes[j]->points->Sumw() + _nodes[i]->points->Sumw() << " sumw for i " << _nodes[i]->points->Sumw() << " sumw for j " << _nodes[j]->points->Sumw() << endl;
		}
		else{
			comp = DNN->NearestNeighbourDistance(i);
			j = DNN->NearestNeighbourIndex(i);
			if(_verb > 1) cout << std::setprecision(20) << "adding entry " << i << " with best distance " << comp << " pair " << j << std::setprecision(5) << endl;
		}
		inmap.insert(CompEntry(comp,verts(i,j)));
}



