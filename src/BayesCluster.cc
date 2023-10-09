#include "BayesCluster.hh"
#include "DynamicNearestNeighbours.hh"

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
void BayesCluster::_cluster(){
	//the 2D Delauney triangulation used in FastJet will be used to seed the 3D clustering
	int n = 10;//(int)_jets.size();
	vector<PointCollection> points(n);
	for (int i = 0; i < n; i++) {
		PointCollection pc = PointCollection();
		Point pt = Point(3);
		pt.SetValue(_jets[i].rap(), 0);
		pt.SetValue(_jets[i].phi_02pi(), 1);
		pt.SetValue(_jets[i].time(), 2);
  		//make sure phi is in the right range
		//if (pt.at(1) <  0)     pt.SetValue(pt.at(1) + twopi, 1); 
  		//if (pt.at(1) >= twopi) pt.SetValue(pt.at(1) - twopi, 1); 
		sanitize(pt);
		pc.AddPoint(pt);
		points[i] = pc;
	}
	const bool verbose = true;
	const bool ignore_nearest_is_mirror = true; //based on _Rparam < twopi, should always be true for this 
	Dnn2piCylinder* DNN = new Dnn2piCylinder(points, ignore_nearest_is_mirror, verbose);

	//need to make a distance map like in FastJet, but instead of clustering
	//based on geometric distance, we are using merge probability (posterior) from BHC
	//all three dimensions will go into calculating the probabilityu
	//but the map will be built only in 2D space
	//structure is typdef'ed in header
	CompareMap ProbMap, DistMap;
	//fill map with initial potential clusterings
	for(int i = 0; i < n; i++){
		_add_entry_to_maps(i, ProbMap, DNN);
		_add_entry_to_maps(i, DistMap, DNN, false);	
	}
	std::pair<std::multimap<double,verts>::iterator, std::multimap<double,verts>::iterator> ret;	
	//run the clustering
	for(int i = 0; i < n; i++){
		// find largest rk value in map (last entry)
		double BestRk;
		verts BestRkPair;
		std::multimap<double,verts>::iterator map_it;
		int jet_i, jet_j;
		bool Valid2;
		double oldrk = 0;
		do{
			map_it = ProbMap.end();
			map_it--;
			BestRk = map_it->first;
			BestRkPair = map_it->second;
			//check for equal rks, break tie with 3d distance
			//right now - only equal rks are for equivalent merges
			ret = ProbMap.equal_range(BestRk);
				for(std::multimap<double,verts>::iterator it = ret.first; it != ret.second; it++){
					cout << "rk: " << it->first << " with points: " << it->second.first << ", " << it->second.second << endl;
			
				}
			if(oldrk == BestRk) cout << "previous was equal to this" << endl;
			oldrk = BestRk;	
			jet_i = BestRkPair.first;
			jet_j = BestRkPair.second;
			if (verbose) cout << "BayesCluster found recombination candidate: " << jet_i << " " << jet_j << " " << BestRk << endl; // GPS debugging
 			ProbMap.erase(map_it);
			Valid2 = DNN->Valid(jet_j);
			if (verbose) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << endl;
		} while(!DNN->Valid(jet_i) || !Valid2);
		//do_ij_recomb - this should be the same as in the OG code (except rk instead of dij)
		int nn;
		if (verbose) cout << "BayesCluster call _do_ij_recomb: " << jet_i << " " << jet_j << " " << BestRk << endl; // GPS debug
      		_do_ij_recombination_step(jet_i, jet_j, BestRk, nn);
		//add new jet to vector of _jets
		_jets.push_back(_jets[nn]);
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
		//update map
		vector<int>::iterator it = updated_neighbors.begin();
		for(; it != updated_neighbors.end(); ++it){
			int ii = *it;
			_add_entry_to_maps(ii, ProbMap, DNN);
			_add_entry_to_maps(ii, DistMap, DNN, false);	
		}
	}
}

void BayesCluster::_add_entry_to_maps(const int i, CompareMap& inmap, const Dnn2piCylinder* DNN, bool prob){
		double comp;
		int j;
		if(prob){
			comp = DNN->NearestNeighbourProb(i);
			j = DNN->NearestNeighbourProbIndex(i);
		}
		else{
			comp = DNN->NearestNeighbourDistance(i);
			j = DNN->NearestNeighbourIndex(i);
		}
		inmap.insert(CompEntry(comp,verts(i,j)));
}
