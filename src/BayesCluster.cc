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
	int n = (int)_jets.size();
	vector<PointCollection> points(n);
  	double twopi = 6.28318530717;
	for (int i = 0; i < n; i++) {
		PointCollection pc = PointCollection();
		Point pt = Point(3);
		pt.SetValue(_jets[i].rap(), 0);
		pt.SetValue(_jets[i].phi_02pi(), 1);
		pt.SetValue(_jets[i].time(), 2);
  		//make sure phi is in the right range
		if (pt.at(1) <  0)     pt.SetValue(pt.at(1) + twopi, 1); 
  		if (pt.at(1) >= twopi) pt.SetValue(pt.at(1) - twopi, 1); 
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
	ProbMap RkMap;
/*
	//fill map with initial potential clusterings
	double rk;
	int j;	
	for(int i = 0; i < n; i++){
		rk = DNN->NearestNeighbourProb(i);
		j = DNN->NearestNeighbourProbIndex(i);
		RkMap.insert(RkEntry(rk,verts(i,j)));
	}

	//run the clustering
	for(int i = 0; i < n; i++){
		// find largest rk value in map (last entry)
		double BestRk;
		verts BestRkPair;
		std::multimap<double,verts>::iterator map_it;
		int jet_i, jet_j;
		bool Valid2;
		do{
			map_it = RkMap.end();
			map_it--;
			BestRk = map_it->first;
			BestRkPair = map_it->second;
			jet_i = BestRkPair.first;
			jet_j = BestRkPair.second;
			if (verbose) cout << "BayesCluster found recombination candidate: " << jet_i << " " << jet_j << " " << BestRk << endl; // GPS debugging
 			RkMap.erase(map_it);
			Valid2 = DNN->Valid(jet_j);
			if (verbose) cout << "BayesCluster validities i & j: " << DNN->Valid(jet_i) << " " << Valid2 << endl;
		} while(!DNN->Valid(jet_i) || !Valid2);
		//do_ij_recomb - this should be the same as in the OG code (except rk instead of dij)
		int nn;
		if (verbose) cout << "BayesCluster call _do_ij_recomb: " << jet_i << " " << jet_j << " " << BestRk << endl; // GPS debug
      		_do_ij_recombination_step(jet_i, jet_j, BestRk, nn);
		//get eta phi (and time!) of new point - centroid of points in combined vertices
		
		//update DNN with RemoveCombinedAddCombination
		//this should also update the merge tree	
	}
*/
}


