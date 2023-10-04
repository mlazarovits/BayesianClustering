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
	vector<EtaPhi> points(n);
	for (int i = 0; i < n; i++) {
		points[i] = EtaPhi(_jets[i].rap(),_jets[i].phi_02pi());
		points[i].sanitize(); // make sure things are in the right range
	}
	const bool verbose = false;
	const bool ignore_nearest_is_mirror = true; //based on _Rparam < twopi, should always be true for this 
	Dnn2piCylinder* DNN = new Dnn2piCylinder(points, ignore_nearest_is_mirror, verbose);

	//need to make a distance map like in FastJet, but instead of clustering
	//based on geometric distance, we are using merge probability (posterior) from BHC
	//all three dimensions will go into calculating the probabilityu
	//but the map will be built only in 2D space
	//structure is typdef'ed in header
	ProbMap RkMap;

	//fill map with initial potential clusterings
	double rk;
	int j;	
	for(int i = 0; i < n; i++){
		rk = DNN->NearestNeighbourProb(i);
		j = DNN->NearestNeighbourProbIndex(i);
		RkMap.insert(RkEntry(rk,verts(i,j)));
	}

	//run the clustering
		// find nearest vertices in map
		//do_ij_recomb
		//get eta phi (and time!) of new point - centroid of points in combined vertices
		//update DNN with RemoveCombinedAddCombination
			//this should also update the merge tree	


}


