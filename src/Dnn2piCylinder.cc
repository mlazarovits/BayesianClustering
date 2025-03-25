// This work was modified from its original form by Margaret Lazarovits on October 2, 2023. 
// The original version of this work was released
// under version 2 of the GNU General Public License. As of v3 of GNU GPL,
// any conditions added in Section 7 also apply. 
//----------------------------------------------------------------------
// Copyright (c) 2005-2021, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER


//#ifndef DROP_CGAL // in case we do not have the code for CGAL
#include <set>
#include "Dnn2piCylinder.hh"
using namespace std;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// initialiser...
Dnn2piCylinder::Dnn2piCylinder(
	const vector<EtaPhi> & input_points, 
	const bool & ignore_nearest_is_mirror,
	const bool & verbose) {
  
  _verbose = verbose;
  _ignore_nearest_is_mirror = ignore_nearest_is_mirror;
  vector<EtaPhi> plane_points;
  vector<int>    plane_point_indices(input_points.size());
  //plane_points.reserve(2*input_points.size());

  for (unsigned int i=0; i < input_points.size(); i++) {
    _RegisterCylinderPoint(input_points[i], plane_points);
    plane_point_indices[i] = i;
  }
  
  if (_verbose) cout << "============== Preparing _DNN" << endl;
  _DNN = new DnnPlane(plane_points, verbose);


  vector<int> updated_point_indices; // we'll not use information from this
  _CreateNecessaryMirrorPoints(plane_point_indices,updated_point_indices);
}
//----------------------------------------------------------------------
/// initialiser for probability merge...
Dnn2piCylinder::Dnn2piCylinder(
	const std::vector<PointCollection>& input_points, 
	const bool & ignore_nearest_is_mirror, MergeTree* mt,
	const bool & verbose) {
 SetMergeTree(mt); 
  _verbose = verbose;
  _ignore_nearest_is_mirror = ignore_nearest_is_mirror;
  vector<PointCollection> plane_points;
  vector<int>    plane_point_indices(input_points.size());
  //plane_points.reserve(2*input_points.size());

  //save 2D points for triangulation
  for (int i=0; i < input_points.size(); i++) {
    _RegisterCylinderPoint(input_points[i], plane_points);
    plane_point_indices[i] = i;
  }
  if (_verbose) cout << "============== Preparing _DNN" << endl;
  _DNN = new DnnPlane(plane_points, mt, verbose);


  vector<int> updated_point_indices; // we'll not use information from this
  _CreateNecessaryMirrorPoints(plane_point_indices,updated_point_indices);

}


//----------------------------------------------------------------------
/// Actions here are similar to those in the
/// Dnn3piCylinder::_RegisterCylinderPoint case, however here we do
/// NOT create the mirror point -- instead we initialise the structure
/// as if there were no need for the mirror point.
///
/// ADDITIONALLY push the cylinder_point onto the vector plane_points.
void Dnn2piCylinder::_RegisterCylinderPoint (const EtaPhi & cylinder_point,
					     vector<EtaPhi> & plane_points) {
  double phi = cylinder_point.second;
  assert(phi >= 0.0 && phi < 2*pi);
  
  // do main point
  MirrorVertexInfo mvi;
  mvi.main_index = _cylinder_index_of_plane_vertex.size();
  _cylinder_index_of_plane_vertex.push_back(_mirror_info.size());
  plane_points.push_back(cylinder_point);
  mvi.mirror_index = INEXISTENT_VERTEX;
  
  // 
  _mirror_info.push_back(mvi);
}

void Dnn2piCylinder::_RegisterCylinderPoint (const PointCollection & cylinder_point,
					     vector<PointCollection> & plane_points) {
  double eta = cylinder_point.mean().at(0);
  double phi = cylinder_point.mean().at(1);
  if(!(phi >= 0.0 && phi < 2*pi)){ cout << "bad phi: " << phi << endl; cylinder_point.Print(); }
  assert(phi >= 0.0 && phi < 2*pi);

 //cout << "RegisterCylinderPoint - rrent cyl to plane vertices: " << _cylinder_index_of_plane_vertex.size() << " points registered so far" << " " << _mirror_info.size() << " mirror_info.size() so far" << endl;
if(_verbose){for(int i = 0; i < _cylinder_index_of_plane_vertex.size(); i++)
	cout << "plane index: " << i << " cyl index: " << _cylinder_index_of_plane_vertex[i] << endl;
 }
  // do main point
  MirrorVertexInfo mvi;
  mvi.main_index = _cylinder_index_of_plane_vertex.size();
  _cylinder_index_of_plane_vertex.push_back(_mirror_info.size());
  plane_points.push_back(cylinder_point);
  mvi.mirror_index = INEXISTENT_VERTEX;
  // 
  _mirror_info.push_back(mvi);
  //cout << "RegisterCylinderPoint end" << endl;
}


//----------------------------------------------------------------------
/// For each plane point specified in the vector plane_indices,
/// establish whether there is a need to create a mirror point
/// according to the following criteria:
///
/// . phi < pi
/// . mirror does not already exist
/// . phi < NearestNeighbourDistance 
///   (if this is not true then there is no way that its mirror point
///   could have a nearer neighbour).
///
/// If conditions all hold, then create the mirror point, insert it
/// into the _DNN structure, adjusting any nearest neighbours, and
/// return the list of plane points whose nearest neighbours have
/// changed (this will include the new neighbours that have just been
/// added)
void Dnn2piCylinder::_CreateNecessaryMirrorPoints(
			  const vector<int> & plane_indices,
			  vector<int> & updated_plane_points) {
if(_verbose) cout << "Dnn2pi - CreateNecesssaryMirrorPoints - start" << endl;



  vector<PointCollection> new_plane_points;
  //manually add node to _merge_tree here because the mirror point (cluster) needs to have
  //the same BHC values as it's originator
  for (size_t i = 0; i < plane_indices.size(); i++) {
    int ip = plane_indices[i]; // plane index
    PointCollection pts = _DNN->points(ip);
   // double phi = pts.mean().at(1);
    double phi = pts.CircularMean(1);


    // require absence of mirror - inexistent says a mirror point needs to be created (if necessary)
    int ic = _cylinder_index_of_plane_vertex[ip];
    if (_mirror_info[ic].mirror_index != INEXISTENT_VERTEX) {continue;}


    // check that we are sufficiently close to the border --
    // i.e. closer than nearest neighbour distance. But RECALL:
    // nearest neighbourDistance actually returns a squared distance
    // (this was really stupid on our part -- caused considerable loss
    // of time ... )
    double nndist = _DNN->NearestNeighbourDistance(ip);
	double offset = 0.1; 
    //if ((phi-offset)*(phi-offset) >= nndist && ((twopi+offset)-(phi))*((twopi+offset)-(phi)) >= nndist) {continue;}
	//do not create mirror points if the nearest neighbor is closer than the 0 point or 2pi point
	//only create mirror point if nearest neighbor is further than 0 and 2pi point
    if (phi*phi >= nndist || (twopi-phi)*(twopi-phi) >= nndist) {continue;}
    if(_verbose) cout << "creating mirror point for " << ip << " with nndist: " << nndist << " phiphi " << phi*phi << " (2pi-phi)^2 " << (twopi-phi)*(twopi-phi) << endl;
    //cout << "creating mirror point for " << ip << " with nndist: " << nndist << " and closest neighbor " << _DNN->NearestNeighbourIndex(ip) << " with phi distance to 0: " << phi << " and phi distance to offset " << phi-offset << endl;
//cout << "point to mirror" << endl;
//pts.Print();

    //copy node
    node* x = new node(*_merge_tree->Get(ip));
    PointCollection* newpts = new PointCollection(_remap_phi(pts));
//cout << "original points" << endl; pts.Print();
//    cout << "mirrored points" << endl;
//	newpts->Print();
    // now proceed to prepare the point for addition
    new_plane_points.push_back(*newpts);
    x->points = newpts;
    x->ismirror = true;
    //make sure this node knows that its mirror exists (and vice versa)
    x->mirror = _merge_tree->Get(ip);
    _merge_tree->Get(ip)->mirror = x;

    //updated mirror index with newly created mirror point
    _mirror_info[ic].mirror_index = _cylinder_index_of_plane_vertex.size();
    _cylinder_index_of_plane_vertex.push_back(ic);
    //add new mirror node to be clustered
    _merge_tree->Insert(x);
  }

  vector<int> indices_to_remove; // empty
  vector<int> indices_added;     // will be filled as result of call
  _DNN->RemoveAndAddPoints(indices_to_remove,new_plane_points,indices_added, 
			   updated_plane_points);

if(_verbose) cout << "Dnn2pi - CreateNecesssaryMirrorPoints - end" << endl;
}



//----------------------------------------------------------------------
/// insertion and removal of points
void Dnn2piCylinder::RemoveAndAddPoints(const vector<int> & indices_to_remove,
				const vector<EtaPhi> & points_to_add,
				vector<int> & indices_added,
				vector<int> & indices_of_updated_neighbours) {

  // translate from "cylinder" indices of points to remove to the
  // plane indices of points to remove, bearing in mind that sometimes
  // there are multple plane points to remove.
  vector<int> plane_indices_to_remove;
  for (unsigned int i=0; i < indices_to_remove.size(); i++) {
    MirrorVertexInfo * mvi;
    mvi = & _mirror_info[indices_to_remove[i]];
    plane_indices_to_remove.push_back(mvi->main_index);
    if (mvi->mirror_index != INEXISTENT_VERTEX) {
      plane_indices_to_remove.push_back(mvi->mirror_index);
    }
  }

  // given "cylinder" points to add get hold of the list of
  // plane-points to add.
  vector<EtaPhi> plane_points_to_add;
  indices_added.clear();
  for (unsigned int i=0; i < points_to_add.size(); i++) {
    indices_added.push_back(_mirror_info.size());
    _RegisterCylinderPoint(points_to_add[i], plane_points_to_add);
  }

  // now get the hard work done (note that we need to supply the
  // plane_indices_added vector but that we will not actually check
  // its contents in any way -- the indices_added that is actually
  // returned has been calculated above).
  vector<int> updated_plane_neighbours, plane_indices_added;
  _DNN->RemoveAndAddPoints(plane_indices_to_remove, plane_points_to_add,
			     plane_indices_added, updated_plane_neighbours);

  vector<int> extra_updated_plane_neighbours;
  _CreateNecessaryMirrorPoints(updated_plane_neighbours,
			       extra_updated_plane_neighbours);

  // extract, from the updated_plane_neighbours, and
  // extra_updated_plane_neighbours, the set of cylinder neighbours
  // that have changed
  set<int> index_set;
  unsigned int i;
  for (i=0; i < updated_plane_neighbours.size(); i++) {
    index_set.insert(
       _cylinder_index_of_plane_vertex[updated_plane_neighbours[i]]);}
  for (i=0; i < extra_updated_plane_neighbours.size(); i++) {
    index_set.insert(
       _cylinder_index_of_plane_vertex[extra_updated_plane_neighbours[i]]);}

  // decant the set into the vector that needs to be returned
  indices_of_updated_neighbours.clear();
  for (set<int>::iterator iter = index_set.begin(); 
       iter != index_set.end(); iter++) {
    indices_of_updated_neighbours.push_back(*iter);
  }
}
//----------------------------------------------------------------------
/// insertion and removal of points
void Dnn2piCylinder::RemoveAndAddPoints(const vector<int> & indices_to_remove,
				const vector<PointCollection>& points_to_add,
				vector<int> & indices_added,
				vector<int> & indices_of_updated_neighbours) {
  // translate from "cylinder" indices of points to remove to the
  // plane indices of points to remove, bearing in mind that sometimes
  // there are multple plane points to remove.
if(_verbose) cout << "Dnn2pi - RemoveAndAddPoints start" << endl;  

vector<int> plane_indices_to_remove;
  for (unsigned int i=0; i < indices_to_remove.size(); i++) {
    MirrorVertexInfo * mvi;
    mvi = & _mirror_info[indices_to_remove[i]];
if(_verbose) cout << "Dnn2pi to remove: " << indices_to_remove[i] << " main idx: " << mvi->main_index << " mirror idx: " << mvi->mirror_index << endl;
    plane_indices_to_remove.push_back(mvi->main_index);
    if (mvi->mirror_index != INEXISTENT_VERTEX) {
      plane_indices_to_remove.push_back(mvi->mirror_index);
    }
  }

if(_verbose){
for(int i = 0; i < plane_indices_to_remove.size(); i++) cout << "remove plane index: " << plane_indices_to_remove[i] << endl;
}

  // given "cylinder" points to add get hold of the list of
  // plane-points to add.
  //vector<EtaPhi> plane_points_to_add;
  vector<PointCollection> plane_points_to_add;
  indices_added.clear();
  
if(_verbose)
cout << "adding points" << endl;
  for (unsigned int i=0; i < points_to_add.size(); i++) {
  if(_verbose)   cout << "adding cyl index: "  << _mirror_info.size() << " for pt #" << i << endl;
	indices_added.push_back(_mirror_info.size());
    _RegisterCylinderPoint(points_to_add[i], plane_points_to_add);
  }

  // now get the hard work done (note that we need to supply the
  // plane_indices_added vector but that we will not actually check
  // its contents in any way -- the indices_added that is actually
  // returned has been calculated above).
  vector<int> updated_plane_neighbours, plane_indices_added;
  _DNN->RemoveAndAddPoints(plane_indices_to_remove, plane_points_to_add,
			     plane_indices_added, updated_plane_neighbours);
if(_verbose){ 
  cout << "added points: ";
  std::vector<int>::iterator it;
  for(it = plane_indices_added.begin(); it != plane_indices_added.end(); ++it)
    cout << *it << " mirror? " << _cylinder_index_of_plane_vertex[*it] << endl;
  }
  vector<int> extra_updated_plane_neighbours;
  //if this merge clusters all points (ie is the last clustering) do not create mirror points
  if(_merge_tree->GetNClusters() == 1) return;


  _CreateNecessaryMirrorPoints(updated_plane_neighbours,
			       extra_updated_plane_neighbours);




  // extract, from the updated_plane_neighbours, and
  // extra_updated_plane_neighbours, the set of cylinder neighbours
  // that have changed
  set<int> index_set;
  unsigned int i;
  for (i=0; i < updated_plane_neighbours.size(); i++) {
  if(_verbose)     cout << "updated main idx: " << _cylinder_index_of_plane_vertex[updated_plane_neighbours[i]] << " plane idx: " << updated_plane_neighbours[i]<< endl;
    index_set.insert(
       _cylinder_index_of_plane_vertex[updated_plane_neighbours[i]]);}
  for (i=0; i < extra_updated_plane_neighbours.size(); i++) {
   if(_verbose)    cout << "updated mirror idx: " << _cylinder_index_of_plane_vertex[extra_updated_plane_neighbours[i]] << " plane idx: " << extra_updated_plane_neighbours[i] << endl;
    index_set.insert(
       _cylinder_index_of_plane_vertex[extra_updated_plane_neighbours[i]]);}
if(_verbose) cout << "number of pts to update: " << index_set.size() << endl;
  // decant the set into the vector that needs to be returned
  indices_of_updated_neighbours.clear();
  for (set<int>::iterator iter = index_set.begin(); 
       iter != index_set.end(); iter++) {
	if(_verbose)cout << "updating idx: "<< *iter << endl;
    indices_of_updated_neighbours.push_back(*iter);
  }
if(_verbose){ 
for(int i = 0; i < _cylinder_index_of_plane_vertex.size(); i++)
	cout << "plane index: " << i << " cyl index: " << _cylinder_index_of_plane_vertex[i] << " valid? " << Valid(_cylinder_index_of_plane_vertex[i]) << endl;}
if(_verbose) cout << "Dnn2pi - RemoveAndAddPoints end" << endl;  
}

//FASTJET_END_NAMESPACE

//#endif //  DROP_CGAL 
