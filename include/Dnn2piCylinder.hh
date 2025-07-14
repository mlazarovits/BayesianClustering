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
#ifndef FASTJET_DNN2PICYLINDER_HH
#define FASTJET_DNN2PICYLINDER_HH

#include "DynamicNearestNeighbours.hh"
#include "DnnPlane.hh"
//#include "fastjet/internal/numconsts.hh"

#include <iostream>
using std::cerr;

//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


/// \if internal_doc
/// @ingroup internal
/// \class Dnn2piCylinder
/// class derived from DynamicNearestNeighbours that provides an
/// implementation for the surface of cylinder (using one 
/// DnnPlane object spanning 0--2pi).
/// \endif
class Dnn2piCylinder : public DynamicNearestNeighbours {
 public:
  /// empty initaliser
  Dnn2piCylinder() {}

  /// Initialiser from a set of points on an Eta-Phi plane, where
  /// eta can have an arbitrary ranges and phi must be in range
  /// 0 <= phi < 2pi;
  /// 
  /// NB: this class is more efficient than the plain Dnn4piCylinder
  /// class, but can give wrong answers when the nearest neighbour is
  /// further away than 2pi (in this case a point's nearest neighbour
  /// becomes itself, because it is considered to be a distance 2pi
  /// away). For the kt-algorithm (e.g.) this is actually not a
  /// problem (the distance need only be accurate when it is less than
  /// R, assuming R<2pi [not necessarily always the case as of
  /// 2010-11-19, when we've removed the requirement R<pi/2 in the
  /// JetDefinition constructor]), so we can tell the routine to
  /// ignore this problem -- alternatively the routine will crash if
  /// it detects it occurring (only when finding the nearest neighbour
  /// index, not its distance).
  Dnn2piCylinder(const std::vector<EtaPhi> &,
		 const bool & ignore_nearest_is_mirror = false,
		 const bool & verbose = false );

  /// initialiser for probability merge...
  /// includes alpha for BHC (a) and GMM (suba)
  Dnn2piCylinder(
	const std::vector<PointCollection>& input_points, 
	const bool & ignore_nearest_is_mirror, MergeTree* mt,
	const bool & verbose);
  

  /// Returns the index of  the nearest neighbour of point labelled
  /// by ii (assumes ii is valid)
  int NearestNeighbourIndex(const int ii) const ;

  /// Returns the distance to the nearest neighbour of point labelled
  /// by index ii (assumes ii is valid)
  double NearestNeighbourDistance(const int ii) const ;

  /// Returns the index of  the nearest neighbour of point labelled
  /// by ii (assumes ii is valid)
  int NearestNeighbourProbIndex(const int ii, int& merge_index) const ;

  /// Returns the index of  the nearest neighbour of point labelled
  /// by ii (assumes ii is valid)
  node* NearestNeighbourProbNode(const int ii) const ;
  

  /// Returns the highest probability of merging point ii
  /// with any of its neighbors
  double NearestNeighbourProb(const int ii) const ;

  /// Returns true iff the given index corresponds to a point that
  /// exists in the DNN structure (meaning that it has been added, and
  /// not removed in the meantime)
  bool Valid(const int index) const;

  void RemoveAndAddPoints(const std::vector<int> & indices_to_remove,
			  const std::vector<EtaPhi> & points_to_add,
			  std::vector<int> & indices_added,
			  std::vector<int> & indices_of_updated_neighbours);

  void RemoveAndAddPoints(const std::vector<int> & indices_to_remove,
			  const std::vector<PointCollection> & points_to_add,
			  std::vector<int> & indices_added,
			  std::vector<int> & indices_of_updated_neighbours);
  ~Dnn2piCylinder();



 private:
  double pi    = 3.14159265358;
  double twopi = 6.28318530717;

  // our extras to help us navigate, find distance, etc.
  const static int INEXISTENT_VERTEX=-3;

  bool _verbose;

  bool _ignore_nearest_is_mirror;

  /// Picture of how things will work... Copy 0--pi part of the 0--2pi
  /// cylinder into a region 2pi--2pi+ a bit of a Euclidean plane. Below we
  /// show points labelled by + that have a mirror image in this
  /// manner, while points labelled by * do not have a mirror image.
  ///      	 			  
  ///      |           .     |		  
  ///      |	       .     |		  
  ///      |           .     |		  
  ///      |           .     |		  
  ///      |        2  .     |		  
  ///      |        *  .     |		  
  ///      | +         . +   |		  
  ///      | 0         . 1   |
  ///      |	       .     |
  ///      0          2pi   2pi + a bit
  ///   	     
  /// Each "true" point has its true "cylinder" index (the index that
  /// is known externally to this class) as well as euclidean plane
  /// indices (main_index and mirror index in the MirrorVertexInfo
  /// structure), which are private concepts of this class.
  /// 
  /// In above picture our structures would hold the following info
  /// (the picture shows the euclidean-plane numbering)
  ///
  /// _mirror_info[cylinder_index = 0] = (0, 1)
  /// _mirror_info[cylinder_index = 1] = (2, INEXISTENT_VERTEX)
  ///
  /// We also need to be able to go from the euclidean plane indices
  /// back to the "true" cylinder index, and for this purpose we use
  /// the std::vector _cylinder_index_of_plane_vertex[...], which in the above example has
  /// the following contents
  ///
  /// _cylinder_index_of_plane_vertex[0] = 0
  /// _cylinder_index_of_plane_vertex[1] = 0
  /// _cylinder_index_of_plane_vertex[2] = 1
  ///

  /// 
  struct MirrorVertexInfo {
    /// index of the given point (appearing in the range 0--2pi) in the 
    /// 0--2pi euclidean plane structure (position will coincide with
    /// that on the 0--2pi cylinder, but index labelling it will be
    /// different)
    int main_index; 
    /// index of the mirror point (appearing in the range 2pi--3pi) in the
    /// 0--3pi euclidean plane structure
    int mirror_index; 
  };

  // for each "true" vertex we have reference to indices in the euclidean
  // plane structure
  std::vector<MirrorVertexInfo> _mirror_info;
  // for each index in the euclidean 0--2pi plane structure we want to
  // be able to get back to the "true" vertex index on the overall
  // 0--2pi cylinder structure
  std::vector<int> _cylinder_index_of_plane_vertex;

  // NB: we define POINTERS here because the initialisation gave
  //     us problems (things crashed!), perhaps because in practice
  //     we were making a copy without being careful and defining
  //     a proper copy constructor.
  DnnPlane * _DNN;

  /// given a phi value in the 0--pi range return one 
  /// in the 2pi--3pi range; whereas if it is in the pi-2pi range then
  /// remap it to be inthe range (-pi)--0.
  inline EtaPhi _remap_phi(const EtaPhi & point) {
    double phi = point.second;
    if (phi < pi) { phi += twopi ;} else {phi -= twopi;}
    return EtaPhi(point.first, phi);}

  inline PointCollection _remap_phi(const PointCollection& points) {
    double phi = points.CircularMean(1);
    double shift;
    if (phi < pi) { shift = twopi ;} else {shift = -twopi;}
//cout << "phi mean for shift " << phi << " circular mean " << points.CircularMean(1) << " shift " << shift << endl;
    BayesPoint transl = BayesPoint({0., -shift, 0.});
    PointCollection pc = PointCollection(points);
    pc.Translate(transl); 
    return pc;}

  //----------------------------------------------------------------------
  /// Actions here are similar to those in the
  /// Dnn3piCylinder::_RegisterCylinderPoint case, however here we do
  /// NOT create the mirror point -- instead we initialise the structure
  /// as if there were no need for the mirror point.
  ///
  /// ADDITIONALLY push the cylinder_point onto the vector plane_points.
  void _RegisterCylinderPoint (const EtaPhi & cylinder_point,
			       std::vector<EtaPhi> & plane_points);

  void _RegisterCylinderPoint (const PointCollection& cylinder_points,
			       std::vector<PointCollection> & plane_points);
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
 // /// return the list of plane points whose nearest neighbours have
  /// changed (this will include the new neighbours that have just been
  /// added)
  void _CreateNecessaryMirrorPoints(
			  const std::vector<int> & plane_indices,
			  std::vector<int> & updated_plane_points);

};


// here follow some inline implementations of the simpler of the
// functions defined above

//----------------------------------------------------------------------
/// Note: one of the difficulties of the 0--2pi mapping is that
/// a point may have its mirror copy as its own nearest neighbour
/// (if no other point is within a distance of 2pi). This does
/// not matter for the kt_algorithm with
/// reasonable values of radius, but might matter for a general use
/// of this algorithm -- depending on whether or not the user has
/// initialised the class with instructions to ignore this problem the
/// program will detect and ignore it, or crash.
inline int Dnn2piCylinder::NearestNeighbourIndex(const int current) const {
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
  int plane_index;
  if (mirror_index == INEXISTENT_VERTEX ) {
    plane_index = _DNN->NearestNeighbourIndex(main_index);
  } else {
    plane_index = (
	_DNN->NearestNeighbourDistance(main_index) < 
	_DNN->NearestNeighbourDistance(mirror_index)) ? 
      _DNN->NearestNeighbourIndex(main_index) : 
      _DNN->NearestNeighbourIndex(mirror_index) ; 
  }
  int this_cylinder_index = _cylinder_index_of_plane_vertex[plane_index];
  // either the user has acknowledged the fact that they may get the
  // mirror copy as the closest point, or crash if it should occur
  // that mirror copy is the closest point.
  assert(_ignore_nearest_is_mirror || this_cylinder_index != current);
  //if (this_cylinder_index == current) {
  //  cerr << "WARNING point "<<current<<
  //    " has its mirror copy as its own nearest neighbour"<<endl;
  //}
  return this_cylinder_index;
}

inline double Dnn2piCylinder::NearestNeighbourDistance(const int current) const {
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
  if (mirror_index == INEXISTENT_VERTEX ) {
    return _DNN->NearestNeighbourDistance(main_index);
  } else {
    return (
	_DNN->NearestNeighbourDistance(main_index) < 
	_DNN->NearestNeighbourDistance(mirror_index)) ? 
      _DNN->NearestNeighbourDistance(main_index) : 
      _DNN->NearestNeighbourDistance(mirror_index) ; 
  }
 
}

inline int Dnn2piCylinder::NearestNeighbourProbIndex(const int current, int& merge_index) const {
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
 //cout << "getting nearest neighbor prob index for current: " << current << " main idx: " << main_index << " mirror_index: " << mirror_index << endl;
  int plane_index;
  if (mirror_index == INEXISTENT_VERTEX ) {
    plane_index = _DNN->NearestNeighbourProbIndex(main_index);
merge_index = main_index;
  } else {
    plane_index = (
	_DNN->NearestNeighbourProb(main_index) > 
	_DNN->NearestNeighbourProb(mirror_index)) ? 
      _DNN->NearestNeighbourProbIndex(main_index) : 
      _DNN->NearestNeighbourProbIndex(mirror_index) ; 
    merge_index = (
	_DNN->NearestNeighbourProb(main_index) > 
	_DNN->NearestNeighbourProb(mirror_index)) ? 
      (main_index) : 
      (mirror_index) ; 
  }

  //infinite vertex - defined in Triangulation
  //no match
  if(plane_index == -1){
	return current;
  }
if(_verbose) cout << "IDX prob main idx: " << _DNN->NearestNeighbourProbIndex(main_index) << " prob mirror idx: " << _DNN->NearestNeighbourProbIndex(mirror_index) << endl; 
if(_verbose)cout << "current phi: ";
if(_verbose){ if( _merge_tree->Get(current) != nullptr) cout << _merge_tree->Get(current)->points->CircularMean(1);
else cout << " current null pts";} 
if(_verbose)cout << " main phi: ";
if(_verbose) if( _merge_tree->Get(main_index)->points != nullptr) cout << _merge_tree->Get(main_index)->points->CircularMean(1);

  int this_cylinder_index = _cylinder_index_of_plane_vertex[plane_index];
  // either the user has acknowledged the fact that they may get the
  // mirror copy as the closest point, or crash if it should occur
  // that mirror copy is the closest point.
  assert(_ignore_nearest_is_mirror || this_cylinder_index != current);
  if (this_cylinder_index == current) {
     cerr << "WARNING point "<<current<<
      " has its mirror copy as its own nearest neighbour"<<endl;
  }
  return this_cylinder_index;
}


//does the same as above but calculates the probability of merging in 3D space (eta, phi, time)
inline double Dnn2piCylinder::NearestNeighbourProb(const int current) const{
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
 if (mirror_index == INEXISTENT_VERTEX) {
    return _DNN->NearestNeighbourProb(main_index);
  } else {
    return (
	_DNN->NearestNeighbourProb(main_index) > 
	_DNN->NearestNeighbourProb(mirror_index)) ? 
      _DNN->NearestNeighbourProb(main_index) : 
      _DNN->NearestNeighbourProb(mirror_index) ; 
  }
}

inline node* Dnn2piCylinder::NearestNeighbourProbNode(const int current) const {
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
if(_verbose) cout << "current: " << current << " main idx: " << main_index << " cyl index: " << _cylinder_index_of_plane_vertex[main_index] <<  " mirror_index: " << mirror_index << " # clusters in merge tree: " << _merge_tree->GetNClusters() << endl;
if(_verbose)cout << "validity: main index - " << _DNN->Valid(main_index) << " mirror index: " << _DNN->Valid(mirror_index) << " cyl index: " << Valid(_cylinder_index_of_plane_vertex[main_index]) << endl;
  int plane_index;
  node* x;
  if (mirror_index == INEXISTENT_VERTEX ) {
    //plane_index = _DNN->NearestNeighbourProbIndex(main_index);
    x = _DNN->NearestNeighbourProbNode(main_index);
//cout << "nearestneighbourprobnode for index " << main_index << " has points " << endl; x->points->Print();
//cout << "node " << main_index << " has nearestneighbourprobenode index " << _DNN->NearestNeighbourProbIndex(main_index) << endl; 
  } else {
    //plane_index = (
    //    _DNN->NearestNeighbourProb(main_index) > 
    //    _DNN->NearestNeighbourProb(mirror_index)) ? 
    //  _DNN->NearestNeighbourProbIndex(main_index) : 
    //  _DNN->NearestNeighbourProbIndex(mirror_index) ; 
    x = (
	_DNN->NearestNeighbourProb(main_index) > 
	_DNN->NearestNeighbourProb(mirror_index)) ? 
        _DNN->NearestNeighbourProbNode(main_index) :
        _DNN->NearestNeighbourProbNode(mirror_index);
  }
	//if(_verbose) cout << "main idx prob: " << _DNN->NearestNeighbourProb(main_index) << " mirror idx prob: " << _DNN->NearestNeighbourProb(mirror_index) << endl;
 //if mirrored points, may need to translate back or set points to points of main idx
 if(x == nullptr) cout << "x null" << endl;
if(_verbose) cout << "phi of selected merge node: " << x->points->CircularMean(1) << endl;
if(x->points->CircularMean(1) > twopi) x->points = _merge_tree->Get(main_index)->points;
  return x;
/*
  //infinite vertex - defined in Triangulation
  //no match
  if(plane_index == -1){
	return current;
  }
cout << "prob main idx: " << _DNN->NearestNeighbourProb(main_index) << " prob mirror idx: " << _DNN->NearestNeighbourProb(mirror_index) << endl; 
cout << "IDX prob main idx: " << _DNN->NearestNeighbourProbIndex(main_index) << " prob mirror idx: " << _DNN->NearestNeighbourProbIndex(mirror_index) << endl; 
cout << "# cyl indices: " << _cylinder_index_of_plane_vertex.size() << endl;
for(int i = 0; i < _cylinder_index_of_plane_vertex.size(); i++)
	cout << "plane index: " << i << " cyl index: " << _cylinder_index_of_plane_vertex[i] << endl;


cout << "plane index: " << plane_index << endl;
  int this_cylinder_index = _cylinder_index_of_plane_vertex[plane_index];
cout << "this cyl index: " << this_cylinder_index << endl;
  // either the user has acknowledged the fact that they may get the
  // mirror copy as the closest point, or crash if it should occur
  // that mirror copy is the closest point.
  assert(_ignore_nearest_is_mirror || this_cylinder_index != current);
  if (this_cylinder_index == current) {
    cerr << "WARNING point "<<current<<
      " has its mirror copy as its own nearest neighbour"<<endl;
  }
  return this_cylinder_index;
*/

}

inline bool Dnn2piCylinder::Valid(const int index) const {
if(_verbose) cout << "Valid - index: " << index << " main index: " << _mirror_info[index].main_index << endl;  
  return (_DNN->Valid(_mirror_info[index].main_index));
}


inline Dnn2piCylinder::~Dnn2piCylinder() {
  delete _DNN; 
}

inline void PointCollection_to_EtaPhi(PointCollection& pc, vector<EtaPhi>& eps){
  for (unsigned int i=0; i < pc.GetNPoints(); i++) {
    EtaPhi ep = EtaPhi(pc.at(i).at(0), pc.at(i).at(1));
    eps.push_back(ep);
  }
  }

//FASTJET_END_NAMESPACE

#endif //  __FASTJET_DNN2PICYLINDER_HH__
//#endif //DROP_CGAL 
