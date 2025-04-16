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


//#ifndef DROP_CGAL // in case we do not have the code for CGAL

#ifndef DNNPLANE_HH
#define DNNPLANE_HH

#include "Triangulation.hh"
#include "DynamicNearestNeighbours.hh"
#include <vector>
#include <iostream>
#include "BaseTree.hh"
using std::cout;
using std::endl;
using std::vector;
//FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


/// \if internal_doc
/// @ingroup internal
/// \class DnnPlane
/// class derived from DynamicNearestNeighbours that provides an
/// implementation for the Euclidean plane
///
/// This class that uses CGAL Delaunay triangulation for most of the
/// work (it allows for easy and efficient removal and addition of
/// points and circulation over a point's neighbours). The treatment
/// of coincident points is not supported by CGAL and is implemented
/// according to the method specified in
/// issue-tracker/2012-02-CGAL-coincident/METHOD
/// \endif
class DnnPlane : public DynamicNearestNeighbours {
 public:
  /// empty initaliser
  DnnPlane() {}

  /// Initialiser from a set of points on an Eta-Phi plane, where both
  /// eta and phi can have arbitrary ranges
  DnnPlane(const std::vector<EtaPhi> &, const bool & verbose = false );

  /// Initialiser from a set of points on an Eta-Phi plane, where both
  /// eta and phi can have arbitrary ranges, with time included for
  /// probabilistic merging
  /// includes alpha for BHC (a) and GMM (suba)
  DnnPlane(const std::vector<PointCollection>& pc, MergeTree* mt, const bool & verbose = false);

  /// Returns the index of neighbour jj of point labelled
  /// by ii (assumes ii is valid)
  int NearestNeighbourIndex(const int ii) const ;

  /// Returns the distance to neighbour jj of point labelled
  /// by index ii (assumes ii is valid)
  double NearestNeighbourDistance(const int ii) const ;
  
  /// Returns the index of neighbour jj of point labelled
  /// by ii (assumes ii is valid)
  int NearestNeighbourProbIndex(const int ii) const ;

  /// Returns the highest probability of merging point ii
  /// with any of its neighbors
  double NearestNeighbourProb(const int ii) const;

  node* NearestNeighbourProbNode(const int ii) const;

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
  
  void AddMirrorNodes(const std::vector<int>& current_indices_to_mirror,
                  const std::vector<node*> & nodes_to_add,
		  std::vector<int> & indices_added,
		  std::vector<int> & indices_of_updated_neighbours);


  /// returns the EtaPhi of point with index i.
  EtaPhi etaphi(const int i) const;
  /// returns the points of point with index i.
  PointCollection points(const int i) const;
  /// returns the eta point with index i.
  double eta(const int i) const;
  /// returns the phi point with index i.
  double phi(const int i) const;

private:
  /// Structure containing a vertex_handle and cached information on
  /// the nearest neighbour.
  struct SuperVertex {
    Vertex_handle vertex; // NULL indicates inexistence...
    double NNdistance;
    int NNindex;
    int coincidence;  // ==vertex->info.val() if no coincidence
                      // points to the coinciding SV in case of coincidence
    // later on for cylinder put a second vertex?
    double MaxRk = -1; //highest probability of merging, defaults to -1 (useful for mirror points)
    int MaxRkindex = -3; //index of vertex to merge with, defaults to inexistent vertex (defined as -3 in Dnn2piCylinder)
    node* n; //node containing point collection, rk value (should match MaxRk), parents, etc.
    node* bestmerge; //node containing bestmerge pair
  };

  //map vertex (via vertex_handle) to 3D points at vertex
  //TODO: update in _setnearestandupdate


  std::vector<SuperVertex> _supervertex;
  //set<Vertex_handle> _vertex_set;
  bool _verbose;

  //static const bool _crash_on_coincidence = true;
  static const bool _crash_on_coincidence = false;

  Triangulation _TR; /// CGAL object for dealing with triangulations

  /// calculates and returns the euclidean distance between points p1
  /// and p2
  inline double _euclid_distance(const CPoint& p1, const CPoint& p2) const {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
//cout << "distx " << distx << " disty " << disty << " for points p1: (" << p1.x() << ", " << p1.y() << ") and p2 (" << p2.x() << ", " << p2.y() << ")" << endl; 
    return distx*distx+disty*disty;
  }
 

  //---------------------------------------------------------------------- 
  /// Determines the index and distance of the nearest neighbour to 
  /// point j and puts the information into the _supervertex entry for j
  void _SetNearest(const int j);

  //----------------------------------------------------------------------
  /// Determines and stores the nearest neighbour of j.
  ///
  /// For each voronoi neighbour D of j if the distance between j and D
  /// is less than D's own nearest neighbour, then update the
  /// nearest-neighbour info in D; push D's index onto 
  /// indices_of_updated_neighbours
  ///
  /// Note that j is NOT pushed onto indices_of_updated_neighbours --
  /// if you want it there, put it there yourself.
  void _SetAndUpdateNearest(const int j, 
			    std::vector<int> & indices_of_updated_neighbours);

  /// given a vertex_handle returned by CGAL on insertion of a new
  /// points, returns the coinciding vertex's value if it turns out
  /// that it corresponds to a vertex that we already knew about
  /// (usually because two points coincide)
  int _CheckIfVertexPresent(const Vertex_handle & vertex, 
			    const int its_index);

  //----------------------------------------------------------------------
  /// if the distance between 'pref' and 'candidate' is smaller (or
  /// equal) than the one between 'pref' and 'near', return true and
  /// set 'mindist' to that distance. Note that it is assumed that
  /// 'mindist' is the euclidian distance between 'pref' and 'near'
  ///
  /// Note that the 'near' point is passed through its vertex rather
  /// than as a point. This allows us to handle cases where we have no min
  /// yet (near is the infinite vertex)
  inline bool _is_closer_to(const CPoint &pref, 
			    const CPoint &candidate,
			    const Vertex_handle &near,
			    double & dist,
			    double & mindist){
    dist = _euclid_distance(pref, candidate);
    return _is_closer_to_with_hint(pref, candidate, near, dist, mindist);
  }

  /// same as '_is_closer_to' except that 'dist' already contains the
  /// distance between 'pref' and 'candidate'
  inline bool _is_closer_to_with_hint(const CPoint &pref, 
				      const CPoint &candidate,
				      const Vertex_handle &near,
				      const double & dist,
				      double & mindist){
    
    // check if 'dist', the pre-computed distance between 'candidate'
    // and 'pref' is smaller than the distance between 'pref' and its
    // currently registered nearest neighbour 'near' (and update
    // things if it is)
    //
    // Interestingly enough, it has to be pointed out that the use of
    // 'abs' instead of 'std::abs' returns wrong results (apparently
    // ints without any compiler warning)
    //
    // The (near != NULL) test is there for one single reason: when
    // checking that a newly inserted point is not closer than a
    // previous NN, if that distance comparison involves a "nearly
    // degenerate" distance we need to access near->point. But
    // sometimes, in the course of RemoveAndAddPoints, its previous NN
    // has been deleted and its vertex (corresponding to 'near') set
    // to NULL. This is not a problem as all points having a deleted
    // point as NN will have their NN explicitly recomputed at the end
    // of RemoveAndAddPoints so here we should just make sure there is
    // no crash... that's done by checking (near != NULL)
    if ((std::abs(dist-mindist)<DISTANCE_FOR_CGAL_CHECKS) &&
	_is_not_null(near) && // GPS cf. note below about != NULL with clang/CGAL
	(_euclid_distance(candidate, near->point())<DISTANCE_FOR_CGAL_CHECKS)){
      // we're in a situation where there might be a rounding issue,
      // use CGAL's distance computation to get it right
      //
      // Note that in the test right above,
      // (abs(dist-mindist)<1e-12) guarantees that the current
      // nearest point is not the infinite vertex and thus
      // nearest->point() is not ill-defined
      if (_verbose) std::cout << "using CGAL's distance ordering" << std::endl;
      if (CGAL::compare_distance_to_point(pref, candidate, near->point())!=CGAL::LARGER){
	mindist = dist;
	return true;
      }
    } else if (dist <= mindist) {
      // Note that the use of a <= in the above expression (instead of
      // a strict ordering <) is important in one case: when checking
      // if a new point is the new NN of one of the points in its
      // neighbourhood, in case of distances being ==, we are sure
      // that 'candidate' is in a cell adjacent to 'pref' while it may
      // no longer be the case for 'near'
      mindist = dist;
      return true;
    } 
    
    return false;
  }


  inline node* _merge_prob(const SuperVertex& v1, const SuperVertex& v2) const{
	node* n1 = v1.n;
	node* n2 = v2.n;
	node* x = _merge_tree->CalculateMerge(n1, n2);
	//double rk = x->val; 
cout << "n1 has " << n1->points->GetNPoints() << " pts" << endl; n1->points->Print();
cout << "n2 has " << n2->points->GetNPoints() << " pts" << endl; n2->points->Print();
cout << "x has " << x->points->GetNPoints() << " pts" << endl; x->points->Print();
	return x;
  }
  /// calculates merge probabilities for neighbor (candidate) of point (pref)
  /// compares to best current merge (pref and best)
  inline bool _best_merge_prob(const SuperVertex& pref,
			       const SuperVertex& candidate,
			       const Vertex_handle &best){
			       //node& bestmerge){
			       //double& rk,
			       //double& maxrk){
//compare x + bestmerge
    node* x = _merge_prob(pref, candidate);
    node* bestnode = _supervertex[best->info().val()].n;
    if (x->log_h1_prior+bestnode->log_didj > x->log_didj+bestnode->log_h1_prior){
	return true;
    }
    else return false;

//    bool ret = _best_merge_prob_with_hint(pref, candidate, best, *x, bestmerge);
//	return ret;
  }
  
  inline bool _best_merge_prob_with_hint(const SuperVertex &pref,
			       const SuperVertex& candidate,
			       const Vertex_handle &best,
			       node& x,
      			       node& bestmerge){
			       //const double& rk,
      			       //double& maxrk){
      //if bestmerge hasn't been set yet, return true
      if(bestmerge.val == -999){
	cout << "bestmerge hasnt been set" << endl;
	bestmerge = node(x);
       cout << "set bestmerge to node with pts " << endl; bestmerge.points->Print();
      cout << "bestmerge log_h1_prior " << bestmerge.log_h1_prior<< " log_didj " <<  bestmerge.log_didj << endl;
	 return true;
      }
      //strictly greater than
      //if (rk > maxrk){
      if (x.log_h1_prior+bestmerge.log_didj > x.log_didj+bestmerge.log_h1_prior){
      cout << "do comparison for bestmerge" << endl;
	cout << "x log_h1_prior " << x.log_h1_prior<< " log_didj " <<  x.log_didj << endl;
      cout << "bestmerge log_h1_prior " << bestmerge.log_h1_prior<< " log_didj " <<  bestmerge.log_didj << endl;
      cout << "comparing " << x.log_h1_prior+bestmerge.log_didj << " to " <<  x.log_didj+bestmerge.log_h1_prior << endl;
      cout << "(rearr) comparing " << x.log_h1_prior - x.log_didj << " to " <<  -bestmerge.log_didj+bestmerge.log_h1_prior << endl;
	cout << "current bestmerge pts" << endl; bestmerge.points->Print();
	//maxrk = rk;
	bestmerge = node(x);
cout << "best merge for node updated " << endl;
	cout << "new bestmerge pts" << endl; bestmerge.points->Print();
	return true;
      }
    //don't update maxrk
    return false;
  }

  /// if a distance between a point and 2 others is smaller than this
  /// and the distance between the two points is also smaller than this
  /// then use CGAL to compare the distances. 
  static const double DISTANCE_FOR_CGAL_CHECKS;  
  static const int max_ndigits; //maximum number of digits for rounding distances

  /// small routine to check if if a vertex handle is not null.
  ///
  /// This is centralised here because of the need for a not-very
  /// readable workaround for certain CGAL/clang++ compiler
  /// combinations.
  inline static bool _is_not_null(const Vertex_handle & vh) {
    // as found in issue 2015-07-clang61+CGAL, the form
    //   vh != NULL
    // appears to be broken with CGAL-3.6 and clang-3.6.0/.1
    //
    // The not very "readable"
    //   vh.operator->()
    // does the job instead
    return vh.operator->();
  }
  
};


// here follow some inline implementations of the simpler of the
// functions defined above
// returns info for vertex ii and its neighbor jj
inline int DnnPlane::NearestNeighbourIndex(const int ii) const {
  return _supervertex[ii].NNindex;}

inline double DnnPlane::NearestNeighbourDistance(const int ii) const {
  return _supervertex[ii].NNdistance;}

inline int DnnPlane::NearestNeighbourProbIndex(const int ii) const{
  return _supervertex[ii].MaxRkindex;}

inline double DnnPlane::NearestNeighbourProb(const int ii) const{
//cout << "node " << ii << " has max rk " << _supervertex[ii].MaxRk << " and points " << endl; _supervertex[ii].n->model->GetData()->Print();
  return _supervertex[ii].MaxRk;}

inline node* DnnPlane::NearestNeighbourProbNode(const int ii) const{
cout << "nearest neighbourprobnode for index " << ii << " is " << endl; _supervertex[ii].n->points->Print();
	return _supervertex[ii].n;}

inline bool DnnPlane::Valid(const int index) const {
  if (index >= 0 && index < static_cast<int>(_supervertex.size())) {
    return _is_not_null(_supervertex[index].vertex);
  }
  else {
    return false;
  }
}


inline EtaPhi DnnPlane::etaphi(const int i) const {
  CPoint * p = & (_supervertex[i].vertex->point());
  return EtaPhi(p->x(),p->y()); }

inline PointCollection DnnPlane::points(const int i) const {
  PointCollection  p =  *(_supervertex[i].n->points);
  return p; }

inline double DnnPlane::eta(const int i) const {
  return _supervertex[i].vertex->point().x(); }

inline double DnnPlane::phi(const int i) const {
  return _supervertex[i].vertex->point().y(); }


//FASTJET_END_NAMESPACE

#endif //  __FASTJET_DNNPLANE_HH__

//#endif // DROP_CGAL
