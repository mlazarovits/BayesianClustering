#ifndef JET_HH
#define JET_HH

// The structure of this method is respectfully repurposed from ClusterSequence_Delaunay in FastJet (Cacciari, Salam, Soyez).
// This work was modified from its original form by Margaret Lazarovits on October 2, 2023. 
// The original version of this work was released
// under version 2 of the GNU General Public License. As of v3 of GNU GPL,
// any conditions added in Section 7 also apply. 

//----------------------------------------------------------------------
// Copyright (c) 2005-2021, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//----------------------------------------------------------------------
#include <math.h>
#include <cmath>
#include "Matrix.hh"
#include "Point.hh"
#include "JetPoint.hh"

//pt, eta, phi, m
//px, py, pz, E
//x, y, z, t
//point with more physics aspects - ecal cell
class Jet{
	public:
		Jet();
		Jet(double px, double py, double pz, double E);
		Jet(JetPoint rh);
		Jet(JetPoint rh, Point vtx);
		Jet(vector<JetPoint> rhs, Point vtx);
		Jet(vector<JetPoint> rhs);
		Jet(const Jet& j); //copy constructor
		virtual ~Jet();		

		bool operator==(Jet& j) const;
		bool operator!=(Jet& j) const;
		void operator =(const Jet& j){
				_px = j.px();
				_py = j.py();
				_pz = j.pz();
				_E = j.E();
				_mom = Point({_px, _py, _pz, _E});
				_eta = j.eta();
				_phi = j.phi();
				_t = j.time();
				
				_kt2 = j.kt2();
				_mass = j.mass();
				
				_child = j._child;
				_parent1 = j._parent1;
				_parent2 = j._parent2;
				_idx = j.GetUserIdx();
				_vtx = j.GetVertex();
				_nRHs = j.GetNConstituents();
				j.GetConstituents(_rhs);
		}

		//return four vector for clustering
		Point four_mom(){ return _mom; }

		void SetVertex(Point vtx){
			if(vtx.Dim() != 3){
				cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
				return;
			}
			_vtx = vtx;
		}
		//AK4 PF jets don't have x, y, z (only eta, phi)
		//return element i in four vector
		double px() const{ return _px; }
		double py() const{ return _py; }
		double pz() const{ return _pz; }
		double E() const{ return _E; }
		double e() const{ return _E; }
		double t() const{ return _t; }
		double time() const{ return _t; }

		//kinematic quantities
		//eta
		double eta() const{
			_ensure_valid_rap_phi();
			return _eta;
		}
		double rap() const{
			_ensure_valid_rap_phi();
			return _eta;
		}
		//phi
		double phi() const{
			_ensure_valid_rap_phi();
			return phi_02pi();
		}


		//transverse energy
		double Et() const{ return (_kt2==0) ? 0.0 : _E/sqrt(1.0+_pz*_pz/_kt2); } 
		//transverse mass
		double mt() const {return sqrt(std::abs(mperp2()));}
  		double mperp() const {return sqrt(std::abs(mperp2()));}
		//squared transverse mass = kt^2+m^2
  		double mperp2() const {return (_E+_pz)*(_E-_pz);}
		//invariant mass squared: m^2 = E^2 - p^2
		double m2() const{ return (_E+_pz)*(_E-_pz)-_kt2; }
		//invariant mass
		double mass() const{return sqrt(m2()); }

	
		//squared transverse momentum
  		double pt2() const {return _kt2;}
  		//the scalar transverse momentum
  		double pt() const {return sqrt(_kt2);} 
	 	//the squared transverse momentum
  		double kt2() const {return _kt2;} 	
		//deltaR between this and another jet pt
		double deltaR(Jet& jet) const{ return sqrt( (_eta - jet.eta())*(_eta - jet.eta()) + (_phi - jet.phi())*(_phi - jet.phi())); }
		
		void SetWeight(double w){ for(int i = 0; i < _nRHs; i++) _rhs[i].SetWeight(w); }
		void SetWeight(vector<double> w){ if(w.size() != _nRHs) return;
			for(int i = 0; i < _nRHs; i++) _rhs[i].SetWeight(w[i]); }
		void GetWeights(vector<double>& ws){
			ws.clear();
			for(int i = 0; i < _nRHs; i++) ws.push_back(_rhs[i].GetWeight()); 
		}

		vector<JetPoint> GetJetPoints() const{return _rhs;}
		
		PointCollection GetJetPoints_mom() const{
			PointCollection pc;
			//for(int i = 0; i < _nRHs; i++) pc += _rhs[i].four_mom();
			return pc;
		}

		//add subjets/pixels to jet
		void add(const Jet& jt);
		void add(const JetPoint& rh);
		
		//constituents (jet points) in jet (clustered or unclustered)
		void GetConstituents(vector<JetPoint>& rhs) const { rhs.clear(); rhs = _rhs; }
		void GetEnergies(vector<double>& energies) const{ energies.clear(); for(int j = 0; j < (int)_rhs.size(); j++) energies.push_back(_rhs[j].E()); }
		void GetEtaPhiConstituents(PointCollection& pc) const{
			pc.Clear();
			for(int i = 0; i < (int)_rhs.size(); i++){
				pc += Point({_rhs[i].eta(), _rhs[i].phi(), _rhs[i].time()});
			}
		}
		
		int GetNConstituents() const{return (int)_rhs.size(); }	
		
		//parents in cluster
		void GetParents(Jet& p1, Jet& p2) const;
		void SetParents(Jet* p1, Jet* p2){ _parent1 = p1; _parent2 = p2; p1->SetBaby(this); p2->SetBaby(this); }

		//children in cluster
		Jet GetBaby() const;

		void AddRecHit(JetPoint rh){
			_rhs.push_back(rh); 
			_nRHs = (int)_rhs.size();
			//recalculate time
			_set_time();
		}
		

		//subjets (jets) in jet (clustered or unclustered)
		void GetSubJets(vector<Jet>& subjets, int depth = 0) const;

		//this has jet?
		bool Has(Jet& jet) const;
		//this has rh?
		bool Has(JetPoint& rh) const;

		//set user idx info
		void SetUserIdx(int i){ _idx = i; }
		int GetUserIdx() const{ return _idx; }
		//history index info
		void SetHistoryIndex(int i){ _hist_idx = i; }
		int GetHistoryIndex() const{ return _hist_idx; }
		
		void GetClusterParams(Matrix& mu, Matrix& cov){ mu = _mu; cov = _cov; }
	
		//define jet time from cluster parameters
		double GetJetTime() const{ return _t; }
		void SetJetTime(double t){ _t = t; }
	
		Point GetVertex() const{return _vtx; }

		//check IR + collinear safety?

		//scale momentum
		void scaleMom(double s){
			_px *= s;
			_py *= s;
			_pz *= s;
			_E *= s;
			_kt2 *= s*s;
		 }
		//sets phi [0,2pi]
		double phi_02pi() const{
			if(_phi < 0) return _phi + 2*acos(-1);
			else if(_phi > 2*acos(-1)) return _phi - 2*acos(-1);
			else return _phi; 
		}

		//wraps phi around pi, [-pi,pi]
		double phi_negPiToPi() const{
 			double pi = acos(-1);
			double o2pi = 1./(2*pi);
			if(fabs(_phi) <= pi)
				return _phi;
			double n = std::round(_phi * o2pi);
			return _phi - n * double(2.* pi);

		}



	
	protected:
		void _ensure_valid_rap_phi() const{
			if(_phi == _invalid_phi) _set_rap_phi();
		}
		//as done in FastJet - https://fastjet.fr/repo/doxygen-3.4.2/PseudoJet_8cc_source.html
		void _set_rap_phi() const{
			if (_kt2 == 0.0) {
				_phi = 0.0; } 
			else {
				_phi = atan2(this->_py,this->_px);
			}
			if (_phi < 0.0) {_phi += twopi;}
			if (_phi >= twopi) {_phi -= twopi;} // can happen if phi=-|eps<1e-15|?
			if (this->_E == abs(this->_pz) && _kt2 == 0) {
				// Point has infinite rapidity -- convert that into a very large
				// number, but in such a way that different 0-pt momenta will have
				// different rapidities (so as to lift the degeneracy between
				// them) [this can be relevant at parton-level]
				double MaxRapHere = _maxRap + abs(this->_pz);
				if (this->_pz >= 0.0) {this->_eta = MaxRapHere;} else {this->_eta = -MaxRapHere;}
			} else {
				// get the rapidity in a way that's modestly insensitive to roundoff
				// error when things pz,E are large (actually the best we can do without
				// explicit knowledge of mass)
				double effective_m2 = std::max(0.0,m2()); // force non tachyonic mass
				double E_plus_pz    = _E + abs(_pz); // the safer of p+, p-
				// p+/p- = (p+ p-) / (p-)^2 = (kt^2+m^2)/(p-)^2
				_eta = 0.5*log((_kt2 + effective_m2)/(E_plus_pz*E_plus_pz));
				if (_pz > 0) {_eta = - _eta;}
			}
		}
		void SetBaby(Jet* child){ _child = child; }

		//set time to be energy weighted time
		void _set_time(){
			_t = 0;
			double norm = 0;
			for(int i = 0; i < _nRHs; i++){
				_t += _rhs[i].t()*_rhs[i].E();
				norm += _rhs[i].E();
			}
			_t /= norm;
		}

	private:
		//momentum four vector
		double _px;
		double _py;
		double _pz;
		mutable double _kt2;
		mutable double _eta;
		mutable double _phi;
		double _E;
		double _mass;
		double _t;
		Point _mom;
		//double _E;
		Point _vtx;

		//mutable double _eta, _phi, _theta;
		//rec hits (ids) in jet
		vector<JetPoint> _rhs;

		int _nRHs;

		//cluster params - modeling jet as gaussian in time, space + energy
		Matrix _mu;
		Matrix _cov;

		//parents + child info
		Jet* _parent1;
		Jet* _parent2;

		Jet* _child;
	
		//user index	
		int _idx;
		//history index
		int _hist_idx;
		

		double _maxRap = 1e5;
		//default values - have yet to be calculated or set
		double _invalid_phi = -100.0;
		double _invalid_eta = -1e200;
		double twopi = 6.28318530717;

};

#endif
