#ifndef JET_HH
#define JET_HH
//
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
#include "BasePDFMixture.hh"
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
		Jet(JetPoint rh, BayesPoint vtx);
		Jet(const vector<JetPoint>& rhs, BayesPoint vtx);
		Jet(const vector<Jet>& jets);
		Jet(const Matrix& mu, const Matrix& cov, double E, double _pi = 1, BayesPoint vtx = BayesPoint({0., 0., 0.})); //constructor from subcluster information
		Jet(BasePDF* pdf, double E, double pi = 1, BayesPoint vtx = BayesPoint({0., 0., 0.})); //constructor from subcluster information
		Jet(BasePDFMixture* model, BayesPoint vtx, double gev, double detR); //need detector radius to convert eta, phi to x, y, z
		Jet(const Jet& j); //copy constructor
		virtual ~Jet();		

		bool operator==(const Jet& j) const;
		bool operator!=(const Jet& j) const;
		void operator =(const Jet& j){
				_px = j.px();
				_py = j.py();
				_pz = j.pz();
				_E = j.E();
				_mom = BayesPoint({_px, _py, _pz, _E});
				_eta = j.eta();
				_phi = j.phi();
				_t = j.time();
				
				_kt2 = j.kt2();
				_mass = j.mass();
				
				_idx = j.GetUserIdx();
				_vtx = j.GetVertex();
				_rhs = j.GetJetPoints();
				_nRHs = (int)_rhs.size();
				_constituents = j._constituents;
				_mu = j._mu;
				_cov = j._cov;
				_pi = j._pi;
		}

		//return four vector for clustering
		BayesPoint four_mom() const{ return _mom; }

		void SetVertex(BayesPoint vtx){
			if(vtx.Dim() != 3){
				cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
				return;
			}
			_vtx = vtx;
		}

		//setting the momentum of eg subclusters with track information
		void SetP(double px, double py, double pz){
			_px = px;
			_py = py;
			_pz = pz;
		
			_kt2 = px*px + py*py;
			_mass = mass();
		}
		//setting momentum to sum of constituents
		void SetP(){
			_px = 0;
			_py = 0;
			_pz = 0;
			for(int c = 0; c < _constituents.size(); c++){
				_px += _constituents[c].px();
				_py += _constituents[c].py();
				_pz += _constituents[c].pz();
			}
			_kt2 = _px*_px + _py*_py;
			_mass = mass();
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
		//[0,2pi]
		double phi() const{
			_ensure_valid_rap_phi();
			return phi_02pi();
		}

		//[-pi, pi]
		double phi_std() const{
			_ensure_valid_rap_phi();
			return phi_negPiToPi();
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
		//invariant mass - https://fastjet.fr/repo/doxygen-3.4.0/PseudoJet_8hh_source.html L1059
		//https://gitlab.cern.ch/CLHEP/CLHEP/-/blob/develop/Vector/Vector/LorentzVector.icc#L150
		double mass() const{return m2() < 0.0 ? -sqrt(-m2()) : sqrt(m2()); }

	
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

		//returns rec hits as jet points
		const vector<JetPoint>& GetJetPoints() const{return _rhs;}
		//returns rec hits as jets
		const void GetJets(vector<Jet>& rhs) const{
			rhs.clear();
			Jet rh;
			for(int r = 0; r < _rhs.size(); r++){
				rh = Jet(_rhs[r], _vtx);
				rh.SetVertex(_vtx);
				rh.SetUserIdx(_rhs[r].rhId());
				rhs.push_back(rh);	
			}
		}

		//add subjets/pixels to jet
		void add(const Jet& jt);
		void add(const JetPoint& rh);
		
		void GetEnergies(vector<double>& energies) const{ energies.clear(); for(int j = 0; j < (int)_rhs.size(); j++) energies.push_back(_rhs[j].E()); }
		void GetEtaPhiTimePoints(PointCollection& pc) const{
			pc.Clear();
			for(int i = 0; i < (int)_rhs.size(); i++){
				pc += BayesPoint({_rhs[i].eta(), _rhs[i].phi_02pi(), _rhs[i].time()});
			}
		}
		
		int GetNPoints() const{return (int)_rhs.size(); }	
		int GetNRecHits() const{return (int)_rhs.size(); }	
		

		void AddRecHit(const JetPoint& rh){
			_rhs.push_back(rh); 
			_nRHs = (int)_rhs.size();
			//recalculate time
			_set_time();
		}
		
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
		
		void GetClusterParams(Matrix& mu, Matrix& cov) const{ mu = _mu; cov = _cov; }
	
		//define jet time from cluster parameters
		double GetJetTime() const{ return _t; }
		void SetJetTime(double t){ _t = t; }
	
		BayesPoint GetVertex() const{return _vtx; }

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
			if(_phi > pi) return _phi - 2*pi;
			else return _phi;
			/*
			double o2pi = 1./(2*pi);
			if(fabs(_phi) <= pi)
				return _phi;
			double n = std::round(_phi * o2pi);
			return _phi - n * double(2.* pi);
			*/
		}

		void Print() const{
			cout << "px: " << _px << " py: " << _py << " pz: " << _pz << " E: " << _E << " mass: " << _mass << endl;
			//for(int i = 0; i < _nRHs; i++) _rhs[i].Print();
		}
		//constituents can be subclusters (from GMM) defined by eta, phi, time center, MM coefficient, and covariance matrix
		//also makes sure constituent and overall jet have same vertex
		//if jet is made of subclusters, set the jet's four vector to be weighted avg of subcluster four-vectors
		void AddConstituent(Jet& jt){
			//check to make sure jet isn't in constituents list
			auto it = find(_constituents.begin(), _constituents.end(), jt);
			if(it != _constituents.end()) return;	
		
			jt.SetVertex(_vtx);
			_constituents.push_back(jt);
			_constituents[_constituents.size()-1].SetVertex(_vtx);
			/*
			//reset overall cluster parameters
			double norm = 0; //should be 1
			if(_constituents.size() == 1){
				_cov = _constituents[0]._cov;
				_mu = _constituents[0]._mu;
				_pi = _constituents[0]._pi; //should be 1
				norm = 1;
			}
			else{ //set mean as probability-weighted center and covariance from that center
				for(int k = 0; k < _constituents.size(); k++){
					Matrix mu;
					mu.mult(_constituents[k]._mu, _constituents[k]._pi);
					_mu.add(mu);
					norm += _constituents[k]._pi;
				}
				_mu.mult(_mu,1./norm);		 
				//for covariance, take into account covariance of constituents with extra term x'_k = x_k +- sig_k
				_cov = Matrix(3,3);
				double xi, xj;
				for(int i = 0; i < 3; i++){
					for(int j = 0; j < 3; j++){
						double entry = 0;
						for(int k = 0; k < _constituents.size(); k++){
							double w = _constituents[k]._pi;
							xi = _constituents[k]._mu.at(i,0);
							xj = _constituents[k]._mu.at(j,0);
							entry += w*(xi - _mu.at(i,0))*(xj - _mu.at(j,0));
						}
						_cov.SetEntry(entry/norm,i,j);	
					}
				}
			}
			//phi wraparound
			//cout << "mu pre phi wraparound " << endl; _mu.Print();
			if(_mu.at(1,0) > 8*atan(1)) _mu.SetEntry(acos(cos(_mu.at(1,0))),1,0);
			//set fourvector by subcluster when hyperparameters are tuned s.t. there are enough subclusters in a jet
			_px = 0;
			_py = 0;
			_pz = 0;
			_E = 0;
			for(auto jet : _constituents){
				cout << "subcl px " << jet.px() << " py " << jet.py() << " pz " << jet.pz() << " E " << jet.E() << " coeff " << jet.GetCoefficient() << endl;
				_px += jet.GetCoefficient()*jet.px();
				_py += jet.GetCoefficient()*jet.py();
				_pz += jet.GetCoefficient()*jet.pz();
				_E  += jet.E();
				
			}
			_px /= norm;
			_py /= norm;
			_pz /= norm;
			_kt2 = _px*_px + _py*_py;
			_mass = mass();
			_update_mom();

		
			//set jet center as weighted sum of subclusters
			_eta = _mu.at(0,0);
			_phi = _mu.at(1,0);
			_t = _mu.at(2,0);
			//cout << "norm " << norm << " mu" << endl; _mu.Print();
			//cout << "const mu" << endl; _constituents[0]._mu.Print();

			_ensure_valid_rap_phi();
//cout << "MAKING jet subcl kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << pt() << " m2 " << m2() <<  endl;

		*/
		}
		//since the GMM has probabilistic assignment of points, these jets will be defined by their center and cov
		vector<Jet>& GetConstituents(){
			return _constituents;
		}
	
		Jet& GetConstituent(int c){
			if(c > _constituents.size()){
				cout << "Error: index " << c << " out of bounds for # of constituents " << _constituents.size() << endl;
				return *this;
			}
			return _constituents[c];
		}
	

		double GetCoefficient() const{ return _pi; }

		//jet_var = sum_i pi_i*var_i / sum_i pi_i
		Matrix GetCovariance() const{
			return _cov;
		}			

		//jet_center = sum_i pi_i*center_i / sum_i pi_i
		Matrix GetCenter() const{
			return _mu;
		}			


		int GetNConstituents() const{
			//if only 1 subcluster, jet is its own constituent
			if(_constituents.size() == 0 && !_cov.empty()) return 1;
			else return (int)_constituents.size();
		}
		void SetCenter(Matrix& mu){
			_eta = mu.at(0,0);
			_phi = mu.at(1,0);
			_t = mu.at(2,0);
			_mu = mu;
		}
		void SetCovariance(Matrix& cov){
			_cov = cov;
		}
		void SetCoefficient(double p){
			_pi = p;
		}
		
		//calculate invariant mass with jet
		double invMass(Jet jet){
			jet.add(*this);
			return jet.mass();	
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


		void _update_mom(){
			_mom.SetValue(_px,0);
			_mom.SetValue(_py,1);
			_mom.SetValue(_pz,2);
			_mom.SetValue(_E,3);
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
		BayesPoint _mom;
		BayesPoint _vtx;

		//mutable double _eta, _phi, _theta;
		//rec hits (ids) in jet
		vector<JetPoint> _rhs;

		int _nRHs;

		//cluster params - modeling jet as gaussian in time, space + energy
		Matrix _mu;
		Matrix _cov;
		double _pi;

		//vector of subjets (NOT rhs)
		vector<Jet> _constituents;	

		
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
