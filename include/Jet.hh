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
		Jet(const Matrix& mu, const Matrix& cov, double E, double _pi = 1, BayesPoint vtx = BayesPoint({0., 0., 0.}), double detR = 129); //constructor from subcluster information - detR in cm
		Jet(BasePDF* pdf, double E, double pi = 1, BayesPoint vtx = BayesPoint({0., 0., 0.}), double detR = 129); //constructor from subcluster information - detR in cm
		Jet(BasePDFMixture* model, BayesPoint vtx, double gev, double detR = 129); //need detector radius to convert eta, phi to x, y, z - detR in cm
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
				_eta = j._eta;
				_phi = j._phi;
				_t = j._t;
				
				_kt2 = j._kt2;
				_mass = j._mass;
				
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
		
		//recalculate momentum from PV
		void CalculatePfromPV(double detR = 129){	
			//get x, y, z from eta, phi
			double x = detR*cos(_phi);
			double y = detR*sin(_phi);
			double theta = 2*atan2(1,exp(_eta));
			double z = detR/tan(theta);

			//calculate momentum vector from PV
			//centered at PV
			double dx = x - _vtx.at(0);
			double dy = y - _vtx.at(1);
			double dz = z - _vtx.at(2);
			//theta is calculated between beamline (z-dir) and x-y vector	
			double p_theta = atan2( sqrt(dx*dx + dy*dy), dz );
			double p_eta = -log(tan(p_theta/2));
			double p_phi = atan2(dy, dx);
			//double pt = _E*sin(theta); //mass = 0
			double pt = _E/cosh(p_eta);
			_px = pt*cos(p_phi);
			_py = pt*sin(p_phi);
			_pz = pt*sinh(p_eta);

			_kt2 = _px*_px + _py*_py;
			_mass = _calc_mass();
		}

		//setting the momentum of eg subclusters with track information
		void SetP(double px, double py, double pz){
			_px = px;
			_py = py;
			_pz = pz;
		
			_kt2 = px*px + py*py;
			_mass = _calc_mass();
		}
		
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
		double _calc_mass() const{return m2() < 0.0 ? -sqrt(-m2()) : sqrt(m2()); }
		double mass() const{ return _mass; }
		double m() const{return _mass; }
		
		//different mass calculation types for comparison
		double mass_rhs() const{
			double px, py, pz, E;
			E = 0;
			px = 0;
			py = 0;
			pz = 0;
			for(int i = 0; i < _nRHs; i++){		
				//theta is calculated between beamline (z-dir) and vector in x-y plane	
				//centered at (0,0,0)
				double x = _rhs[i].x() - _vtx.at(0);
				double y = _rhs[i].y() - _vtx.at(1);
				double z = _rhs[i].z() - _vtx.at(2);
				double theta = atan2( sqrt(x*x + y*y), z );
				double phi = atan2(y, x);
				double eta = -log(tan(theta/2.)); 
				//see https://cmssdt.cern.ch/lxr/source/DataFormats/CaloTowers/src/CaloTower.cc L145
				//pt = E*sin(theta); //mass = 0, equivalent to below
				double pt = _rhs[i].E()/cosh(eta);
				px += pt*cos(phi);
				py += pt*sin(phi);
				pz += pt*sinh(eta);
			//cout << "i " << i << " px " << pt*cos(phi) << " py " << pt*sin(phi) << " pz " << pt*cosh(eta) << " " << pt*sinh(eta) << endl;
				
				E += _rhs[i].E();

			}
			double kt2 = px*px + py*py;
			double m2 = (E+pz)*(E-pz)-kt2;
			return m2 < 0.0 ? -sqrt(-m2) : sqrt(m2); 

		}
		double mass_subcls() const{
			double px, py, pz, E;
			E = 0;
			px = 0;
			py = 0;
			pz = 0;
			for(int k = 0; k < GetNConstituents(); k++){
				Jet subcl = _constituents[k];
	
				//set momentum from subclusters
				px += subcl.px();
				py += subcl.py();
				pz += subcl.pz();
				E += subcl.E(); //set from rhs
			
			}
			double kt2 = px*px + py*py;
			double m2 = (E+pz)*(E-pz)-kt2;
			return m2 < 0.0 ? -sqrt(-m2) : sqrt(m2); 

		}



	
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
				rh.CalculatePfromPV();
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
			else if(_phi >= 2*acos(-1)) return _phi - 2*acos(-1);
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
		
		//constituents are subclusters (from GMM) defined by eta, phi, time center, MM coefficient, and covariance matrix
		//also makes sure constituent and overall jet have same vertex
		//if jet is made of subclusters, set the jet's four vector to be weighted avg of subcluster four-vectors
		//this does not affect the rhs within a jet, when a jet is created, the list of rhs or GMM sets the rhs
		//these points set the jet's four-vector + overall _mu + _cov
		void AddConstituent(Jet& jt){
			//check to make sure jet isn't in constituents list
			auto it = find(_constituents.begin(), _constituents.end(), jt);
			if(it != _constituents.end()) return;	
		
			jt.SetVertex(_vtx);
			//jt.CalculatePfromPV();
			_constituents.push_back(jt);
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
		
		//calculate invariant mass with jet (without adding all the jet components - ie rechits, recalc time, etc)
		double invMass(const Jet& jet){
			double px = _px + jet.px();
			double py = _py + jet.py();
			double pz = _pz + jet.pz();
			double E = _E + jet.E();
			double kt2 = px*px + py*py;

			double m2 = (E+pz)*(E-pz)-kt2;
			return m2 < 0.0 ? -sqrt(-m2) : sqrt(m2);
			
	
		}
		
		void CalculateCenter(){
			double norm = 0;
			PointCollection phipts;
			_mu = Matrix(3,1);
			_eta = 0;
			_phi = 0;
			_t = 0;
			for(int i = 0; i < _nRHs; i++){		
				//eta, phi centered at (0,0,0)
				_eta += _rhs[i].eta()*_rhs[i].E();
				_t += _rhs[i].t()*_rhs[i].E();
				norm += _rhs[i].E();
				
				BayesPoint phipt(1);
				phipt.SetValue(_rhs[i].phi_02pi(),0);
				phipt.SetWeight(_rhs[i].E());
				
				phipts += phipt;

			}
			_eta /= norm;
			_t /= norm;
			_phi = phipts.CircularCentroid(0);
			if (_phi < 0.0) {_phi += twopi;}
			if (_phi >= twopi) {_phi -= twopi;} // can happen if phi=-|eps<1e-15|?
			_ensure_valid_rap_phi();
			
			//set mean from member variables
			_mu.SetEntry(_eta,0,0);
			_mu.SetEntry(_phi,1,0);
			_mu.SetEntry(_t,2,0);
		}


		void CalculateCovariance(){
			double deta, dphi, dtime, eta_phi, eta_time, phi_time;
			//start from scratch - clear cov
			_cov = Matrix(3,3);
			double norm = 0;
			//energy-weighted rh cov
			for(int i = 0; i < _nRHs; i++){
				deta = _rhs[i].eta() - _mu.at(0,0);	
				dphi = _rhs[i].phi() - _mu.at(1,0);	
				dphi = acos(cos(dphi));
				dtime = _rhs[i].t() - _mu.at(2,0);
			
				norm += _rhs[i].E();
				Matrix cov_entry = Matrix(3,3);
				cov_entry.SetEntry(_rhs[i].E()*deta*deta,0,0);
				cov_entry.SetEntry(_rhs[i].E()*deta*dphi,1,0);
				cov_entry.SetEntry(_rhs[i].E()*deta*dtime,2,0);
				cov_entry.SetEntry(_rhs[i].E()*dphi*deta,0,1);
				cov_entry.SetEntry(_rhs[i].E()*dphi*dphi,1,1);
				cov_entry.SetEntry(_rhs[i].E()*dtime*dphi,2,1);
				cov_entry.SetEntry(_rhs[i].E()*deta*dtime,0,2);
				cov_entry.SetEntry(_rhs[i].E()*dphi*dtime,1,2);
				cov_entry.SetEntry(_rhs[i].E()*dtime*dtime,2,2);
				

				_cov.add(cov_entry);	

			}
			_cov.mult(_cov,1/norm);	
		}


		//add mixture model to associated jet that is already made (ie from AK4, etc) with correct transfer factor
		//jet needs to have rechits set already to correctly add them to respective subclusters
		void SetModel(BasePDFMixture* model, double gev){
			Matrix r_nk = model->GetPosterior();
			vector<double> norms;
			model->GetNorms(norms);
			for(int k = 0; k < model->GetNClusters(); k++){
				Jet subcl(model->GetModel(k), norms[k]/gev, model->GetPi(k), _vtx); 
				//add rechits as "effective" crystals (ie weighted by their responsibility to this cluster)
				//their associated responsibility will be saved as the energy of that crystal for this subcluster
				double subcl_px = 0;
				double subcl_py = 0;
				double subcl_pz = 0;

				//posterior is r_nk*w_n for each pt n s.t. sum_k r_nk*w_n = w_n -> need just 
				for(int n = 0; n < _nRHs; n++){
					JetPoint effRh = _rhs[n];
					effRh.SetWeight(r_nk.at(n,k)/(_rhs[n].E()*gev));
					effRh.SetEnergy(_rhs[n].E()*effRh.GetWeight());
					subcl.AddRecHit(effRh);
					//recalculate momentum vector accordingly
					//calculate momentum vector from PV
					//centered at PV
					double dx = effRh.x() - _vtx.at(0);
					double dy = effRh.y() - _vtx.at(1);
					double dz = effRh.z() - _vtx.at(2);
					//theta is calculated between beamline (z-dir) and x-y vector	
					double p_theta = atan2( sqrt(dx*dx + dy*dy), dz );
					double p_eta = -log(tan(p_theta/2));
					double p_phi = atan2(dy, dx);

					//cout << "eta " << effRh.eta() << " p_eta " << p_eta << " phi " << effRh.phi() << " p_phi " << p_phi << " x " << effRh.x() << " dx " << dx << " y " << effRh.y() << " dy " << dy << " z " << effRh.z() << " dz " << dz << " PV x " << _vtx.at(0) << " PV y " << _vtx.at(1) << " PV z " << _vtx.at(2) << endl;
					//double pt = _E*sin(theta); //mass = 0
					double pt = effRh.E()/cosh(p_eta); 
					subcl_px += pt*cos(p_phi);
					subcl_py += pt*sin(p_phi);
					subcl_pz += pt*sinh(p_eta);
				}
				//set subcluster momentum three-vector and mass
				subcl.SetP(subcl_px, subcl_py, subcl_pz);
				_constituents.push_back(subcl);
			}

		}

		void Get2DMat(const Matrix& inmat, Matrix& outmat){
			if(!outmat.square()) return;
			if(outmat.GetDims()[0] != 2) return;
			outmat.reset();
			outmat.SetEntry(inmat.at(0,0),0,0);	
			outmat.SetEntry(inmat.at(0,1),0,1);	
			outmat.SetEntry(inmat.at(1,0),1,0);	
			outmat.SetEntry(inmat.at(1,1),1,1);
		}
		double CalcSpatialSize(){
			if(_cov.GetDims()[0] != 3 || _cov.GetDims()[1] != 3){
				cout << "Error: can't calculate size for matrix of size " << _cov.GetDims()[0] << " x " << _cov.GetDims()[1] << endl;
				return -1;
			}
			vector<double> eigvals;
			vector<Matrix> eigvecs;
			Matrix cov2D(2,2);
			Get2DMat(_cov,cov2D);	
			cov2D.eigenCalc(eigvals, eigvecs);
			//define jet size as length of major axis
			return sqrt(eigvals[1]);
		}

		//add PU cleaning method
		//if remove == false, rechits are downweighted by 1 - r_nk for each subcluster k that doesnt pass PU cleaning reqs
		Jet CleanOutPU(bool remove = false){
			if(_constituents.size() < 2) return *this; //if no subclusters or only 1, return current jet

			Matrix cov = GetCovariance();
			Jet cleanedJet;
			vector<JetPoint> cleanedRhs;
			//loop through constituents and reset rh weights based on above
			//loop through rhs - for rechit n
			//DOWNWEIGHTING
			// if max(r_nk) belongs to subcluster k that FAILS PU cleaning criteria, weight energy by 1 - max(r_nk)
			// if max(r_nk) belongs to subcluster k that PASSES PU cleaning criteria, weight energy by max(r_nk)
			//REMOVING
			// if max(r_nk) belongs to subcluster k that FAILS PU cleaning criteria, weight energy by 0.
			// if max(r_nk) belongs to subcluster k that PASSES PU cleaning criteria, weight energy by 1.
			//(unweighted) r_nk's are saved as rh weights
			vector<bool> pass; //true = pass, false = fail
			for(int k = 0; k < _constituents.size(); k++){
				double relE = _constituents[k].e() / this->e();
				Matrix subcl_cov = _constituents[k].GetCovariance();
				double relTimeVar = sqrt(subcl_cov.at(2,2)) / sqrt(cov.at(2,2));
				double relSize = _constituents[k].CalcSpatialSize() / this->CalcSpatialSize();

				if(relTimeVar - relE <= 0 && relTimeVar + relSize <= 1)
					pass.push_back(true);
				else
					pass.push_back(false);
			cout <<	"subcluster #" << k << " rel time var " << relTimeVar << " rel E " << relE << " rel size " << relSize << " pass? " << pass[k] << endl;	
			}
			//return empty jet if no subclusters pass criteria
			if(find(pass.begin(), pass.end(), true) == pass.end()){
				Jet ret;
				return ret;
			}

			double totE = 0;
			for(int n = 0; n < _nRHs; n++){
				double maxRnk = 0;
				int assignedK = -1;
				JetPoint effRh;
				double totR = 0;
				//cout << "rh #" << n;
				for(int k = 0; k < _constituents.size(); k++){
					effRh = _constituents[k]._rhs[n];
					if(effRh.GetWeight() > maxRnk){
						maxRnk = effRh.GetWeight();
						assignedK = k;
					}
					if(pass[k]) //responsibility of good subclusters sum to overall weight
						totR += effRh.GetWeight();
				}
			
				if(remove){
					if(pass[assignedK]){
						effRh.SetWeight(1.);
					}
					else
						effRh.SetWeight(0.);

				}
				else{ //downweight
					effRh.SetWeight(totR);
				}
				effRh.SetEnergy(_rhs[n].E()*effRh.GetWeight());
				totE += effRh.e();

				cleanedRhs.push_back(effRh);
			}
			if(totE == 0){
				Jet ret;
				return ret;
			}
			cleanedJet = Jet(cleanedRhs, _vtx);
			for(int k = 0; k < pass.size(); k++){
				if(remove){
					if(pass[k])
						cleanedJet.AddConstituent(_constituents[k]);
				}
				else{
					cleanedJet.AddConstituent(_constituents[k]);
				}
			}	
			return cleanedJet;
		}
		/*
		Jet CleanOutPU(double maxRelTimeVar = 1, double minRelPt = 0.2, bool remove = false){
			if(_constituents.size() < 2) return *this; //if no subclusters or only 1, return current jet

			Matrix cov = GetCovariance();
			Jet cleanedJet;
			vector<JetPoint> cleanedRhs;
			//loop through constituents and reset rh weights based on above
			//loop through rhs - for rechit n
			//DOWNWEIGHTING
			// if max(r_nk) belongs to subcluster k that FAILS PU cleaning criteria, weight energy by 1 - max(r_nk)
			// if max(r_nk) belongs to subcluster k that PASSES PU cleaning criteria, weight energy by max(r_nk)
			//REMOVING
			// if max(r_nk) belongs to subcluster k that FAILS PU cleaning criteria, weight energy by 0.
			// if max(r_nk) belongs to subcluster k that PASSES PU cleaning criteria, weight energy by 1.
			//(unweighted) r_nk's are saved as rh weights
			vector<bool> pass; //true = pass, false = fail
			for(int k = 0; k < _constituents.size(); k++){
				double rel_subcl_pt = _constituents[k].pt() / this->pt();
				Matrix subcl_cov = _constituents[k].GetCovariance();
				
				double rel_subcl_size = sqrt(subcl_cov.at(2,2)) / sqrt(cov.at(2,2));

				if(rel_subcl_pt > minRelPt && rel_subcl_size < maxRelTimeVar)
					pass.push_back(true);
				else
					pass.push_back(false);
			}
			//return empty jet if no subclusters pass criteria
			if(find(pass.begin(), pass.end(), true) == pass.end()){
				Jet ret;
				return ret;
			}

			double totE = 0;
			for(int n = 0; n < _nRHs; n++){
				double maxRnk = 0;
				int assignedK = -1;
				JetPoint effRh;
				double totR = 0;
				for(int k = 0; k < _constituents.size(); k++){
					effRh = _constituents[k]._rhs[n];
					if(effRh.GetWeight() > maxRnk){
						maxRnk = effRh.GetWeight();
						assignedK = k;
					}
					if(pass[k]) //responsibility of good subclusters sum to overall weight
						totR += effRh.GetWeight();
				}
			
				if(remove){
					if(pass[assignedK]){
						effRh.SetWeight(1.);
					}
					else
						effRh.SetWeight(0.);

				}
				else{ //downweight
					effRh.SetWeight(totR);
				}
				effRh.SetEnergy(_rhs[n].E()*effRh.GetWeight());
				totE += effRh.e();

				cleanedRhs.push_back(effRh);
			}
			if(totE == 0){
				Jet ret;
				return ret;
			}
			cleanedJet = Jet(cleanedRhs, _vtx);
			for(int k = 0; k < pass.size(); k++){
				if(remove){
					if(pass[k])
						cleanedJet.AddConstituent(_constituents[k]);
				}
				else{
					cleanedJet.AddConstituent(_constituents[k]);
				}
			}	
			return cleanedJet;
		}
		*/

	
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
