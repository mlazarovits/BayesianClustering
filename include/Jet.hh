#ifndef JET_HH
#define JET_HH

#include <math.h>
#include "Matrix.hh"
#include "JetPoint.hh"
#include "Point.hh"

//point with more physics aspects - ecal cell
class Jet{
	public:
		Jet();
		Jet(JetPoint rh);
		Jet(JetPoint rh, Point vtx);
		Jet(vector<JetPoint> rhs, Point vtx);
		Jet(vector<JetPoint> rhs);
		virtual ~Jet();		

		bool operator==(Jet& j) const;
		bool operator!=(Jet& j) const;

		//return four vector for clustering
		Point four_mom(){ return m_p; }
		Point four_space(){ return m_space; }

		void SetFourMom(Point pt);
		void SetFourPos(Point pt);
		void SetVertex(Point vtx){
			if(vtx.Dim() != 3){
				cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
				return;
			}
			m_vtx = vtx;
		}

		//return element i in four vector
		double p(int i){ return m_p.at(i); }
		double x(int i){ return m_space.at(i); }
	
		//scale momentum
		void scaleMom(double s){
			m_px *= s;
			m_py *= s;
			m_pz *= s;
			m_E *= s;
			m_kt2 *= s*s;
		 };


		vector<JetPoint> GetJetPoints() const{return m_rhs;}
		
		PointCollection GetJetPoints_mom() const{
			PointCollection pc;
			//for(int i = 0; i < m_nRHs; i++) pc += m_rhs[i].four_mom();
			return pc;
		}

		//add subjets/pixels to jet
		void add(Jet& jt);
		void add(JetPoint& rh);
		
		//constituents (jet points) in jet (clustered or unclustered)
		void GetConstituents(vector<JetPoint>& rhs) const { rhs.clear(); rhs = m_rhs; }
		void GetEnergies(vector<double>& energies) const{ energies.clear(); for(int j = 0; j < (int)m_rhs.size(); j++) energies.push_back(m_rhs[j].E()); }
		void GetXYZConstituents(PointCollection& pc) const{
			pc.Clear();
			for(int i = 0; i < (int)m_rhs.size(); i++){
				pc += m_rhs[i].four_space();
			}
		}
		void GetEtaPhiConstituents(PointCollection& pc) const{
			pc.Clear();
			for(int i = 0; i < (int)m_rhs.size(); i++){
				pc += Point({m_rhs[i].eta(), m_rhs[i].phi(), m_rhs[i].time()});
			}
		}
		
		int GetNConstituents() const{return (int)m_rhs.size(); }	
		
		//parents in cluster
		void GetParents(Jet& p1, Jet& p2) const;
		void SetParents(Jet* p1, Jet* p2){ m_parent1 = p1; m_parent2 = p2; p1->SetBaby(this); p2->SetBaby(this); }

		//children in cluster
		void GetBaby(Jet& child) const;

		
		//subjets (jets) in jet (clustered or unclustered)
		void GetSubJets(vector<Jet>& subjets, int depth = 0) const;

		//this has jet?
		bool Has(Jet& jet) const;
		//this has rh?
		bool Has(JetPoint& rh) const;
		
		void SetEta(double eta){ m_eta = eta; }
		void SetPhi(double phi){ 
			m_phi = setPhi_negPiToPi(phi);
		}

		//sets phi [0,2pi]
		double setPhi_02pi(double phi){
			if(phi < 0.0) return phi + 2*acos(-1);
			else if(phi >= 2*acos(-1)) return phi - 2*acos(-1);
			return phi; 
		}

		//wraps phi around pi, [-pi,pi]
		double setPhi_negPiToPi(double phi){
 			double pi = acos(-1);
			double o2pi = 1./(2*pi);
			if(fabs(phi) <= pi)
				return phi;
			double n = std::round(phi * o2pi);
			return phi - n * double(2.* pi);

		}

		void SetEnergy(double enr){ m_E = enr; }
		
		void SetPt(double pt){ m_kt2 = pt*pt; }

		//set user idx info
		void SetUserIdx(int i){ m_idx = i; }
		int GetUserIdx(){ return m_idx; }
		

		//boost function
		//transverse energy
		double Et() const{ return (m_kt2==0) ? 0.0 : m_E/sqrt(1.0+m_pz*m_pz/m_kt2); } 
		//transverse mass
		double mt() const {return sqrt(std::abs(mperp2()));}
  		double mperp() const {return sqrt(std::abs(mperp2()));}
		//squared transverse mass = kt^2+m^2
  		double mperp2() const {return (m_E+m_pz)*(m_E-m_pz);}
		//invariant mass squared: m^2 = E^2 - p^2
		double m2() const{ return (m_E+m_pz)*(m_E-m_pz)-m_kt2; }
		//invariant mass
		double mass() const{return sqrt(m2()); }

	
		//squared transverse momentum
  		double pt2() const {return m_kt2;}
  		//the scalar transverse momentum
  		double pt() const {return sqrt(m_kt2);} 
	 	//the squared transverse momentum
  		double kt2() const {return m_kt2;} 	
		//deltaR between this and another jet pt
		double deltaR(Jet& jet) const{ return sqrt( (m_eta - jet.eta())*(m_eta - jet.eta()) + (m_phi - jet.phi())*(m_phi - jet.phi())); }
		//eta
		double eta() const{return m_eta;}
		double phi() const{return m_phi;}
		
		void GetClusterParams(Matrix& mu, Matrix& cov){ mu = m_mu; cov = m_cov; }
	
		//define jet time from cluster parameters
		double GetJetTime() const{ return 0.; }
	
		Point GetVertex() const{return m_vtx; }

		//check IR + collinear safety?
	
	protected:
		void SetBaby(Jet* child){ m_child = child; }

		void RecalcKT2(){ m_kt2 = m_px*m_px + m_py*m_py; }
		void RecalcPhi(){ if (m_kt2 == 0.0) {
			m_phi = 0.0; } 
			else {
			  m_phi = atan2(m_py,m_px);
			}
			m_phi = setPhi_negPiToPi(m_phi);
		}
		
		void RecalcEta(){
			if(m_pz == 0){
				m_theta = 0;
				m_eta = m_maxRap;
			}
			else{
				m_theta = atan2(m_kt2 , m_pz);
				if(m_theta < 0) m_theta += acos(-1);
				m_eta = -log(tan(m_theta/2.));
			}

		}
	private:
		//momentum four vector
		double m_E;
		double m_px;
		double m_py;
		double m_pz;
		double m_kt2;
		Point m_p;

		//spatial four vector
		double m_t;
		double m_x;
		double m_y;
		double m_z;
		Point m_space;

		//vertex from which momemtum direction is calculated
		Point m_vtx;

		//rec hits (JetPoints) in jet
		vector<JetPoint> m_rhs;

		int m_nRHs;

		double m_eta;
		double m_phi;
		double m_theta;

		double m_maxRap = 1e5;
		//default values - have yet to be calculated or set
		double m_invalid_phi = -100.0;
		double m_invalid_eta = -1e200;


		//cluster params - modeling jet as gaussian in time, space + energy
		Matrix m_mu;
		Matrix m_cov;

		//parents + child info
		Jet* m_parent1;
		Jet* m_parent2;

		Jet* m_child;
	
		//user index	
		int m_idx;


};

#endif
