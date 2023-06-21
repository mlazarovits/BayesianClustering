#ifndef JET_HH
#define JET_HH

#include "Point.hh"
#include "Matrix.hh"

//point with more physics aspects - ecal cell
class Jet{

	public:
		Jet();
		//Jet(double E, double px, double py, double pz);
		Jet(Point pt, bool mom = true);
		virtual ~Jet();

		//return four vector for clustering
		Point four_mom(){ return m_momvec; }
		Point four_space(){ return m_spacevec; }

		void SetFourMom(Point pt);
		void SetFourPos(Point pt);

		//return element i in four vector
		double mom_at(int i){ return m_momvec.at(i); }
		double space_at(int i){ return m_spacevec.at(i); }
		//have scaling function
		void scale(double s);
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

		void SetEta(double eta){ m_eta = eta; }
		void SetPhi(double phi){ 
			if(phi < 0.0) m_phi = phi + 2*acos(-1);
			else if(phi >= 2*acos(-1)) m_phi = phi - 2*acos(-1);
			else m_phi = phi; }
		
		void SetRecHitId(unsigned int id){ m_rhId = id; }
		unsigned int rhId(){ return m_rhId; }

		//parents in cluster
		void GetParents(Jet& p1, Jet& p2);

		//children in cluster
		void GetBaby(Jet& child);

		//constituents in jet (clustered or unclustered)
		void GetSubjets(vector<Jet>& subjets, int depth = 0);

		//this has jet?
		bool Has(Jet& jet);

		void GetClusterParams(Matrix& mu, Matrix& cov){ mu = m_mu; cov = m_cov; }
	
		//define jet time from cluster parameters
		double GetJetTime(){ return 0.; }

		//check IR + collinear safety?

	private:
		//momentum four vector
		double m_E;
		double m_px;
		double m_py;
		double m_pz;
		double m_kt2;
		Point m_momvec;

		//spatial four vector
		double m_t;
		double m_x;
		double m_y;
		double m_z;
		Point m_spacevec;

		double m_eta;
		double m_phi;
		double m_theta;

		double m_maxRap = 1e5;
		//default values - have yet to be calculated or set
		double m_invalid_phi = -100.0;
		double m_invalid_eta = -1e200;

		//rechit info
		unsigned int m_rhId;


		//cluster params - modeling jet as gaussian in time, space + energy
		Matrix m_mu;
		Matrix m_cov;

		//parents + child info
		Jet* m_parent1;
		Jet* m_parent2;

		Jet* m_child;









};

#endif
