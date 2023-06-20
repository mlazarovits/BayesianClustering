#ifndef JETPOINT_HH
#define JETPOINT_HH

#include "Point.hh"

//point with more physics aspects
class JetPoint{

	public:
		JetPoint();
		JetPoint(double E, double px, double py, double pz);
		JetPoint(Point pt, bool mom = true);
		JetPoint(double t, double x, double y, double z);
		virtual ~JetPoint(){ };

		//return four vector for clustering
		Point four_mom(){ return m_momvec; }
		Point four_space(){ return m_spacevec; }

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
  		double mperp2() const {return (_E+_pz)*(_E-_pz);}
		
		//squared transverse momentum
  		double pt2() const {return _kt2;}
  		//the scalar transverse momentum
  		double  pt() const {return sqrt(_kt2);} 
	 	//the squared transverse momentum
  		double kt2() const {return _kt2;} 	
		//deltaR between this and another jet pt
		double deltaR(JetPoint& jt) const{ return sqrt( (m_eta - jet.eta())*(m_eta - jet.eta()) + (m_phi - jet.phi())*(m_phi - jet.phi())); }
		//eta
		double eta() const{return m_eta;}
		doube phi() const{return m_phi;}

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















};

#endif
