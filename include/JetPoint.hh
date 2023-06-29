#ifndef JETPOINT_HH
#define JETPOINT_HH

#include "Point.hh"

//point that is physics-specific - ecal cell (rechit)
class JetPoint{

	public:
		JetPoint();
		JetPoint(double x, double y, double z, double t);
		JetPoint(Point pt);
		virtual ~JetPoint();

		bool operator==(JetPoint& j) const;
		bool operator!=(JetPoint& j) const;

		//return four vector for clustering
		Point four_space() const{ return m_vec; }

		void SetFourSpace(Point pt);

		//return element i in position (space-time) four vector
		double x(int i) const{ return m_vec.at(i); }


		//eta
		double eta() const{return m_eta;}
		double phi() const{return m_phi;}
		double Eta() const{return m_eta;}
		double Phi() const{return m_phi;}

		double e() const{return m_E;}
		double E() const{return m_E;}
		double energy() const{return m_E;}
		double Energy() const{return m_E;}

		double theta() const{return m_theta;}
		double Theta() const{return m_theta;}

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
		
		void SetRecHitId(unsigned int id){ m_rhId = id; }
		unsigned int rhId(){ return m_rhId; }

		//set user idx info
		void SetUserIdx(int i){ m_idx = i; }
		int userIdx(){ return m_idx; }


	private:
		//spatial four vector
		double m_t;
		double m_x;
		double m_y;
		double m_z;
		Point m_vec;

		double m_eta;
		double m_phi;
		double m_theta;
		double m_E;


		//rechit info
		unsigned int m_rhId;
		
		//user idx info
		int m_idx;

		double m_maxRap = 1e5;
		//default values - have yet to be calculated or set
		double m_invalid_phi = -100.0;
		double m_invalid_eta = -1e200;





};

#endif
