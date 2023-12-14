#ifndef JETPOINT_HH
#define JETPOINT_HH

#include "Point.hh"
#include <math.h>

class JetPoint{

	public:
		JetPoint();
		JetPoint(double x, double y, double z, double t);
		JetPoint(const JetPoint& jp);
		virtual ~JetPoint();

		bool operator==(JetPoint& j) const;
		bool operator!=(JetPoint& j) const;


		void SetEnergy(double e){ _E = e;}
		void SetPhi(double p){ _phi = p; _ensure_valid_rap_phi(); }
		void SetEta(double e){ _eta = e; _ensure_valid_rap_phi(); }
		void SetWeight(double w){ _w = w; } 
		double GetWeight(){ return _w; } 


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
			return _phi;
		}


		double e() const{return _E;}
		double E() const{return _E;}
		double energy() const{return _E;}
		double Energy() const{return _E;}

		double theta() const{return _theta;}
		double Theta() const{return _theta;}

		double time() const{return _t;}
		double t() const{return _t;}
		double x() const{return _x;}
		double y() const{return _y;}
		double z() const{return _z;}


		double perp2() const{return _x*_x + _y*_y;}	

		//calculations taken from ROOT
		// Authors: W. Brown, M. Fischler, L. Moneta    2005
		 
		 /**********************************************************************
		  *                                                                    *
		  * Copyright (c) 2005 , FNAL MathLib Team                             *
		  *                                                                    *
		  *                                                                    *
		  **********************************************************************/
		 
		 
		// Header source file for function calculating eta
		//
		// Created by: Lorenzo Moneta  at 14 Jun 2007
		/**
		Calculate eta given rho and zeta.
		This formula is faster than the standard calculation (below) from log(tan(theta/2)
		but one has to be careful when rho is much smaller than z (large eta values)
		Formula is  eta = log( zs + sqrt(zs^2 + 1) )  where zs = z/rho
		
		For large value of z_scaled (tan(theta) ) one can appoximate the sqrt via a Taylor expansion
		We do the approximation of the sqrt if the numerical error is of the same order of second term of
		the sqrt.expansion:
		eps > 1/zs^4   =>   zs > 1/(eps^0.25)
		
		When rho == 0 we use etaMax (see definition in etaMax.h)
		*/
		/**
		etaMax: Function providing the maximum possible value of pseudorapidity for
		a non-zero rho, in the Scalar type with the largest dynamic range.
		*/
		double _set_rap() const{
			double rho = sqrt(perp2());
			//from ROOT
			double etaMax = 22756.0;
			if (rho > 0) {
				// value to control Taylor expansion of sqrt
				static const double big_z_scaled = pow(std::numeric_limits<double>::epsilon(), -.25);
				double z_scaled = _z/rho;
				if (std::fabs(z_scaled) < big_z_scaled) {
				   return log(z_scaled + std::sqrt(z_scaled * z_scaled + 1.0));
				} else {
				   // apply correction using first order Taylor expansion of sqrt
				   return _z > 0 ? log(2.0 * z_scaled + 0.5 / z_scaled) : -log(-2.0 * z_scaled);
				}
			}
			// case vector has rho = 0
			else if (_z==0) {
			   return 0;
			}
			else if (_z>0) {
			   return _z + etaMax;
			}
			else {
			   return _z - etaMax;
			}
 
		}
		void _set_rap_phi() const{
			_eta = _set_rap();
			_phi = atan2(_y, _x);
		}
		void _ensure_valid_rap_phi() const{
			if(_phi == _invalid_phi) _set_rap_phi();
		}


		void Print() const{
                        std::string out = "(";
                        out += std::to_string(eta())+",";
                        out += std::to_string(phi_02pi())+",";
                        out += std::to_string(time())+")";
                        cout << out << " w = " << _w << endl;
		
		}

		
		void SetRecHitId(unsigned int id){ _rhId = id; }
		unsigned int rhId(){ return _rhId; }

		//set user idx info
		void SetUserIdx(int i){ _idx = i; }
		int userIdx(){ return _idx; }
		

		double _maxRap = 1e5;
		//default values - have yet to be calculated or set
		double _invalid_phi = -100.0;
		double _invalid_eta = -1e200;
		double twopi = 6.28318530717;


	protected:
		//spatial four vector
		double _t;
		double _x;
		double _y;
		double _z;
		
		mutable double _eta, _phi, _theta;
		double _E;


		//rechit info
		unsigned int _rhId;
		
		//user idx info
		int _idx;


		//weight of point
		double _w;


};

#endif
