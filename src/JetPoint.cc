#include "JetPoint.hh"


JetPoint::JetPoint(){
	_x = -999;
	_y = -999;
	_z = -999;
	_t = -999;

	_eta = _invalid_eta;
	_phi = _invalid_phi;
	_theta = -999;
	_E = -999;
	
	_idx = -999;
	_rhId = -999;
	
	_w = 1;

}




JetPoint::JetPoint(double x, double y, double z, double t){
	_x = x;
	_y = y;
	_z = z;
	_t = t;

	_phi = atan2(y,x);//_invalid_phi;
	double rho = sqrt(x*x + y*y);
	_eta = atan2(rho,z);//_invalid_eta;
	_theta = acos( z / sqrt(x*x + y*y + z*z) );	
	_E = 0;

	_idx = -999;
	_rhId = -999;
	_w = 1;

}


JetPoint::JetPoint(const JetPoint& jp){
	_x = jp._x;
	_y = jp._y;
	_z = jp._z;
	_t = jp._t;
	
	_phi = jp._phi;
	_eta = jp._eta;
	_theta = jp._theta;
	_E = jp._E;
	

	_idx = jp._idx;
	_rhId = jp._rhId;
	_w = jp._w;
}



JetPoint::~JetPoint(){
}


bool JetPoint::operator ==(JetPoint& jet) const{
	if(_rhId != -999) 
		if(_rhId == jet.rhId()) return true;
	return ((_x == jet.x()) && (_y == jet.y()) && (_z == jet.z()) && (_t == jet.t()));
}

bool JetPoint::operator !=(JetPoint& jet) const{
	return !(*this == jet);
}






