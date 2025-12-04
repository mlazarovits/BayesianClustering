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

	_skipTime = false;
}




JetPoint::JetPoint(double x, double y, double z, double t){
	_x = x;
	_y = y;
	_z = z;
	_t = t;

	_set_rap_phi();
	_theta = atan2( sqrt(_x*_x + _y*_y) , z );	
	_E = 0;

	_idx = -999;
	_rhId = -999;
	_w = 1;

	_skipTime = false;
}



bool JetPoint::operator ==(JetPoint& jet) const{
	if(_rhId != -999) 
		if(_rhId == jet.rhId()) return true;
	return ((_x == jet.x()) && (_y == jet.y()) && (_z == jet.z()) && (_t == jet.t()));
}

bool JetPoint::operator !=(JetPoint& jet) const{
	return !(*this == jet);
}






