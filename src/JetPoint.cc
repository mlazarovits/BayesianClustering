#include "JetPoint.hh"


JetPoint::JetPoint(){
	_x = -999;
	_y = -999;
	_z = -999;
	_t = -999;

	_space = Point({_x, _y, _z, _t});

	_eta = _invalid_eta;
	_phi = _invalid_phi;
	
	_idx = -999;
	_rhId = -999;
}



JetPoint::JetPoint(double x, double y, double z, double t){
	_x = x;
	_y = y;
	_z = z;
	_t = t;

	_space = Point({_x, _y, _z, _t});
	
	_eta = _invalid_eta;
	_phi = _invalid_phi;
	
	_idx = -999;
	_rhId = -999;

}

JetPoint::JetPoint(Point pt){
	//check that Point is valid spacial four vector
	if(pt.Dim() != 4){
		cout << "Error: Point for JetPoint ctor must be of dimension 4. Dimension is " << pt.Dim() << endl;
		return;
	}

	_space = pt;
	_x = _space.at(0);
	_y = _space.at(1);
	_z = _space.at(2);
	_t = _space.at(3);
	
	_eta = _invalid_eta;
	_phi = _invalid_phi;

	_idx = -999;
	_rhId = -999;
}



JetPoint::~JetPoint(){
}


bool JetPoint::operator ==(JetPoint& jet) const{
	if(_rhId != -999) 
		if(_rhId == jet.rhId()) return true;
	return _space == jet.four_space();
}

bool JetPoint::operator !=(JetPoint& jet) const{
	return !(*this == jet);
}


void JetPoint::SetFourSpace(Point pt){
	if(pt.Dim() != 4){
		cout << "Error: spatial four vector for JetPoint must have dimension 4. Dimension is " << pt.Dim() << endl;
		return;
	}
	_space = pt;
}




