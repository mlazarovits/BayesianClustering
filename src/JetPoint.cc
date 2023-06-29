#include "JetPoint.hh"


JetPoint::JetPoint(){
	m_t = -999;
	m_x = -999;
	m_y = -999;
	m_z = -999;

	m_vec = Point({m_x, m_y, m_z, m_t});

	m_eta = m_invalid_eta;
	m_phi = m_invalid_phi;
	
	m_idx = -999;
	m_rhId = -999;
}



JetPoint::JetPoint(double x, double y, double z, double t){
	m_t = t;
	m_x = x;
	m_y = y;
	m_z = z;

	m_vec = Point({m_x, m_y, m_z, m_t});
	
	m_eta = m_invalid_eta;
	m_phi = m_invalid_phi;
	
	m_idx = -999;
	m_rhId = -999;

}

JetPoint::JetPoint(Point pt){
	//check that Point is valid spacial four vector
	if(pt.Dim() != 4){
		cout << "Error: Point for JetPoint ctor must be of dimension 4. Dimension is " << pt.Dim() << endl;
		return;
	}

	m_vec = pt;
	m_t = m_vec.at(0);
	m_x = m_vec.at(1);
	m_y = m_vec.at(2);
	m_z = m_vec.at(3);
	
	m_eta = m_invalid_eta;
	m_phi = m_invalid_phi;

	m_idx = -999;
	m_rhId = -999;
}



JetPoint::~JetPoint(){
}


bool JetPoint::operator ==(JetPoint& jet) const{
	if(m_rhId != -999) 
		if(m_rhId == jet.rhId()) return true;
	return m_vec == jet.four_space();
}

bool JetPoint::operator !=(JetPoint& jet) const{
	return !(*this == jet);
}


void JetPoint::SetFourSpace(Point pt){
	if(pt.Dim() != 4){
		cout << "Error: spatial four vector for JetPoint must have dimension 4. Dimension is " << pt.Dim() << endl;
		return;
	}
	m_vec = pt;
}




