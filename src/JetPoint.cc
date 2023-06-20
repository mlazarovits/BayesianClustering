#include "JetPoint.hh"


JetPoint::JetPoint(){
	m_E = 0;
	m_px = 0;
	m_py = 0;
	m_pz = 0;

	m_t = 0;
	m_x = 0;
	m_y = 0;
	m_z = 0;

	m_eta = -999;
	m_phi = -999;


}

JetPoint::JetPoint(double E, double px, double py, double pz){
	m_E = E;
	m_px = px;
	m_py = py;
	m_pz = pz;

	m_kt2 = m_px*m_px + m_py*m_py;
	
	if(m_pz == 0){
		m_theta = 0;
		m_eta = m_maxRap;
	}
	else{
		m_theta = atan(m_kt2 / pz);
		if(m_theta < 0) m_theta += acos(-1);
		m_eta = -log(tan(m_theta/2.));
	}
	m_momvec = Point({E, px, py, pz});
	
	m_t = 0;
	m_x = 0;
	m_y = 0;
	m_z = 0;
	
	m_spacevec = Point({t, x, y, z});

	m_eta = -999;
	m_phi = -999;
}


JetPoint::JetPoint(Point pt, bool mom){
	if(mom){
		m_momvec = pt;
		m_E = m_momvec.at(0);
		m_px = m_momvec.at(1);
		m_py = m_momvec.at(2);
		m_pz = m_momvec.at(3);
		m_kt2 = m_px*m_px + m_py*m_py;
	}
	else{
		m_spacevec = pt;
		m_t = m_spacevec.at(0);
		m_x = m_spacevec.at(1);
		m_y = m_spacevec.at(2);
		m_z = m_spacevec.at(3);
	}

}

JetPoint::JetPoint(double t, double x, double y, double z){
	m_E = 0;
	m_px = 0;
	m_py = 0;
	m_pz = 0;
	
	m_momvec = Point({E, px, py, pz});
	m_t = t;
	m_x = x;
	m_y = y;
	m_z = z;
	m_spacevec = Point({t, x, y, z});

	m_eta = 0;
	m_phi = 0;
}
