#include "Jet.hh"


Jet::Jet(){
	m_E = 0;
	m_px = 0;
	m_py = 0;
	m_pz = 0;

	m_t = 0;
	m_x = 0;
	m_y = 0;
	m_z = 0;

	m_eta = m_invalid_eta;
	m_phi = m_invalid_phi;
	
	m_parent1 = nullptr;
	m_parent2 = nullptr;
	m_child = nullptr; 
}

//Jet::Jet(double E, double px, double py, double pz){
//	m_E = E;
//	m_px = px;
//	m_py = py;
//	m_pz = pz;
//
//	m_kt2 = m_px*m_px + m_py*m_py;
//	
//	if(m_pz == 0){
//		m_theta = 0;
//		m_eta = m_maxRap;
//	}
//	else{
//		m_theta = atan(m_kt2 / pz);
//		if(m_theta < 0) m_theta += acos(-1);
//		m_eta = -log(tan(m_theta/2.));
//	}
//	m_momvec = Point({m_E, m_px, m_py, m_pz});
//	
//	m_t = 0;
//	m_x = 0;
//	m_y = 0;
//	m_z = 0;
//	
//	m_spacevec = Point({m_t, m_x, m_y, m_z});
//
//	m_eta = -999;
//	m_phi = -999;
//	m_parent1 = nullptr;
//	m_parent2 = nullptr;
//	m_child = nullptr; 
//}


Jet::Jet(Point pt, bool mom){
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
	m_parent1 = nullptr;
	m_parent2 = nullptr;
	m_child = nullptr; 

}



Jet::~Jet(){
}

void Jet::SetFourMom(Point pt){
	m_momvec = pt;
}

void Jet::SetFourPos(Point pt){
	m_spacevec = pt;
}



void Jet::GetParents(Jet& p1, Jet& p2){
	//check if jet has parents that combined to make this
	//if no parents found, fills with nullptr

	p1 = *m_parent1;
	p2 = *m_parent2;
}


void Jet::GetBaby(Jet& child){
	child = *m_child;
}

//set after clustering
void Jet::GetSubjets(vector<Jet>& subjets, int depth){


}


//set after clustering
bool Jet::Has(Jet& jet){
	return true;

}
