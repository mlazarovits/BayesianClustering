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
	
	m_nRHs = 0;
	
	m_parent1 = nullptr;
	m_parent2 = nullptr;
	m_child = nullptr; 

}


Jet::Jet(JetPoint rh, Point vtx){
	m_t = rh.x(0);
	m_x = rh.x(1);
	m_y = rh.x(2);
	m_z = rh.x(3);
	
	m_rhs.push_back(rh);
	m_nRHs = (int)m_rhs.size();
	m_space = rh.four_space();

	m_E = rh.E();
	m_eta = rh.eta();
	m_phi = rh.eta();
	
	if(vtx.Dim() != 3){
		cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
		return;
	}
	m_vtx = vtx;
	Point dir = Point({m_x - vtx.at(0), m_y - vtx.at(1), m_z - vtx.at(3)});

	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(m_x*m_x + m_y*m_y), m_z );
	SetPt(m_E*sin(theta)); //consistent with mass = 0
	SetFourMom(Point({pt()*cos(m_phi), pt()*sin(m_phi), pt()*cosh(m_eta),m_E}));

		
	m_parent1 = nullptr;
	m_parent2 = nullptr;
	m_child = nullptr; 

	m_idx = 999;
}



Jet::Jet(vector<JetPoint> rhs, Point vtx){
	m_nRHs = (int)rhs.size();	
	for(int i = 0; i < m_nRHs; i++) m_rhs.push_back(rhs[i]);

	if(vtx.Dim() != 3){
		cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
		return;
	}
	m_vtx = vtx;

	double theta, pt, x, y, z;
	for(int i = 0; i < m_nRHs; i++){
		Point dir = Point({m_x - vtx.at(0), m_y - vtx.at(1), m_z - vtx.at(3)});
		
		//theta is calculated between beamline (z-dir) and x-y vector	
		x = rhs[i].x(0);
		y = rhs[i].x(1);
		z = rhs[i].x(2);
		theta = atan2( sqrt(x*x + y*y), z );
		pt = rhs[i].E()*sin(theta); //consistent with mass = 0
		m_px += pt*cos(rhs[i].phi());
		m_py += pt*sin(rhs[i].phi());
		m_pz += pt*cosh(rhs[i].eta());
		
		m_E += rhs[i].E();


	}
	RecalcKT2();
	RecalcPhi();
	RecalcEta();
//TODO: need to recalculate space four vector from clustering algo
/*	
	m_t = rh.at(0);
	m_x = rh.at(1);
	m_y = rh.at(2);
	m_z = rh.at(3);
	m_s = rh.four_space();
*/

}


Jet::~Jet(){
	m_rhs.clear();
}


bool Jet::operator ==(Jet& jet) const{
	return m_space == jet.four_space();
}

bool Jet::operator !=(Jet& jet) const{
	return m_space != jet.four_space();
}

//TODO: make sure it's consistent with the rec hits (JetPoints)
void Jet::SetFourMom(Point pt){
	if(pt.Dim() != 4){
		cout << "Error: must provide 4 vector for momentum." << endl;
		return;
	}
	
	m_p = pt;
}

//TODO: make sure it's consistent with the rec hits (JetPoints)
void Jet::SetFourPos(Point pt){
	if(pt.Dim() != 4){
		cout << "Error: must provide 4 vector for space." << endl;
		return;
	}
	m_space = pt;
}

//add jet jt to this
//adding four vectors - recalculate invariant mass and other kinematic quantities
void Jet::add(Jet& jt){
	//add rhs from jt
	vector<JetPoint> rhs;
	jt.GetConstituents(rhs);

	for(int i = 0; i < (int)rhs.size(); i++) m_rhs.push_back(rhs[i]);
	m_nRHs += (int)rhs.size();
	
	//add momentum
	m_px += jt.p(0);
	m_py += jt.p(1);
	m_pz += jt.p(2);
	m_E  += jt.p(3);

	//recalculate kt2 of cluster
	RecalcKT2();
	//recalculate phi of cluster
	RecalcPhi();
	//recalculate pseudorap of cluster - from FastJet
	RecalcEta();

//TODO: need to recalculate space four vector from clustering algo
/*
	m_x = 0;
	m_y = 0;
	m_z = 0;
	m_t = 0;


	for(int i = 0; i < m_nRHs; i++){
		m_x += m_space_vec_coll.at(i).at(0);
		m_y += m_space_vec_coll.at(i).at(1);
		m_z += m_space_vec_coll.at(i).at(2);
		m_t += m_space_vec_coll.at(i).at(3);

	}
	m_x /= m_nRHs;
	m_y /= m_nRHs;
	m_z /= m_nRHs;
	m_t /= m_nRHs;
*/
}


void Jet::add(JetPoint& rh){
	m_rhs.push_back(rh);
	m_nRHs += 1;

	if(m_vtx.Dim() != 3){
		cout << "Error: must provide 3 dimensional spatial coordinates for vertex for momentum direction." << endl;
		return;
	}
	//direction between ECAL position and vertex:q

	Point dir = Point({m_x - m_vtx.at(0), m_y - m_vtx.at(1), m_z - m_vtx.at(3)});

	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(m_x*m_x + m_y*m_y), m_z );
	double pt = rh.E()*sin(theta); //consistent with mass = 0

	//add momentum
	//rh is massless => pT = Esin(theta)
	m_px += pt*cos(rh.phi());
	m_py += pt*sin(rh.phi());
	m_pz += pt*sinh(rh.eta());
	m_E  += rh.E();

	//recalculate kt2
	RecalcKT2();
	RecalcPhi();
	//recalculate pseudorap
	RecalcEta();
//TODO: need to recalculate space four vector from clustering algo
/*
	m_x = 0;
	m_y = 0;
	m_z = 0;
	m_t = 0;


	for(int i = 0; i < m_nRHs; i++){
		m_x += m_space_vec_coll.at(i).at(0);
		m_y += m_space_vec_coll.at(i).at(1);
		m_z += m_space_vec_coll.at(i).at(2);
		m_t += m_space_vec_coll.at(i).at(3);

	}
	m_x /= m_nRHs;
	m_y /= m_nRHs;
	m_z /= m_nRHs;
	m_t /= m_nRHs;
*/





}

void Jet::GetParents(Jet& p1, Jet& p2) const{
	//check if jet has parents that combined to make this
	//if no parents found, fills with nullptr
	p1 = *m_parent1;
	p2 = *m_parent2;
}


void Jet::GetBaby(Jet& child) const{
	child = *m_child;
}


//TODO: set after clustering
void Jet::GetSubJets(vector<Jet>& subjets, int depth) const{


}


//set after clustering
bool Jet::Has(Jet& jet) const{
	return true;

}

bool Jet::Has(JetPoint& rh) const{
	for(int i = 0; i < m_nRHs; i++)
		if(m_rhs[i] == rh)
			return true;
	return false;

}

