#include "Jet.hh"
//TODO: change default space coords to eta phi

Jet::Jet(){
	_E = 0;
	_px = 0;
	_py = 0;
	_pz = 0;
	
	_t = 0;
	_x = 0;
	_y = 0;
	_z = 0;
	
	_nRHs = 0;
	
	_parent1 = nullptr;
	_parent2 = nullptr;
	_child = nullptr; 

}


Jet::Jet(JetPoint rh, Point vtx){
	_t = rh.x(0);
	_x = rh.x(1);
	_y = rh.x(2);
	_z = rh.x(3);
	
	_rhs.push_back(rh);
	_nRHs = (int)_rhs.size();
	_space = rh.four_space();

	_E = rh.E();
	_eta = rh.eta();
	_phi = rh.eta();
	
	if(vtx.Dim() != 3){
		cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
		return;
	}
	_vtx = vtx;
	Point dir = Point({_x - vtx.at(0), _y - vtx.at(1), _z - vtx.at(3)});

	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(_x*_x + _y*_y), _z );
	SetPt();
	SetFourMom(Point({pt()*cos(_phi), pt()*sin(_phi), pt()*cosh(_eta),_E}));

		
	_parent1 = nullptr;
	_parent2 = nullptr;
	_child = nullptr; 

	_idx = 999;
}


Jet::Jet(vector<JetPoint> rhs){
	_nRHs = (int)rhs.size();	
	for(int i = 0; i < _nRHs; i++) _rhs.push_back(rhs[i]);

	double theta, pt, x, y, z;
	for(int i = 0; i < _nRHs; i++){
		
		//theta is calculated between beamline (z-dir) and x-y vector	
		x = rhs[i].x(0);
		y = rhs[i].x(1);
		z = rhs[i].x(2);
		theta = atan2( sqrt(x*x + y*y), z );
		pt = rhs[i].E()*sin(theta); //consistent with mass = 0
		_px += pt*cos(rhs[i].phi());
		_py += pt*sin(rhs[i].phi());
		_pz += pt*cosh(rhs[i].eta());
		
		_E += rhs[i].E();


	}
	RecalcKT2();
	RecalcPhi();
	RecalcEta();
//TODO: need to recalculate space four vector from clustering algo
/*	
	_t = rh.at(0);
	_x = rh.at(1);
	_y = rh.at(2);
	_z = rh.at(3);
	m_s = rh.four_space();
*/

}

Jet::Jet(vector<JetPoint> rhs, Point vtx){
	_nRHs = (int)rhs.size();	
	for(int i = 0; i < _nRHs; i++) _rhs.push_back(rhs[i]);

	if(vtx.Dim() != 3){
		cout << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
		return;
	}
	_vtx = vtx;

	double theta, pt, x, y, z;
	for(int i = 0; i < _nRHs; i++){
		Point dir = Point({_x - vtx.at(0), _y - vtx.at(1), _z - vtx.at(3)});
		
		//theta is calculated between beamline (z-dir) and x-y vector	
		x = rhs[i].x(0);
		y = rhs[i].x(1);
		z = rhs[i].x(2);
		theta = atan2( sqrt(x*x + y*y), z );
		pt = rhs[i].E()*sin(theta); //consistent with mass = 0
		_px += pt*cos(rhs[i].phi());
		_py += pt*sin(rhs[i].phi());
		_pz += pt*cosh(rhs[i].eta());
		
		_E += rhs[i].E();


	}
	RecalcKT2();
	RecalcPhi();
	RecalcEta();
//TODO: need to recalculate space four vector from clustering algo
/*	
	_t = rh.at(0);
	_x = rh.at(1);
	_y = rh.at(2);
	_z = rh.at(3);
	m_s = rh.four_space();
*/

}


Jet::~Jet(){
	_rhs.clear();
}


bool Jet::operator ==(Jet& jet) const{
	return _space == jet.four_space();
}

bool Jet::operator !=(Jet& jet) const{
	return _space != jet.four_space();
}

//TODO: make sure it's consistent with the rec hits (JetPoints)
void Jet::SetFourMom(Point pt){
	if(pt.Dim() != 4){
		cout << "Error: must provide 4 vector for momentum." << endl;
		return;
	}
	
	_mom = pt;
}

//TODO: make sure it's consistent with the rec hits (JetPoints)
void Jet::SetFourPos(Point pt){
	if(pt.Dim() != 4){
		cout << "Error: must provide 4 vector for space." << endl;
		return;
	}
	_space = pt;
}

//add jet jt to this
//adding four vectors - recalculate invariant mass and other kinematic quantities
void Jet::add(Jet& jt){
	//add rhs from jt
	vector<JetPoint> rhs;
	jt.GetConstituents(rhs);

	for(int i = 0; i < (int)rhs.size(); i++) _rhs.push_back(rhs[i]);
	_nRHs += (int)rhs.size();
	
	//add momentum
	_px += jt.p(0);
	_py += jt.p(1);
	_pz += jt.p(2);
	_E  += jt.p(3);

	//recalculate kt2 of cluster
	RecalcKT2();
	//recalculate phi of cluster
	RecalcPhi();
	//recalculate pseudorap of cluster - from FastJet
	RecalcEta();

//TODO: need to recalculate space four vector from clustering algo
/*
	_x = 0;
	_y = 0;
	_z = 0;
	_t = 0;


	for(int i = 0; i < _nRHs; i++){
		_x += _space_vec_coll.at(i).at(0);
		_y += _space_vec_coll.at(i).at(1);
		_z += _space_vec_coll.at(i).at(2);
		_t += _space_vec_coll.at(i).at(3);

	}
	_x /= _nRHs;
	_y /= _nRHs;
	_z /= _nRHs;
	_t /= _nRHs;
*/
}


void Jet::add(JetPoint& rh){
	_rhs.push_back(rh);
	_nRHs += 1;

	if(_vtx.Dim() != 3){
		cout << "Error: must provide 3 dimensional spatial coordinates for vertex for momentum direction." << endl;
		return;
	}
	//direction between ECAL position and vertex:q

	Point dir = Point({_x - _vtx.at(0), _y - _vtx.at(1), _z - _vtx.at(3)});

	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(_x*_x + _y*_y), _z );
	double pt = rh.E()*sin(theta); //consistent with mass = 0

	//add momentum
	//rh is massless => pT = Esin(theta)
	_px += pt*cos(rh.phi());
	_py += pt*sin(rh.phi());
	_pz += pt*sinh(rh.eta());
	_E  += rh.E();

	//recalculate kt2
	RecalcKT2();
	RecalcPhi();
	//recalculate pseudorap
	RecalcEta();
//TODO: need to recalculate space four vector from clustering algo
/*
	_x = 0;
	_y = 0;
	_z = 0;
	_t = 0;


	for(int i = 0; i < _nRHs; i++){
		_x += _space_vec_coll.at(i).at(0);
		_y += _space_vec_coll.at(i).at(1);
		_z += _space_vec_coll.at(i).at(2);
		_t += _space_vec_coll.at(i).at(3);

	}
	_x /= _nRHs;
	_y /= _nRHs;
	_z /= _nRHs;
	_t /= _nRHs;
*/





}

void Jet::GetParents(Jet& p1, Jet& p2) const{
	//check if jet has parents that combined to make this
	//if no parents found, fills with nullptr
	p1 = *_parent1;
	p2 = *_parent2;
}


void Jet::GetBaby(Jet& child) const{
	child = *_child;
}


//TODO: set after clustering
void Jet::GetSubJets(vector<Jet>& subjets, int depth) const{


}


//set after clustering
bool Jet::Has(Jet& jet) const{
	return true;

}

bool Jet::Has(JetPoint& rh) const{
	for(int i = 0; i < _nRHs; i++)
		if(_rhs[i] == rh)
			return true;
	return false;

}

