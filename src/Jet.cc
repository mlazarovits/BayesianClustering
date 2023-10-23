#include "Jet.hh"
//TODO: change default space coords to eta phi
#include <iostream>
using std::cerr;

Jet::Jet(){
	_E = 0;
	_px = 0;
	_py = 0;
	_pz = 0;
	_mass = 0;
	_phi = _invalid_phi;
	_eta = _invalid_eta;

	//will set from GMM	
	_t = 0;
	
	_nRHs = 0;
	
	_parent1 = nullptr;
	_parent2 = nullptr;
	_child = nullptr; 

}


Jet::Jet(double px, double py, double pz, double E){
	_E = E;
	_px = px;
	_py = py;
	_pz = pz;
	_kt2 = px*px + py*py;
	_mass = mass();
	_phi = _invalid_phi;
	_eta = _invalid_eta;

	//sets eta + phi
	_ensure_valid_rap_phi();

	_t = 0;
	
	_nRHs = 0;
	
	_parent1 = nullptr;
	_parent2 = nullptr;
	_child = nullptr; 


}

Jet::Jet(JetPoint rh, Point vtx){
	
	_rhs.push_back(rh);
	_nRHs = (int)_rhs.size();

	_E = rh.E();
	_eta = rh.eta();
	_phi = rh.phi();
	_t = rh.t();
	_mass = mass();
	
	if(vtx.Dim() != 3){
		cerr << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
		return;
	}
	_vtx = vtx;

	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(rh.x()*rh.x() + rh.y()*rh.y()), rh.z() );
	double pt = _E*sin(theta); //mass = 0
	
	_px = pt*cos(_phi);
	_py = pt*sin(_phi);
	_pz = pt*sinh(_eta);
	_kt2 = sqrt(pt); 		
	_parent1 = nullptr;
	_parent2 = nullptr;
	_child = nullptr; 

	_idx = 999;
	_set_time();
}


Jet::Jet(vector<JetPoint> rhs){
	_nRHs = (int)rhs.size();	
	for(int i = 0; i < _nRHs; i++) _rhs.push_back(rhs[i]);

	double theta, pt, x, y, z;
	for(int i = 0; i < _nRHs; i++){
		
		//theta is calculated between beamline (z-dir) and x-y vector	
		x = rhs[i].x();
		y = rhs[i].y();
		z = rhs[i].z();
		theta = atan2( sqrt(x*x + y*y), z );
		pt = rhs[i].E()*sin(theta); //consistent with mass = 0
		_px += pt*cos(rhs[i].phi());
		_py += pt*sin(rhs[i].phi());
		_pz += pt*cosh(rhs[i].eta());
		
		_E += rhs[i].E();


	}
	_kt2 = _px*_px + _py*_py;
	_ensure_valid_rap_phi();
	_set_time();
}

Jet::Jet(vector<Jet> jets){
	vector<JetPoint> rhs;
	for(int i = 0; i < (int)jets.size(); i++){
		jets[i].GetConstituents(rhs);
		for(int j = 0; j < (int)rhs.size(); j++)
		_rhs.push_back(rhs[j]);
	}
	_nRHs = (int)_rhs.size();	
	double theta, pt, x, y, z;
	for(int i = 0; i < _nRHs; i++){
		//theta is calculated between beamline (z-dir) and x-y vector	
		x = rhs[i].x();
		y = rhs[i].y();
		z = rhs[i].z();
		theta = atan2( sqrt(x*x + y*y), z );
		pt = rhs[i].E()*sin(theta); //consistent with mass = 0
		_px += pt*cos(rhs[i].phi());
		_py += pt*sin(rhs[i].phi());
		_pz += pt*cosh(rhs[i].eta());
		
		_E += rhs[i].E();


	}
	_kt2 = _px*_px + _py*_py;
	_ensure_valid_rap_phi();
	_set_time();
}

Jet::Jet(vector<JetPoint> rhs, Point vtx){
	_nRHs = (int)rhs.size();	
	for(int i = 0; i < _nRHs; i++) _rhs.push_back(rhs[i]);

	if(vtx.Dim() != 3){
		cerr << "Error: must provide 3 dimensional spacial coordinates for vertex for momentum direction." << endl;
		return;
	}
	_vtx = vtx;

	double theta, pt, x, y, z;
	for(int i = 0; i < _nRHs; i++){
		//theta is calculated between beamline (z-dir) and x-y vector	
		x = rhs[i].x();
		y = rhs[i].y();
		z = rhs[i].z();
		theta = atan2( sqrt(x*x + y*y), z );
		pt = rhs[i].E()*sin(theta); //consistent with mass = 0
		_px += pt*cos(rhs[i].phi());
		_py += pt*sin(rhs[i].phi());
		_pz += pt*cosh(rhs[i].eta());
		
		_E += rhs[i].E();


	}
	_kt2 = _px*_px + _py*_py;
	_ensure_valid_rap_phi();
	_set_time();
}


Jet::Jet(const Jet& j){
	_px = j.px();
	_py = j.py();
	_pz = j.pz();
	_E = j.E();
	_mom = Point({_px, _py, _pz, _E});
	_eta = j.eta();
	_phi = j.phi();
	_t = j.time();

	_kt2 = j.kt2();
	_mass = j.mass();
	
	_child = j._child;
	_parent1 = j._parent1;
	_parent2 = j._parent2;
	_idx = j.GetUserIdx();
	_vtx = j.GetVertex();
	_rhs = j.GetJetPoints();
	_nRHs = (int)_rhs.size();

}

Jet::~Jet(){
	_rhs.clear();
}


bool Jet::operator ==(Jet& jet) const{
	return _mom == jet.four_mom();
}

bool Jet::operator !=(Jet& jet) const{
	return _mom != jet.four_mom();
}

//add jet jt to this
//adding four vectors - recalculate invariant mass and other kinematic quantities
void Jet::add(const Jet& jt){
	//add rhs from jt
	vector<JetPoint> rhs;
	jt.GetConstituents(rhs);

	for(int i = 0; i < (int)rhs.size(); i++) _rhs.push_back(rhs[i]);
	_nRHs += (int)rhs.size();
	
	//add momentum
	_px += jt.px();
	_py += jt.py();
	_pz += jt.pz();
	_E  += jt.E();
	//set time to be energy-weighted average of rec hit times
	_set_time();

	//recalculate kt2 of cluster
	_kt2 = _px*_px + _py*_py;
	//recalculate eta and phi of cluster
	_set_rap_phi();
}


void Jet::add(const JetPoint& rh){
	_rhs.push_back(rh);
	_nRHs += 1;

	if(_vtx.Dim() != 3){
		cerr << "Error: must provide 3 dimensional spatial coordinates for vertex for momentum direction." << endl;
		return;
	}
	//direction between ECAL position and vertex:q

	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(rh.x()*rh.x() + rh.y()*rh.y()), rh.z() );
	double pt = rh.E()*sin(theta); //consistent with mass = 0

	//add momentum
	//rh is massless => pT = Esin(theta)
	_px += pt*cos(rh.phi());
	_py += pt*sin(rh.phi());
	_pz += pt*sinh(rh.eta());
	_E  += rh.E();
	//set time to be energy-weighted average of rec hit times
	_set_time();
	//recalculate kt2 of cluster
	_kt2 = sqrt(pt);
	//recalculate eta and phi of cluster
	_set_rap_phi();


}

void Jet::GetParents(Jet& p1, Jet& p2) const{
	//check if jet has parents that combined to make this
	//if no parents found, fills with nullptr
	p1 = *_parent1;
	p2 = *_parent2;
}


Jet Jet::GetBaby() const{
	return *_child;
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

