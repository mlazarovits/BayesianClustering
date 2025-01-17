#include "Jet.hh"
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

	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	//will set from GMM	
	_t = 0;
	
	_nRHs = 0;
	
	_cov = Matrix();
	_mu = Matrix();
	_pi = 0;
	_update_mom();
}


Jet::Jet(double px, double py, double pz, double E){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
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
	_cov = Matrix();
	_mu = Matrix();
	_pi = 0;
	
	_update_mom();

}

Jet::Jet(JetPoint rh, BayesPoint vtx){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	_rhs.push_back(rh);
	_nRHs = (int)_rhs.size();
	
	_vtx = vtx;

	_E = rh.E();
	_eta = rh.eta();
	_phi = rh.phi();
	_t = rh.t();

	//theta is calculated between beamline (z-dir) and x-y vector	
	//centered at (0,0,0)
	double theta = atan2( sqrt(rh.x()*rh.x() + rh.y()*rh.y()), rh.z() );
	//double pt = _E*sin(theta); //mass = 0
	double pt = _E/cosh(_eta);

	_px = pt*cos(_phi);
	_py = pt*sin(_phi);
	_pz = pt*sinh(_eta);
	_kt2 = _px*_px + _py*_py; 		
	_mass = mass();

	
	_idx = 999;
	_set_time();
	_ensure_valid_rap_phi();
	
	_cov = Matrix();
	_mu = Matrix();
	_pi = 0;
	_update_mom();
}


Jet::Jet(const vector<JetPoint>& rhs, BayesPoint vtx){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	for(int i = 0; i < (int)rhs.size(); i++) _rhs.push_back(rhs[i]);
	_nRHs = (int)_rhs.size();	
	double theta, pt, x, y, z;
	_phi = _invalid_phi;
	_eta = _invalid_eta;
	
	_E = 0;
	_px = 0;
	_py = 0;
	_pz = 0;

	_eta = 0;
	_phi = 0;

	_vtx = vtx;
	double phi, eta;
	double norm = 0;
	for(int i = 0; i < _nRHs; i++){		
		//theta is calculated between beamline (z-dir) and vector in x-y plane	
		//centered at (0,0,0)
		x = rhs[i].x();// - _vtx.at(0);
		y = rhs[i].y();// - _vtx.at(1);
		z = rhs[i].z();// - _vtx.at(2);
		theta = atan2( sqrt(x*x + y*y), z );
		phi = atan2(y, x);
		eta = -log(tan(theta/2.)); 
		//see https://cmssdt.cern.ch/lxr/source/DataFormats/CaloTowers/src/CaloTower.cc L145
		//pt = _E*sin(theta); //mass = 0, equivalent to below
		pt = rhs[i].E()/cosh(eta);
		_px += pt*cos(phi);
		_py += pt*sin(phi);
		_pz += pt*sinh(eta);
	//cout << "i " << i << " px " << pt*cos(phi) << " py " << pt*sin(phi) << " pz " << pt*cosh(eta) << " " << pt*sinh(eta) << endl;
		
		_E += rhs[i].E();

		_eta += rhs[i].eta();
		//if(rhs[i].phi() < 0) 
		//	_phi += rhs[i].phi()+2*acos(-1);
		//else
		_phi += rhs[i].phi();
		norm += 1.;//rhs[i].GetWeight();
	}
	
	_eta /= norm;
	_phi /= norm;
	
	_kt2 = _px*_px + _py*_py;
	//wraparound
	//_phi = acos(cos(_phi));
	_mass = mass();
//cout << "MAKING jet pts kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << sqrt(_kt2) << " mass " << _mass << " n pts " << rhs.size() << endl;

	_idx = 999;
	_ensure_valid_rap_phi();
	_set_time();
	
	_cov = Matrix();
	_mu = Matrix();
	_pi = 0;
	_update_mom();

}

Jet::Jet(const vector<Jet>& jets){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	for(int i = 0; i < (int)jets.size(); i++){
		vector<JetPoint> rhs = jets[i].GetJetPoints();
		for(int j = 0; j < (int)rhs.size(); j++)
		_rhs.push_back(rhs[j]);
	}
	_nRHs = (int)_rhs.size();	
	double pt;
	_E = 0;
	_px = 0;
	_py = 0;
	_pz = 0;
	for(int i = 0; i < _nRHs; i++){
		pt = _rhs[i].E()*cosh(_rhs[i].eta()); //consistent with mass = 0
		_px += pt*cos(_rhs[i].phi());
		_py += pt*sin(_rhs[i].phi());
		_pz += pt*sinh(_rhs[i].eta());
		
		_E += _rhs[i].E();


	}
	_kt2 = _px*_px + _py*_py;

	_idx = 999;
	_ensure_valid_rap_phi();
	_set_time();
	
	_cov = Matrix();
	_mu = Matrix();
	_pi = 0;
	_update_mom();
}

//_pi = 1 ==> 1 subcluster in jet, _pi != 1 ==> >1 subcluster in jet
//no weighting yet because if this is the only component, weight (ie pi) should be 1
Jet::Jet(const Matrix& mu, const Matrix& cov, double E, double pi, BayesPoint vtx){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	_E = E;
	_eta =  mu.at(0,0);
	_phi =  mu.at(1,0);
	_t = mu.at(2,0);
	//0 mass hypothesis default sets four vector
	double pt = E/cosh(_eta);
	_px = pt*cos(_phi);
	_py = pt*sin(_phi);
	_pz = pt*sinh(_eta);

	_kt2 = _px*_px + _py*_py;
	_mass = mass();
//cout << "jet subcl kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << pt << " m2 " << m2() << endl;

	_idx = 999;
	_ensure_valid_rap_phi();

	_mu = mu;
	_cov = cov;
	_pi = pi;

	_vtx = vtx;

	_update_mom();
}

Jet::Jet(BasePDF* pdf, double E, double pi, BayesPoint vtx){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	_E = E;

	map<string, Matrix> params = pdf->GetParameters();
	_eta = params["mean"].at(0,0);
	_phi = params["mean"].at(1,0);
	_t =   params["mean"].at(2,0);
	//0 mass hypothesis default sets four vector
	double pt = E/cosh(_eta);
	_px = pt*cos(_phi);
	_py = pt*sin(_phi);
	_pz = pt*sinh(_eta);

	_kt2 = _px*_px + _py*_py;
	_mass = mass();
//cout << "jet subcl kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << pt << " m2 " << m2() << endl;

	_idx = 999;
	_ensure_valid_rap_phi();

	_mu = params["mean"];
	_cov = params["cov"];
	_pi = pi;

	_vtx = vtx;

	_update_mom();
}


Jet::Jet(BasePDFMixture* model, BayesPoint vtx, double gev, double detR = 129){
	_vtx = vtx;
	_mom = BayesPoint(4);
	
	Matrix r_nk = model->GetPosterior();
	double Ek;
	vector<double> norms;
	model->GetNorms(norms);
	int nsubcl = model->GetNClusters();
	_nRHs = model->GetData()->GetNPoints();
	double pt = 0;
	_E = 0;
	_px = 0;
	_py = 0;
	_pz = 0;
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
	_pi = 0;
	double x, y, z, t, eta, phi, theta, w;
	map<string, Matrix> params;
	//set mean
	model->GetJetParameters(params);
	_mu = params["mean"];
	_cov = params["cov"];	
	_eta = _mu.at(0,0);
	_phi = _mu.at(1,0);
	_t = _mu.at(2,0);


	for(int i = 0; i < _nRHs; i++){
		//add rhs to jet
		BayesPoint rh = model->GetData()->at(i);
		//transform eta, phi to x, y, z
		eta = rh.at(0);
		phi = rh.at(1);
		t = rh.at(2);
		x = detR*cos(phi);
		y = detR*sin(phi);
		theta = 2*atan2(1,exp(eta));
		z = detR/tan(theta);
		//push back as jet point to _rhs	
		_rhs.push_back(JetPoint(x,y,z,t));
		_rhs[i].SetEnergy(rh.w()/gev);	
	
		pt = _rhs[i].E()/cosh(_rhs[i].eta()); //consistent with mass = 0
		_px += pt*cos(_rhs[i].phi());
		_py += pt*sin(_rhs[i].phi());
		_pz += pt*sinh(_rhs[i].eta());
		//cout << "total pts r " << i << " pt " << pt << " E " << _rhs[i].E() << " px " << pt*cos(_rhs[i].phi()) << " py " << pt*sin(_rhs[i].phi()) << " pz " << pt*sinh(_rhs[i].eta()) << " eta " << _rhs[i].eta() << " phi " << _rhs[i].phi() << endl;
		
		_E += _rhs[i].E();

		//weights are sum of weights over all subclusters
		//without additional weights, this w_rh = sum_k r_rh,k = 1
		//if subcluster is downweighted this w_rh != 1
		w = 0;
		for(int k = 0; k < nsubcl; k++)
			w += r_nk.at(i,k);	 
		//set probabilistic coefficient
		//pi = sum_rh sum_k r_nk
		_pi += w;	
	}	
	//set subcluster (ie constituent) parameters + 4vectors
	double pxt = 0;
	double pyt = 0;
	double pzt = 0;
	double Et = 0;
	double mt = 0;	
	for(int k = 0; k < nsubcl; k++){
		params = model->GetPriorParameters(k);
		Ek = norms[k]/gev;
		Jet subcl(model->GetModel(k), Ek, model->GetPi(k), _vtx);
		//subcluster momentom 3vector is r_nk weighted average over rechits for cluster k
		double pxk = 0;
		double pyk = 0;
		double pzk = 0;
		for(int r = 0; r < _nRHs; r++){
			Jet rh(_rhs[r], _vtx);
			pxk += r_nk.at(r,k)*rh.px();	
			pyk += r_nk.at(r,k)*rh.py();	
			pzk += r_nk.at(r,k)*rh.pz();	
			//cout << "subcl pts k " << k << " r " << r << " pt " << rh.pt() << " E " << rh.E() << " px " << rh.px() << " py " << rh.py() << " pz " << rh.pz() << " eta " << rh.eta() << " phi " << rh.phi() << endl;
		}

		subcl.SetP(pxk, pyk, pzk);
	//cout << "subcluster k " << k << " px " << subcl.px() << " py " << subcl.py() << " pz " << subcl.pz() << " E " << subcl.E() << " m2 " << subcl.m2() << " mass " << subcl.mass() << endl;
		pxt += pxk;
		pyt += pyk;
		pzt += pzk;
		Et += Ek;
		_constituents.push_back(subcl);
	}
	//cout << "jet from subcls px " << pxt << " py " << pyt << " pz " << pzt << " E " << Et << " m2 " << (Et+pzt)*(Et-pzt)-(pxt*pxt + pyt*pyt) << endl;
	//set momentum from subclusters
	_px = pxt;
	_py = pyt;
	_pz = pzt;


	_kt2 = _px*_px + _py*_py;

	_idx = 999;
	_ensure_valid_rap_phi();
	_set_time();


	_update_mom();
//cout << "jet from GMM px " << _px << " py " << _py << " pz " << _pz << " m2 " << m2() << " mass " << mass() << " E " << _E << endl;

}

Jet::Jet(const Jet& j){
	_px = j.px();
	_py = j.py();
	_pz = j.pz();
	_E = j.E();
	_mom = BayesPoint({_px, _py, _pz, _E});
	_eta = j.eta();
	_phi = j.phi();
	_t = j.time();

	_kt2 = j.kt2();
	_mass = j.mass();
	
	_idx = j.GetUserIdx();
	_vtx = j.GetVertex();
	_rhs = j.GetJetPoints();
	_nRHs = (int)_rhs.size();

	_constituents = j._constituents;
	_mu = j._mu;
	_cov = j._cov;
	_pi = j._pi;
}

Jet::~Jet(){
	_rhs.clear();
}


bool Jet::operator ==(const Jet& jet) const{
	return _mom == jet.four_mom();
}

bool Jet::operator !=(const Jet& jet) const{
	return _mom != jet.four_mom();
}

//add jet jt to this
//adding four vectors - recalculate invariant mass and other kinematic quantities
void Jet::add(const Jet& jt){
	//add rhs from jt
	vector<JetPoint> rhs = jt.GetJetPoints();

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

	//add subclusters (aka constituents)
	for(auto j : jt._constituents)
		_constituents.push_back(j);
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

