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
	
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
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
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
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
	
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
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
	
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
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
	
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
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

	_phi = params["mean"].at(1,0);
	_mu = params["mean"];
	//put phi on 02pi
	//if pt is negative
	double math_pi = acos(-1);
	if(_phi < 0){
		_phi = _phi + 2*math_pi;
	}
	//if pt is geq 2*pi
	else if(_phi >= 2*math_pi) _phi = _phi - 2*math_pi;
	_mu.SetEntry(_phi,1,0);
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


	//set subcluster (ie constituent) parameters + 4vectors
	double pxt = 0;
	double pyt = 0;
	double pzt = 0;
	double Et = 0;
	double mt = 0;
	_mu = Matrix(3,1);
	_cov = Matrix(3,3);	
	_eta = 0;
	_phi = 0;
	_t = 0;

	PointCollection means;
	for(int k = 0; k < nsubcl; k++){
		auto params = model->GetLikelihoodParameters(k);
		Ek = norms[k]/gev;
		Jet subcl(model->GetModel(k), Ek, model->GetPi(k), _vtx);
		_constituents.push_back(subcl);
		
		//set momentum from subclusters
		_px += subcl.px();
		_py += subcl.py();
		_pz += subcl.pz();
	//	_E += subcl.E(); //set from rhs
	
		//_mu.add(params["mean"]);
		BayesPoint pt(3);
		pt.SetValue(subcl.eta(),0);	
		pt.SetValue(subcl.phi(),1);	
		pt.SetValue(subcl.time(),2);	
cout << "subcl #" << k << " center" << endl; pt.Print();
		means += pt;
		//_eta += subcl.eta();
		//_phi += subcl.phi();
		//_t += subcl.time();

		//cout << "subcl k " << k << " center " << subcl.eta() << " " << subcl.phi() << endl; params["mean"].Print();
	}
	BayesPoint mean = means.mean();
	mean.SetValue(means.CircularMean(1),1);
	_mu = Matrix(mean);
	_eta = mean.at(0);
	_phi = mean.at(1);
	//put phi on 02pi
	//if pt is negative
	double pi = acos(-1);
	if(_phi < 0){
		_phi = _phi + 2*pi;
	}
	//if pt is geq 2*pi
	else if(_phi >= 2*pi) _phi = _phi - 2*pi;
	_mu.SetEntry(_phi,1,0);
	_t = mean.at(2);

	double deta, dphi, dtime, eta_phi, eta_time, phi_time;
	for(int k = 0; k < nsubcl; k++){
		auto params = model->GetLikelihoodParameters(k);
		//do point-wise covariance with mean set by subclusters
		deta = params["mean"].at(0,0) - _mu.at(0,0);	
		dphi = params["mean"].at(1,0) - _mu.at(1,0);	
		dphi = acos(cos(dphi));
		dtime = params["mean"].at(2,0) - _mu.at(2,0);
	
		Matrix cov_entry = Matrix(3,3);
		cov_entry.SetEntry(deta*deta,0,0);
		cov_entry.SetEntry(deta*dphi,1,0);
		cov_entry.SetEntry(deta*dtime,2,0);
		cov_entry.SetEntry(dphi*deta,0,1);
		cov_entry.SetEntry(dphi*dphi,1,1);
		cov_entry.SetEntry(dtime*dphi,2,1);
		cov_entry.SetEntry(deta*dtime,0,2);
		cov_entry.SetEntry(dphi*dtime,1,2);
		cov_entry.SetEntry(dtime*dtime,2,2);
		cov_entry.SetEntry(dtime*dtime,2,2);
		

		_cov.add(cov_entry);	

	}
	_cov.mult(_cov,1/double(nsubcl));	
	//cout << "jet from subcls px " << pxt << " py " << pyt << " pz " << pzt << " E " << Et << " m2 " << (Et+pzt)*(Et-pzt)-(pxt*pxt + pyt*pyt) << endl;

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
	
	//	pt = _rhs[i].E()/cosh(_rhs[i].eta()); //consistent with mass = 0
	//	_px += pt*cos(_rhs[i].phi());
	//	_py += pt*sin(_rhs[i].phi());
	//	_pz += pt*sinh(_rhs[i].eta());
	//cout << "rh #" << i << " time " << t << " phi " << phi << " weight " << rh.w() << endl;	
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


	_kt2 = _px*_px + _py*_py;

	_idx = 999;
	_ensure_valid_rap_phi();
	_set_time();


	_update_mom();
	_mass = mass();

	//double kt2 = pxt*pxt + pyt*pyt;
	//double m2t = (Et+pzt)*(Et-pzt)-kt2; 
	//mt = m2t < 0.0 ? -sqrt(-m2t) : sqrt(m2t); 
	//cout << "jet from rh px " << _px << " py " << _py << " pz " << _pz << " E " << _E << " m2 " << m2() << " mass " << _mass << endl;
	//cout << "jet from GMM px " << pxt << " py " << pyt << " pz " << pzt << " E " << Et << " m2 " << m2t << " mass " << mt << "\n" << endl;

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

