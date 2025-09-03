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
	_idx = 999;	
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
	_mass = _calc_mass();
	_phi = _invalid_phi;
	_eta = _invalid_eta;

	//sets eta + phi from p vector
	_ensure_valid_rap_phi();

	_t = 0;
	
	_nRHs = 0;
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
	_pi = 0;
	
	_idx = 999;	
	_update_mom();

}

Jet::Jet(JetPoint rh, BayesPoint vtx){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	_rhs.push_back(rh);
	_nRHs = (int)_rhs.size();
	
	_vtx = vtx; //vtx = (x,y,z)

	_E = rh.E();
	_eta = rh.eta();
	_phi = rh.phi();
	_t = rh.t();

	//calculate momentum vector from PV
	//centered at PV
	double dx = rh.x() - _vtx.at(0);
	double dy = rh.y() - _vtx.at(1);
	double dz = rh.z() - _vtx.at(2);
	//theta is calculated between beamline (z-dir) and x-y vector	
	double theta = atan2( sqrt(dx*dx + dy*dy), dz );
	double p_eta = -log(tan(theta/2));
	double p_phi = atan2(dy, dx);
	//double pt = _E*sin(theta); //mass = 0
	double pt = _E/cosh(p_eta);
//cout << "PV x " << _vtx.at(0) << " PV y " << _vtx.at(1) << " PV z " << _vtx.at(2) << " _eta " << _eta << " p_eta " << p_eta << " _phi " << _phi << " p_phi " << p_phi << " test_p_eta " << test_p_eta << endl;
	_px = pt*cos(p_phi);
	_py = pt*sin(p_phi);
	_pz = pt*sinh(p_eta);
	//_px = pt*cos(_phi);
	//_py = pt*sin(_phi);
	//_pz = pt*sinh(_eta);
	_kt2 = _px*_px + _py*_py; 		
	_mass = _calc_mass();

	
	_idx = 999;
	_ensure_valid_rap_phi();
	
	_cov = Matrix(3,3);
	_mu = Matrix(3,1);
	_mu.SetEntry(rh.eta(),0,0);
	_mu.SetEntry(rh.phi(),1,0);
	_mu.SetEntry(rh.t(),2,0);
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
	_t = 0;

	_vtx = vtx;
	double phi, eta;
	for(int i = 0; i < _nRHs; i++){		
		//for momentum vector centered at PV
		//theta is calculated between beamline (z-dir) and vector in x-y plane	
		x = rhs[i].x() - _vtx.at(0);
		y = rhs[i].y() - _vtx.at(1);
		z = rhs[i].z() - _vtx.at(2);
		theta = atan2( sqrt(x*x + y*y), z );
		phi = atan2(y, x);
		eta = -log(tan(theta/2.)); 
		//see https://cmssdt.cern.ch/lxr/source/DataFormats/CaloTowers/src/CaloTower.cc L145
		//pt = _E*sin(theta); //mass = 0, equivalent to below
		pt = rhs[i].E()/cosh(eta);
		_px += pt*cos(phi);
		_py += pt*sin(phi);
		_pz += pt*sinh(eta);
	//cout << "rhs ctor - i " << i << " px " << pt*cos(phi) << " py " << pt*sin(phi) << " pz " << pt*sinh(eta) << " E " << _rhs[i].E() << " pt " << pt << " eta " << eta << " phi " << phi << endl;
		
		_E += rhs[i].E();

	}
	CalculateCenter();
	
	
	_kt2 = _px*_px + _py*_py;
	//wraparound
	//_phi = acos(cos(_phi));
	_mass = _calc_mass();
//cout << "MAKING jet pts kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << sqrt(_kt2) << " mass " << _mass << " n pts " << rhs.size() << endl;

	_idx = 999;
	_ensure_valid_rap_phi();
	
	_pi = 0;
	_update_mom();
	
	CalculateCovariance();
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
	for(int i = 0; i < jets.size(); i++){
		//pt = _rhs[i].E()*cosh(_rhs[i].eta()); //consistent with mass = 0
		_px += jets[i].px();//pt*cos(_rhs[i].phi());
		_py += jets[i].py();//pt*sin(_rhs[i].phi());
		_pz += jets[i].pz();//pt*sinh(_rhs[i].eta());
		
		_E +=  jets[i].E();


	}
	_kt2 = _px*_px + _py*_py;
	_ensure_valid_rap_phi();

	_mass = _calc_mass();

	_idx = 999;

	CalculateCenter();	
	_pi = 0;
	_update_mom();
	CalculateCovariance();	
}

//_pi = 1 ==> 1 subcluster in jet, _pi != 1 ==> >1 subcluster in jet
//no weighting yet because if this is the only component, weight (ie pi) should be 1
Jet::Jet(const Matrix& mu, const Matrix& cov, double E, double pi, BayesPoint vtx, double detR){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	_E = E;
	_eta =  mu.at(0,0);
	_phi =  mu.at(1,0);
	_t = mu.at(2,0);
	
	//get x, y, z from eta, phi
	double x = detR*cos(_phi);
	double y = detR*sin(_phi);
	double theta = 2*atan2(1,exp(_eta));
	double z = detR/tan(theta);

	//calculate momentum vector from PV
	//centered at PV
	double dx = x - _vtx.at(0);
	double dy = y - _vtx.at(1);
	double dz = z - _vtx.at(2);
	//theta is calculated between beamline (z-dir) and x-y vector	
	double p_theta = atan2( sqrt(dx*dx + dy*dy), dz );
	double p_eta = -log(tan(p_theta/2));
	double p_phi = atan2(dy, dx);
	//double pt = _E*sin(theta); //mass = 0
	double pt = _E/cosh(p_eta);
	_px = pt*cos(p_phi);
	_py = pt*sin(p_phi);
	_pz = pt*sinh(p_eta);

	_kt2 = _px*_px + _py*_py;
	_mass = _calc_mass();
//cout << "jet subcl kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << pt << " m2 " << m2() << endl;

	_idx = 999;
	_ensure_valid_rap_phi();

	_mu = mu;
	_cov = cov;
	_pi = pi;

	_vtx = vtx;

	_update_mom();
}

Jet::Jet(BasePDF* pdf, double E, double pi, BayesPoint vtx, double detR){
	_vtx = BayesPoint(3);
	_mom = BayesPoint(4);
	_E = E;

	map<string, Matrix> params = pdf->GetParameters();
	_eta = params["mean"].at(0,0);
	_phi = params["mean"].at(1,0);
	_t =   params["mean"].at(2,0);
	//get x, y, z from eta, phi
	double x = detR*cos(_phi);
	double y = detR*sin(_phi);
	double theta = 2*atan2(1,exp(_eta));
	double z = detR/tan(theta);

	//calculate momentum vector from PV
	//centered at PV
	double dx = x - _vtx.at(0);
	double dy = y - _vtx.at(1);
	double dz = z - _vtx.at(2);
	//theta is calculated between beamline (z-dir) and x-y vector	
	double p_theta = atan2( sqrt(dx*dx + dy*dy), dz );
	double p_eta = -log(tan(p_theta/2));
	double p_phi = atan2(dy, dx);
	//double pt = _E*sin(theta); //mass = 0
	double pt = _E/cosh(p_eta);
	_px = pt*cos(p_phi);
	_py = pt*sin(p_phi);
	_pz = pt*sinh(p_eta);

	_kt2 = _px*_px + _py*_py;
	_mass = _calc_mass();
//cout << "jet subcl kt2 " << _kt2 << " px " << _px << " py " << _py << " pz " << _pz << " eta " << _eta << " phi " << _phi << " mass " << _mass << " energy " << _E << " pt " << pt << " m2 " << m2() << endl;

	_idx = 999;
	_ensure_valid_rap_phi();

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


Jet::Jet(BasePDFMixture* model, BayesPoint vtx, double gev, double detR){
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
	for(int i = 0; i < _nRHs; i++){
		//add rhs to jet
		BayesPoint rh = model->GetData()->at(i);
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
		_rhs[i].SetWeight(rh.w());	
	
		//calculate momentum vector from PV
		//centered at PV
		double dx = x - _vtx.at(0);
		double dy = y - _vtx.at(1);
		double dz = z - _vtx.at(2);
		//theta is calculated between beamline (z-dir) and x-y vector	
		double p_theta = atan2( sqrt(dx*dx + dy*dy), dz );
		double p_eta = -log(tan(p_theta/2));
		double p_phi = atan2(dy, dx);

		//cout << "eta " << rh.at(0) << " p_eta " << p_eta << " phi " << rh.at(1) << " p_phi " << p_phi << " x " << x << " dx " << dx << " y " << y << " dy " << dy << " z " << z << " dz " << dz << " PV x " << _vtx.at(0) << " PV y " << _vtx.at(1) << " PV z " << _vtx.at(2) << endl;
		//double pt = _E*sin(theta); //mass = 0
		pt = _rhs[i].E()/cosh(p_eta); 
		_px += pt*cos(p_phi);
		_py += pt*sin(p_phi);
		_pz += pt*sinh(p_eta);
	//cout << "basepdfmix ctor - i " << i << " px " << pt*cos(p_phi) << " py " << pt*sin(p_phi) << " pz " << pt*sinh(p_eta) << " E " << _rhs[i].E() << " pt " << pt << " p_eta " << p_eta << " eta " << eta << " p_phi " << p_phi << " phi " << phi << endl;
		
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
		
		_t += _rhs[i].t()*_rhs[i].E();

	}
	//_eta = model->GetData()->Centroid(0);
	//_phi = model->GetData()->CircularCentroid(1);
	////CircularCentroid returns via atan2 which has range [-pi,pi]
	//if (_phi < 0.0) {_phi += twopi;}
	//if (_phi >= twopi) {_phi -= twopi;} // can happen if phi=-|eps<1e-15|?
	//_t /= _E;
	//_mu = Matrix(3,1);
	//_ensure_valid_rap_phi();
	//_mu.SetEntry(_eta,0,0);
	//_mu.SetEntry(_phi,1,0);
	//_mu.SetEntry(_t,2,0);
	CalculateCenter();
	double ret = 0;
	for(auto rh : _rhs) ret += rh.phi_02pi();
	ret /= (double)_rhs.size();
	cout << "phi centroid from CalculateCenter " << _phi << " euclidean phi mean " << ret << " circular phi mean " << model->GetData()->CircularMean(1) << endl;
	//if(fabs(_phi - model->GetData()->CircularCentroid(1)) > 1e-10){
	//	cout << "Jet::WARNING - CalculateCenter() _phi " << _phi << " doesn't match model->GetData()->CircularCentroid() phi " << model->GetData()->CircularCentroid(1) << endl;
	//	//check if inputs to CalculateCenter() _phi and CircularCentroid phi are the same (ie phi coord and weights)
	//	cout << "CalculateCenter inputs" << endl;
	//	for(int i = 0; i < _nRHs; i++){
	//		cout << "phi " << _rhs[i].phi_02pi() << " w " << _rhs[i].GetWeight() << endl;
	//	}
	//	cout << "model data" << endl; model->GetData()->Print();
	//}
	
	_kt2 = _px*_px + _py*_py;

	_idx = 999;

	_update_mom();
	_mass = _calc_mass();
	CalculateCovariance();	

	//set constituents (subclusters)
	for(int k = 0; k < nsubcl; k++){
		auto params = model->GetLHPosteriorParameters(k);
		Ek = norms[k]/gev;
		Jet subcl(model->GetModel(k), Ek, model->GetPi(k), _vtx);
	
		//add rechits as "effective" crystals (ie weighted by their responsibility to this cluster)
		//their associated responsibility will be saved as the energy of that crystal for this subcluster
		double subcl_px = 0;
		double subcl_py = 0;
		double subcl_pz = 0;
		double totw = 0;
		double subcl_norm = 0;

		Matrix cov(3,3);
		double deta = 0;
		double dphi = 0;
		double dtime = 0;
		for(int n = 0; n < _nRHs; n++){
			JetPoint effRh = _rhs[n];
			subcl_norm += r_nk.at(n,k);	

			effRh.SetWeight(r_nk.at(n,k)/(_rhs[n].E()*gev));
			effRh.SetEnergy(_rhs[n].E()*effRh.GetWeight());

			BayesPoint effRh_pt({effRh.eta(), effRh.phi(), effRh.t()});
			effRh_pt.SetWeight(_rhs[n].E()*effRh.GetWeight());
			subcl.AddRecHit(effRh);
			//recalculate momentum vector accordingly
			//calculate momentum vector from PV
			//centered at PV
			double dx = effRh.x() - _vtx.at(0);
			double dy = effRh.y() - _vtx.at(1);
			double dz = effRh.z() - _vtx.at(2);
			//theta is calculated between beamline (z-dir) and x-y vector	
			double p_theta = atan2( sqrt(dx*dx + dy*dy), dz );
			double p_eta = -log(tan(p_theta/2));
			double p_phi = atan2(dy, dx);

			//cout << "eta " << rh.at(0) << " p_eta " << p_eta << " phi " << rh.at(1) << " p_phi " << p_phi << " x " << x << " dx " << dx << " y " << y << " dy " << dy << " z " << z << " dz " << dz << " PV x " << _vtx.at(0) << " PV y " << _vtx.at(1) << " PV z " << _vtx.at(2) << endl;
			//double pt = _E*sin(theta); //mass = 0
			pt = effRh.E()/cosh(p_eta); 
			subcl_px += pt*cos(p_phi);
			subcl_py += pt*sin(p_phi);
			subcl_pz += pt*sinh(p_eta);
		}
		//set subcluster momentum three-vector and mass
		//cov.mult(cov,1/totw);
		//subcl.SetCovariance(cov);
		subcl.CalculateCenter();
		subcl.CalculateCovariance();
//if(nsubcl == 1){
//cout << "subcl norm " << subcl_norm << " " << norms[k] << " with total w for cluster #" << k << ": " << totw << " for jet with " << _nRHs << " rechits" << endl;
//cout << "subcl center eta " << subcl.eta() << " phi " << subcl.phi() << " jet eta " << _eta << " phi " << _phi << endl;
//cout << "subcl cov from CalcCov" << endl; subcl._cov.Print();
//cout << "jet cov" << endl; _cov.Print();
//}
		subcl.SetP(subcl_px, subcl_py, subcl_pz);
		_constituents.push_back(subcl);
	}

}

//copy ctor
Jet::Jet(const Jet& j){
	_px = j.px();
	_py = j.py();
	_pz = j.pz();
	_E = j.E();
	_mom = BayesPoint({_px, _py, _pz, _E});
	_eta = j._eta;
	_phi = j._phi;
	_t = j._t;

	_kt2 = j._kt2;
	_mass = j._mass;
	
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
	_constituents.clear();
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
	//recalculate kt2 of cluster
	_kt2 = _px*_px + _py*_py;
	_mass = _calc_mass();

	//set time to be energy-weighted average of rec hit times
	//_set_time();

	////recalculate eta and phi of cluster
	//_set_rap_phi();

	CalculateCenter();

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
	//_set_time();
	//recalculate kt2 of cluster
	_kt2 = sqrt(pt);
	//recalculate eta and phi of cluster
	//_set_rap_phi();
	CalculateCenter();

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

