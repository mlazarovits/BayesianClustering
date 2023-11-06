#include "BasicDetectorSim.hh"
#include "TMath.h"
#include <algorithm>

BasicDetectorSim::BasicDetectorSim(){
	//init parameters to CMS ECAL geometry
	//see CMS TDR (ISBN 978-92-9083-268-3)
	_rmax = 1.29; //inner radius of ECAL barrel (m)
	_b = 3.8; //field of CMS solenoid (T)
	_netacal = 170; //number of cells in eta in ECAL barrel (2*85)
	_nphical = 360; //number of cells in phi in ECAL barrel (neta*nphi = 61200 total cells in barrel)
	_deta = 0.0174; //eta component of cell cross-section (2.2 cm - Moliere radius)
	_dphi = 0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_calEres = 0.00445; //energy resolution approximated as radiation length/2. (rad length = 0.89 cm) 
	_calTres = 0.2; //time resolution for CMS ECAL (ns) (200 ps)
	_sagres = 0.000013; //value from LHC parameters in PGS (examples/par/lhc.par)
	_rs = RandomSample(); //random sampler
	_gev = 0; //default
	_default_transfer = true;
	_nevts = 1000;
	//initialize cal
	for(int i = 0; i < _netacal; i++)
		_cal.push_back(vector<double>(_nphical));
}

//ctor with input pythia cmnd file
BasicDetectorSim::BasicDetectorSim(string infile){
	//init parameters to CMS ECAL geometry
	//see CMS TDR (ISBN 978-92-9083-268-3)
	_rmax = 1.29; //inner radius of ECAL barrel (m)
	_b = 3.8; //field of CMS solenoid (T)
	_netacal = 170; //number of cells in eta in ECAL barrel (2*85)
	_nphical = 360; //number of cells in phi in ECAL barrel (neta*nphi = 61200 total cells in barrel)
	_deta = 0.0174; //eta component of cell cross-section (2.2 cm - Moliere radius)
	_dphi = 0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_calEres = 0.00445; //energy resolution approximated as radiation length/2. (rad length = 0.89 cm) 
	_calTres = 0.2; //time resolution for CMS ECAL (ns) (200 ps)
	_sagres = 0.000013; //value from LHC parameters in PGS (examples/par/lhc.par)
	_rs = RandomSample(); //random sampler
	_gev = 999; //default
	_default_transfer = true;
	//initialize cal
	for(int i = 0; i < _netacal; i++)
		_cal.push_back(vector<double>(_nphical));

	//sets pythia settings by given .cmnd file
	_pythia.readFile(infile);
	_nevts = _pythia.mode("Main:numberOfEvents");

}

//use pythia8 to simulate events
//TODO: find syntax to simulate min bias events for PU
void BasicDetectorSim::SimQCD(){
	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 13 TeV.
	_pythia.readString("HardQCD:all = on");
	_pythia.readString("PhaseSpace:pTHatMin = 20.");
	_pythia.readString("Beams:eCM = 13000.");
	_pythia.init();
  
}

void BasicDetectorSim::SimTTbar(){
	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 13 TeV.
	_pythia.readString("Top:all = on");
	_pythia.readString("PhaseSpace:pTHatMin = 20.");
	_pythia.readString("Beams:eCM = 13000.");
	_pythia.init();
  
}



//default arg is nevts = 1
void BasicDetectorSim::SimulateEvents(int evt){
	double maxeta = 1.749;
	//(nominal values of) constants for energy resolution
	//these are from CMS TDR Fig. 1.7
	//the more conservative values are taken (s.t. sig/E is larger than other values)
	_etamin = -(_netacal/2.)*_deta;
	_etamax = (_netacal/2.)*_deta;

	for(int i = 0; i < _nevts; i++){
		if(!_pythia.next() || i != evt) continue;
		cout << "getting event " << evt << " with " << _pythia.event.size() << " particles" << endl;
		//loop through all particles
		//make sure to only record those that would
		//leave RecHits in ECAL (ie EM particles (ie ie photons and electrons))
		for(int p = 0; p < _pythia.event.size(); p++){
			//reset reco particle four momentum
			Particle particle = _pythia.event[p];
			//make sure particle is final-state and (probably) stable
			if(particle.statusHepMC() != 1) continue;
			//make sure particle is in detector acceptance
			//since this is a CMS ECAL sim, use CMS ECAL geometry
			if(fabs(particle.eta()) > maxeta) continue;
			
			//min pt cut here
			//2 GeV particle cut
			//if(particle.pT() < 2) continue;
			
			//create new particle for reco one
			RecoParticle rp(particle); 
			//calculate new pt (does full pvec but same pz)
			CalcTrajectory(rp);
			//check if in call cell crack
			if(_in_cell_crack(rp))
				continue;
			if(rp.particle.e() < 0) continue;
			//fill ecal cell with reco particle
	cout << "particle " << p << endl;
			FillCal(rp);
			
			//save gen particle four vector
		//	_genps.push_back(_pythia.event[p]);
			//save reco particle four vector
		//	_recops.push_back(rp);	
	
		}

	}

}

//updates pvec of p
//this code is based on ParticlePropagator class in Delphes
//https://github.com/delphes/delphes/blob/master/modules/ParticlePropagator.cc
void BasicDetectorSim::CalcTrajectory(RecoParticle& rp){
	//calculate halflength from max eta
	double theta = 2*atan(exp(-_etamax));
	double halfLength = _rmax/tan(theta);

	//pythia units are in mm (or mm/c for time, natural units)
	//convert to cm (or cm/c where c is in cm/ns)
	double x = rp.particle.xProd() * 1e-1;
	double y = rp.particle.yProd() * 1e-1;
	double z = rp.particle.zProd() * 1e-1;
	//conver to ns: mm/c * 1cm/10 mm = cm/c * cm/ns = ns
	double t = rp.particle.tProd() * 1e-1 * _sol;

	double q = rp.particle.charge();

	//eta check is done in SimulateEvents but double check with halflength
	if(fabs(z) > halfLength) return;

	double px = rp.particle.px();
	double py = rp.particle.py();
	double pz = rp.particle.pz();
	double pt = rp.particle.pT();
	double pt2 = rp.particle.pT2();
	double e = rp.particle.e();

	double tmp, tr, tz, x_t, y_t, z_t, r_t;
	double gammam, omega, r, alpha;
	double x_c, y_c, r_c, vz;
	double phi0, phid, phit, pio;
	double xd, yd, zd, td;
	//uncharged trajectory or no mag field
	//TODO: check math
	if(fabs(q) < 1e-9 || fabs(_b) < 1e-9){
		//calculate time to detector
		// solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
		tmp = px*y - py*x;
		tr = (sqrt(pt2 * _rmax - tmp*tmp) - px * x - py * y)/pt2;
		tz = (TMath::Sign(halfLength, pz) - z) / pz;
		t = fmin(tr, tz);

		//calculate new x, y, z based on time to detector
		x_t = x + px*t;
		y_t = y + py*t;
		z_t = z + pz*t;	
	
		rp.Position.SetCoordinates(x_t, y_t, z_t, rp.Position.T() + t * e);
		//set to gen values
		rp.Momentum.SetCoordinates(rp.Momentum.pt(), rp.Momentum.eta(), rp.Momentum.phi(), rp.Momentum.e());

	}
	//charged particles in magnetic field
	else{
		//do helix calculation
		// 1. initial transverse momentum p_{T0}: Part->pt
		//    initial transverse momentum direction phi_0 = -atan(p_{X0} / p_{Y0})
      		//    relativistic gamma: gamma = E / mc^2; gammam = gamma * m
      		//    gyration frequency omega = q * Bz / (gammam)
      		//    helix radius r = p_{T0} / (omega * gammam)

		gammam = e * 1e9 / (_sol*_sol); //gammam in [eV/c^2]
		omega = q * _b / gammam;
		r = pt / (q * _b) * 1e9 / _sol; //in [m]

		phi0 = atan2(py, px); // [rad] in [-pi, pi]
		
		//2. helix axis coordinates
		x_c = x + r * sin(phi0);
		y_c = y - r * cos(phi0);
		r_c = sqrt(x_c*x_c + y_c*y_c);
		
		// time of closest approach
		td = (phi0 + TMath::ATan2(x_c, y_c)) / omega;
		
		// remove all the modulo pi that might have come from the atan
		pio = TMath::Abs(TMath::Pi() / omega);
		while(TMath::Abs(td) > 0.5 * pio)
		{
		  td -= TMath::Sign(1.0, td) * pio;
		}


		vz = pz * _sol/e;

		// calculate coordinates of closest approach to z axis
		//new phi
		phid = phi0 - omega * td;
		xd = x_c - r * TMath::Sin(phid);
		yd = y_c + r * TMath::Cos(phid);
		zd = z + vz * td;
		
		// momentum at closest approach
		px = pt * TMath::Cos(phid);
		py = pt * TMath::Sin(phid);

		//reset momentum
		rp.Momentum.SetCoordinates(pt, rp.Momentum.Eta(), phid, rp.Momentum.E());
		
		// 3. time evaluation t = TMath::Min(t_r, t_z)
		//    t_r : time to exit from the sides
		//    t_z : time to exit from the front or the back
		tz = (vz == 0.0) ? 1.0E99 : (TMath::Sign(halfLength, pz) - z) / vz;
		
		if(r_c + fabs(r) < _rmax)
		{
			// helix does not cross the cylinder sides
			t = tz;
		}
		else
		{
			alpha = acos((r * r + r_c * r_c - _rmax * _rmax) / (2 * fabs(r) * r_c));
			tr = td + fabs(alpha / omega);
			
			t = fmin(tr, tz);
		}
		
		// 4. position in terms of x(t), y(t), z(t)
		phit = phi0 - omega * t;
		x_t = x_c - r * TMath::Sin(phit);
		y_t = y_c + r * TMath::Cos(phit);
		z_t = z + vz * t;
		r_t = TMath::Hypot(x_t, y_t);
		
		//reset spatial coordinates	
		if(r_t > 0.0)
		{
			rp.Position.SetXYZT(x_t * 1.0E3, y_t * 1.0E3, z_t * 1.0E3, rp.Position.T() + t * _sol * 1.0E3);
		}
			
		if(rp.Momentum.e() < 1.05)
		{
			cout << "x: " << x_t * 1.0E3 << " y: " <<  y_t * 1.0E3 << " z: " <<  z_t * 1.0E3 << " t: " << rp.Position.T() + t * _sol * 1.0E3 << " e: " << rp.Momentum.e() << endl;
		}
			

	}
	
}





//fill ecal cells (and create rechits)
void BasicDetectorSim::FillCal(const RecoParticle& rp){
	Particle p = rp.particle;
	ROOT::Math::PtEtaPhiEVector Momentum = rp.Momentum;
	ROOT::Math::XYZTVector Position = rp.Position;
	//set energy resolution based on particle type
	double e, e_cell, e_sig, t_cell;
	e = Momentum.e();
	//if particle is a photon or electron
	//use approx CMS ECAL energy resolution
	if(fabs(p.id()) != 11 && fabs(p.id()) != 22){
		//do energy resolution smearing with a gaussian
		//centered on nominal energy and with an energy-dependent spread 
		//that follows equation 1.2 in CMS TDR
		//(sig/E)^2 = (s/sqrt(E))^2 + (n/E)^2 + c^2
		e_sig = sqrt( _s*_s*sqrt(e) + _n*_n + _c*_c*e*e );
	}
	else{
		//else if it's a charged hadron
		//use approx CMS HCAL energy resolution
		//using formula from Thomson pg. 21
		e_sig = sqrt(e)*0.5;

	}

	double eta, phi;
	int ieta, iphi;
	double x, y, z, t, r, theta; //coordinates of JetPoint (rechit)
	//assume all energy gets deposited in ecal
	//get eta,phi indices from fcn
	eta = Momentum.eta();
	phi = Momentum.phi();
	//get particle time of arrival
	t = Position.t();
	//length in xy plane
	if(p.isNeutral() || _b < 1e-9)
		r = Momentum.pt() / _sol;
	else
		r = Momentum.pt() / (p.charge() * _b) / _sol;
		//TODO: check units	

	_get_etaphi_idx(eta, phi, ieta, iphi);	
	
	//simulate shower cell radius
	//or just take it to be constant :)
	//the Moliere radius of the CMS ECAL is 2.2 cm = 0.1 is deta or dphi 
	//two Moliere radii contain 95% of shower's energy deposition (on average)
	//for constant shower radius, take 2*2.2 to be 2*0.1 (0.1 corresponds to 2.2 cm or 1 deg in deltaPhi or deltaEta in barrel)
	//assume circular (orthogonal in eta-phi) shower
	double showrad = 0.2;

	//energy threshold for reconstruction
	double e_thresh = 0.1; //(GeV)

	//check all cells in a (generous) 5x5 grid around the central ieta, iphi cell
	//ieta-2 < eta < ieta+2
	double aeta, beta, aphi, bphi, ceta, cphi; //integral bounds
	int iieta, iiphi;
	cout << "doing 5x5 grid loop in eta, phi" << endl;
	for(int i = -2; i < 2; i++){
		for(int j = -2; j < 2; j++){
			//get ieta, iphi for cell in grid
			iieta = ieta+i;
			iiphi = iphi+j;

			//get integral bounds
			_get_etaphi(iieta, iiphi, ceta, cphi);
	 		aeta = ceta - _deta/2.;
	 		beta = ceta + _deta/2.;
		
			aphi = cphi - _dphi/2.;
			bphi = cphi + _dphi/2.;
	
			//do (independent) 2D (eta-phi) gaussian integral
			//with mean = particle eta, phi
			//and sigma = shower radius
			//assume normalized gaussian distribution 
			//eta integral * phi integral
			//percent of particle energy in this cell
			e_cell = e*0.5*( erf((beta - eta)/(showrad*sqrt(2))) - erf((aeta - eta)/(showrad*sqrt(2))) ) * 0.5*( erf((bphi - phi)/(showrad*sqrt(2))) - erf((aphi - phi)/(showrad*sqrt(2))) );
	
			if(e_cell < 1e-3) continue;
cout << std::setprecision(15) << endl;
cout << "cell energy unsmeared: " << e_cell << " esig: " << e_sig << " e: " << e << " pt " << Momentum.pt() << endl;
			//smear energy in each cell
			//make sure e can't be negative
			do{	
				cout << "original energy: " << e_cell << endl;
				e_cell = _rs.SampleGaussian(e_cell, e_sig, 1).at(0); //returns a vector, take first (and only) element
				cout << "smeared energy: " << e_cell << endl;
			}while(e_cell < 0);
	cout << "done while loop" << endl;	
			//in each cell to find energy deposited
			//if integral is below some threshold e_cell = 0
			if(e_cell < e_thresh) continue;

			//add energy to right ieta, iphi cell	
			//energy has already been smeared
			_cal[ieta][iphi] += e_cell;
cout << "added energy to cal" << endl;
cout << "t: " << t << " calTres: " << _calTres << endl;
			//update rec hit information
			//rec hit time is based on nominal particle time
			//with spread 
			//t can be negative (early times)
			t_cell = _rs.SampleGaussian(t, _calTres, 1).at(0);
cout << "done t sampling" << endl;
			//get x, y, z based on cell eta phi
			x = r*cos(cphi);
			y = r*sin(cphi);
			theta = 2*tan(exp(-ceta));
			z = r*cos(theta); 

			_cal_rhs.push_back( JetPoint(x, y, z, t_cell) );
			_cal_rhs[i].SetEnergy(e_cell);	
			if(_default_transfer) _gev += e_cell;
		cout << "done particle loop" << endl;
		}
	}

}



//get "Jets" for clustering
void BasicDetectorSim::GetRecHits(vector<Jet>& rhs){
	rhs.clear();
	if(_default_transfer) _gev /= (double)_cal_rhs.size();
	for(int j = 0; j < _cal_rhs.size(); j++){
		_cal_rhs[j].SetWeight(_cal_rhs[j].e()/_gev);
		rhs.push_back(Jet(_cal_rhs[j]));
	}

}




bool BasicDetectorSim::_in_cell_crack(const RecoParticle& rp){
	Particle p = rp.particle;
	double eta = p.eta();
	double phi = p.phi();
	
	if(fabs(eta/_deta) - double(int(eta/_deta)) < _crack_frac){
		return true;
	}
	if(fabs(phi/_dphi) - double(int(phi/_dphi)) < _crack_frac){
		return true;
	}

	return false;
}


void BasicDetectorSim::_get_etaphi_idx(double eta, double phi, int& ieta, int& iphi){
	//get phi index
	//offset for indexing starting at 0
	iphi = int((phi/_dphi)+0.5) - 1;
	if(iphi > _nphical) iphi -= _nphical;
	else if(iphi < 0) iphi += _nphical; 

	//get eta index
	ieta = 0;
	if(eta > _etamin && eta < _etamax){
		//offset for indexing starting at 0
		ieta = int(((eta - _etamin)/_deta) + 0.5) - 1; 
      		if(ieta < _netacal) ieta = 0;
		else if(ieta < 0) ieta = 0;
      
	}
}

//get eta-phi center given ieta, iphi cell indices
void BasicDetectorSim::_get_etaphi(int ieta, int iphi, double& eta, double& phi){
	//eta      
	eta = _etamin + _deta/2. + double(ieta-1)*_deta;      
	//phi 
	phi = _dphi/2. + double(iphi-1)*_dphi;
}
