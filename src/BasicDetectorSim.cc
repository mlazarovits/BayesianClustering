#include "BasicDetectorSim.hh"
#include "TMath.h"
#include "fastjet/ClusterSequence.hh"

#include <TFile.h>
#include <algorithm>

///////////////////////////
//all lengths are in meters (saving to JetPoints in cm)
//all times are in seconds (saving to JetPoints in ns)
//all energies, masses and momenta are in GeV (same for JetPoints)
///////////////////////////


BasicDetectorSim::BasicDetectorSim(){
	//init parameters to CMS ECAL geometry
	//see CMS TDR (ISBN 978-92-9083-268-3)
	_rmax = 1.29; //inner radius of ECAL barrel (m)
	_b = 3.8; //field of CMS solenoid (T)
	_netacal = 170; //number of cells in eta in ECAL barrel (2*85)
	_nphical = 360; //number of cells in phi in ECAL barrel (neta*nphi = 61200 total cells in barrel)
	_deta = 2*acos(-1)/360.; //0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_dphi = 2*acos(-1)/360.; //0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_calEres = 0.00445; //energy resolution approximated as radiation length/2. (rad length = 0.89 cm) 
	_calTres = 0.2 * 1e-9; //time resolution for CMS ECAL (s) (200 ps)
	_sagres = 0.000013; //value from LHC parameters in PGS (examples/par/lhc.par)
	_rs = RandomSample(); //random sampler
	_gev = 0; //default
	_default_transfer = true;
	_nevts = 1000;
	//initialize cal - save e, t, n emissions
	_etamax = 1.479 + _deta/2.; //puts outermost corner at true etamax
	_etamin = -_etamax;
	_phimin = -acos(-1);
	//energy threshold for reconstruction
	_ethresh = 0.7; //(GeV)
	_ncell = 1;
	_pu = false; //pileup switch
	//create container for cal cell info
	//E, t, nEmissions
	for(int i = 0; i < _netacal; i++){
		_cal.push_back({});
		for(int j = 0; j < _nphical; j++)
			_cal[i].push_back(Point({0.,0.,0.}));
	}
}

//ctor with input pythia cmnd file
BasicDetectorSim::BasicDetectorSim(string infile){
	//init parameters to CMS ECAL geometry
	//see CMS TDR (ISBN 978-92-9083-268-3)
	_rmax = 1.29; //inner radius of ECAL barrel (m)
	_b = 3.8; //field of CMS solenoid (T)
	_netacal = 170; //number of cells in eta in ECAL barrel (2*85)
	_nphical = 360; //number of cells in phi in ECAL barrel (neta*nphi = 61200 total cells in barrel)
	_deta = 2*acos(-1)/360.; //0.0174; //eta component of cell cross-section (2.2 cm - Moliere radius)
	_dphi = 2*acos(-1)/360.; //0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_calEres = 0.00445; //energy resolution approximated as radiation length/2. (rad length = 0.89 cm) 
	_calTres = 0.2 * 1e-9; //time resolution for CMS ECAL (s) (200 ps)
	_sagres = 0.000013; //value from LHC parameters in PGS (examples/par/lhc.par)
	_rs = RandomSample(); //random sampler
	_gev = 999; //default
	_default_transfer = true;
	_etamax = 1.479;
	_etamin = -_etamax;
	_phimin = -acos(-1);
	//energy threshold for reconstruction
	_ethresh = 0.1; //(GeV)
	_ncell = 1;
	_pu = false; //pileup switch
	//initialize cal - save e, t, n emissions
	for(int i = 0; i < _netacal; i++){
		_cal.push_back({});
		for(int j = 0; j < _nphical; j++)
			_cal[i].push_back(Point({0.,0.,0.}));
	}

	//sets pythia settings by given .cmnd file
	_pythia.readFile(infile);
	_nevts = _pythia.mode("Main:numberOfEvents");

}


//use pythia8 to simulate events
void BasicDetectorSim::_simQCD(){
	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 13 TeV.
	_pythia.readString("HardQCD:all = on");
	_pythia.readString("PhaseSpace:pTHatMin = 20.");
	_pythia.readString("Beams:eCM = 13000.");
	if(_verb > 1) cout << "Simulating QCD" << endl;
}

void BasicDetectorSim::_simTTbar(){
	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 13 TeV.
	_pythia.readString("Top:all = on");
	_pythia.readString("PhaseSpace:pTHatMin = 20.");
	_pythia.readString("Beams:eCM = 13000.");
	if(_verb > 1) cout << "Simulating ttbar" << endl;
}



//default arg is nevts = 1
void BasicDetectorSim::SimulateEvents(int evt){
	Pythia8::Pythia pileup(_pythia.settings, _pythia.particleData);
	if(find(_procs_to_sim.begin(), _procs_to_sim.end(), ttbar) != _procs_to_sim.end()){
		_simTTbar();
	}	
	if(find(_procs_to_sim.begin(), _procs_to_sim.end(), qcd) != _procs_to_sim.end()){
		_simQCD();
	}	
	if(_pu){
  		if(_verb == 0) pileup.readString("Print:quiet = on");
		pileup.readString("Random:setSeed = on");
  		pileup.readString("Random:seed = 10000002");
		pileup.readString("SoftQCD:nonDiffractive = on"); //minbias events only
  		pileup.settings.parm("Beams:eCM = 13000.");
		pileup.init();			
		if(_verb > 1) cout << "Simulating pileup" << endl;
	}

	_pythia.init();
	//sigma for z-smearing
	//beamspot spread is ~5 cm = 0.05 m to use with _sol
	double zig = 0.05/2.;
	double znew, tnew;
	//calculate halflength from max eta
	double theta = 2*atan(exp(-_etamax));
	double zmax = _rmax/tan(theta); //[mm]

	Pythia8::Event sumEvent; //one object where individual events are collected
	
	//declare fjinput/output containers
	vector<fastjet::PseudoJet> fjinputs, fjoutputs;
	//fastjet objects - "AK4" jets
	double Rparam = 0.4;
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recomb = fastjet::E_scheme;
	fastjet::JetDefinition jetdef = fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam, recomb, strategy); 
	
	for(int i = 0; i < _nevts; i++){
		if(!_pythia.next() || i != evt) continue;
		//store event info if pileup is on
		sumEvent = _pythia.event;
		if(_pu){
			//simulate n pileup events
			int nPU = _rs.SamplePoisson(_nPUavg,1).at(0);
			for(int p = 0; p < nPU; p++){
				pileup.next();	
				sumEvent += pileup.event;
			}
		}
		// Reset Fastjet input
    		fjinputs.clear();
		fjinputs.resize(0);
		
		//loop through all particles
		//make sure to only record those that would
		//leave RecHits in ECAL (ie EM particles (ie ie photons and electrons))
		for(int p = 0; p < sumEvent.size(); p++){
			//reset reco particle four momentum
			Pythia8::Particle particle = sumEvent[p];
			//make sure particle is final-state and (probably) stable
			if(particle.statusHepMC() != 1) continue;
			// No neutrinos
      			if (_pythia.event[i].idAbs() == 12 || _pythia.event[i].idAbs() == 14 ||
      			    _pythia.event[i].idAbs() == 16)     continue;
			
			
			//set production vertex from z-smearing
			_rs.SetRange(particle.zProd()*1e-3 - zig, particle.zProd()*1e-3 + zig);
			znew = _rs.SampleGaussian(particle.zProd()*1e-3, zig, 1).at(0);	
			//set time from linear model with slope = { z > 0 ? 1/sol : -1/sol }
			tnew = particle.zProd()*1e-3 >= 0 ? znew*1./(_sol) : -1./(_sol)*znew;
			//original pythia coords are in m, convert to mm
			particle.vProd(particle.xProd(), particle.yProd(), znew*1e3, tnew*1e-3);
			
			//create new particle for reco one
			RecoParticle rp(particle); 
			//make sure particle is in detector acceptance
			//since this is a CMS ECAL sim, use CMS ECAL geometry
			//this is for gen particles
			//include z offset
			//z and eta are one to one
			if(fabs(particle.eta() + znew*(_etamax/(zmax))) > _etamax) continue;
			//if(fabs(rp.Position.eta() + znew*(_etamax/(zmax))) > _etamax) continue;
			//reset phi for reco particle to include zshift
		
			//if particle energy is less than single rechit thresh, skip	
			if(particle.e() < _ethresh) continue;
			
			//calculate new pt (does full pvec but same pz)
			CalcTrajectory(rp);
			//check if in cal cell crack
			if(_in_cell_crack(rp))
				continue;
			//if track curls up/exceeds zmax
			if(fabs(rp.Position.z()) >= zmax || fabs(rp.Position.eta()) > _etamax) continue;
			//fill ecal cell with reco particle
			FillCal(rp);
		

			//add particle to fastjet
			//running fastjet on gen particles, no shower, etc.
      			fjinputs.push_back( fastjet::PseudoJet( _pythia.event[p].px(),
      			  _pythia.event[p].py(), _pythia.event[p].pz(), _pythia.event[p].e() ) );
			

			//save reco particle four vector (with corresponding gen info)
			_recops.push_back(rp);	
			
		}
		//run fastjet
		fastjet::ClusterSequence cs(fjinputs, jetdef);
		//get jets - min 5 pt
		fjoutputs = cs.inclusive_jets(5.);
		for(int j = 0; j < fjoutputs.size(); j++) _jets.push_back(fjoutputs[j]);

	}
	MakeRecHits();
	ReconstructEnergy();
	//sort jets by pt
	_jets = sorted_by_pt(_jets);
}

//updates pvec of p
//this code is based on ParticlePropagator class in Delphes
//https://github.com/delphes/delphes/blob/master/modules/ParticlePropagator.cc
void BasicDetectorSim::CalcTrajectory(RecoParticle& rp){
	ROOT::Math::PtEtaPhiEVector Momentum = rp.Momentum;
	ROOT::Math::XYZTVector Position = rp.Position;
	//calculate halflength from max eta
	double theta = 2*atan(exp(-_etamax));
	double halfLength = _rmax/tan(theta); //[m]

	//pythia units are in mm (or mm/c for time, natural units)
	//convert to m (or m/c where c is in m/s)
	//initialized to production coordinates
	double x = Position.x()*1e-3;
	double y = Position.y()*1e-3;
	double z = Position.z()*1e-3;
                           
	double q = rp.Particle.charge();

	//eta check is done in SimulateEvents but double check with halflength
	if(fabs(z) > halfLength) return;

	double px = Momentum.px();
	double py = Momentum.py();
	double pz = Momentum.pz();
	double pt = Momentum.pt();
	double pt2 = Momentum.Perp2();
	double e = Momentum.e();

	double tmp, tr, tz, x_t, y_t, z_t, r_t, t;
	double gammam, omega, r, alpha;
	double x_c, y_c, r_c, vz;
	double phi0, phid, phit, pio, etad;
	double xd, yd, zd, td, dpv;
	//uncharged trajectory or no mag field
	if(fabs(q) < 1e-9 || fabs(_b) < 1e-9){
		//calculate time to detector
		// solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
		tmp = px * y - py * x;
      		tr = (TMath::Sqrt(pt2 * _rmax*_rmax - tmp * tmp) - px * x - py * y) / pt2;
		tz = (TMath::Sign(halfLength, pz) - z) / pz;
		t = fmin(tr, tz)*(e/_sol); //t*e/c ~ [m/GeV]*[GeV*s/m] = s
		//calculate new x, y, z based on time to detector
		x_t = x + (px/e)*_sol*t;
		y_t = y + (py/e)*_sol*t;
		z_t = z + (pz/e)*_sol*t;	
		
		//time = r/(beta*c) = r/(p*c/E) = r*E/c*p	
		rp.Position.SetCoordinates(x_t, y_t, z_t, Position.T() + t);
		//keep momentum at gen values - update eta, phi
		rp.Momentum.SetCoordinates(pt, rp.Momentum.eta(), rp.Momentum.phi(), e); 	
	}
	//charged particles in magnetic field
	else{
		//do helix calculation
		// 1. initial transverse momentum p_{T0}: Part->pt
		//    initial transverse momentum direction phi_0 = -atan(p_{X0} / p_{Y0})
      		//    relativistic gamma: gamma = E / mc^2; gammam = gamma * m
      		//    gyration frequency omega = q * Bz / (gammam)
      		//    helix radius r = p_{T0} / (omega * gammam)

		gammam = e * 1e9 / (_sol * _sol); //gammam in [eV/c^2], c needs to be in [m/s]
		omega = q * _b / (gammam); //in [89875518/s] - should be [rad/s]?
		r = pt * 1e9 / (q * _b * _sol); //in [m]

		phi0 = atan2(py, px); // [rad] in [-pi, pi]
		//2. helix axis coordinates
		x_c = x + r * sin(phi0);
		y_c = y - r * cos(phi0);
		r_c = sqrt(x_c*x_c + y_c*y_c);
		
		// time of closest approach
		td = (phi0 + atan2(x_c, y_c)) / omega; //original delphes

		// remove all the modulo pi that might have come from the atan
		pio = fabs(TMath::Pi()/omega);
		while(fabs(td) > 0.5 * pio)
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
		rp.Momentum.SetCoordinates(pt, rp.Momentum.Eta(), phid, e);
		// 3. time evaluation t = TMath::Min(t_r, t_z)
		//    t_r : time to exit from the sides
		//    t_z : time to exit from the front or the back
		tz = (vz == 0.0) ? 1.0E99 : (TMath::Sign(halfLength, pz) - z) / vz;
		
		if(r_c + fabs(r) < _rmax)
		{
			// helix does not cross the cylinder sides
			//convert to [s] from [m/Gev]
			t = tz;
		}
		else
		{
			alpha = acos((r * r + r_c * r_c - _rmax * _rmax) / (2 * fabs(r) * r_c));
			tr = td + fabs(alpha / omega);
			//convert to [s] from [m/Gev]
			t = fmin(tr, tz);
		}
		// 4. position in terms of x(t), y(t), z(t)
		phit = phi0 - omega * t;
		x_t = x_c - r * TMath::Sin(phit);
		y_t = y_c + r * TMath::Cos(phit);
		z_t = z + vz * t;
		r_t = TMath::Hypot(x_t, y_t);
		//position in m, t in s
		if(r_t > 0.0)
			rp.Position.SetCoordinates(x_t, y_t, z_t, (Position.T() + t));
	}

}





//fill ecal cells (and create rechits)
void BasicDetectorSim::FillCal(RecoParticle& rp){
	Pythia8::Particle p = rp.Particle;
	PtEtaPhiEVector Momentum = rp.Momentum;
	XYZTVector Position = rp.Position;
	//set energy resolution based on particle type
	double e, e_cell, e_sig, t_cell, pt, q;
	e = Momentum.e();
	pt = Momentum.pt();
	q = p.charge();



	double eta, phi;
	int ieta, iphi;
	double x, y, z, t; //coordinates of JetPoint (rechit)
	double dpv, r, theta;
	double pz = rp.Momentum.pz();
	//assume all energy gets deposited in ecal
	//get eta,phi indices from fcn
	eta = Position.eta();
	phi = Position.phi();
	//get ieta, iphi
	_get_etaphi_idx(eta, phi, ieta, iphi);
	//get particle time of arrival
	//in seconds
	t = Position.t();
	//simulate shower cell radius
	//or just take it to be constant :)
	//the Moliere radius of the CMS ECAL is 2.2 cm = 0.01 is deta or dphi 
	//this should have ~90% of the shower's energy on average which would correspond to ~2 sigma
	//assume circular (orthogonal in eta-phi) shower
	double showrad = _deta/2.; //same for phi

	double aeta, beta, aphi, bphi, ceta, cphi; //integral bounds
	int iieta, iiphi;
	for(int i = -_ncell; i < _ncell+1; i++){
		for(int j = -_ncell; j < _ncell+1; j++){
			//get ieta, iphi for cell in grid
			iieta = ieta+i;
			iiphi = iphi+j;

			//if these indices are out of bounds for detector indices, skip
			if(iieta >= _netacal || iieta < 0) continue;
			if(iiphi >= _nphical || iiphi < 0) continue;

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
			e_cell = e*(0.5*( erf((beta - eta)/(showrad*sqrt(2))) - erf((aeta - eta)/(showrad*sqrt(2))) ) * 0.5*( erf((bphi - phi)/(showrad*sqrt(2))) - erf((aphi - phi)/(showrad*sqrt(2))) ));
			
			//add energy to right ieta, iphi cell	
			_cal[iieta][iiphi].SetValue(_cal[iieta][iiphi].at(0)+e_cell, 0);	
			//add time to right ieta, iphi cell	
			_cal[iieta][iiphi].SetValue(_cal[iieta][iiphi].at(1)+t, 1);	
			//add number of emissions to right ieta, iphi cell	
			_cal[iieta][iiphi].SetValue(_cal[iieta][iiphi].at(2)+1,2);


		}
	}
}



void BasicDetectorSim::MakeRecHits(){
	double e, t, x, y, z, theta;
	double e_cell, t_cell;
	double e_sig, t_sig;
	double eta, phi;
	int nrhs = 0;
	int etot = 0;
	//do energy first to get transfer factor
	for(int i = 0; i < _netacal; i++){
		for(int j = 0; j < _nphical; j++){
			e = _cal[i][j].at(0);
			if(e == 0.) continue;
			//time in cell ieta, jphi is energy weighted average
			//do amplitude dependent energy smearing
			//that follows equation 1.2 in CMS TDR
			//(sig/E)^2 = (s/sqrt(E))^2 + (n/E)^2 + c^2
			e_sig = _s*_s/e + _n*_n/(e*e) + _c*_c;
			e_sig = e*sqrt(e_sig);
			//smear total energy in cell
			//update random sampling range to match
			//energy in this cell
			//out to 5 sigma
			_rs.SetRange(e - 5*e_sig, e + 5*e_sig);
			//smear energy in each cell
			e_cell = _rs.SampleGaussian(e, e_sig, 1).at(0); //returns a vector, take first (and only) element
			//make sure e can't be negative
			if(e_cell < 0) continue;	
			//in each cell to find energy deposited
			//if integral is below some threshold e_cell = 0
			if(e_cell < _ethresh) continue;

			t = _cal[i][j].at(1)/((double)_cal[i][j].at(2));
			
			//do amplitude dependent time smearing
			//with constant c = 200 ps
			t_sig = _calTres*_calTres;// + _n*_n/(e*e);
			t_sig = sqrt(t_sig);
			//smear time in cell
			//t can be negative (early times)
			//update range to be centered on t, up to 5 sigma (calTres)
			_rs.SetRange(t - 5*t_sig, t + 5*t_sig);
			t_cell = _rs.SampleGaussian(t, t_sig, 1).at(0);
			
			//reset e and t for cal cells
			etot += e_cell;
			nrhs++;			
			_cal[i][j].SetValue(e_cell, 0);
			_cal[i][j].SetValue(t_cell, 1);		
	
		}
	}
	if(_default_transfer) _gev = etot/(double)nrhs;
}



//loop through all the particles and reconstruct energy
//from the filled rec hits
//this needs to be done after all particles have been showered
//to account for overlap
void BasicDetectorSim::ReconstructEnergy(){
	XYZTVector Position;
	PtEtaPhiEVector Momentum;
	Pythia8::Particle Particle;
	double eta, phi, theta, ceta, cphi;
	int ieta, iphi, iieta, iiphi;	
	double reco_e, reco_t, reco_nrh;
	double x, y, z, t, r, q;
	for(int p = 0; p < _recops.size(); p++){
		Position = _recops[p].Position;
		Momentum = _recops[p].Momentum;
		Particle = _recops[p].Particle;
		eta = Position.eta();
		phi = Position.phi();
		q = _recops[p].Particle.charge();
		
		reco_e = 0;
		reco_t = 0;
		reco_nrh = 0;
		//get ieta, iphi
		_get_etaphi_idx(eta, phi, ieta, iphi);
		//loop through corresponding nxn grind
		for(int i = -_ncell; i < _ncell+1; i++){
			for(int j = -_ncell; j < _ncell+1; j++){
				//get ieta, iphi for cell in grid
				iieta = ieta+i;
				iiphi = iphi+j;

				//if these indices are out of bounds for detector indices, skip
				if(iieta >= _netacal || iieta < 0) continue;
				if(iiphi >= _nphical || iiphi < 0) continue;
				
				//also do zero suppression here so
				//this rh doesn't get counted in time average
				if(_cal[iieta][iiphi].at(0) < _ethresh) continue;
			
				//get integral bounds
				_get_etaphi(iieta, iiphi, ceta, cphi);
	
				//get x, y, z based on cell eta phi
				x = _rmax*cos(cphi);
				y = _rmax*sin(cphi);
				theta = 2*tan(exp(-ceta));
				z = sqrt(x*x + y*y)/tan(theta);
				t = _cal[iieta][iiphi].at(1); 

				//get cal energies for iieta and iiphi;
				reco_e += _cal[iieta][iiphi].at(0);		
				//get cal time for iieta and iiphi (this time is itself an E-weighted average)
				reco_t += _cal[iieta][iiphi].at(1)*_cal[iieta][iiphi].at(0);
				reco_nrh++;
		
				//add emission to particle
				//doesn't include tof correction to time
				//save time in ns, space in cm
				JetPoint jet(x*1e2, y*1e2, z*1e2, t*1e9);
				jet.SetEnergy(_cal[iieta][iiphi].at(0));
				jet.SetEta(ceta);
				jet.SetPhi(cphi);
				_recops[p].AddEmission(jet);
				_cal_rhs.push_back(jet);
	

			}
		}
		//reset e and t
		_recops[p].Momentum.SetE(reco_e);
		_recops[p].Position.SetCoordinates(_recops[p].Position.x()*1e2, _recops[p].Position.y()*1e2, _recops[p].Position.z()*1e2, reco_t/((double)reco_nrh)*1e9);


	}	


}



//get "Jets" for clustering
void BasicDetectorSim::GetRecHits(vector<Jet>& rhs){
	rhs.clear();
	for(int j = 0; j < _cal_rhs.size(); j++){
		_cal_rhs[j].SetWeight(_cal_rhs[j].E()/_gev);
		rhs.push_back(Jet(_cal_rhs[j]));
	}
}

void BasicDetectorSim::GetEmissions(vector<vector<JetPoint>>& ems){
	ems.clear();
	for(int i = 0; i < _recops.size(); i++){
		ems.push_back({});
		for(int j = 0; j < _recops[i].ems.size(); j++)
			ems[i].push_back(_recops[i].ems[j]);
	}

}

void BasicDetectorSim::GetParticlesMom(vector<PtEtaPhiEVector>& genps, vector<PtEtaPhiEVector>& recops){
	genps.clear();
	recops.clear();
	PtEtaPhiEVector genvec;
	for(int i = 0; i < _recops.size(); i++){
		recops.push_back(_recops[i].Momentum);
		genvec.SetCoordinates(_recops[i].Particle.pT(), _recops[i].Particle.eta(), _recops[i].Particle.phi(), _recops[i].Particle.e());
		genps.push_back(genvec);
	}
}

void BasicDetectorSim::GetParticlesPos(vector<XYZTVector>& genps, vector<XYZTVector>& recops){
	genps.clear();
	recops.clear();
	XYZTVector genvec;
	for(int i = 0; i < _recops.size(); i++){
		recops.push_back(_recops[i].Position);
		//genvec.SetCoordinates(_recops[i].particle.xProd(), _recops[i].particle.yProd(), _recops[i].particle.zProd(), _recops[i].particle.tProd());
		genvec.SetCoordinates(_recops[i].Particle.px(), _recops[i].Particle.py(), _recops[i].Particle.pz(), _recops[i].Particle.e());
		genps.push_back(genvec);
	}
}



void BasicDetectorSim::GetTrueJets(vector<Jet>& jets){
	jets.clear();
	for(int i = 0; i < _jets.size(); i++)
		jets.push_back(Jet( _jets[i].px(), _jets[i].py(), _jets[i].pz(), _jets[i].e()) );

}

bool BasicDetectorSim::_in_cell_crack(const RecoParticle& rp){
	Pythia8::Particle p = rp.Particle;
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
	//offset index for indexing starting at 0
	//offset phi for [-pi, pi] or [0, 2pi]
	//iphi = int(((phi - _phimin)/_dphi)+0.5) - 1;
	iphi = floor(((phi - _phimin)/_dphi));
	if(iphi > _nphical) iphi -= _nphical;
	else if(iphi < 0) iphi += _nphical; 

	//get eta index
	ieta = 0;
	if(eta > _etamin && eta < _etamax){
		//offset for indexing starting at 0
		//ieta = int(((eta - _etamin)/_deta) + 0.5) - 1; 
		ieta = floor(((eta - _etamin)/_deta)); 
		if(ieta > _netacal) ieta = 0;
		else if(ieta < 0) ieta = 0;
      
	}
}

//get eta-phi center given ieta, iphi cell indices
void BasicDetectorSim::_get_etaphi(int ieta, int iphi, double& eta, double& phi){
 	double pi = acos(-1);
	//eta      
	eta = _etamin + _deta/2. + double(ieta)*_deta;      
	//phi 
	phi = _phimin + _dphi/2. + double(iphi)*_dphi;
}



