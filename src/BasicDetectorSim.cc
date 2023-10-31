#include "BasicDetectorSim.hh"



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

	//sets pythia settings by given .cmnd file
	_pythia.readFile(infile);

}

//use pythia8 to simulate events
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




void BasicDetectorSim::SimulateEvents(){
	double maxeta = 1.749;
	double e = 0;
	//(nominal values of) constants for energy resolution
	//these are from CMS TDR Fig. 1.7
	//the more conservative values are taken (s.t. sig/E is larger than other values)
	double s = 3.63; //stochastic term
	double n = 124; //(MeV), noise term
	double c = 0.26; //constant term
	double e_sig = 1; //sigma of energy resolution smearing

	for(int i = 0; i < _nevts; i++){
		if(!_pythia.next()) continue;
		//loop through all particles
		//make sure to only record those that would
		//leave RecHits in ECAL (ie EM particles (ie ie photons and electrons))
		for(int p = 0; p < _pythia.event.size(); p++){
			Particle p = _pythia.event[p];
			//make sure particle is a photon or electron
			if(fabs(p.id()) != 11 && fabs(p.id()) != 22) continue;
			//make sure particle is final-state and (probably) stable
			if(p.statusHepMC() != 1) continue;
			//make sure particle is in detector acceptance
			//since this is a CMS ECAL sim, use CMS ECAL geometry
			if(fabs(p.eta()) > maxeta) continue;
			
			//could do min pt cut here
			//if i HAD one

			//do energy resolution smearing with a gaussian
			//centered on nominal energy and with an energy-dependent spread 
			//that follows equation 1.2 in CMS TDR
			//(sig/E)^2 = (s/sqrt(E))^2 + (n/E)^2 + c^2
			e_sig = sqrt( s*s*sqrt(p.eCalc()) + n*n + c*c*p.eCalc()*p.eCalc() );
			e = _rs.SampleGaussian(p.eCalc(), e_sig, 1).at(0); //returns a vector, take first (and only) element
			
			//save particle four vector
			
		
	
		}

	}



}
