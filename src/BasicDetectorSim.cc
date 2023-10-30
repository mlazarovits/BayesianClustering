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
