//----------------------------------------------------------------------
// This class is based on the PGS detector simulation by John Conway.
// https://conway.physics.ucdavis.edu/research/software/pgs/pgs4-general.htm
// A basic detector simulation with calorimeter only. Calorimeter default parameters based on
// CMS ECAL geometry.
//----------------------------------------------------------------------

#ifndef BasicDetectorSim_HH
#define BasicDetectorSim_HH

#include "RandomSample.hh"
#include "Pythia8/Pythia.h"
#include <string>

using std::string;
//using Pythia = Pythia8::Pythia;

class BasicDetectorSim{
	public:
		BasicDetectorSim();
		BasicDetectorSim(string infile);
		virtual ~BasicDetectorSim(){ };

		void SimQCD(); //use pythia to simulate QCD events
		void SimTTbar(); //use pythia to simulate ttbar events

		void calcTrajectories(); //calculate trajectories/tracks for 
					// particles from PV to detector (cal) face
		void fillCal(); // for energy depositions
	
		// see pgs_find* for finding/reconstructing specific physics objects


	private:
		double _rmax; //max radius of detector (m)
		double _b; //magnetic field (T)
		double _netacal; //number of cells in eta in calorimeter
		double _nphical; //number of cells in phi in calorimeter
		double _deta; //width of cell in eta	
		double _dphi; //width of cell in phi	
		double _calEres; //calorimeter energy resolution
		double _calTres; //calorimeter time resolution
		double _sagres; //sagitta resolution
		RandomSample _rs; //random sampler for smearing

		Pythia8::Pythia _pythia; //pythia object for event generation






};
#endif
