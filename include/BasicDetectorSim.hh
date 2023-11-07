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
#include "Jet.hh"
#include <Math/Vector4D.h>
#include <TH1D.h>
#include <string>

using std::string;
using Pythia = Pythia8::Pythia;
using Particle = Pythia8::Particle; 

using PtEtaPhiEVector = ROOT::Math::PtEtaPhiEVector;
using XYZTVector = ROOT::Math::XYZTVector;
class BasicDetectorSim{
struct RecoParticle;
	public:
		BasicDetectorSim();
		BasicDetectorSim(string infile);
		virtual ~BasicDetectorSim(){ };

		void SimQCD(); //use pythia to simulate QCD events
		void SimTTbar(); //use pythia to simulate ttbar events

		//this is what does the detector effects on the tracks
		void CalcTrajectory(RecoParticle& rp); //calculate trajectories/tracks for 
					// particles from PV to detector (cal) face
					// see pgs_track_res in PGS

		//this is what creates the showers from the reco particles
		void FillCal(RecoParticle& rp); // for energy depositions
				// see pgs_fill_cal in PGS

		void MakeRecHits(); //loops through filled ecal cells
		//does time and energy smearings
		
	
		void SimulateEvents(int evt = -1); //loop through pythia events to get gen level info

		//get cal rec hits
		void GetRecHits(vector<Jet>& rhs); 
		//get gen + reco momentum
		void GetParticlesMom(vector<PtEtaPhiEVector>& genps, vector<PtEtaPhiEVector>& recops);
		//get gen + reco position
		void GetParticlesPos(vector<XYZTVector>& genps, vector<XYZTVector>& recops);

		//sets transfer factor for rec hit weights
		void SetEnergyTransferFactor(double gev){ _gev = gev; _default_transfer = false; }
		void SetNEvents(int e){ _nevts = e; _pythia.readString("Main:numberOfEvents = "+std::to_string(_nevts)); }
		void SetVerbosity(int v = 0){
			if(v == 0)
				_pythia.readString("Print:quiet = on");
			else if(v == 1){
				//_pythia.readString("Print:next = off");
				_pythia.readString("Next:numberShowEvent = 0");
				_pythia.readString("Next:numberShowLHA = 0");
				_pythia.readString("Next:numberShowInfo = 0");
				_pythia.readString("Next:numberShowProcess = 0");
			}
			_pythia.readString("Print:errors = "+std::to_string(v));
		}

	private:
		double _rmax; //max radius of detector (cm)
		double _b; //magnetic field (T)
		double _netacal; //number of cells in eta in calorimeter
		double _nphical; //number of cells in phi in calorimeter
		double _deta; //width of cell in eta	
		double _dphi; //width of cell in phi
		double _etamin; //min eta	
		double _etamax; //max eta
		double _calEres; //calorimeter energy resolution
		double _calTres; //calorimeter time resolution
		double _sagres; //sagitta resolution
		double _crack_frac; //calorimeter cell edge crack fraction
		
		//terms for ECAL energy resolution
		double _s = 0.0363; //% stochastic term
		double _n = 124*1e-3; //(124 MeV), noise term
		double _c = 0.0026; //& constant term
	
		//speed of light in (cm/ns)
		constexpr static double _sol = 29.9792458;
	
		RandomSample _rs; //random sampler for smearing
		int _nevts; //number of events to simulate
		vector<vector<double>> _calE; //energy in ecal cells in [eta][phi]
		vector<vector<double>> _calT; //time in ecal cells in [eta][phi]
		vector<vector<int>>    _calN; //number of emissions in an [eta][phi] cell
		vector<vector<Point>>  _cal; //3-dim point where each point is (e, t, n) for individual emissions in [eta][phi] cell
		vector<JetPoint> _cal_rhs; //ecal rec hits

		vector<Particle> _genps; //gen particles
		vector<RecoParticle> _recops; //reco particles
		Pythia8::Pythia _pythia; //pythia object for event generation


		double _gev; //transfer factor for weighting rechits
		bool _default_transfer;

		bool _in_cell_crack(const RecoParticle& rp); //check if particle's current four vector means it is in between cells
		void _get_etaphi_idx(double eta, double phi, int& ieta, int& iphi);
		void _get_etaphi(int ieta, int iphi, double& eta, double& phi);


		struct RecoParticle{
			//associated gen particle
			Particle particle;
			//associated emissions
			vector<JetPoint> ems;
			//reco position and momentum vectors
			PtEtaPhiEVector Momentum;
			XYZTVector Position;
			
			//ctor from gen particle
			RecoParticle(Particle& p){ 
				particle = p;
				//init momentum to gen momentum
				Momentum.SetCoordinates(p.pT(), p.eta(), p.phi(), p.e());
				//init to production coordinates
				//convert units to cm and ns
				//originally in mm
				Position.SetCoordinates(p.xProd() * 1e-1, p.yProd() * 1e-1, p.zProd() * 1e-1, p.tProd() * 1e-1 / _sol );
				
			}
			void AddEmission(JetPoint& j){ ems.push_back(j); }				
		
		};

};
#endif
