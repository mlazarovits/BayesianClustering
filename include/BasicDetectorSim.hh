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
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include <Math/Vector4D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <string>
#include <memory>

using std::string;
using Pythia = Pythia8::Pythia;

using PtEtaPhiEVector = ROOT::Math::PtEtaPhiEVector;
using XYZTVector = ROOT::Math::XYZTVector;
class BasicDetectorSim{
struct RecoParticle;
	public:
		BasicDetectorSim();
		BasicDetectorSim(string infile);
		virtual ~BasicDetectorSim(){ 

		};

		void SimTTbar(){ _procs_to_sim.push_back(ttbar); };
		void SimQCD(){ _procs_to_sim.push_back(qcd); };

		//this is what does the detector effects on the tracks
		void CalcTrajectory(RecoParticle& rp); //calculate trajectories/tracks for 
					// particles from PV to detector (cal) face
					// see pgs_track_res in PGS

		//this add "track" information to ntuples
		//not real tracks because there is no tracker in this sim :)
		//this is just the gen momentum information at the detector face
		//you could also smear this information here OR in CalcTrajectory where you would smear the radius of curvature (see PGS)
		void SaveTracks(RecoParticle& rp);

		//this is what creates the showers from the reco particles
		void FillCal(RecoParticle& rp); // for energy depositions
				// see pgs_fill_cal in PGS

		
		void MakeRecHits(); //loops through filled ecal cells
		//does time and energy smearings
	

		void ReconstructEnergy(); //reconstruct e and t from smeared and summed (overlapping showers) rhs	
	
		void SimulateEvents(int evt = -1); //loop through pythia events to get gen level info


		void TurnOnPileup(int npuavg = 5){
			_pu = true;
			_nPUavg = npuavg; //number of pu events on average (used for poisson sampling)
		}

		void TurnOnSpikes(double sprob = 0.01){ _spikes = true; _spikeprob = sprob; }


		//set resolution constants
		void SetTimeResCts(double cte, double rate){ _calTresCte = cte; _calTresRate = rate; }
		//set energy threshold for zero suppression
		void SetEnergyThreshold(double e){ _ethresh = e; }

		//get cal rec hits
		void GetRecHits(vector<Jet>& rhs); 
		//get gen + reco momentum
		void GetParticlesMom(vector<PtEtaPhiEVector>& genps, vector<PtEtaPhiEVector>& recops);
		//get gen + reco position
		void GetParticlesPos(vector<XYZTVector>& genps, vector<XYZTVector>& recops);
		//get emissions per reco particle
		void GetEmissions(vector<vector<JetPoint>>& ems);
		//get fastjet jets - run on gen particles
		void GetTrueJets(vector<Jet>& jets);

		void SetNEvents(int e){ _nevts = e; _pythia.settings.readString("Main:numberOfEvents = "+std::to_string(_nevts)); }
		void SetVerbosity(int v = 0){
			_verb = v;
			if(v == 0)
				_pythia.settings.readString("Print:quiet = on");
			else if(v == 1){
				//_pythia.settings.readString("Print:next = off");
				_pythia.settings.readString("Next:numberShowEvent = 0");
				_pythia.settings.readString("Next:numberShowLHA = 0");
				_pythia.settings.readString("Next:numberShowInfo = 0");
				_pythia.settings.readString("Next:numberShowProcess = 0");
			}
		}


		void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj; }	

		//init tree
		void InitTree(string fname);
		void WriteTree();

		//write gen info
		void FillGenJets();
		//write reco info
		void FillRecoJets();

		//set energy smear constant
		void SetEnergySmear(double c){ _c = c; }
	private:
		double _rmax; //max radius of detector (m)
		double _b; //magnetic field (T)
		double _netacal; //number of cells in eta in calorimeter
		double _nphical; //number of cells in phi in calorimeter
		double _deta; //width of cell in eta	
		double _dphi; //width of cell in phi
		double _etamin; //min eta	
		double _etamax; //max eta
		double _phimin; //sets [-pi, pi] or [0, 2pi] for phi indexing
		double _calEres; //calorimeter energy resolution
		double _calTresCte; //calorimeter time resolution - constant term
		double _calTresRate; //calorimeter time resolution - energy dependent term
		double _sagres; //sagitta resolution
		double _crack_frac; //calorimeter cell edge crack fraction
		
		//terms for ECAL energy resolution
		//(nominal values of) constants for energy resolution
		//these are from CMS TDR Fig. 1.7
		//the more conservative values are taken (s.t. sig/E is larger than other values)
		double _s = 3.63; //% stochastic term
		double _n = 124*1e-3; //(124 MeV), noise term
		double _c = 0.26; //& constant term
	
		//2*_ncell + 1 = n
		//nxn grid of reconstructing energy
		int _ncell;
		//energy threshold for reconstruction
		double _ethresh; //(GeV)
		
		//speed of light in (m/s)
		constexpr static double _sol = 2.99792458e8; //cm/ns = 29.9792458;
	
		RandomSample _rs; //random sampler for smearing
		int _nevts; //number of events to simulate
		//this needs to be separate from JetPoint bc there is no field in JetPoint to track how many emissions are in a cell
		vector<vector<BayesPoint>>  _cal; //3-dim point where each point is (e, t, n) for individual emissions in [eta][phi] cell
		//vector<JetPoint> _cal_rhs; //ecal rec hits
		vector<Jet> _cal_rhs; //ecal rec hits
		vector<fastjet::PseudoJet>  _jets; //outputs from fastjet
		vector<fastjet::PseudoJet>  _jetsReco; //outputs from fastjet
		fastjet::JetDefinition _jetdef; //fastjet clustering definition 
		double _Rparam;
		fastjet::Strategy _strategy; //fastjet clustering strategy
		fastjet::RecombinationScheme _recomb; //fastjet recombination strategy

		vector<RecoParticle> _recops; //reco particles
		Pythia8::Pythia _pythia; //pythia object for main event generation
		bool _pu; //pileup switch
		int _nPUavg; //number of pu events on average
		bool _spikes; //spikes switch
		double _spikeprob; //probability of spiking in rec hit

		//process enums
		enum _proc {ttbar, qcd};
		vector<int> _procs_to_sim;
		void _simQCD(); //use pythia to simulate QCD events
		void _simTTbar(); //use pythia to simulate ttbar events


		bool _in_cell_crack(const RecoParticle& rp); //check if particle's current four vector means it is in between cells
		void _get_etaphi_idx(double eta, double phi, int& ieta, int& iphi);
		void _get_etaphi(int ieta, int iphi, double& eta, double& phi);

		int _verb;

		//for writing
		TTree* _tree = nullptr;
		std::unique_ptr<TFile> _file;
		void _reset();
		vector<double> _rhE, _rhx, _rhy, _rhz, _rht, _rheta, _rhphi;
		//gen jets 
		vector<double> _jgeta, _jgphi, _jgenergy, _jgpt, _jgmass;
		//reco jets
		vector<double> _jeta, _jphi, _jenergy, _jpt, _jmass;
		//pv info
		double _pvx, _pvy, _pvz;
		//track info
		vector<double> _trackpx, _trackpy, _trackpz, _tracketa, _trackphi;
		int _npredjets, _ntruejets;
		vector<double> _predjeteta, _predjetphi, _predjetpt, _predjetmass, _predjetnparts;
		vector<double> _truejeteta, _truejetphi, _truejetpt, _truejetmass, _truejetnparts;
		vector<double> _spikeE;
		int _evt, _nRhs, _nSpikes, _nRecoParticles;
		BayesPoint _PV;

		int _evti, _evtj;

		struct RecoParticle{
			//associated gen particle
			Pythia8::Particle Particle;
			//associated emissions
			vector<JetPoint> ems;
			//reco position and momentum vectors
			PtEtaPhiEVector Momentum;
			XYZTVector Position;
			
			//ctor from gen particle
			RecoParticle(Pythia8::Particle& p){ 
				Particle = p;
				//init momentum to gen momentum
				Momentum.SetCoordinates(p.pT(), p.eta(), p.phi(), p.e());
				//init to production coordinates
				//originally in mm
				//convert to m and s
				Position.SetCoordinates(p.xProd()*1e-3, p.yProd()*1e-3, p.zProd()*1e-3, p.tProd() / (_sol*1e-3) );
			}
			void AddEmission(JetPoint& j){ ems.push_back(j); }				
		
		};

};
#endif
