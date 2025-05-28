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
using std::set;
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
		void FillTracks(RecoParticle& rp);

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


		//set resolution constants - 
		void SetTimeResCts(double cte, double stoch, double noise){ _calTresCte = cte; _calTresStoch = stoch; _calTresNoise = noise; }
		//set energy threshold for zero suppression
		void SetEnergyThreshold(double e){ _ethresh = e; }
		//set min gen particle pt to be included in gen AK4 jet
		double _genpart_minpt;
		void SetMinGenPartPt(double p){_genpart_minpt = p; cout << "Minimum gen particle pt for subcluster analysis " << _genpart_minpt << " GeV" << endl; }



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

		//write gen jet info
		void FillGenJets();
		//write reco jet info
		void FillRecoJets();
		//write gen particle info
		void FillGenParticles();

		//set energy smear constant
		void SetEnergySmear(double c){ _c = c; }



	protected:
		void RecordMomInfo(Pythia8::Particle particle){
			//mother idx in gen particle list - assuming normal mother case where mother1 > 0 && mother2 == 0
			//'top level' gen particles are tops - if top do not save mother info
			int momidx_pythia = particle.mother1();
			if(fabs(particle.id()) != 6)
				_genmoms.insert(momidx_pythia);
		}
		void SaveGenInfo(Pythia8::Particle particle, int genmom){
			RecordMomInfo(particle);

			//assume that gen particle has not been propagated to detector face
			RecoParticle genpart(particle);
			CalcTrajectory(genpart);
			fastjet::PseudoJet fj_genpart( genpart.Momentum.px(), genpart.Momentum.py(), genpart.Momentum.pz(), genpart.Momentum.e() );
			fj_genpart.set_user_index(_genparts.size());
			_genparts.push_back(fj_genpart);
			_genpartIdx.push_back(_genparts.size()-1);
			_genpartids.push_back((int)genpart.Particle.id());
			//save top specific info
			if(fabs(genpart.Particle.id()) == 6){
				_topPt.push_back(genpart.Momentum.pt());
			}
			//cout << "saving gen mom " << genmom << endl;
			_genpartMomIdx.push_back(genmom);
		}

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
		double _calTresNoise; //calorimeter time resolution - energy squared dependent term
		double _calTresStoch; //calorimeter time resolution - energy dependent term
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
	
		vector<fastjet::PseudoJet>  _genAK4jets; //gen AK4 outputs from fastjet
		int _ngenAK4jets;
		vector<fastjet::PseudoJet>  _genAK8jets; //gen AK8 outputs from fastjet
		int _ngenAK8jets;
		vector<fastjet::PseudoJet>  _genAK15jets; //gen AK15 outputs from fastjet
		int _ngenAK15jets;
			
	
		vector<fastjet::PseudoJet>  _recoAK4jets; //reco outputs from fastjet
		int _nrecoAK4jets;
		vector<fastjet::PseudoJet>  _recoAK8jets; //reco outputs from fastjet
		int _nrecoAK8jets;
		vector<fastjet::PseudoJet>  _recoAK15jets; //reco outputs from fastjet
		int _nrecoAK15jets;

		vector<fastjet::PseudoJet>  _genparts; //gen particles
		vector<int> _genpartids; //genpart ids
		int _ngenparts; //# of objects
		fastjet::JetDefinition _jetdef_AK4; //fastjet clustering definition 
		fastjet::JetDefinition _jetdef_AK8; //fastjet clustering definition 
		fastjet::JetDefinition _jetdef_AK15; //fastjet clustering definition 
		double _Rparam; //default Rparam for default jets (not fat jets)
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

		//cluster sequences for gen and reco jets
		fastjet::ClusterSequence _gencsAK4;
		fastjet::ClusterSequence _gencsAK8;
		fastjet::ClusterSequence _gencsAK15;

		fastjet::ClusterSequence _recocsAK4;
		fastjet::ClusterSequence _recocsAK8;
		fastjet::ClusterSequence _recocsAK15;
		

		//for writing
		TTree* _tree = nullptr;
		std::unique_ptr<TFile> _file;
		void _reset();
		vector<double> _rhE, _rhx, _rhy, _rhz, _rht, _rheta, _rhphi;
		vector<unsigned int> _rhids;
	
		//gen AK4 jets 
		vector<double> _jgAK4eta, _jgAK4phi, _jgAK4energy, _jgAK4pt, _jgAK4mass, _jgAK4pz;
		//# particles in gen AK4 jets
		vector<int> _jgAK4nparts;
		//indices of gen AK4 particles in gen AK4 jets
		vector<vector<int>> _jgAK4partIdxs;
		
		//gen AK8 jets 
		vector<double> _jgAK8eta, _jgAK8phi, _jgAK8energy, _jgAK8pt, _jgAK8mass, _jgAK8pz;
		//# particles in gen AK8 jets
		vector<int> _jgAK8nparts;
		//indices of gen AK8 particles in gen AK8 jets
		vector<vector<int>> _jgAK8partIdxs;
		
		//gen AK15 jets 
		vector<double> _jgAK15eta, _jgAK15phi, _jgAK15energy, _jgAK15pt, _jgAK15mass, _jgAK15pz;
		//# particles in gen AK15 jets
		vector<int> _jgAK15nparts;
		//indices of gen AK15 particles in gen AK15 jets
		vector<vector<int>> _jgAK15partIdxs;
		
		
		//gen top info
		vector<double> _topPt_had, _topPt_hadlep, _topPt_lep, _topPt;
		vector<int> _topDecayId;
		//gen particle info (even intermediate particles)
		vector<double> _genparteta, _genpartphi, _genpartenergy, _genpartpt, _genpartmass, _genpartpz;
		vector<int> _genpartMomIdx;
		set<double> _genmoms;
		vector<int> _genpartIdx,_genpartEvtIdx;
	
		//reco AK4 jets
		vector<double> _jAK4eta, _jAK4phi, _jAK4energy, _jAK4pt, _jAK4mass;
		vector<vector<unsigned int>> _jAK4rhids;
		
		//reco AK8 jets
		vector<double> _jAK8eta, _jAK8phi, _jAK8energy, _jAK8pt, _jAK8mass;
		vector<vector<unsigned int>> _jAK8rhids;
		
		//reco AK15 jets
		vector<double> _jAK15eta, _jAK15phi, _jAK15energy, _jAK15pt, _jAK15mass;
		vector<vector<unsigned int>> _jAK15rhids;
	
		//pv info
		double _pvx, _pvy, _pvz;
		//track info
		vector<double> _trackpx, _trackpy, _trackpz, _tracketa, _trackphi;
		//spike info	
		vector<double> _spikeE;
		//event info
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
				Position.SetCoordinates(p.xProd()*1e-3, p.yProd()*1e-3, p.zProd()*1e-3, p.tProd() / (_sol*1e3) );
			}
			void AddEmission(JetPoint& j){ ems.push_back(j); }				
		
		};

};
#endif
