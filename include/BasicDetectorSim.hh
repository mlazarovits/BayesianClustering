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
using PxPyPzEVector = ROOT::Math::PxPyPzEVector;
using XYZTVector = ROOT::Math::XYZTVector;
class BasicDetectorSim{
struct RecoParticle;
	public:
		BasicDetectorSim();
		BasicDetectorSim(string infile);
		virtual ~BasicDetectorSim(){ 

		};

		void SimTTbar(){ _procs_to_sim.push_back(ttbar); _simttbar = true;};
		void SimQCD(){ _procs_to_sim.push_back(qcd); _simqcd = true;};
		void SimSingleW(){ _procs_to_sim.push_back(w); _simw = true;};

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


		void TurnOnPileup(int npuavg = 5,bool oot = false){
			_pu = true;
			_nPUavg = npuavg; //number of pu events on average (used for poisson sampling)
			_oot = oot;
		}

		void TurnOnSpikes(double sprob = 0.01){ _spikes = true; _spikeprob = sprob; }


		//set resolution constants - 
		void SetTimeResCts(double cte, double stoch, double noise){ _calTresCte = cte; _calTresStoch = stoch; _calTresNoise = noise; }
		//set energy threshold for zero suppression
		void SetEnergyThreshold(double e){ _ethresh = e; }
		//set min gen particle pt to be included in gen AK4 jet
		double _genpart_minpt;
		void SetMinGenPartPt(double p){_genpart_minpt = p; cout << "Minimum gen particle pt for subcluster analysis " << _genpart_minpt << " GeV" << endl; }

		void SetPtHatMin(double p){ _ptHatMin = p; cout << "ptHat minimum = " << _ptHatMin << endl;}

	
		//get cal rec hits
		void GetRecHits(vector<Jet>& rhs); 
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

		void TurnOffShower(){
			_noShower = true;
		}


	protected:
		//void SaveGenInfo(Pythia8::Particle particle, int genmom){
		void SaveGenInfo(int evtidx, int genmom){
			Pythia8::Particle particle = _sumEvent[evtidx];
			_genpartEvtIdx.push_back(evtidx);

			//assume that gen particle has not been propagated to detector face
			RecoParticle genpart(particle);
			//cout << "pre CalcTraj - SaveGenInfo - saving gen part with eta " << genpart.Momentum.eta() << " " << particle.eta() << " phi " << genpart.Momentum.phi() << " " << particle.phi() << " energy " << genpart.Momentum.e() << " " << particle.e() << " pz " << genpart.Momentum.pz() << " " << particle.pz() << endl;
			//cout << "position - pre CalcTraj - SaveGenInfo - saving gen part with eta " << genpart.Position.eta() << " " << particle.eta() << " phi " << genpart.Position.phi() << " " << particle.phi() << " energy " << genpart.Position.e() << " " << particle.e() << " pz " << genpart.Position.pz() << " " << particle.pz() << endl;
			CalcTrajectory(genpart);
			fastjet::PseudoJet fj_genpart( genpart.Momentum.px(), genpart.Momentum.py(), genpart.Momentum.pz(), genpart.Momentum.e() );
			
			//cout << "position - post CalcTraj - SaveGenInfo - saving gen part with eta " << genpart.Position.eta() << " " << particle.eta() << " phi " << genpart.Position.phi() << " " << particle.phi() << " energy " << genpart.Position.e() << " " << particle.e() << " pz " << genpart.Position.pz() << " " << particle.pz() << endl;

			//cout << "post CalcTraj - SaveGenInfo - saving gen part with eta " << fj_genpart.eta() << " " <<  genpart.Momentum.eta() << " " << particle.eta() << " phi " << fj_genpart.phi() << " " << genpart.Momentum.phi() << " " << particle.phi() << " energy " << fj_genpart.e() << " " << genpart.Momentum.e() << " " << particle.e() << " pz " << fj_genpart.pz() << " " << genpart.Momentum.pz() << " " << particle.pz() << endl;


			fj_genpart.set_user_index(_genparts.size());
			//_genparts.push_back(fj_genpart);
			_genparts.push_back(genpart);
			_genpartIdx.push_back(_genparts.size()-1);
			_genpartids.push_back((int)genpart.Particle.id());
			//cout << "saving gen mom " << genmom << endl;
			_genpartMomIdx.push_back(genmom);
		}

		void FindMom(vector<int> mothers_idx, vector<int> mothers_id, int id, set<int>& idxs){
			vector<int>::iterator momit = find(mothers_id.begin(), mothers_id.end(), id);
			if(momit != mothers_id.end()){
				int idx = mothers_idx[momit - mothers_id.begin()];
				int momidx = _sumEvent[idx].mother1();
				while(_sumEvent[idx].mother1() == _sumEvent[idx].mother2() || fabs(_sumEvent[momidx].id()) == id){
					idx = _sumEvent[idx].mother1();
					momidx = _sumEvent[idx].mother1();
				}
				//cout << "mother of particle " << idx << ": " << momidx << " " << _sumEvent[idx].mother2() << " mother1 id " << _sumEvent[momidx].id() << " id " << _sumEvent[idx].id() << " grandmother " << _sumEvent[_sumEvent[momidx].mother1()].id() << " status " << _sumEvent[idx].status() << endl;
				idxs.insert(idx);
			}

		}


	private:
		Pythia8::Event _sumEvent; //one object where individual events are collected

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

		//vector<fastjet::PseudoJet>  _genparts; //gen particles
		vector<RecoParticle>  _genparts; //gen particles
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
		enum _proc {ttbar, qcd, w};
		vector<int> _procs_to_sim;
		void _simQCD(); //use pythia to simulate QCD events
		void _simTTbar(); //use pythia to simulate ttbar events
		void _simSingleW();
		bool _simttbar, _simqcd, _simw;

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
		vector<int> _topDecayId, _wDecayId;
		//gen particle info (even intermediate particles)
		vector<double> _genparteta, _genpartphi, _genpartenergy, _genpartpt, _genpartmass, _genpartpz;
		vector<int> _genpartMomIdx;
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
		double _pvx, _pvy, _pvz, _pvt;
		//pu pv info
		vector<double> _pu_pvx, _pu_pvy, _pu_pvz, _pu_pvt;
		bool _oot;
		//track info
		vector<double> _trackpx, _trackpy, _trackpz, _tracketa, _trackphi;
		//spike info	
		vector<double> _spikeE;
		//event info
		int _evt, _nRhs, _nSpikes, _nRecoParticles;
		BayesPoint _PV;

		//flag for turning off showering - saves final state gen particles 
		//(ie not the initial ones of the hard process saved in _genparts) 
		//as ECALRecHits to be used in clustering
		bool _noShower;
		int _evti, _evtj;

		double _ptHatMin;

		struct RecoParticle{
			//associated gen particle
			Pythia8::Particle Particle;
			//associated emissions
			vector<JetPoint> ems;
			//reco position and momentum vectors
			PtEtaPhiEVector Momentum;
			//PxPyPzEVector Momentum;
			XYZTVector Position;
			
			//ctor from gen particle
			RecoParticle(Pythia8::Particle& p){ 
				Particle = p;
				//init momentum to gen momentum
				Momentum.SetCoordinates(p.pT(), p.eta(), p.phi(), p.e());
				//Momentum.SetCoordinates(p.px(), p.py(), p.pz(), p.e());
				//init to production coordinates
				//originally in mm
				//convert to m and s
				Position.SetCoordinates(p.xProd()*1e-3, p.yProd()*1e-3, p.zProd()*1e-3, p.tProd() / (_sol*1e3) );
			}
			void AddEmission(JetPoint& j){ ems.push_back(j); }				
		
		};

};
#endif
