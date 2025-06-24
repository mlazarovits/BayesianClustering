#include "BasicDetectorSim.hh"
#include "TMath.h"
#include "Matrix.hh"
#include "fastjet/ClusterSequence.hh"
#include <TSystem.h>
#include <TFile.h>
#include <algorithm>

///////////////////////////
//all lengths are in meters (saving to JetPoints in cm)
//all times are in seconds (saving to JetPoints in ns)
//all energies, masses and momenta are in GeV (same for JetPoints)
///////////////////////////

BasicDetectorSim::BasicDetectorSim(){
	gSystem->Load("lib/libvecDict.so");
	//init parameters to CMS ECAL geometry
	//see CMS TDR (ISBN 978-92-9083-268-3)
	_rmax = 1.29; //inner radius of ECAL barrel (m)
	_b = 3.8; //field of CMS solenoid (T)
	_netacal = 344; //number of cells in eta in ECAL barrel (2*72)
	_nphical = 360; //number of cells in phi in ECAL barrel (neta*nphi = 61200 total cells in barrel)
	_deta = 2*acos(-1)/360.; //0.0174; //eta component of cell cross-section (2.2 cm - Moliere radius)
	_dphi = 2*acos(-1)/360.; //0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_calEres = 0.00445; //energy resolution approximated as radiation length/2. (rad length = 0.89 cm) 
	_calTresCte = 0.1727 * 1e-9;//0.2*1e-9; 
	_calTresNoise = 2.106 * 1e-9;//0.34641*1e-9; 
	_calTresStoch = 0.5109 * 1e-9;//1.60666*1e-9;
	_sagres = 0.000013; //value from LHC parameters in PGS (examples/par/lhc.par)
	_rs = RandomSample(); //random sampler
	_nevts = 1000;
	//initialize cal - save e, t, n emissions
	_etamax = 3.;//1.479 + _deta/2.; //puts outermost corner at true etamax
	_etamin = -_etamax;
	_phimin = -acos(-1);
	//energy threshold for reconstruction
	_ethresh = 0.1;//0.7; //(GeV)
	_ncell = 2;
	_pu = false; //pileup switch
	_spikes = false; //spikes switch
	_spikeprob = 0.01; //event-by-event probability of spike occurring
	//create container for cal cell info
	//E, t, nEmissions
	for(int i = 0; i < _netacal; i++){
		_cal.push_back({});
		for(int j = 0; j < _nphical; j++)
			_cal[i].push_back(BayesPoint({0.,0.,0.}));
	}
	_nSpikes = 0;
	_evti = 0;
	_evtj = 0;
	//fastjet objects - "AK4" jets
	_Rparam = 0.4;
	_strategy = fastjet::Best;
	_recomb = fastjet::E_scheme;
	_jetdef_AK4 = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4, _recomb, _strategy); 
	_jetdef_AK8 = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8, _recomb, _strategy); 
	_jetdef_AK15 = fastjet::JetDefinition(fastjet::antikt_algorithm, 1.5, _recomb, _strategy); 
	_genpart_minpt = 0;

	//default PV is detector center
	_PV = BayesPoint({0.,0.,0.});
	_pvx = _PV.at(0);
	_pvy = _PV.at(1);
	_pvz = _PV.at(2);
	_pvt = 0;
	
	_oot_pvx = -999;
	_oot_pvy = -999;
	_oot_pvz = -999;
	_oot_pvt = -999;


	//set beam spot spread in z (mm) and time (mm/c)
	//t spread = 100 ps => 0.1 ns * 30 cm/ns * 1e1 mm/cm = 30 mm z spread
	_pythia.settings.readString("Beams:allowVertexSpread = on");
	_pythia.settings.readString("Beams:sigmaVertexZ = 30"); //given in mm
	_pythia.settings.readString("Beams:sigmaTime = "+std::to_string(30)); //sigmaTime is given in mm/c
	_pythia.settings.readString("Beams:maxDevVertex = 1");
	_pythia.settings.readString("Beams:maxDevTime = 1");	

	_simttbar = false;
	_simqcd = false;
	_simw = false;


	_ptHatMin = 200;
}

//ctor with input pythia cmnd file
BasicDetectorSim::BasicDetectorSim(string infile){
	gSystem->Load("lib/libvecDict.so");
	//init parameters to CMS ECAL geometry
	//see CMS TDR (ISBN 978-92-9083-268-3)
	_rmax = 1.29; //inner radius of ECAL barrel (m)
	_b = 3.8; //field of CMS solenoid (T)
	_netacal = 344; //number of cells in eta in ECAL barrel (2*172)
	_nphical = 360; //number of cells in phi in ECAL barrel (neta*nphi = 61200 total cells in barrel)
	_deta = 2*acos(-1)/360.; //0.0174; //eta component of cell cross-section (2.2 cm - Moliere radius)
	_dphi = 2*acos(-1)/360.; //0.0174; //phi component of cell cross-section (2.2 cm - Moliere radius)
	_calEres = 0.00445; //energy resolution approximated as radiation length/2. (rad length = 0.89 cm) 
	_calTresCte = 0.1727 * 1e-9;//0.2*1e-9; 
	_calTresNoise = 2.106 * 1e-9;//0.34641*1e-9; 
	_calTresStoch = 0.5109 * 1e-9;//1.60666*1e-9;
	_sagres = 0.000013; //value from LHC parameters in PGS (examples/par/lhc.par)
	_rs = RandomSample(); //random sampler
	_etamax = 3.;//1.479;
	_etamin = -_etamax;
	_phimin = -acos(-1);
	//energy threshold for reconstruction
	_ethresh = 0.1; //(GeV)
	_ncell = 2;
	_pu = false; //pileup switch
	_spikes = false; //spikes switch
	_spikeprob = 0.01; //event-by-event probability of spike occurring
	//initialize cal - save e, t, n emissions
	for(int i = 0; i < _netacal; i++){
		_cal.push_back({});
		for(int j = 0; j < _nphical; j++)
			_cal[i].push_back(BayesPoint({0.,0.,0.}));
	}

	//sets pythia settings by given .cmnd file
	_pythia.readFile(infile);
	_nevts = _pythia.mode("Main:numberOfEvents");
	_nSpikes = 0;
	_evti = 0;
	_evtj = 0;
	//fastjet objects - "AK4" jets
	_Rparam = 0.4;
	_strategy = fastjet::Best;
	_recomb = fastjet::E_scheme;
	_jetdef_AK4 = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4, _recomb, _strategy); 
	_jetdef_AK8 = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8, _recomb, _strategy); 
	_jetdef_AK15 = fastjet::JetDefinition(fastjet::antikt_algorithm, 1.5, _recomb, _strategy); 
	_genpart_minpt = 0;

	//default PV is detector center
	_PV = BayesPoint({0.,0.,0.});
	_pvx = _PV.at(0);
	_pvy = _PV.at(1);
	_pvz = _PV.at(2);
	_pvt = 0;
	
	_oot_pvx = -999;
	_oot_pvy = -999;
	_oot_pvz = -999;
	_oot_pvt = -999;
	
	_simttbar = false;
	_simqcd = false;
	_simw = false;
	
	_ptHatMin = 200;
}


//use pythia8 to simulate events
void BasicDetectorSim::_simQCD(){
	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 13 TeV.
	//_pythia.settings.readString("HardQCD:all = on");
	_pythia.settings.readString("HardQCD:gg2qqbar = on");
	_pythia.settings.readString("PhaseSpace:pTHatMin = "+std::to_string(_ptHatMin));
	_pythia.settings.readString("Beams:eCM = 13000.");
	if(_verb > 1) cout << "Simulating QCD" << endl;
}

//https://pythia.org/latest-manual/TopProcesses.html
void BasicDetectorSim::_simTTbar(){
	// Create Pythia instance and set it up to generate hard QCD processes
	// above pTHat = 20 GeV for pp collisions at 13 TeV.
	//ttbar specific (not tqbar production)
	_pythia.settings.readString("Top:gg2ttbar = on");
	_pythia.settings.readString("Top:qqbar2ttbar = on");
	_pythia.settings.readString("Top:ffbar2ttbar(s:gmZ) = on");
	_pythia.settings.readString("Top:ffbar2ttbar(s:gmZ) = on");
	_pythia.settings.readString("PhaseSpace:pTHatMin = "+std::to_string(_ptHatMin));
	_pythia.settings.readString("Beams:eCM = 13000.");
	if(_verb > 1) cout << "Simulating ttbar" << endl;

}

void BasicDetectorSim::_simSingleW(){
	//_pythia.settings.readString("WeakSingleBoson:ffbar2W = on");
	_pythia.settings.readString("WeakBosonAndParton:ffbar2Wgm = on");
	_pythia.settings.readString("PhaseSpace:pTHatMin = "+std::to_string(_ptHatMin));
	_pythia.settings.readString("Beams:eCM = 13000.");
	if(_verb > 1) cout << "Simulating single W" << endl;


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
	if(find(_procs_to_sim.begin(), _procs_to_sim.end(), w) != _procs_to_sim.end()){
		_simSingleW();
	}	
	if(_pu){
  		if(_verb == 0) pileup.settings.readString("Print:quiet = on");
		pileup.settings.readString("Random:setSeed = on");
  		pileup.settings.readString("Random:seed = 10"+std::to_string(_evti));
		pileup.settings.readString("SoftQCD:nonDiffractive = on"); //minbias events only
  		pileup.settings.readString("Beams:eCM = 13000.");
		//have PU be out-of-time by 25 ns
		//25 ns * 29.9792458 cm/ns * 1e1 mm/cm = 7500 mm
		pileup.settings.readString("Beams:offsetTime = -7494.8114");
		//for doing z-spread of PU beamspot
		pileup.settings.readString("Beams:allowVertexSpread = on");
		pileup.settings.readString("Beams:sigmaVertexZ = 30"); //given in mm
		pileup.settings.readString("Beams:sigmaTime = "+std::to_string(30)); //sigmaTime is given in mm/c
		pileup.settings.readString("Beams:maxDevVertex = 1");
		pileup.settings.readString("Beams:maxDevTime = 1");
		pileup.init();			
		cout << "Simulating pileup" << endl;
	}
	cout << "Using tres_cte = " << _calTresCte*1e9 << " ns and tres_stoch = " << _calTresStoch*1e9 << " ns and tres_noise " << _calTresNoise*1e9 << " ns" << endl;
	//set random number seed - 
	//The seed to be used, if setSeed is on.
	//A negative value gives the default seed,
	//a value 0 gives a random seed based on the time, and
	//a value between 1 and 900,000,000 a unique different random number sequence.
	_pythia.settings.readString("Random:setSeed = on");
	_pythia.settings.readString("Random:seed = "+std::to_string(_evti+1));	

	//list changed parameters
	//_pythia.settings.listChanged();
	_pythia.init();

	//sigma for z-smearing
	//beamspot spread is ~5 cm = 0.05 m to use with _sol
	double zig = 0.05/2.;
	double zshift, tnew;
	//calculate halflength from max eta
	double theta = 2*atan2(1,exp(_etamax));
	double zmax = _rmax/tan(theta); //[mm]

	
	//declare fjinput/output containers
	vector<fastjet::PseudoJet> fjinputs;
	
	if(_evti == _evtj){
		_evti = 0;
		_evtj = _nevts;
	}
	double evt_Etot = 0;
	int nhad = 0;
	int nlep = 0;
	int nsemilep = 0;
	for(int i = _evti; i < _evtj; i++){
		if(evt != -1){
			if(i != evt) continue;
		}
		if(!_pythia.next()) continue;
		//store event info if pileup is on
		//cout << std::setprecision(5) << "event #" << i << endl;
		_sumEvent = _pythia.event;
		cout << endl;
		cout << std::setprecision(5) << "event #" << i << " has " << _sumEvent.size() << " particles" << endl;
		_evt = i;		
		set<int> w_idxs, top_idxs, d_idxs, u_idxs, s_idxs;
		
		//set PV for event - look at first particle in record
		Pythia8::Particle evtRec = _sumEvent[1];//_sumEvent.back();
		_PV = BayesPoint({evtRec.xProd(), evtRec.yProd(), evtRec.zProd()});
		_pvx = _PV.at(0)*1e-1; //put in cm from mm
		_pvy = _PV.at(1)*1e-1; //put in cm from mm
		_pvz = _PV.at(2)*1e-1; //put in cm from mm
		_pvt = evtRec.tProd()*1e-3*(1/_sol)*1e9;	

		if(_pu){
			//simulate n pileup events
			int nPU = _rs.SamplePoisson(_nPUavg,1).at(0);
			for(int p = 0; p < nPU; p++){
				pileup.next();
				Pythia8::Event pu_event = pileup.event;
				//save PV info of OOT pu event
				_oot_pvx = pu_event[1].xProd()*1e-1;
				_oot_pvy = pu_event[1].yProd()*1e-1;
				_oot_pvz = pu_event[1].zProd()*1e-1;
				_oot_pvt = pu_event[1].tProd()*1e-3*(1/_sol)*1e9;
				_sumEvent += pu_event;
			}
		}
		//loop through all particles
		//make sure to only record those that would
		//leave RecHits in ECAL (ie EM particles (ie ie photons and electrons))
		//cout << "event size: " << _sumEvent.size() << endl;
		

		//set production vertex for this event from z-smearing
		//simulate z-shift from Gaussian
		//zig is nominal beam spot spread - should be 3 sigma for distribution
		//_rs.SetRange(-zig/3., zig/3.);
		//zshift = _rs.SampleGaussian(0., zig/3., 1).at(0);
		for(int p = 0; p < _sumEvent.size(); p++){
			//reset reco particle four momentum
			Pythia8::Particle particle = _sumEvent[p];
			//if simulating QCD, want to save the hard partons (status 23) produced
			//may want to separate partons from hard QCD dijets and pileup (diffractive) - in which case save particles in each event before PU is added
			//may want to save with a bool (pu, !pu)
			if(_simqcd){
				if(fabs(particle.status()) == 23){
					SaveGenInfo(p, -1);
				}
			}
			if(_simw){
				if(fabs(particle.id()) == 24 && fabs(particle.status()) == 22){
cout << "particle eta " << particle.eta() << " phi " << particle.phi() << endl;
					SaveGenInfo(p, -1);
				}
					
			}
			//make sure particle is final-state and (probably) stable
			if(particle.statusHepMC() != 1) continue;
			// No neutrinos
      			if (_sumEvent[p].idAbs() == 12 || _sumEvent[p].idAbs() == 14 ||
      			    _sumEvent[p].idAbs() == 16)     continue;
		

			//extreme gen momentum eta cut - for CMS tracker (not necessary here)	
			//if(fabs(particle.eta()) > 3.) continue;	

			//puT in pT cut hehe (charged particles only)
			//muon would need ~3.5 GeV to get to muon chambers so this should be the ceiling for the cut
			//shower as many gen particles as possible
			//if(particle.pT() < 3. && particle.charge() != 0) continue;
		
			//set time from linear model with slope = { z > 0 ? 1/sol : -1/sol }
			//tnew = (particle.zProd()*1e-3 +zshift) >= 0 ? (particle.zProd()*1e-3+zshift)*1./(_sol) : -1./(_sol)*(particle.zProd()*1e-3+zshift);
			//original pythia coords are in m, convert to mm
			//particle.vProd(particle.xProd(), particle.yProd(), particle.zProd()+(zshift*1e3), tnew);
			
			//make sure particle is in detector acceptance
			//since this is a CMS ECAL sim, use CMS ECAL geometry
			//this is for gen particles
			//include z offset
			//z and eta are one to one
			//if(fabs(particle.eta() + particle.zProd()*(_etamax/(zmax))) > _etamax) continue;
			//if(fabs(rp.Position.eta() + znew*(_etamax/(zmax))) > _etamax) continue;
		
			//zero suppression threshold
			//if(particle.e() < _ethresh) continue;

			if(particle.mother1() > 0 && particle.mother2() == 0){
				vector<int> mothers_id;
				vector<int> mothers_idx;
				int momidx = 999;
				int thisp = p;
				//cout << "mother search for particle " << p << ": " << _sumEvent[p].id() << endl;
				while(thisp > 0){
					momidx = _sumEvent[thisp].mother1();
					mothers_id.push_back(fabs(_sumEvent[momidx].id()));
					mothers_idx.push_back(momidx);
					//cout << "mother of " << thisp << " (id: " << _sumEvent[thisp].id() << ") is " << momidx << " (id: " << _sumEvent[momidx].id() << ")" << endl;
					thisp = momidx;
				}
				
				//get top from this particle's history (if it exists)
				//need to catch W's from tops - check if top is in mother chain - starting point (need to go all the way up mother chain to get original particle that hasn't recoiled, etc)
				FindMom(mothers_idx, mothers_id, 6, top_idxs);
				//find mothers that are Ws
				FindMom(mothers_idx, mothers_id, 24, w_idxs);
				//for(auto t = top_idxs.begin(); t != top_idxs.end(); t++)
				//	cout << "top quark idx test " << t->second << " status " << _sumEvent[t->second].status() << endl;
				//for(auto d = d_idxs.begin(); d != d_idxs.end(); d++)
				//	cout << "quark idx test " << *d << " status " << _sumEvent[*d].status() << " name " << _sumEvent[*d].name() << endl;
				//cout << "id of particle at idx 5 " << _sumEvent[5].id() << " at idx 6 " << _sumEvent[6].id() << " status " << _sumEvent[6].status() << endl;
				//cout << endl;	
			}

			//dont reconstruct muons - they would only mildly interact with an EM cal anyway
			if(particle.idAbs() == 13)
				continue;
			//double test_t = particle.tProd()/(_sol*1e3)*1e9;	
			//if(test_t > 1e2) cout << "particle energy " << particle.e() << " time (ns) " << test_t << endl;
			//create new particle for reco one
			RecoParticle rp(particle); 
			//calculate new pt (does full pvec but same pz)
			CalcTrajectory(rp);
			//check if in cal cell crack
			if(_in_cell_crack(rp))
				continue;
			//if track curls up/exceeds zmax
			if(fabs(rp.Position.z()) >= zmax || fabs(rp.Position.eta()) > _etamax) continue;
		
			//if gen particle doesn't exceed min pt, skip
			//if(rp.Momentum.pt() < _genpart_minpt) continue;
			//add gen particle to be clustered for gen jet
			//don't include electrons in gen jet clustering (or save to gen particles collection), muons are skipped above bc they are not showered
			//add particle to fastjet
			if(rp.Particle.idAbs() != 11){
				fastjet::PseudoJet fjinput( rp.Momentum.px(), rp.Momentum.py(), rp.Momentum.pz(), rp.Momentum.e() );
				fjinput.set_user_index(_genparts.size());
				fjinputs.push_back(fjinput);
			}

			//fill ecal cell with reco particle
			FillCal(rp);
	
			//save tracks (gen information of momentum propagated to detector)
			//don't save if particle doesn't make it to detector face (with if statement above)
			FillTracks(rp);

			//int tieta, tiphi;
			//evt_Etot += rp.Momentum.e();
			//_get_etaphi_idx(rp.Position.eta(), rp.Position.phi(), tieta, tiphi);
			//cout << "gen input px " << rp.Momentum.px() <<  " py " << rp.Momentum.py() << " pz " << rp.Momentum.pz() << " energy " << rp.Momentum.e() << " eta " << rp.Position.eta() << " tieta " << tieta << " phi " << rp.Position.phi() << " tiphi " << tiphi << " q " << rp.Particle.charge() << " pt " << rp.Momentum.pt() << endl; 
      			
			//save reco particle four vector (with corresponding gen info)
			_recops.push_back(rp);	
		}
		//cout << "event total energy " << evt_Etot << endl;
		evt_Etot = 0;
	

	
		//running fastjet on gen particles, no shower, etc.
		_gencsAK4 = fastjet::ClusterSequence(fjinputs, _jetdef_AK4);
		//get jets - min 5 pt
		vector<fastjet::PseudoJet> fjoutputs_AK4 = _gencsAK4.inclusive_jets(5.);
		for(int j = 0; j < fjoutputs_AK4.size(); j++) _genAK4jets.push_back(fjoutputs_AK4[j]);
		//sort jets by pt
		_genAK4jets = sorted_by_pt(_genAK4jets);
		
		_gencsAK8 = fastjet::ClusterSequence(fjinputs, _jetdef_AK8);
		//get jets - min 5 pt
		vector<fastjet::PseudoJet> fjoutputs_AK8 = _gencsAK8.inclusive_jets(5.);
		for(int j = 0; j < fjoutputs_AK8.size(); j++) _genAK8jets.push_back(fjoutputs_AK8[j]);
		//sort jets by pt
		_genAK8jets = sorted_by_pt(_genAK8jets);
		
		_gencsAK15 = fastjet::ClusterSequence(fjinputs, _jetdef_AK15);
		//get jets - min 5 pt
		vector<fastjet::PseudoJet> fjoutputs_AK15 = _gencsAK15.inclusive_jets(5.);
		for(int j = 0; j < fjoutputs_AK15.size(); j++) _genAK15jets.push_back(fjoutputs_AK15[j]);
		//sort jets by pt
		_genAK15jets = sorted_by_pt(_genAK15jets);

		//cout << fjoutputs.size() << " " << _jets.size() << " gen jets from " << fjinputs.size() << " inputs " <<  endl;
		//Fill gen jet information
		FillGenJets();	

	
		//make rhs and reconstruct time + energy for particles in evt	
		MakeRecHits();
		ReconstructEnergy();

		//fill reco jets after cells have been reconstructed
		FillRecoJets();


		//save info on top decays!
		vector<int> had = {1, 2, 3, 4, 5};
		vector<int> lep = {11, 12, 13, 14, 15, 16};
		if(_simw)
			cout << "w_idxs size " << w_idxs.size() << endl;
		//save gen w + direct daughters
		for(auto w = w_idxs.begin(); w != w_idxs.end(); w++){
			cout << "w status " << _sumEvent[*w].status() << " daughter1 " << _sumEvent[*w].daughter1() << " daughter2 " << _sumEvent[*w].daughter2() << endl;
			//get direct daughters of W
			//if daughter(s) is copy, reset W idx and daughters
			//continue until daughters are not W copy
			int cur_idx = *w;
			int d1 = _sumEvent[cur_idx].daughter1();
			int d2 = _sumEvent[cur_idx].daughter2();
			int cur_id = _sumEvent[cur_idx].id();
			while(fabs(_sumEvent[d1].id()) == 24){
				cur_idx = d1;
				d1 = _sumEvent[cur_idx].daughter1();
				d2 = _sumEvent[cur_idx].daughter2();
				cur_id = _sumEvent[cur_idx].id();
			}
			//save gen info of daughters
			//find w idx for mom in _genpartEvtIdx
			int widx = *w;
			vector<int>::iterator w_it = find(_genpartEvtIdx.begin(), _genpartEvtIdx.end(), widx);
			int momidx = std::distance(_genpartEvtIdx.begin(), w_it);
			SaveGenInfo(d1,momidx);
			SaveGenInfo(d2,momidx);
			int id_d1 = fabs(_sumEvent[d1].id());
			int id_d2 = fabs(_sumEvent[d2].id());
			bool lep_d1 = (find(lep.begin(), lep.end(), id_d1) != lep.end());
			bool lep_d2 = (find(lep.begin(), lep.end(), id_d2) != lep.end());

			bool had_d1 = (find(had.begin(), had.end(), id_d1) != had.end());
			bool had_d2 = (find(had.begin(), had.end(), id_d2) != had.end());

			//save w decay info
			if(lep_d1 && lep_d2){
				_wDecayId.push_back(1);
			}
			else if(had_d1 && had_d2){
				_wDecayId.push_back(0);
			}
			else{
				_wDecayId.push_back(-1);
			}
			cout << "final idx " << cur_idx << " d1 " << d1 << " d2 " << d2 << " final_id " << cur_id << " d1 id " << _sumEvent[d1].id() << " d2 id " << _sumEvent[d2].id() << " w decay id " << _wDecayId[_wDecayId.size()-1] << endl;
		}



		//get top decay gen info
		if(_simttbar)
			cout << "top_idxs size " << top_idxs.size() << endl;
		//sort top_idxs by top energy
		map<double, int, std::greater<double>> topE_idxs;
		for(auto t = top_idxs.begin(); t != top_idxs.end(); t++){
			topE_idxs[_sumEvent[*t].e()] = *t;
		}
		int nW = 0;
		//top gen info for tops + direct daughters + direct grand-daughters
		for(auto t = topE_idxs.begin(); t != topE_idxs.end(); t++){
			//save gen info for top
			SaveGenInfo(t->second,-1);
			//find children - going down decay chain
			vector<int> kids_id;
			vector<int> kids_idx = _sumEvent[t->second].daughterListRecursive();
			//cout << "top idx " << *t << " has " << kids_idx.size() << " daughters" << endl;
			for(auto k : kids_idx){ kids_id.push_back(fabs(_sumEvent[k].id()));}
			vector<int>::iterator w_it = find(kids_id.begin(), kids_id.end(), 24);
			vector<int>::iterator b_it = find(kids_id.begin(), kids_id.end(), 5);
		
			//save gen W info - W, idx to its mom, and daughter info
			int Widx = -999;
			int bidx = -999;
			if(w_it != kids_id.end()){
				Widx = std::distance(kids_id.begin(), w_it);
				Widx = kids_idx[Widx];
				//make sure W doesnt decay into a copy of itself, or the next decay isn't just a radiation 
				while(_sumEvent[Widx].daughter1() == _sumEvent[Widx].daughter2() || fabs(_sumEvent[_sumEvent[Widx].daughter1()].id()) == 24 || fabs(_sumEvent[_sumEvent[Widx].daughter2()].id()) == 24){
					//get daughters of W copy decay
					Widx = _sumEvent[Widx].daughter1();
				}
				//get W index
				//make sure W doesnt decay into a copy of itself, or the next decay isn't just a radiation 
				if(!(find(_procs_to_sim.begin(), _procs_to_sim.end(), qcd) != _procs_to_sim.end()))
					cout << "W idx " << Widx << " id " << _sumEvent[Widx].id() << " daughter1 " << _sumEvent[Widx].daughter1() << " daughter2 " << _sumEvent[Widx].daughter2() << endl;
				//save gen info of W at detector
				vector<int>::iterator t_it = find(_genpartEvtIdx.begin(),_genpartEvtIdx.end(),t->second);
				int genmomidx = std::distance(_genpartEvtIdx.begin(),t_it);
				SaveGenInfo(Widx,genmomidx);
				//cout << "W mom evt idx " << t->second << " gen mom idx " << genmomidx << " for particle idx " << _genpartIdx[_genpartIdx.size()-1] << " particle id " << _genpartids[_genpartids.size()-1] << endl;
				
			}
			//save gen b info
			if(b_it != kids_id.end()){
				//get b index
				bidx = std::distance(kids_id.begin(), b_it);
				bidx = kids_idx[bidx];
				//make sure b doesnt decay into a copy of itself, or the next decay isn't just a radiation 
				while(_sumEvent[bidx].daughter1() == _sumEvent[bidx].daughter2() || fabs(_sumEvent[_sumEvent[bidx].daughter1()].id()) == 5 || fabs(_sumEvent[_sumEvent[bidx].daughter2()].id()) == 5){
					//get daughters of b copy decay
					bidx = _sumEvent[bidx].daughter1();
				}
				if(!(find(_procs_to_sim.begin(), _procs_to_sim.end(), qcd) != _procs_to_sim.end()))
					cout << "b quark idx " << bidx << " id " << _sumEvent[bidx].id() << " daughter1 " << _sumEvent[bidx].daughter1() << " daughter2 " << _sumEvent[bidx].daughter2() << endl;
				//save gen info of b at detector
				vector<int>::iterator t_it = find(_genpartEvtIdx.begin(),_genpartEvtIdx.end(),t->second);
				int genmomidx = std::distance(_genpartEvtIdx.begin(),t_it);
				SaveGenInfo(bidx,genmomidx);
				//cout << "b mom evt idx " << t->second << " gen mom idx " << genmomidx << " for particle idx " << _genpartIdx[_genpartIdx.size()-1] << " particle id " << _genpartids[_genpartids.size()-1] << endl;
			}
			
			//save gen W daughters info
			if(Widx != -999){
				//cout << " daughters of W " << Widx << ": " << _sumEvent[_sumEvent[Widx].daughter1()].id() << " " << _sumEvent[_sumEvent[Widx].daughter2()].id() << endl;

				vector<int>::iterator Wmom_it = find(_genpartEvtIdx.begin(),_genpartEvtIdx.end(),Widx);
				int momidx = std::distance(_genpartEvtIdx.begin(),Wmom_it);
				//and save gen info of W daughters at detector
				//daughter 1
				SaveGenInfo(_sumEvent[Widx].daughter1(),momidx);	
				//cout << " W daughter 1 mom evt idx " << Widx << " gen mom idx " << momidx  << " for particle idx " << _genpartIdx[_genpartIdx.size()-1] << " particle id " << _genpartids[_genpartids.size()-1] << endl;
				//vector<int>::iterator d1_it_had = find(had.begin(),had.end(),fabs(_genpartids[_genpartids.size()-1]));	
				//vector<int>::iterator d1_it_lep = find(lep.begin(),lep.end(),fabs(_genpartids[_genpartids.size()-1]));
				//cout << "d1 id " << fabs(_genpartids[_genpartids.size()-1]) << " d1 had " << (d1_it_had != had.end()) << " d1 lep " << (d1_it_lep != lep.end()) << endl;
				//daughter 2
				SaveGenInfo(_sumEvent[Widx].daughter2(),momidx);	
				//cout << " W daughter 2 mom evt idx " << Widx << " gen mom idx " << momidx << " for particle idx " << _genpartIdx[_genpartIdx.size()-1] << " particle id " << _genpartids[_genpartids.size()-1] << " " << _sumEvent[_sumEvent[Widx].daughter2()].id() << endl;
				//vector<int>::iterator d2_it_had = find(had.begin(),had.end(),fabs(_genpartids[_genpartids.size()-1]));		
				//vector<int>::iterator d2_it_lep = find(lep.begin(),lep.end(),fabs(_genpartids[_genpartids.size()-1]));
				//cout << "d2 id " << fabs(_genpartids[_genpartids.size()-1]) << " d2 had " << (d2_it_had != had.end()) << " d2 lep " << (d2_it_lep != lep.end()) << endl;
				int w_d1 = fabs(_sumEvent[_sumEvent[Widx].daughter1()].id());	
				int w_d2 = fabs(_sumEvent[_sumEvent[Widx].daughter2()].id());	
		
				bool lep_d1 = (find(lep.begin(), lep.end(), w_d1) != lep.end());
				bool lep_d2 = (find(lep.begin(), lep.end(), w_d2) != lep.end());

				bool had_d1 = (find(had.begin(), had.end(), w_d1) != had.end());
				bool had_d2 = (find(had.begin(), had.end(), w_d2) != had.end());
				
				//save info on gen top decays
				//0  : fully hadronic (to light quarks)
				//1  : fully leptonic (no taus, including neutrinos)
				if(lep_d1 && lep_d2){
					_topDecayId.push_back(1);
					_wDecayId.push_back(1);
				}
				else if(had_d1 && had_d2){
					_topDecayId.push_back(0);
					_wDecayId.push_back(0);
				}
				else{
					_topDecayId.push_back(-1);
					_wDecayId.push_back(-1);
				}
				//cout << "top decay id " << _topDecayId[_topDecayId.size()-1] << " d1 id " << w_d1 << " d2 id " << w_d2 << " had1 " << had_d1 << " had_d2 " << had_d2 << " lep_d1 " << lep_d1 << " lep_d2 " << lep_d2 << endl;	
				nW++;

			}
			else{ cout << "W not found" << endl; }
			//save gen b daughter info 
			vector<int> bkids_idx;
			if(bidx != -999){
				cout << " daughters of b " << bidx << ": " << _sumEvent[_sumEvent[bidx].daughter1()].id() << " " << _sumEvent[_sumEvent[bidx].daughter2()].id() << endl;
		
				//and save gen info of b daughters at detector
				vector<int>::iterator bmom_it = find(_genpartEvtIdx.begin(),_genpartEvtIdx.end(),bidx);
				int momidx = std::distance(_genpartEvtIdx.begin(),bmom_it);
				//and save gen info of W daughters at detector
				SaveGenInfo(_sumEvent[bidx].daughter1(),momidx);	
				//cout << " b daughter 1 mom evt idx " << bidx << " gen mom idx " << momidx << " for particle idx " << _genpartIdx[_genpartIdx.size()-1] << " particle id " << _genpartids[_genpartids.size()-1] << endl;
				SaveGenInfo(_sumEvent[bidx].daughter2(),momidx);	
				//cout << " b daughter 2 mom evt idx " << bidx << " gen mom idx " << momidx << " for particle idx " << _genpartIdx[_genpartIdx.size()-1] << " particle id " << _genpartids[_genpartids.size()-1] << endl;
	
	
				

			}
			else{ cout << "b not found" << endl; }
			//for(int g = 0; g < _genparts.size(); g++) cout << "gen part evt idx " << _genpartEvtIdx[g] << " gen part idx " << _genpartIdx[g] << " gen part id " << _genpartids[g] << " genpartmomidx #" << g << ": " << _genpartMomIdx[g] << endl;

			//cout << endl;
		}
		//fill gen particles
		FillGenParticles();
		
		if(_tree != nullptr) _tree->Fill();
		//reset event
		_sumEvent.clear();
		_reset();
		//if only simulating one event, save for GetRecHits
		if(evt != -1) _cal_rhs.clear();
		// Reset Fastjet input
	
		fjinputs.resize(0);
		fjoutputs_AK4.clear();
		fjoutputs_AK4.resize(0);
		fjoutputs_AK8.clear();
		fjoutputs_AK8.resize(0);
		fjoutputs_AK15.clear();
		fjoutputs_AK15.resize(0);
		//cout << endl;
	}
//cout << "nhad " << nhad << " nlep " << nlep << " nsemilep " << nsemilep << endl;
}

//updates pvec of p
//this code is based on ParticlePropagator class in Delphes
//https://github.com/delphes/delphes/blob/master/modules/ParticlePropagator.cc
void BasicDetectorSim::CalcTrajectory(RecoParticle& rp){
	ROOT::Math::PtEtaPhiEVector Momentum = rp.Momentum;
	ROOT::Math::XYZTVector Position = rp.Position;
	//calculate halflength from max eta
	double theta = 2*atan2(1,exp(_etamax));
	double halfLength = _rmax/tan(theta); //[m]

	//pythia units are in mm (or mm/c for time, natural units)
	//already converted to m in RecoParticle ctor (or m/c where c is in m/s)
	//initialized to production coordinates
	double x = Position.x();
	double y = Position.y();
	double z = Position.z();
                           
	double q = rp.Particle.charge();

	//eta check is done in SimulateEvents but double check with halflength
	if(fabs(z) > halfLength) return;

	//all momenta and energy are in GeV
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
	//if(rp.Particle.tProd()/(_sol*1e3)*1e9 > 1e2) cout << "LARGE TIME - CalcTrajectory - original time " << rp.Position.T()*1e9 << " " << Position.T()*1e9 << " " << rp.Particle.tProd()/(_sol*1e3)*1e9 << endl;
	//cout << "original position x: " << rp.Position.x() << " y: " << rp.Position.y() << " z: " << rp.Position.z() << endl;
	//cout << "original momentum px: " << rp.Momentum.px() << " py: " << rp.Momentum.py() << " pz: " << rp.Momentum.pz() << " eta: " << rp.Momentum.eta() << " phi: " << rp.Momentum.phi() << " pt: " << pt << endl; 
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
		//    initial transverse momentum direction phi_0 = -atan(p_{X0} / p_{Y0})                omega = q * _b / (gammam); //in [89875518/s] - should be [rad/s]?

      		//    relativistic gamma: gamma = E / mc^2; gammam = gamma * m
      		//    gyration frequency omega = q * Bz / (gammam)
      		//    helix radius r = p_{T0} / (omega * gammam)

		gammam = e * 1e9 / (_sol * _sol); //gammam in [eV/c^2], c needs to be in [m/s]
		omega = q * _b / (gammam); //in [89875518/s] - should be [rad/s]?
		r = pt * 1e9 / (q * _b * _sol); //in [m]
//cout << "radius of curvature " << r << " energy " << rp.Momentum.E() << endl;
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
		//cout << "r " << r << " phi0 " << phi0 << " og phi " << rp.Momentum.phi() << " phid " << phid << " omega " << omega << " td " << td << " atan(x_c, y_c) " << atan2(x_c,y_c) << " id " << rp.Particle.id() << " pt " << rp.Particle.pT() << endl;//" " << atan2(y_c, x_c) << endl;
		//reset momentum
		rp.Momentum.SetCoordinates(pt, rp.Momentum.Eta(), phid, e);
		//rp.Momentum.SetCoordinates(px, py, pz, e);
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
	//if(rp.Particle.tProd()/(_sol*1e3)*1e9 > 1e2) cout << " energy " << e << " pt " << rp.Momentum.pt() << " new time " << rp.Position.T()*1e9 << " charge " << q << endl;
	//cout << "new position x: " << rp.Position.x() << " y: " << rp.Position.y() << " z: " << rp.Position.z() << " eta: " << rp.Position.eta() << " phi: " << rp.Position.phi() << endl;
	//cout << "new momentum px: " << rp.Momentum.px() << " py: " << rp.Momentum.py() << " pz: " << rp.Momentum.pz() << " eta: " << rp.Momentum.eta() << " phi: " << rp.Momentum.phi() << endl; 

	
}


//this also fills the "track" information which is just gen because there is no tracker in this sim :)
//this is used to get the momentum of the resulting groups of rhs s.t. quantities like masses can be calculated
void BasicDetectorSim::FillTracks(RecoParticle& rp){
	//particle has already been propagated to detector face
	_trackpx.push_back(rp.Momentum.px());
	_trackpy.push_back(rp.Momentum.py());
	_trackpz.push_back(rp.Momentum.pz());

	_tracketa.push_back(rp.Momentum.eta());
	//put phi on [0, 2pi] domain
	if(rp.Momentum.phi() < 0)
		_trackphi.push_back(rp.Momentum.phi()+2*acos(-1));
	else if(rp.Momentum.phi() >= 2*acos(-1))
		_trackphi.push_back(rp.Momentum.phi()-2*acos(-1));
	else
		_trackphi.push_back(rp.Momentum.phi());
}



//fill ecal cells - create showers
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
	double e_check = 0;
	//add energy to right ieta, iphi cell
	//set as max energy
//	if(_cal[ieta][iphi].at(0) < e){	
		//_cal[ieta][iphi].SetValue(e+_cal[ieta][iphi].at(0), 0);	
		//add time to right ieta, iphi cell	
		//_cal[ieta][iphi].SetValue(t, 1);	
		//add number of emissions to right ieta, iphi cell	
		//_cal[ieta][iphi].SetValue(1+_cal[ieta][iphi].at(2),2);
//	}
	double teta, tphi;
	//check transformation
	_get_etaphi(ieta, iphi, teta, tphi);
	//cout << "FillCal check transformation - eta " << eta << " ieta " << ieta << " phi " << phi << " iphi " << iphi << " teta " << teta << " tphi " << tphi << " energy " << _cal[ieta][iphi].at(0) << " emissions " << _cal[ieta][iphi].at(2) << endl;
	
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
			e_check += e_cell;
			//cout << "filling cell ieta " << iieta << " iphi " << iiphi << endl;
			//cout << "eta " << eta << " phi " << phi << " iieta " << iieta << " iiphi " << iiphi << " e_cell " << e_cell << " e " << _cal[iieta][iiphi].at(0) << endl;		

		}
	}
	//cout << "original energy " << e << " showered energy " << e_check << " ratio " << e_check/e << " _ncell " << _ncell << endl;
}


//add up emissions in rec hits
//some may be overlapping so need to fill call
//with all particles first
void BasicDetectorSim::MakeRecHits(){
	double e, t, x, y, z, theta;
	double e_cell, t_cell;
	double e_sig, t_sig;
	double eta, phi;
	int nrhs = 0;
	double etot = 0;
	double etot_og = 0;
	//declare fjinput/output containers
	vector<fastjet::PseudoJet> fjinputs;
	//do energy first to get transfer factor
	for(int i = 0; i < _netacal; i++){
		for(int j = 0; j < _nphical; j++){
			e = _cal[i][j].at(0);
			if(e == 0.) continue;
			//get cell bounds
			_get_etaphi(i, j, eta, phi);
			
			//time in cell ieta, jphi is energy weighted average
			//do amplitude dependent energy smearing
			//that follows equation 1.2 in CMS TDR
			//(sig/E)^2 = (s/sqrt(E))^2 + (n/E)^2 + c^2
			e_sig = _s*_s/e + _n*_n/(e*e) + _c*_c;
			e_sig = sqrt(e_sig);
			//e_sig is in units %, need to convert to units energy
			e_sig = e_sig/100. * e;
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

			//cout << "# emissions " << _cal[i][j].at(2) << " original e " << e << " smeared e " << e_cell << " e sig " << e_sig << endl;

			t = _cal[i][j].at(1)/((double)_cal[i][j].at(2));
			
			//do amplitude dependent time smearing
			t_sig = _calTresCte*_calTresCte + _calTresNoise*_calTresNoise/(e_cell*e_cell) + (_calTresStoch*_calTresStoch)/e_cell;
			t_sig = sqrt(t_sig/2.); //divide by 2 bc params were taken from measurements of 2 rechits
			//smear time in cell
			//t can be negative (early times)
			//update range to be centered on t, up to 5 sigma (calTres)
			_rs.SetRange(t - 5*t_sig, t + 5*t_sig);
			t_cell = t;//_rs.SampleGaussian(t, t_sig, 1).at(0);
			//if(e_cell > 1) cout << "t " << t*1e9 << " e " << e_cell << " tsig " << t_sig*1e9 << " t_cell " << t_cell*1e9 << endl;
			//cout << "filling cell ieta " << i << " iphi " << j << " og e " << e << " ecell " << e_cell << " esig " << e_sig << " e_sig % " << e_sig/e << endl;	
			//reset e and t for cal cells
			etot += e_cell;
			etot_og += e;
			_cal[i][j].SetValue(e_cell, 0);
			_cal[i][j].SetValue(t_cell, 1);		
			
			//get x, y, z based on cell eta phi
			x = _rmax*cos(phi);
			y = _rmax*sin(phi);
			theta = 2*atan2(1,exp(eta));
			z = _rmax/tan(theta);
			t = _cal[i][j].at(1); 
			
			JetPoint jp(x*1e2, y*1e2, z*1e2, t*1e9);
			jp.SetEnergy(_cal[i][j].at(0));
			Jet jet(jp, _PV);
			_cal_rhs.push_back(jet);
			//add particle to fastjet
			//running fastjet on reco cells
			//TURN ON HERE TO RUN FASTJET ON RECHITS
			fastjet::PseudoJet input(jet.px(), jet.py(), jet.pz(), jet.e());
//cout << "rh eta " << eta << " rh jet eta " << jet.eta() << " input eta " << input.eta() << endl;
			input.set_user_index(i*1000 + j); //i = idx / 1000, j = idx % 1000
			nrhs++;			
			fjinputs.push_back(input); 
			
			//save rec hits to tree
			_rhE.push_back(_cal[i][j].at(0));			
			_rhx.push_back(x*1e2); //in cm
			_rhy.push_back(y*1e2); //in cm
			_rhz.push_back(z*1e2); //in cm
			_rht.push_back(t*1e9); //in ns
			_rheta.push_back(eta);
			_rhphi.push_back(phi);
			_rhids.push_back(i*1000 + j);
			_nRhs++;

		}
	}
	

	//cout << "reco event total energy " << etot << endl;
	//cluster reco particles with fastjet	
	//run fastjet
	_recocsAK4 = fastjet::ClusterSequence(fjinputs, _jetdef_AK4);
	//get jets - min 5 pt
	vector<fastjet::PseudoJet> fjoutputs_AK4 = _recocsAK4.inclusive_jets(5.);
	for(int j = 0; j < fjoutputs_AK4.size(); j++) _recoAK4jets.push_back(fjoutputs_AK4[j]);
	//sort jets by pt
	_recoAK4jets = sorted_by_pt(_recoAK4jets);
	
	_recocsAK8 = fastjet::ClusterSequence(fjinputs, _jetdef_AK8);
	//get jets - min 5 pt
	vector<fastjet::PseudoJet> fjoutputs_AK8 = _recocsAK8.inclusive_jets(5.);
	for(int j = 0; j < fjoutputs_AK8.size(); j++) _recoAK8jets.push_back(fjoutputs_AK8[j]);
	//sort jets by pt
	_recoAK8jets = sorted_by_pt(_recoAK8jets);
	
	_recocsAK15 = fastjet::ClusterSequence(fjinputs, _jetdef_AK15);
	//get jets - min 5 pt
	vector<fastjet::PseudoJet> fjoutputs_AK15 = _recocsAK15.inclusive_jets(5.);
	for(int j = 0; j < fjoutputs_AK15.size(); j++) _recoAK15jets.push_back(fjoutputs_AK15[j]);
	//sort jets by pt
	_recoAK15jets = sorted_by_pt(_recoAK15jets);
	//cout << fjoutputs.size() << " " << " reco jets from " << fjinputs.size() << " inputs " << endl;
	//for(auto j : _jetsReco){
	//	cout << "reco jet e " << j.E() << " pt " << j.pt() << " eta " << j.eta() << " phi " << j.phi_std() << " m " << j.m() << endl;
	//}
	//cout << "gen jets" << endl;
	//for(auto j : _jets){
	//	cout << "gen jet e " << j.E() << " pt " << j.pt() << " eta " << j.eta() << " phi " << j.phi_std() << " m " << j.m() << endl;
	//}

	fjinputs.clear();
	fjinputs.resize(0);
	fjoutputs_AK4.clear();
	fjoutputs_AK4.resize(0);
	fjoutputs_AK8.clear();
	fjoutputs_AK8.resize(0);
	fjoutputs_AK15.clear();
	fjoutputs_AK15.resize(0);

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
		//loop through corresponding nxn grid
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
				//cout << "particle #" << p << " filling cell ieta " << iieta << " iphi " << iiphi << " cell e " << _cal[iieta][iiphi].at(0) << endl;
				if(_cal[iieta][iiphi].at(0) < _ethresh) continue; //_cal[iieta][iiphi].SetValue(0., 0.);
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
				reco_t += _cal[iieta][iiphi].at(1);
				reco_nrh++;
		
				//add emission to particle
				//doesn't include tof correction to time
				//save time in ns, space in cm
				JetPoint jp(x*1e2, y*1e2, z*1e2, t*1e9);
				jp.SetEnergy(_cal[iieta][iiphi].at(0));
				jp.SetEta(ceta);
				jp.SetPhi(cphi);
				//rhid = iieta0iiphi
				jp.SetRecHitId(int(iieta*1e4 + iiphi));
				_recops[p].AddEmission(jp);
				//Jet jet(jp,_PV);
				//_cal_rhs.push_back(jet);
				//double theta = atan2( sqrt(jp.x()*jp.x() + jp.y()*jp.y()), jp.z() );
				//cout << "iieta " << iieta << " iiphi " << iiphi << " e " << _cal[iieta][iiphi].at(0) << " jet e " << jet.e() << " jet pt " << jet.pt() << " jet px " << jet.px() << " py " << jet.py() << " pz " << jet.pz() << " e*sin(theta) " << jp.e()*sin(theta) << " e/cosh(eta) " << jp.e()/cosh(jp.eta()) << endl;	
			
	
				//simulate spikes
				//emissions separate from showers
				//don't want to conflate spike energy/times with shower energy/times
				if(_spikes){
					//roll dice to see if spike occurs
					_rs.SetRange(0.,1.);
					r = _rs.SampleFlat();
					if(r > _spikeprob) continue;				
			
					//if yes, roll dice for energy
					//assume spikes have a characteristic/on average energy of ~80 GeV
					//gain switch is not calibrated for energies above 120 GeV
					//so the time reco gets weird (in CMS)
					//0 is 4sigma away so captures most of the distribution
					_rs.SetRange(0.,120.);
					reco_e = _rs.SampleGaussian(80., 20., 1).at(0);
					//roll dice for time
					//assume spikes have an early time between 5 to 15 ns early
					//peaked at -10 ns
					//5 sigma spread, make sure spike is OUT of time
					//t can be negative (early times)
					_rs.SetRange(-35*1e-9, 0.);
					reco_t = _rs.SampleGaussian(-10.*1e-9,5.*1e-9,1).at(0); 
					//smear eta, phi based on cell dimensions
					//this is to avoid 2D overlap with true rec hits
					_rs.SetRange(ceta - _deta/2., ceta + _deta/2.);
					ceta = _rs.SampleFlat();
					_rs.SetRange(cphi - _dphi/2., cphi + _dphi/2.);
					cphi = _rs.SampleFlat();
					//that are measured in cell center
					//get x, y, z based on cell eta phi
					x = _rmax*cos(cphi);
					y = _rmax*sin(cphi);
					theta = 2*tan(exp(-ceta));
					z = sqrt(x*x + y*y)/tan(theta);

					//doesn't include tof correction to time
					//save time in ns, space in cm
					JetPoint jp(x*1e2, y*1e2, z*1e2, reco_t*1e9);
					jp.SetEnergy(reco_e);
					jp.SetEta(eta);
					jp.SetPhi(phi);
					jp.SetRecHitId(int(iieta*1e4 + iiphi));
					Jet j(jp,_PV);
					_cal_rhs.push_back(j);
				
					//save spikes to tree
					_spikeE.push_back(reco_e);			
					_nSpikes++;

				}
		
	

			}
		}
		//reset e and t for reco particle
		_recops[p].Momentum.SetE(reco_e);
		//set pt wrt to new energy and m = 0
		//p = cosh(eta) = E for m = 0
		_recops[p].Momentum.SetPt(reco_e / cosh(_recops[p].Momentum.eta()) );
		//double pt = reco_e / cosh(_recops[p].Momentum.eta());
		//_recops[p].Momentum.SetPx(pt*cos(_recops[p].Momentum.phi()));
		//_recops[p].Momentum.SetPy(pt*sin(_recops[p].Momentum.phi()));
		//_recops[p].Momentum.SetPz(pt*sinh(_recops[p].Momentum.eta()));
		//_recops[p].Position.SetCoordinates(_recops[p].Position.x()*1e2, _recops[p].Position.y()*1e2, _recops[p].Position.z()*1e2, reco_t/((double)reco_nrh)*1e9);
		//cout << " reco particle " << p << " eta " << _recops[p].Position.eta() << " phi " << _recops[p].Position.phi() << " gen eta " << _recops[p].Particle.eta() << " gen phi " << _recops[p].Particle.phi() << " reco pt " << _recops[p].Momentum.pt() << " gen pt " << _recops[p].Particle.pT() << " reco_e " << reco_e << " gen e " << _recops[p].Particle.e() << " ratio reco E / gen E " << _recops[p].Momentum.E()/_recops[p].Particle.e() << " mom eta " << _recops[p].Momentum.eta() << " pos eta " << _recops[p].Position.eta() << endl;
      		//RUN FASTJET ON RECO PARTICLES (NOT RECHITS)
		//fjinputs.push_back( fastjet::PseudoJet( _recops[p].Momentum.px(),
      		//  _recops[p].Momentum.py(), _recops[p].Momentum.pz(), _recops[p].Momentum.e() ) );
		_nRecoParticles++;

	}

}



void BasicDetectorSim::FillGenParticles(){
	_ngenparts = (int)_genparts.size();
	for(int g = 0; g < _genparts.size(); g++){
		_genparteta.push_back(_genparts[g].Position.eta());
		_genpartphi.push_back(_genparts[g].Position.phi());
		_genpartenergy.push_back(_genparts[g].Momentum.e());
		_genpartpt.push_back(_genparts[g].Momentum.pt());
		_genpartpz.push_back(_genparts[g].Momentum.pz());
		_genpartmass.push_back(_genparts[g].Momentum.mass());
	}
}

void BasicDetectorSim::FillGenJets(){
	//vector<fastjet::PseudoJet> consts;
	_ngenAK4jets = _genAK4jets.size();
	vector<fastjet::PseudoJet> consts;
	for(auto jet : _genAK4jets){
		_jgAK4eta.push_back(jet.eta());
		_jgAK4phi.push_back(jet.phi());
		_jgAK4energy.push_back(jet.e());
		_jgAK4pt.push_back(jet.pt());
		_jgAK4pz.push_back(jet.pz());
		_jgAK4mass.push_back(jet.m());
		consts = jet.constituents();
		_jgAK4nparts.push_back((int)consts.size());
		
		_jgAK4partIdxs.push_back({});
		for(auto c : consts)
			_jgAK4partIdxs[_jgAK4partIdxs.size()-1].push_back(c.user_index());
	}

	_ngenAK8jets = _genAK8jets.size();
	consts.clear();
	for(auto jet : _genAK8jets){
		_jgAK8eta.push_back(jet.eta());
		_jgAK8phi.push_back(jet.phi());
		_jgAK8energy.push_back(jet.e());
		_jgAK8pt.push_back(jet.pt());
		_jgAK8pz.push_back(jet.pz());
		_jgAK8mass.push_back(jet.m());
		consts = jet.constituents();
		_jgAK8nparts.push_back((int)consts.size());
		
		_jgAK8partIdxs.push_back({});
		for(auto c : consts)
			_jgAK8partIdxs[_jgAK8partIdxs.size()-1].push_back(c.user_index());
	}
	
	_ngenAK15jets = _genAK15jets.size();
	consts.clear();
	for(auto jet : _genAK15jets){
		_jgAK15eta.push_back(jet.eta());
		_jgAK15phi.push_back(jet.phi());
		_jgAK15energy.push_back(jet.e());
		_jgAK15pt.push_back(jet.pt());
		_jgAK15pz.push_back(jet.pz());
		_jgAK15mass.push_back(jet.m());
		consts = jet.constituents();
		_jgAK15nparts.push_back((int)consts.size());
		
		_jgAK15partIdxs.push_back({});
		for(auto c : consts)
			_jgAK15partIdxs[_jgAK15partIdxs.size()-1].push_back(c.user_index());
	}

}

void BasicDetectorSim::FillRecoJets(){
	_nrecoAK4jets = _recoAK4jets.size();
	vector<fastjet::PseudoJet> consts;
	for(auto jet : _recoAK4jets){
		consts = jet.constituents();
		_jAK4rhids.push_back({});
		for(auto c : consts){
			_jAK4rhids[_jAK4rhids.size()-1].push_back(c.user_index());
		}	
		_jAK4eta.push_back(jet.eta());
		_jAK4phi.push_back(jet.phi());
		_jAK4energy.push_back(jet.e());
		_jAK4pt.push_back(jet.pt());
		_jAK4mass.push_back(jet.m());
	}
	
	_nrecoAK8jets = _recoAK8jets.size();
	consts.clear();
	for(auto jet : _recoAK8jets){
		consts = jet.constituents();
		_jAK8rhids.push_back({});
		for(auto c : consts){
			_jAK8rhids[_jAK8rhids.size()-1].push_back(c.user_index());
		}	
		_jAK8eta.push_back(jet.eta());
		_jAK8phi.push_back(jet.phi());
		_jAK8energy.push_back(jet.e());
		_jAK8pt.push_back(jet.pt());
		_jAK8mass.push_back(jet.m());
	}
	
	_nrecoAK15jets = _recoAK15jets.size();
	consts.clear();
	for(auto jet : _recoAK15jets){
		consts = jet.constituents();
		_jAK15rhids.push_back({});
		for(auto c : consts){
			_jAK15rhids[_jAK15rhids.size()-1].push_back(c.user_index());
		}	
		_jAK15eta.push_back(jet.eta());
		_jAK15phi.push_back(jet.phi());
		_jAK15energy.push_back(jet.e());
		_jAK15pt.push_back(jet.pt());
		_jAK15mass.push_back(jet.m());
	}
}

//get "Jets" for clustering
void BasicDetectorSim::GetRecHits(vector<Jet>& rhs){
	rhs.clear();
	for(int j = 0; j < _cal_rhs.size(); j++){
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



//returns gen AK4 jets
void BasicDetectorSim::GetTrueJets(vector<Jet>& jets){
	jets.clear();
	for(int i = 0; i < _genAK4jets.size(); i++)
		jets.push_back(Jet( _genAK4jets[i].px(), _genAK4jets[i].py(), _genAK4jets[i].pz(), _genAK4jets[i].e()) );

}


//init tree
void BasicDetectorSim::InitTree(string fname){
	std::unique_ptr<TFile> f(new TFile(fname.c_str(), "RECREATE"));
	_file = std::move(f);	
	TDirectory* dir = _file->mkdir("tree");
	dir->cd();

	_tree = new TTree("llpgtree","llpgtree");
	_tree->Branch("event", &_evt)->SetTitle("Event");
	_tree->Branch("ECALRecHit_energy", &_rhE)->SetTitle("rec hit energy (GeV)");
	_tree->Branch("ECALRecHit_rhx", &_rhx)->SetTitle("rec hit x coord (cm)");
	_tree->Branch("ECALRecHit_rhy", &_rhy)->SetTitle("rec hit y coord (cm)");
	_tree->Branch("ECALRecHit_rhz", &_rhz)->SetTitle("rec hit z coord (cm)");
	_tree->Branch("ECALRecHit_time", &_rht)->SetTitle("rec hit time (ns)");
	_tree->Branch("ECALRecHit_eta", &_rheta)->SetTitle("rec hit eta");
	_tree->Branch("ECALRecHit_phi", &_rhphi)->SetTitle("rec hit phi");
	_tree->Branch("ECALRecHit_ID", &_rhids)->SetTitle("rec hit id");
	_tree->Branch("nRHs",&_nRhs)->SetTitle("Number of rec hits");
	
	_tree->Branch("PV_x",&_pvx)->SetTitle("x coordinate PV (cm)");
	_tree->Branch("PV_y",&_pvy)->SetTitle("y coordinate PV (cm)");
	_tree->Branch("PV_z",&_pvz)->SetTitle("z coordinate PV (cm)");
	_tree->Branch("PV_t",&_pvt)->SetTitle("t coordinate PV (ns)");

	_tree->Branch("ootPV_x",&_oot_pvx)->SetTitle("x coordinate of oot PV (cm)");
	_tree->Branch("ootPV_y",&_oot_pvy)->SetTitle("y coordinate of oot PV (cm)");
	_tree->Branch("ootPV_z",&_oot_pvz)->SetTitle("z coordinate of oot PV (cm)");
	_tree->Branch("ootPV_t",&_oot_pvt)->SetTitle("t coordinate of oot PV (ns)");

	_tree->Branch("ECALSpike_energy", &_spikeE)->SetTitle("spike energy (GeV)");
	_tree->Branch("nSpikes", &_nSpikes)->SetTitle("Number of spikes");
	_tree->Branch("nRecoParticles", &_nRecoParticles)->SetTitle("Number of reco particles");

	//gen AK4 jets - gen particles clustered with FJ AK4
	_tree->Branch("AK4Jet_genEta", &_jgAK4eta)->SetTitle("Gen jet eta - FastJet AK4");
	_tree->Branch("AK4Jet_genPhi", &_jgAK4phi)->SetTitle("Gen jet phi - FastJet AK4");
	_tree->Branch("AK4Jet_genEnergy",&_jgAK4energy)->SetTitle("Gen jet energy - FastJet AK4");
	_tree->Branch("AK4Jet_genPt",&_jgAK4pt)->SetTitle("Gen jet pt - FastJet AK4");
	_tree->Branch("AK4Jet_genPz",&_jgAK4pz)->SetTitle("Gen jet pz - FastJet AK4");
	_tree->Branch("AK4Jet_genMass",&_jgAK4mass)->SetTitle("Gen jet mass - FastJet AK4");
	_tree->Branch("AK4Jet_genNJet",&_ngenAK4jets)->SetTitle("Number of gen jets - FastJet AK4");
	_tree->Branch("AK4Jet_genConstituentIdxs", &_jgAK4partIdxs)->SetTitle("Gen jet constituent indices - FastJet AK4");
	_tree->Branch("AK4Jet_genNConstituents", &_jgAK4nparts)->SetTitle("Gen jet n constituents - FastJet AK4");
	
	//gen AK8 jets - gen particles clustered with FJ AK8
	_tree->Branch("AK8Jet_genEta", &_jgAK8eta)->SetTitle("Gen jet eta - FastJet AK8");
	_tree->Branch("AK8Jet_genPhi", &_jgAK8phi)->SetTitle("Gen jet phi - FastJet AK8");
	_tree->Branch("AK8Jet_genEnergy",&_jgAK8energy)->SetTitle("Gen jet energy - FastJet AK8");
	_tree->Branch("AK8Jet_genPt",&_jgAK8pt)->SetTitle("Gen jet pt - FastJet AK8");
	_tree->Branch("AK8Jet_genPz",&_jgAK8pz)->SetTitle("Gen jet pz - FastJet AK8");
	_tree->Branch("AK8Jet_genMass",&_jgAK8mass)->SetTitle("Gen jet mass - FastJet AK8");
	_tree->Branch("AK8Jet_genNJet",&_ngenAK8jets)->SetTitle("Number of gen jets - FastJet AK8");
	_tree->Branch("AK8Jet_genConstituentIdxs", &_jgAK8partIdxs)->SetTitle("Gen jet constituent indices - FastJet AK8");
	_tree->Branch("AK8Jet_genNConstituents", &_jgAK8nparts)->SetTitle("Gen jet n constituents - FastJet AK8");
	
	//gen AK15 jets - gen particles clustered with FJ AK15
	_tree->Branch("AK15Jet_genEta", &_jgAK15eta)->SetTitle("Gen jet eta - FastJet AK15");
	_tree->Branch("AK15Jet_genPhi", &_jgAK15phi)->SetTitle("Gen jet phi - FastJet AK15");
	_tree->Branch("AK15Jet_genEnergy",&_jgAK15energy)->SetTitle("Gen jet energy - FastJet AK15");
	_tree->Branch("AK15Jet_genPt",&_jgAK15pt)->SetTitle("Gen jet pt - FastJet AK15");
	_tree->Branch("AK15Jet_genPz",&_jgAK15pz)->SetTitle("Gen jet pz - FastJet AK15");
	_tree->Branch("AK15Jet_genMass",&_jgAK15mass)->SetTitle("Gen jet mass - FastJet AK15");
	_tree->Branch("AK15Jet_genNJet",&_ngenAK15jets)->SetTitle("Number of gen jets - FastJet AK15");
	_tree->Branch("AK15Jet_genConstituentIdxs", &_jgAK15partIdxs)->SetTitle("Gen jet constituent indices - FastJet AK15");
	_tree->Branch("AK15Jet_genNConstituents", &_jgAK15nparts)->SetTitle("Gen jet n constituents - FastJet AK15");
	
	
	//gen decay info
	_tree->Branch("Top_decayId",&_topDecayId)->SetTitle("gen top decay type (0 = had, 1 = lep)");
	_tree->Branch("W_decayId",&_wDecayId)->SetTitle("gen W decay type (0 = had, 1 = lep)");


	//reco jets - cells clustered with FJ AK4
	_tree->Branch("AK4Jet_eta",&_jAK4eta)->SetTitle("Jet eta - FastJet AK4, reco");
	_tree->Branch("AK4Jet_phi",&_jAK4phi)->SetTitle("Jet phi - FastJet AK4, reco");
	_tree->Branch("AK4Jet_energy",&_jAK4energy)->SetTitle("Jet energy - FastJet AK4, reco");
	_tree->Branch("AK4Jet_pt",&_jAK4pt)->SetTitle("Jet pt - FastJet AK4, reco");
	_tree->Branch("AK4Jet_mass",&_jAK4mass)->SetTitle("Jet mass - FastJet AK4, reco");
	_tree->Branch("AK4Jet_RhIDs",&_jAK4rhids)->SetTitle("Jet rh ids - FastJet AK4");
	_tree->Branch("AK4Jet_NJet",&_nrecoAK4jets)->SetTitle("Number of jets - FastJet AK4");
	
	//reco jets - cells clustered with FJ AK8
	_tree->Branch("AK8Jet_eta",&_jAK8eta)->SetTitle("Jet eta - FastJet AK8, reco");
	_tree->Branch("AK8Jet_phi",&_jAK8phi)->SetTitle("Jet phi - FastJet AK8, reco");
	_tree->Branch("AK8Jet_energy",&_jAK8energy)->SetTitle("Jet energy - FastJet AK8, reco");
	_tree->Branch("AK8Jet_pt",&_jAK8pt)->SetTitle("Jet pt - FastJet AK8, reco");
	_tree->Branch("AK8Jet_mass",&_jAK8mass)->SetTitle("Jet mass - FastJet AK8, reco");
	_tree->Branch("AK8Jet_RhIDs",&_jAK8rhids)->SetTitle("Jet rh ids - FastJet AK8");
	_tree->Branch("AK8Jet_NJet",&_nrecoAK8jets)->SetTitle("Number of jets - FastJet AK8");
	
	//reco jets - cells clustered with FJ AK15
	_tree->Branch("AK15Jet_eta",&_jAK15eta)->SetTitle("Jet eta - FastJet AK15, reco");
	_tree->Branch("AK15Jet_phi",&_jAK15phi)->SetTitle("Jet phi - FastJet AK15, reco");
	_tree->Branch("AK15Jet_energy",&_jAK15energy)->SetTitle("Jet energy - FastJet AK15, reco");
	_tree->Branch("AK15Jet_pt",&_jAK15pt)->SetTitle("Jet pt - FastJet AK15, reco");
	_tree->Branch("AK15Jet_mass",&_jAK15mass)->SetTitle("Jet mass - FastJet AK15, reco");
	_tree->Branch("AK15Jet_RhIDs",&_jAK15rhids)->SetTitle("Jet rh ids - FastJet AK15");
	_tree->Branch("AK15Jet_NJet",&_nrecoAK15jets)->SetTitle("Number of jets - FastJet AK15");


	_tree->Branch("Track_px", &_trackpx)->SetTitle("Track px");
	_tree->Branch("Track_py", &_trackpy)->SetTitle("Track py");
	_tree->Branch("Track_pz", &_trackpz)->SetTitle("Track pz");
	_tree->Branch("Track_eta", &_tracketa)->SetTitle("Track eta");
	_tree->Branch("Track_phi", &_trackphi)->SetTitle("Track phi");
	
	_tree->Branch("genpart_eta", &_genparteta)->SetTitle("genpart eta");
	_tree->Branch("genpart_phi", &_genpartphi)->SetTitle("genpart phi");
	_tree->Branch("genpart_energy",&_genpartenergy)->SetTitle("genpart energy");
	_tree->Branch("genpart_pt",&_genpartpt)->SetTitle("genpart pt");
	_tree->Branch("genpart_pz",&_genpartpz)->SetTitle("genpart pz");
	_tree->Branch("genpart_mass",&_genpartmass)->SetTitle("genpart mass");
	_tree->Branch("genpart_id",&_genpartids)->SetTitle("genpart pdg id");
	_tree->Branch("genpart_momIdx",&_genpartMomIdx)->SetTitle("genpart momIdx");
	_tree->Branch("genpart_idx",&_genpartIdx)->SetTitle("indices of genparts");
	_tree->Branch("genpart_ngenpart",&_ngenparts)->SetTitle("number of genparts");
}


void BasicDetectorSim::_reset(){
	_rhE.clear();
	_rhx.clear();
	_rhy.clear();
	_rhz.clear();
	_rht.clear();
	_rheta.clear();
	_rhphi.clear();
	_rhids.clear();
	_spikeE.clear();
	_evt = 0;
	_nRhs = 0;
	_nRecoParticles = 0;
	_nSpikes = 0;
	for(int i = 0; i < _netacal; i++)
		for(int j = 0; j < _nphical; j++)
			_cal[i][j] = BayesPoint({0., 0., 0.});
	_cal_rhs.clear();
	_recops.clear();
	_recoAK4jets.clear();
	_recoAK8jets.clear();
	_recoAK15jets.clear();
	
	_genAK4jets.clear();
	_genAK8jets.clear();
	_genAK15jets.clear();

	_jgAK4eta.clear();
	_jgAK4phi.clear();
	_jgAK4energy.clear();
	_jgAK4pt.clear();
	_jgAK4pz.clear();
	_jgAK4mass.clear();
	_jgAK4nparts.clear();
	for(auto j : _jgAK4partIdxs)
		j.clear();
	_jgAK4partIdxs.clear();
	
	_jgAK8eta.clear();
	_jgAK8phi.clear();
	_jgAK8energy.clear();
	_jgAK8pt.clear();
	_jgAK8pz.clear();
	_jgAK8mass.clear();
	_jgAK8nparts.clear();
	for(auto j : _jgAK8partIdxs)
		j.clear();
	_jgAK8partIdxs.clear();
	
	_jgAK15eta.clear();
	_jgAK15phi.clear();
	_jgAK15energy.clear();
	_jgAK15pt.clear();
	_jgAK15pz.clear();
	_jgAK15mass.clear();
	_jgAK15nparts.clear();
	for(auto j : _jgAK15partIdxs)
		j.clear();
	_jgAK15partIdxs.clear();

	_topDecayId.clear();
	_wDecayId.clear();

	_jAK4eta.clear();
	_jAK4phi.clear();
	_jAK4energy.clear();
	_jAK4pt.clear();
	_jAK4mass.clear();
	
	_jAK8eta.clear();
	_jAK8phi.clear();
	_jAK8energy.clear();
	_jAK8pt.clear();
	_jAK8mass.clear();

	_jAK15eta.clear();
	_jAK15phi.clear();
	_jAK15energy.clear();
	_jAK15pt.clear();
	_jAK15mass.clear();

	_genparts.clear();
	_genpartids.clear();
	_genpartMomIdx.clear();
	_genpartEvtIdx.clear();
	_genparteta.clear();
	_genpartphi.clear();
	_genpartenergy.clear();
	_genpartpt.clear();
	_genpartpz.clear();
	_genpartmass.clear();
	_genpartIdx.clear();
	

	for(auto j : _jAK4rhids)
		j.clear();
	_jAK4rhids.clear();
	
	for(auto j : _jAK8rhids)
		j.clear();
	_jAK8rhids.clear();
	
	for(auto j : _jAK15rhids)
		j.clear();
	_jAK15rhids.clear();
	
	_trackpx.clear();
	_trackpy.clear();
	_trackpz.clear();
	_tracketa.clear();
	_trackphi.clear();

	_pvx = 0;
	_pvy = 0;
	_pvz = 0;
	_pvt = 0;
	_oot_pvx = -999;
	_oot_pvy = -999;
	_oot_pvz = -999;
	_oot_pvt = -999;
}

void BasicDetectorSim::WriteTree(){
	cout << "Writing to " << _file->GetName() << endl;
	_file->cd();	
	_file->Write();
	_file->Close();
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
//cout << " eta " << eta << " etamax " << _etamax << " etamin " << _etamin << " ieta " << ieta << endl;
}

//get eta-phi center given ieta, iphi cell indices
void BasicDetectorSim::_get_etaphi(int ieta, int iphi, double& eta, double& phi){
 	double pi = acos(-1);
	//eta      
	eta = _etamin + _deta/2. + double(ieta)*_deta;      
	//phi 
	phi = _phimin + _dphi/2. + double(iphi)*_dphi;
}



