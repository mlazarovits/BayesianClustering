#ifndef BHCJETSKIMMER_HH
#define BHCJETSKIMMER_HH
#include "JetSimProducer.hh"
#include "BaseSkimmer.hh"
#include "BaseTree.hh"
#include "TGraph.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TExec.h"
#include "TLine.h"
#include <set>
#include <Math/Vector4D.h>

using node = BaseTree::node;
using procCat = BaseSkimmer::procCat;

class BHCJetSkimmer{
	public:
		BHCJetSkimmer(){
			_evti = 0;
			_evtj = 0;
			_gev = 1./10.;
			_oname = "";
			_radius = 0;
			_smear = false;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_thresh = 1.;
			_evt2disp = 0;
			_evt2disp_z = 0;

			_tree = new TTree("tree","tree");

			_zoom_window = false;
			
			_minTopPt = 0;
			_minTopE = 0;
			_minWPt = 0;		

			_nGhosts = 0;
			_check_merges = false;
	
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);
						
			_cell = 0;
			_tresCte = 0;
			_tresNoise = 0;
			_tresStoch = 0;

			_strategy = NlnN;
			_sel = def;

			_infile = nullptr;
		}

		virtual ~BHCJetSkimmer(){ }

		//BHCJetSkimmer(TFile* file){
		BHCJetSkimmer(string file){
			if(gSystem->AccessPathName(file.c_str())){ cout << "Error: file " << file << " doesn't exist." << endl; return; }
			InitHists();
			_infile = TFile::Open(file.c_str());
			_prod = std::make_unique<JetSimProducer>(_infile);
			_strategy = NlnN;
			_sel = def;
			_tree = new TTree("tree","tree");

			_oname = "";
			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_gev = 1./10.;
			_radius = 0;
			_smear = false;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_thresh = 1.;
			_minTopPt = 0;
			_minWPt = 0;		
			_minTopE = 0;
			_evt2disp = 0;
			_evt2disp_z = 0;
			
			_zoom_window = false;
			
			_nGhosts = 0;
			_check_merges = false;
			//beta
			_prior_params["scale"] = Matrix(1e-3);
			//nu
			_prior_params["dof"] = Matrix(3);
			//W
			Matrix W(3,3);
			W.InitIdentity();
			W.mult(W,1./3);
			_prior_params["scalemat"] = W;
			//m
			_prior_params["mean"] = Matrix(3,1);
			
			_cell = acos(-1)/180;
			_tresCte = 0.1727;//times given in ns//0.133913 * 1e-9;
			_tresStoch = 0.5109;//1.60666 * 1e-9; 
			_tresNoise = 2.106;//0.00691415 * 1e-9;
			
		}

		//instructions for adding a branch:
		//branch must be declared in InitMapTree()
		//by initializing the appropriate map (ie _obs, _vobs, _vvobs)
		//and then setting the branch address to that branch name (type+"_"+obsname)
		//then in the appropriate FillBranches* function (or on its own)
		//call appropriate *FillBranch(val, type, obsname) (ie vFillBranch(val, type, obsname) for val = jet.eta()) to fill that branch

		void InitMapTree(){
			vector<string> obs = {"Energy","Pt","Mass","EtaCenter","PhiCenter","TimeCenter", "EtaVar", "PhiVar", "TimeVar", "EtaPhiCov", "EtaTimeCov", "PhiTimeCov"};
			_obs["evt"] = -1;
			_tree->Branch("evt",&_obs.at("evt"))->SetTitle("event");
			
			//gen particles
			_vobs["genpart_EtaCenter"] = {};
			_tree->Branch("genpart_EtaCenter",&_vobs.at("genpart_EtaCenter"))->SetTitle("gen particle eta");
			_vobs["genpart_PhiCenter"] = {};
			_tree->Branch("genpart_PhiCenter",&_vobs.at("genpart_PhiCenter"))->SetTitle("gen particle phi");
			_vobs["genpart_Energy"] = {};
			_tree->Branch("genpart_Energy",&_vobs.at("genpart_Energy"))->SetTitle("gen particle energy");
			_vobs["genpart_Pt"] = {};
			_tree->Branch("genpart_Pt",&_vobs.at("genpart_Pt"))->SetTitle("gen particle pt");
			_vobs["genpart_Mass"] = {};
			_tree->Branch("genpart_Mass",&_vobs.at("genpart_Mass"))->SetTitle("gen particle mass");
			_vobs["genpart_ID"] = {};
			_tree->Branch("genpart_ID",&_vobs.at("genpart_ID"))->SetTitle("gen particle id");


			string obsname;
			for(int t = 0; t < _types.size(); t++){
				//avoid doing nJets, EtaCenter, etc plots for all gentypes, just use indexing to match gen-matched jets to their gen particles
				for(int g = 0; g < _genmatches.size(); g++){
					//jet
					string key = _types[t]+"Jet_"+_genmatches[g];
					obsname = "MatchedIdx";
					_vobs[key+obsname] = {};
					_tree->Branch((key+obsname).c_str(),&_vobs.at(key+obsname))->SetTitle(("index of gen "+_genmatches[g]+" this "+_types[t]+" jet is matched to").c_str());
				
					if(_types[t].find("gen") == string::npos){
						//index of parton subcluster is matched to
						obsname = "subclusterPartonMatchedIdx";
						key = _types[t]+"Jet"+_genmatches[g]+"_";
						_vvobs[key+obsname] = {};
						_tree->Branch((key+obsname).c_str(),&_vvobs.at(key+obsname))->SetTitle("index of daughter parton matched to subcluster given that its jet is matched to the parent (or same gen particle, ie jet--q && subcl--q)");
					}
				}
				string key = _types[t]+"Jet_";
				string extraTitle = "";
					
				obsname = "nJets";
				_obs[key+obsname] = -1;
				_tree->Branch((key+obsname).c_str(),&_obs.at(key+obsname))->SetTitle(("number of "+_types[t]+" jets"+extraTitle).c_str());
				
				if(_types[t].find("gen") == string::npos){
					obsname = "nSubclustersJet";
					_vobs[key+obsname] = {};
					_tree->Branch((key+obsname).c_str(),&_vobs.at(key+obsname))->SetTitle(("number of subclusters per"+_types[t]+" jet"+extraTitle).c_str());
					obsname = "JetSize";
					_vobs[key+obsname] = {};
					_tree->Branch((key+obsname).c_str(),&_vobs.at(key+obsname))->SetTitle((_types[t]+" jet size"+extraTitle).c_str());
				}

				
				if(_types[t] == "BHC"){
					obsname = "nGhosts";
					_vobs[key+obsname] = {};
					_tree->Branch((key+obsname).c_str(),&_vobs.at(key+obsname))->SetTitle(("number of ghost subclusters per"+_types[t]+" jet"+extraTitle).c_str());
				}
		
				//	
				for(int o = 0; o < obs.size(); o++){
					//skip shape variables for gen jets
					if(_types[t].find("gen") != string::npos && (obs[o].find("Var") != string::npos || obs[o].find("Cov") != string::npos)) continue;

					string title = "jet "+obs[o]+extraTitle;
					string units = "";
					if(obs[o] == "Energy" || obs[o] == "Mass" || obs[o] == "Pt")
						units = " [GeV]";
					if(obs[o] == "TimeCenter")
						units = " [ns]";
				
					//jet	
					_vobs[key+obs[o]] = {};
					title = "jet "+obs[o]+extraTitle+units;
					_tree->Branch((key+obs[o]).c_str(),&_vobs.at(key+obs[o]))->SetTitle(title.c_str());
				
					if(_types[t].find("gen") == string::npos){
						//subcluster	
						_vvobs[key+"subcluster"+obs[o]] = {};
						title = "subcluster "+obs[o]+" per jet"+extraTitle+units;
						_tree->Branch((key+"subcluster"+obs[o]).c_str(),&_vvobs.at(key+"subcluster"+obs[o]))->SetTitle(title.c_str());
					}
				}

			}
		}


		void InitHists(){
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_W);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_W2);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_gluon);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_q1);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_q2);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_b1);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_b2);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_top1);
			_evtdisps_obj.push_back(EvtDisplay_etaCell_phiCell_top2);
		}
		void SetMinRhE(double r){ _prod->SetMinRhE(r); }
		void SetRecoMinPt(double r){ _prod->SetRecoMinPt(r); }
		void SetRecoMinE(double r){ _prod->SetRecoMinE(r); }
		void SetGenMinPt(double r){ _prod->SetGenMinPt(r); }
		void SetGenMinE(double r){ _prod->SetGenMinE(r); }
		void SetMinNrhs(int r){ _prod->SetMinNrhs(r); }
		void SetMinNGenConsts(int r){ _prod->SetMinNGenConsts(r); }
		void SetGenTopMinPt(double r){ _minTopPt = r; 
			cout << "Minimum gen top pt: " << _minTopPt << endl;}
		void SetGenWMinPt(double r){ _minWPt = r; 
			cout << "Minimum gen W pt: " << _minWPt << endl;}
		bool _check_merges;
		void CheckMerges(bool t){ _check_merges = t; if(_check_merges) cout << "Checking merges" << endl; }
		void Skim();
		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else if(i == 2) _strategy = gmmOnly;
			else if(i == 3) _strategy = NlnNonAK4;
			else return; 
		}
		void SetEventSelection(int i){
			//default selection (generic hadronic)
			if(i == 0) _sel = def;
			else if(i == 1){
				_sel = singW;//single W
				_genmatches = {"W","q"};
			}
			else if(i == 2){
				_sel = boostTop;//ttbar
				_genmatches = {"W","q","b","Top"};
			}
			else if(i == 3){
				_sel = QCDdijets;//QCD
				_genmatches = {"q"};
			}
			else{
			} 
		}
		int _nGhosts;
		void SetNGhosts(int t){ _nGhosts = t; cout << "Adding " << _nGhosts << " ghosts to each BHC merging step." << endl;}
		//set evt display hist name according to event #
		//z is what to put on z axis (0 : energy, 1 : time)
		void SetEvent2Display(int e, int z = 0){
			_evt2disp = e;
			_evt2disp_z = z;	
		}		


		void TreesToJets(){
			_predJets.clear();
			vector<JetPoint> rhs;
			double x, y, z, eta, phi, t, theta, px, py, pz;
			BayesPoint vertex({_pvx, _pvy, _pvz});
			int njets_tot = 0;
			for(int i = 0; i < _trees.size(); i++){
				//get points from tree
				//at least 2 points (rhs)
				if(_trees[i]->points->GetNPoints() < 2) continue;
				rhs.clear();
				njets_tot++;
				//loop over points
				//cout << "TREE " << i << endl;
				//cout << "tree points " << endl; pc->Print();
				//create new Jet
				//Jet predJet(_trees[i]->model, BayesPoint({_pvx, _pvy, _pvz}), _gev, _radius);
				Jet predJet(_trees[i].get(), BayesPoint({_pvx, _pvy, _pvz}), _gev, _radius);
			//cout << "pre pt cut - pred jet px " << predJet.px() << " py " << predJet.py() << " pz " << predJet.pz() << " pt " << predJet.pt() << " E " << predJet.E() << " m2 " << predJet.m2() << " mass " << predJet.mass() << " eta " << predJet.eta() << " phi " << predJet.phi() << endl;
				//put pt cut in for predjets of 5 GeV to match reco AK4 definition
				if(predJet.pt() < 5) continue; 
				//add Jet to jets	
				_predJets.push_back(predJet);	
				vFillBranch(_trees[i]->model->GetNGhosts(),"BHCJet","nGhosts");
				//fill hists for mixture model
			}
			//sort predicted jets
			sort(_predJets.begin(), _predJets.end(), ptsort);
			cout << _predJets.size()  << " pred jets total" << endl;
			//cout << _predJets.size() << " pred jets pt > 20 GeV" << endl;
			for(auto j : _predJets) cout << "pred jet px " << j.px() << " py " << j.py() << " pz " << j.pz() << " E " << j.E() << " m2 " << j.m2() << " mass " << j.mass() << " eta " << j.eta() << " phi " << j.phi() << " pt " << j.pt() << endl;
		}



		void FillBranch(double obs, string type, string obsname){
			string key = type+"_"+obsname;
			if(_obs.find(key) == _obs.end()){
				cout << "vobs key " << key << " not found in map" << endl;
				return;
			}
			_obs.at(key) = obs;
		}

		void vFillBranch(double obs, string type, string obsname){
			string key;
			key = type+"_"+obsname;
			if(_vobs.find(key) == _vobs.end()){
				cout << "vobs key " << key << " not found in map" << endl;
				return;
			}
			_vobs.at(key).push_back(obs);
		}
		
		void vvFillBranch(double obs, string type, string obsname, int idx){
			string key = type+"_"+obsname;
			if(_vvobs.find(key) == _vvobs.end()){
				cout << "vvobs key " << key << " not found in map" << endl;
				return;
			}
			_vvobs.at(key)[idx].push_back(obs);
		}


		void addVector(string type, string obsname){
			string key = type+"_"+obsname;
			if(_vvobs.find(key) == _vvobs.end()){
				cout << "vvobs key " << key << " not found in map" << endl;
				return;
			}
			_vvobs[key].push_back({});

		}

		void FillBranchesGen(const Jet& jet){
			string type = "genpart";
			vFillBranch(jet.e(), type, "Energy");
			vFillBranch(jet.pt(), type, "Pt");
			vFillBranch(jet.m(), type, "Mass");
			vFillBranch(jet.eta(), type, "EtaCenter");
			vFillBranch(jet.phi(), type, "PhiCenter");
			int genidx = jet.GetUserIdx();
			int id = _base->genpart_id->at(genidx);
			vFillBranch(id, type, "ID");
		}


		//parent idx is the idx of the gen particle that the jet was matched to
		//gmit->second[jidx] = parentidx;
		//genWidx = parentidx
		//type = jettype+genparent = type+gmit->first 
		void SubclusterPartonGenMatching(string type, int parentidx, const Jet& jet, int jidx, bool matchDaughters = true){
			//if jet was matched to this type of gen particle, find subcluster matches to W partons
			vector<Jet> consts;
			jet.GetConstituents(consts);
			vector<int> subclGenMatchIdx(consts.size(), -1);	
			if(parentidx != -1){
				//subcluster-parton matching
				vector<Jet> partons;
				map<int, int> partonIdxToGlobalIdx;
				if(matchDaughters){
					//get gen partons from parent decay
					int ggenidx = _genparts[parentidx].GetUserIdx();
					for(int g = 0; g < _genparts.size(); g++){
						int genidx = _genparts[g].GetUserIdx();
						if(_base->genpart_momIdx->at(genidx) != ggenidx) continue;
						partons.push_back(_genparts[g]);
						partonIdxToGlobalIdx[partons.size()-1] = g;
					}
				}
				else{ //match to parent particle (ie q or b or gluon)
					partons.push_back(_genparts[parentidx]);
					partonIdxToGlobalIdx[partons.size()-1] = parentidx;
				}
				GenericMatchJet(consts,partons,subclGenMatchIdx); //match subclusters to W partons
				for(int k = 0; k < consts.size(); k++){
					//W partons need to map back to genparts idx
					int partonIdx = subclGenMatchIdx[k];
					if(partonIdx != -1)
						vvFillBranch(partonIdxToGlobalIdx[partonIdx], type, "subclusterPartonMatchedIdx", jidx);
					else
						vvFillBranch(partonIdx, type, "subclusterPartonMatchedIdx", jidx);
				}
			}
			//else jet was not matched to this type of gen particle, subclusters are not matched to partons at all
			else{
				for(int k = 0; k < consts.size(); k++){
					vvFillBranch(-1, type, "subclusterPartonMatchedIdx", jidx);
				}
			}

		}



		void FillBranchesJet(const Jet& jet, string type, int jidx, const map<string, vector<int>>& genmatches){
			if(find(_types.begin(), _types.end(), type) == _types.end()){
				cout << "Type " << type << " not in predetermined list. Not filling hists." << endl;
				return;
			}
			type += "Jet";
			//jet observables
			vFillBranch(jet.e(), type, "Energy");
			vFillBranch(jet.pt(), type, "Pt");
			vFillBranch(jet.m(), type, "Mass");
			vFillBranch(jet.eta(), type, "EtaCenter");
			vFillBranch(jet.phi(), type, "PhiCenter");
			vFillBranch(jet.t(), type, "TimeCenter");
			if(type.find("gen") == string::npos){
				Matrix cov = jet.GetCovariance();
				vFillBranch(cov.at(0,0), type, "EtaVar");
				vFillBranch(cov.at(1,1), type, "PhiVar");
				vFillBranch(cov.at(2,2), type, "TimeVar");
				vFillBranch(cov.at(0,1), type, "EtaPhiCov");
				vFillBranch(cov.at(0,2), type, "EtaTimeCov");
				vFillBranch(cov.at(2,1), type, "PhiTimeCov");
				double jetsize = CalcSize(cov);
				vFillBranch(jetsize, type, "JetSize");
			}
			//gen matching
			//in genmatches map: key is particle type, value is index matched to (in object's own vector)
			//need to save idx in overall genparticles list
			for(auto gmit = genmatches.begin(); gmit != genmatches.end(); gmit++){
				string obsname = gmit->first+"MatchedIdx";
				//jet gen matching
				vFillBranch(gmit->second[jidx], type, obsname);
				
				//do parton subcluster matching
				//match subclusters to DAUGHTERS of this parent particle the jet was matched to
				bool matchDaughters = true;
				if(gmit->first != "W" && gmit->first != "Top")
					matchDaughters = false; //match subclusters to parent particles (bc qs, bs, gluons do not have daughters that are saved in gen record)
				//don't do subcluster matching for gen jets (they don't have subclusters)
				if(type.find("gen") == string::npos){
					addVector(type+gmit->first, "subclusterPartonMatchedIdx");
					SubclusterPartonGenMatching(type+gmit->first, gmit->second[jidx], jet, jidx, matchDaughters);
				}
			}
			//skip filling subcluster branches for gen jets 
			if(type.find("gen") != string::npos) return;
			vFillBranch(jet.GetNConstituents(), type, "nSubclustersJet");
			//subcluster observables
			addVector(type, "subclusterEnergy");
			addVector(type, "subclusterPt");
			addVector(type, "subclusterMass");
			addVector(type, "subclusterEtaCenter");
			addVector(type, "subclusterPhiCenter");
			addVector(type, "subclusterTimeCenter");
			addVector(type, "subclusterEtaVar");
			addVector(type, "subclusterPhiVar");
			addVector(type, "subclusterTimeVar");
			addVector(type, "subclusterEtaPhiCov");
			addVector(type, "subclusterEtaTimeCov");
			addVector(type, "subclusterPhiTimeCov");
			Jet subcl;
			for(int k = 0; k < jet.GetNConstituents(); k++){
				jet.GetConstituent(k, subcl);
				Matrix subcl_cov = subcl.GetCovariance();

				vvFillBranch(subcl.e(), type, "subclusterEnergy",jidx);
				vvFillBranch(subcl.pt(), type,  "subclusterPt",jidx);
				vvFillBranch(subcl.m(), type,   "subclusterMass",jidx);
				vvFillBranch(subcl.eta(), type, "subclusterEtaCenter",jidx);
				vvFillBranch(subcl.phi(), type, "subclusterPhiCenter",jidx);
				vvFillBranch(subcl.t(), type,   "subclusterTimeCenter",jidx);
				vvFillBranch(subcl_cov.at(0,0), type, "subclusterEtaVar",jidx);
				vvFillBranch(subcl_cov.at(1,1), type, "subclusterPhiVar",jidx);
				vvFillBranch(subcl_cov.at(2,2), type, "subclusterTimeVar",jidx);
				vvFillBranch(subcl_cov.at(0,1), type, "subclusterEtaPhiCov",jidx);
				vvFillBranch(subcl_cov.at(0,2), type, "subclusterEtaTimeCov",jidx);
				vvFillBranch(subcl_cov.at(2,1), type, "subclusterPhiTimeCov",jidx);
				
			}
		}


		void DoGenMatching(const vector<Jet>& jets, map<string, vector<int>>& genmatches){
			genmatches.clear();

			vector<int> genWMatchIdxs_genParts(jets.size(),-1);
			vector<int> genqMatchIdxs_genParts(jets.size(),-1);
			vector<int> genbMatchIdxs_genParts(jets.size(),-1);
			vector<int> genTopMatchIdxs_genParts(jets.size(),-1);
			vector<int> genGluonMatchIdxs_genParts(jets.size(),-1);
			GenericMatchJet(jets,_genparts, genWMatchIdxs_genParts, 24); //match BHC jets to good gen Ws
			GenericMatchJet(jets,_genparts, genqMatchIdxs_genParts, 0); //match BHC jets to good gen qs
			GenericMatchJet(jets,_genparts, genbMatchIdxs_genParts, 5); //match BHC jets to good gen qs
			GenericMatchJet(jets,_genparts, genTopMatchIdxs_genParts, 6); //match BHC jets to good gen qs
			GenericMatchJet(jets,_genparts, genGluonMatchIdxs_genParts, 21); //match BHC jets to good gen qs
			int val = -1;
			for(int g = 0; g < _genmatches.size(); g++){
				if(_genmatches[g] == "") continue;
				genmatches[_genmatches[g]] = {};
				vector<int> vals_global;
				if(_genmatches[g] == "W"){
					vals_global = genWMatchIdxs_genParts;
				}
				else if(_genmatches[g] == "q"){
					vals_global = genqMatchIdxs_genParts;
				}
				else if(_genmatches[g] == "b"){
					vals_global = genbMatchIdxs_genParts;
				}
				else if(_genmatches[g] == "Top"){
					vals_global = genTopMatchIdxs_genParts;
				}
				else if(_genmatches[g] == "gluon"){
					vals_global = genGluonMatchIdxs_genParts;
				}
				else continue;
				for(int j = 0; j < jets.size(); j++){
					genmatches.at(_genmatches[g]).push_back(vals_global[j]);
				}
			}

		}

		void FillPredJetHists(){
			//event observables
			int njets = _predJets.size();
			FillBranch((double)njets, "BHCJet", "nJets");

			//should be able to do in post - BHCJet_WMatched == 1 && BHCJet_subclusterqMatched == 1
			//for e in events:
			//	tree->GetEntry(i);
			//	set branch addresses
			//	njets = BHCJet_nJets;
			//	for j in njets:
			//		if(!BHCJet_WMatched[j]) continue;
			//		nclusters = BHCJet_nSubclustersJet[j];
			//		do subcluster dR matching with daughters of W that BHCJet is matched to (genq_Widx)
			
			map<string, vector<int>> genmatches;
			DoGenMatching(_predJets, genmatches);
			int njets_noPU = 0;
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _predJets.size(); j++){
				FillBranchesJet(_predJets[j], "BHC", j, genmatches);
				//do PU cleaning - downweighted not complete removal
				vector<bool> scores;
				Jet jet_puCleaned = _predJets[j].CleanOutPU(scores,false);
				//if no subclusters in jet pass PU cleaning, "remove" jet (ie don't count it towards nJets_noPU or plot its properties
				if(find(scores.begin(), scores.end(), true) == scores.end())
					continue;
				FillBranchesJet(jet_puCleaned, "BHCnoPU", njets_noPU, genmatches); //use njets_noPU as jet idx since this is incremented in the correct way for noPU jets
				njets_noPU++;
			}
			FillBranch((double)njets_noPU, "BHCnoPUJet", "nJets");

		}


		//fill hists for gen jets and gen particles
		//ONLY GEN PARTICLES CONSIDERED - tops, b's from tops, quarks from W's, direct daughters of b's
		void FillGenParticleHists(){
			for(int g = 0; g < _genparts.size(); g++){
				FillBranchesGen(_genparts[g]); 
			}
		}
	
		void FillGenJetHists(){
			map<string, vector<int>> genmatches;	
			//genAK4
			//event observables
			int njets = _genAK4jets.size();
			FillBranch((double)njets, "genAK4Jet", "nJets");
			DoGenMatching(_genAK4jets, genmatches);
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _genAK4jets.size(); j++){
				FillBranchesJet(_genAK4jets[j], "genAK4", j, genmatches);
			}
			
			//genAK8
			njets = _genAK8jets.size();
			FillBranch((double)njets, "genAK8Jet", "nJets");
			DoGenMatching(_genAK8jets, genmatches);
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _genAK8jets.size(); j++){
				FillBranchesJet(_genAK8jets[j], "genAK8", j, genmatches);
			}
			
			//genAK15
			njets = _genAK15jets.size();
			FillBranch((double)njets, "genAK15Jet", "nJets");
			DoGenMatching(_genAK15jets, genmatches);
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _genAK15jets.size(); j++){
				FillBranchesJet(_genAK15jets[j], "genAK15", j, genmatches);
			}
	
		}


		void FillRecoJetHists(){
			map<string, vector<int>> genmatches;	
			//recoAK4
			//event observables
			int njets = _recoAK4jets.size();
			FillBranch((double)njets, "recoAK4Jet", "nJets");
			DoGenMatching(_recoAK4jets, genmatches);
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _recoAK4jets.size(); j++){
				FillBranchesJet(_recoAK4jets[j], "recoAK4", j, genmatches);
			}
			
			//recoAK8
			njets = _recoAK8jets.size();
			FillBranch((double)njets, "recoAK8Jet", "nJets");
			DoGenMatching(_recoAK8jets, genmatches);
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _recoAK8jets.size(); j++){
				FillBranchesJet(_recoAK8jets[j], "recoAK8", j, genmatches);
			}
			
			//recoAK15
			njets = _recoAK15jets.size();
			FillBranch((double)njets, "recoAK15Jet", "nJets");
			DoGenMatching(_recoAK15jets, genmatches);
			//jet observables (subclusters observables inside)
			for(int j = 0; j < _recoAK15jets.size(); j++){
				FillBranchesJet(_recoAK15jets[j], "recoAK15", j, genmatches);
			}
		}
	
			
		//draw ellispes + tmarkers to canvas and write canvas to file
		//center_coords = center of event display in [eta_c, phi_c] per gen object s.t. center_coords[i] = BayesPoint(eta_c, phi_c) for object i
		//window_width = width of event display in [deta, dphi] per gen object s.t. window_width[i] = BayesPoint(deta, dphi) for object i
		//void WriteEventDisplays(TFile* ofile, vector<BayesPoint> center_coords = {}, vector<BayesPoint> window_width = {}){
		void WriteEventDisplays(TFile* ofile, const map<string,BayesPoint>& center_coords = {}, const map<string,BayesPoint>& window_width = {}){
			if(EvtDisplay_etaCell_phiCell->Integral() == 0){ cout << "skipping all evt disps" << endl; return;} //don't write canvases if this event display isn't filled
			string plot_title;
			if(_sel == singW){
				plot_title = "single W^{#pm}";
				if(_oname.find("Wgluon") != string::npos){
					plot_title = "single W^{#pm}+gluon";
				}
			}
			if(_sel == QCDdijets){
				plot_title = "QCD dijets";
			}
			if(_sel == boostTop){
				plot_title = "t#bar{t}";
			}
			cout << "writing event display hist" << endl;
			ofile->cd();
			//write overall event display
			TCanvas* cv = new TCanvas("evtdisp","evtdisp");
			cv->SetTitle("");
			cv->SetRightMargin(0.15);
			if(_evt2disp_z == 1){ //time
				EvtDisplay_etaCell_phiCell->GetZaxis()->SetTitle("time [ns]");
			}
			else if(_evt2disp_z % 2 == 0 && _evt2disp_z > 0){ //responsibility is event numbers where _evt2disp_z / 2 == k of subcl responsibility to display
				EvtDisplay_etaCell_phiCell->GetZaxis()->SetTitle("responsibility");
			}
			else{
			}
			EvtDisplay_etaCell_phiCell->SetTitle("");
			//do formatting
			EvtDisplay_etaCell_phiCell->GetXaxis()->CenterTitle(true);
			EvtDisplay_etaCell_phiCell->GetYaxis()->CenterTitle(true);
			EvtDisplay_etaCell_phiCell->GetZaxis()->CenterTitle(true);
			EvtDisplay_etaCell_phiCell->GetYaxis()->SetLabelFont(132);
			EvtDisplay_etaCell_phiCell->GetXaxis()->SetLabelFont(132);
			EvtDisplay_etaCell_phiCell->GetZaxis()->SetLabelFont(132);
			EvtDisplay_etaCell_phiCell->GetYaxis()->SetTitleFont(132);
			EvtDisplay_etaCell_phiCell->GetXaxis()->SetTitleFont(132);
			EvtDisplay_etaCell_phiCell->GetZaxis()->SetTitleFont(132);
			EvtDisplay_etaCell_phiCell->GetYaxis()->SetTitleSize(0.04);
			EvtDisplay_etaCell_phiCell->GetXaxis()->SetTitleOffset(1.05);
			EvtDisplay_etaCell_phiCell->GetXaxis()->SetTitleSize(0.04);
			EvtDisplay_etaCell_phiCell->GetZaxis()->SetTitleSize(0.04);
			EvtDisplay_etaCell_phiCell->Draw("colz1");
			//plot jets
			for(int j = 0; j < _predJets.size(); j++){
				_jellipses[j].Draw();
				_jcenters[j].Draw();
				
				for(int k = 0; k < _predJets[j].GetNConstituents(); k++){
					_subclellipses[j][k].Draw();
					_subclcenters[j][k].Draw();
				}
			}
			//plot gen particles
			for(int m = 0; m < _plot_particles.size(); m++){
				_plot_particles[m].Draw();
			}
			//write labels
			string lat_cms = "#font[22]{Pythia 8} #font[132]{event generation, #sqrt{s} = 13 TeV}";
        		TLatex lat;
        		lat.SetNDC();
        		lat.SetTextSize(0.04);
        		lat.SetTextFont(42);
        		lat.DrawLatex(0.03,0.92,lat_cms.c_str());
			TLatex lat2;
        		lat2.SetNDC();
        		lat2.SetTextSize(0.04);
        		lat2.SetTextFont(42);
        		plot_title = "#font[132]{"+plot_title+"}";
			lat2.DrawLatex(0.8,0.92,plot_title.c_str());
			cv->Write();
		
			vector<string> names;
			for(auto it = center_coords.begin(); it != center_coords.end(); it++){
				names.push_back(it->first);
			}
cout << _evtdisps_obj.size() << " # of obj evtdisps" << endl;	
			//write object specific plots
			for(int h = 0; h < _evtdisps_obj.size(); h++){
				//get match string for center + width
				string name = _evtdisps_obj[h]->GetName();
				string objmatch = "EvtDisplay_etaCell_phiCell_";
				name = name.substr(objmatch.size());
				//skip hists that aren't filled
				cout << "hist #" << h << " name of obj evtdisp to be filled " << name << endl;
				if(find(names.begin(), names.end(), name) == names.end()) continue;
		
				BayesPoint center = center_coords.at(name);
				double eta_max = 0;
				double phi_max = 0;
				double eta_min = 999;
				double phi_min = 999;
				BayesPoint width = BayesPoint({eta_max, phi_max});// - set by rhs drawn to get everything in frame = window_width[name];
				BayesPoint max_width = BayesPoint({eta_max, phi_max});// - set by rhs drawn to get everything in frame = window_width[name];
				BayesPoint min_width = BayesPoint({eta_min, phi_min});// - set by rhs drawn to get everything in frame = window_width[name];
cout << "drawing hist #" << h << " of " << _evtdisps_obj.size() << " with name " << name << endl;
				for(int j = 0; j < _predJets.size(); j++){
					//if(find(jet_idxs.begin(), jet_idxs.end(), j) == jet_idxs.end()) continue;
					BayesPoint ell_center({_predJets[j].eta(), _predJets[j].phi()});
					ell_center.Translate(center.at(0),0);
					ell_center.CircularTranslate(center.at(1),1);
					cout << "jet #" << j << " window width eta " << window_width.at(name).at(0) << " phi " << window_width.at(name).at(1) << " this jet center eta " << ell_center.at(0) << " phi " << ell_center.at(1) << endl;
					if(fabs(ell_center.at(0)) > window_width.at(name).at(0)) continue; //out of frame in eta
					if(fabs(ell_center.at(1)) > window_width.at(name).at(1)) continue; //out of frame in phi

cout << "drawing rhs from jet #" << j << endl;

					vector<JetPoint> rhs;
					Jet subcl;
					if(_evt2disp_z % 2 == 0 && _evt2disp_z > 0){
						int k = _evt2disp_z / 2;
						if(k >= _predJets[j].GetNConstituents()){
							_predJets[j].GetJetPoints(rhs);
						}
						else{
							_predJets[j].GetConstituent(k, subcl);
							subcl.GetJetPoints(rhs);
						}
					}
					else{
						_predJets[j].GetJetPoints(rhs);
					}
					for(auto rh : rhs){
						double w;
						if(_evt2disp_z == 0){ //energy
							w = rh.E();
						}
						else if(_evt2disp_z == 1){ //time
							w = rh.t();
						}
						else if(_evt2disp_z % 2 == 0){ //responsibility is event numbers where _evt2disp_z / 2 == k of subcl responsibility to display
							w = rh.GetWeight();
						}
						else{
							w = rh.E(); //default energy weighted
						}
						if(w == 0) continue;
						//center according to main gen particle
						BayesPoint rh_pt({rh.eta(), rh.phi()}); //save as BayesPoint to do correct circular translation to (0,0)
						rh_pt.SetWeight(w);
						//translate into local eta, phi coords
						rh_pt.Translate(center.at(0),0);
						rh_pt.CircularTranslate(center.at(1),1);
						if(rh_pt.at(0) > eta_max)
							eta_max = rh_pt.at(0);
						if(rh_pt.at(1) > phi_max)
							phi_max = rh_pt.at(1);
						if(rh_pt.at(0) < eta_min)
							eta_min = rh_pt.at(0);
						if(rh_pt.at(1) < phi_min)
							phi_min = rh_pt.at(1);
//cout << "at jet #" << j << " rh eta " << rh_pt.at(0) << " rh_phi " << rh_pt.at(1) << " eta_max " << eta_max << " eta_min " << eta_min << " phi_max " << phi_max << " phi_min " << phi_min << endl;
						_evtdisps_obj[h]->Fill(rh_pt.at(0), rh_pt.at(1), rh_pt.w());
					}
				}
//TODO - check eta/phi_max/min values (see drawn canvases)
cout << "eta_max " << eta_max << " eta_min " << eta_min << " phi_max " << phi_max << " phi_min " << phi_min << endl;
				if(_evt2disp_z == 1){ //time
					_evtdisps_obj[h]->GetZaxis()->SetTitle("time [ns]");
				}
				else if(_evt2disp_z % 2 == 0 && _evt2disp_z > 0){ //responsibility is event numbers where _evt2disp_z / 2 == k of subcl responsibility to display
					_evtdisps_obj[h]->GetZaxis()->SetTitle("responsibility");
				}
				else{
				}
				min_width.SetValue(eta_min, 0);
				min_width.SetValue(phi_min, 1);
				max_width.SetValue(eta_max, 0);
				max_width.SetValue(phi_max, 1);
				width.SetValue(max(fabs(eta_max), fabs(eta_min)), 0);
				width.SetValue(max(fabs(phi_max), fabs(phi_min)), 1);
				if(_evtdisps_obj[h]->GetEntries() == 0) continue; //don't draw if not filled for this particle gen obj
				TCanvas* cv_obj = new TCanvas(_evtdisps_obj[h]->GetName(),_evtdisps_obj[h]->GetTitle());
				cv_obj->SetRightMargin(0.15);
				cv_obj->cd();
				_evtdisps_obj[h]->GetXaxis()->SetRangeUser(min_width.at(0), max_width.at(0));
				_evtdisps_obj[h]->GetYaxis()->SetRangeUser(min_width.at(1), max_width.at(1));
				_evtdisps_obj[h]->SetTitle("");
				_evtdisps_obj[h]->GetXaxis()->CenterTitle(true);
				_evtdisps_obj[h]->GetYaxis()->CenterTitle(true);
				_evtdisps_obj[h]->GetZaxis()->CenterTitle(true);
				_evtdisps_obj[h]->GetYaxis()->SetLabelFont(132);
				_evtdisps_obj[h]->GetXaxis()->SetLabelFont(132);
				_evtdisps_obj[h]->GetZaxis()->SetLabelFont(132);
				_evtdisps_obj[h]->GetYaxis()->SetTitleFont(132);
				_evtdisps_obj[h]->GetXaxis()->SetTitleFont(132);
				_evtdisps_obj[h]->GetZaxis()->SetTitleFont(132);
				_evtdisps_obj[h]->GetYaxis()->SetTitleSize(0.04);
				_evtdisps_obj[h]->GetXaxis()->SetTitleOffset(1.05);
				_evtdisps_obj[h]->GetXaxis()->SetTitleSize(0.04);
				_evtdisps_obj[h]->GetZaxis()->SetTitleSize(0.04);
				_evtdisps_obj[h]->Draw("colz1");
cout << "hist for " << name << " integral " << _evtdisps_obj[h]->Integral() << " entries " << _evtdisps_obj[h]->GetEntries() << endl;
				//do for gen particles too 
				for(int m = 0; m < _plot_particles.size(); m++){
					BayesPoint m_center({_plot_particles[m].GetX(), _plot_particles[m].GetY()});
					m_center.Translate(center.at(0),0);
					m_center.CircularTranslate(center.at(1),1);
				
					//double dr = dR(m_center.at(0), m_center.at(1), 0., 0.);
					//if(dr > sqrt(width.at(0)*width.at(0) + width.at(1)*width.at(1))) continue;
					if(m_center.at(0) > max_width.at(0) || m_center.at(0) < min_width.at(0)) continue;
					if(m_center.at(1) > max_width.at(1) || m_center.at(1) < min_width.at(1)) continue;
					_plot_particles[m].DrawMarker(m_center.at(0), m_center.at(1));

				}
				for(int j = 0; j < _predJets.size(); j++){
					//if(find(jet_idxs.begin(), jet_idxs.end(), j) == jet_idxs.end()) continue;
					BayesPoint ell_center({_jellipses[j].GetX1(), _jellipses[j].GetY1()});
					ell_center.Translate(center.at(0),0);
					ell_center.CircularTranslate(center.at(1),1);
					cout << "ellipse (jet) #" << j << " window width eta " << width.at(0) << " phi " << width.at(1) << " this ellipse center eta " << ell_center.at(0) << " phi " << ell_center.at(1) << endl;
					//make sure ellipse center is in window
					if(ell_center.at(0) > max_width.at(0)) continue; //out of frame in eta
					if(ell_center.at(1) > max_width.at(1)) continue; //out of frame in eta
					if(ell_center.at(0) < min_width.at(0)) continue; //out of frame in eta
					if(ell_center.at(1) < min_width.at(1)) continue; //out of frame in eta
					
					double ell_center_eta = ell_center.at(0);//_jellipses[j].GetX1(); 
					double ell_center_phi = ell_center.at(1);//_jellipses[j].GetY1();
					double ell_maj_r = _jellipses[j].GetR1();
					double ell_min_r = _jellipses[j].GetR2();
					double theta = _jellipses[j].GetTheta();
					//put in rad
					theta *= acos(-1)/180;
					double r_eta = ell_maj_r*cos(theta);
					double r_phi = ell_maj_r*sin(theta);
					//if full ellipse cannot be drawn in window, skip
					double dr = dR(ell_center.at(0), ell_center.at(1), 0., 0.);
					if(dr > sqrt(width.at(0)*width.at(0) + width.at(1)*width.at(1))) continue;
					
					//if(ell_center.at(0) + r_eta > max_width.at(0)) continue; 
					//if(ell_center.at(0) - r_eta < min_width.at(0)) continue; 

					//if(ell_center.at(1) + r_phi > max_width.at(1)) continue; 
					//if(ell_center.at(1) - r_phi < min_width.at(1)) continue; 
					//if(fabs(r_eta) > fabs(width.at(0))) continue;
					//if(fabs(r_phi) > fabs(width.at(1))) continue;

				cout << "drawing ellipse (jet) #" << j << endl; 

		
					_jellipses[j].DrawEllipse(ell_center.at(0), ell_center.at(1), ell_maj_r, ell_min_r, 0, 360, _jellipses[j].GetTheta());
				

					
					_jcenters[j].DrawMarker(ell_center.at(0), ell_center.at(1));
				
					for(int k = 0; k < _predJets[j].GetNConstituents(); k++){
						BayesPoint subcl_center({_subclellipses[j][k].GetX1(), _subclellipses[j][k].GetY1()});
						subcl_center.Translate(center.at(0),0);
						subcl_center.CircularTranslate(center.at(1),1);
						ell_center_eta = _subclellipses[j][k].GetX1(); 
						ell_center_phi = _subclellipses[j][k].GetY1();
						ell_maj_r = _subclellipses[j][k].GetR1();
						ell_min_r = _subclellipses[j][k].GetR2();
						theta = _subclellipses[j][k].GetTheta();
						
						_subclellipses[j][k].DrawEllipse(subcl_center.at(0), subcl_center.at(1), ell_maj_r, ell_min_r, 0, 360, theta);
					
						_subclcenters[j][k].DrawMarker(subcl_center.at(0), subcl_center.at(1));
					}

				}

        			TLatex lat1;
        			lat1.SetNDC();
        			lat1.SetTextSize(0.04);
        			lat1.SetTextFont(42);
        			lat1.DrawLatex(0.03,0.92,lat_cms.c_str());
				TLatex lat3;
        			lat3.SetNDC();
        			lat3.SetTextSize(0.04);
        			lat3.SetTextFont(42);
        			lat3.DrawLatex(0.8,0.92,plot_title.c_str());
				cv_obj->SetTitle("");
				cv_obj->Write();

			}


		}

		void WriteOutput(TFile* ofile, const map<string,BayesPoint>& center_coords = {}, const map<string,BayesPoint>& window_width = {}){
			ofile->cd();
			EvtDisplay_etaCell_phiCell->Write();
			WriteEventDisplays(ofile, center_coords, window_width);
			ofile->cd();
			ofile->Close();
		}	
		
		//44 - 129 - eta-phi event display of rechits for specified _evt2disp with cell energy on the z axis (overall event)
		TH2D* EvtDisplay_etaCell_phiCell = new TH2D("EvtDisplay_etaCell_phiCell","EvtDisplay_etaCell_phiCell;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,0,8*atan(1));

		//event display histograms (per gen object)
		vector<TH2D*> _evtdisps_obj;
		//these should be centered on their respective gen particles
		//0 - single W, W+gluon, first W in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_W = new TH2D("EvtDisplay_etaCell_phiCell_W","EvtDisplay_etaCell_phiCell_W;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//1 - second W in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_W2 = new TH2D("EvtDisplay_etaCell_phiCell_W2","EvtDisplay_etaCell_phiCell_W2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//2 - gluon in W+gluon
		TH2D* EvtDisplay_etaCell_phiCell_gluon = new TH2D("EvtDisplay_etaCell_phiCell_gluon","EvtDisplay_etaCell_phiCell_gluon;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//3 - q1 in QCD dijets
		TH2D* EvtDisplay_etaCell_phiCell_q1 = new TH2D("EvtDisplay_etaCell_phiCell_q1","EvtDisplay_etaCell_phiCell_q1;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//4 - q2 in QCD dijets
		TH2D* EvtDisplay_etaCell_phiCell_q2 = new TH2D("EvtDisplay_etaCell_phiCell_q2","EvtDisplay_etaCell_phiCell_q2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//5 - b1 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_b1 = new TH2D("EvtDisplay_etaCell_phiCell_b1","EvtDisplay_etaCell_phiCell_b1;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//6 - b2 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_b2 = new TH2D("EvtDisplay_etaCell_phiCell_b2","EvtDisplay_etaCell_phiCell_b2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//7 - top1 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_top1 = new TH2D("EvtDisplay_etaCell_phiCell_top1","EvtDisplay_etaCell_phiCell_top1;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));
		//8 - top2 in ttbar
		TH2D* EvtDisplay_etaCell_phiCell_top2 = new TH2D("EvtDisplay_etaCell_phiCell_top2","EvtDisplay_etaCell_phiCell_top2;Pseudorapidity (#eta);Azimuthal angle (#phi);Energy [GeV]",344,-3,3,360,-4*atan(1),4*atan(1));


		void SetSmear(bool t){ _smear = t; }
		double _cell, _tresCte, _tresNoise, _tresStoch;
		void SetMeasErrParams(double spatial, double tresCte, double tresStoch, double tresNoise){ _cell = spatial; _tresCte = tresCte; _tresStoch = tresStoch; _tresNoise = tresNoise; 
	cout << "Using tres_cte = " << _tresCte << " ns, tres_stoch = " << _tresStoch << " ns and tres_noise = " << _tresNoise << endl;
 }
		void SetOutfile(string fname){ _oname = fname; }
		void SetTransferFactor(double gev){
			_gev = gev;
			_prod->SetTransferFactor(_gev);
		}
		void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj; }
		void SetVerbosity(int verb){_verb = verb;}
		int _verb;

	//find gen object that most closely is dr matched to jet
	//when gen parts are used, need to specify id to match to (ie W - 24, top - 6, any light quark - 0)
	void GenericMatchJet(const vector<Jet>& injets, const vector<Jet>& matchjets, vector<int>& bestMatchIdxs, int id = -1){
		//loop through gen particles
		//dr match to jet
		double bestDr, dr;
		vector<Jet> matchobjs;
		//map 'local' match obj idxs to 'global' matchjets idxs
		map<int,int> matchIdxToObjIdx;
		//if id is specified, make matchobjs the genparts that satisfy this criteria
		if(id != -1){
			vector<int> qids = {1,2,3,4};
			for(int g = 0; g < matchjets.size(); g++){
				int idx = matchjets[g].GetUserIdx();
				int this_id = fabs(_base->genpart_id->at(idx));
				if(id != 0 && this_id != id) continue;
				if(id == 0){
					if(find(qids.begin(), qids.end(), this_id) == qids.end()) continue;
				}
				matchobjs.push_back(matchjets[g]);
				matchIdxToObjIdx[matchobjs.size()-1] = g;
			}
			//add no match option
			matchIdxToObjIdx[-1] = -1;
		}
		else{ matchobjs = matchjets; }

		int nMatch = matchobjs.size();
		bestMatchIdxs.clear();
		bestMatchIdxs = {};
	
		//no gen particles to match, set all jets to unmatched	
		if(nMatch < 1){
			for(auto j : injets) bestMatchIdxs.push_back(-1);
			return;
		}
		if(injets.size() < 1) return;

		//vector<double> drs;
		//drs[i][j] = dr for jet i and gen particle j
		vector<vector<double>> drs;
		vector<int> idxs;
		for(int j = 0; j < injets.size(); j++){
			drs.push_back({});
			for(int g = 0; g < nMatch; g++){
				drs[j].push_back(999);
				//cout << "jet #" << j << " eta " << injets[j].eta() << " phi " << injets[j].phi() << " matchobj #" << g << " eta " << matchobjs[g].eta() << " phi " << matchobjs[g].phi() << endl;
				dr = dR(matchobjs[g].eta(), matchobjs[g].phi(), injets[j].eta(), injets[j].phi());
				//cout << "jet # " << j << " matchobj # " << g << " dr " << dr << endl;		
				drs[j][g] = dr;
			}
		}
		vector<int> best_idxs; //one per jet
		int otherJet, thismatchidx, thisJet;
		//go back through jets can check to see if there are overlapping matches
		for(int j = 0; j < injets.size(); j++){
			//for(int g = 0; g < nMatch; g++){
			//	cout << "jet " << j << " and match jet " << g << " have dr " << drs[j][g] << endl;
			//}
			//cout << "jet " << j << " has best dr " << *min_element(drs[j].begin(), drs[j].end()) << " at match jet " << find(drs[j].begin(), drs[j].end(), *min_element(drs[j].begin(), drs[j].end())) - drs[j].begin() << endl;
			double mindr = *min_element(drs[j].begin(), drs[j].end());
			int matchidx = find(drs[j].begin(), drs[j].end(), mindr) - drs[j].begin();
			if(mindr == 999) matchidx = -1; //no match found (ie no available match jet for best match)
			best_idxs.push_back(matchidx);
			thismatchidx = matchidx;
			thisJet = j;
			//if other jets have the same genidx matched, go through and disambiguate until there is only 1 instance of genidx
			while(count(best_idxs.begin(), best_idxs.end(), matchidx) > 1 && matchidx != -1){
				otherJet = find(best_idxs.begin(), best_idxs.end(), matchidx) - best_idxs.begin();
				//this happens if the "otherJet" to be analyzed comes before thisjet (ie it gets found first)
				//skip otherJet (ie thisJet) in this case and look at all other jets
				if(otherJet == thisJet){
					otherJet = find(best_idxs.begin()+otherJet+1, best_idxs.end(), matchidx) - best_idxs.begin();
				}
				//for(int b = 0; b < best_idxs.size(); b++) cout << "b " << b << " bestidx " << best_idxs[b] << endl;
				//cout << " found another match at jet " << otherJet << " with other dr " << drs[otherJet][matchidx] << " against this jet " << thisJet << endl;
				//if other dr is less than current mindr
				if(drs[otherJet][matchidx] < mindr){
					//set this dr to 999 (is invalid), find new min for this jet, reset genidx to this index
					drs[thisJet][matchidx] = 999;
					mindr = *min_element(drs[thisJet].begin(), drs[thisJet].end());
					if(mindr == 999) matchidx = -1;
					else matchidx = find(drs[thisJet].begin(), drs[thisJet].end(), mindr) - drs[thisJet].begin();
					best_idxs[thisJet] = matchidx;
					//cout << " reset gen match of this jet " << thisJet << " to gen jet " << matchidx << " with dr " << mindr << endl;
		
				}
				//if this dr is less than (or equal to) current mindr
				else{
					//set other dr to 999 (is invalid), find new min for other jet, reset other genidx to index of new mind
					drs[otherJet][matchidx] = 999;
					thismatchidx = matchidx;
					mindr = *min_element(drs[otherJet].begin(), drs[otherJet].end());
					if(mindr == 999) matchidx = -1;
					else matchidx = find(drs[otherJet].begin(), drs[otherJet].end(), mindr) - drs[otherJet].begin();
					thisJet = otherJet;
					best_idxs[thisJet] = matchidx;
					//cout << " reset match of other jet " << otherJet << " to jet " << matchidx << " with dr " << mindr << endl;
				}	
				//cout << "matchidx is now " << matchidx << " with count " << count(best_idxs.begin(), best_idxs.end(), matchidx) << " for jet " << thisJet << endl;

			}
			//cout << "jet " << j << " has best exclusive match with " << best_idxs[j] << "\n" << endl;
		}
		bestMatchIdxs = best_idxs;
		//map 'local' match obj index to 'global' all obj index
		if(!matchIdxToObjIdx.empty()){
			for(int i = 0; i < bestMatchIdxs.size(); i++){
				//cout << "# " << i << " best match idx " << bestMatchIdxs[i] << " best obj idx " << matchIdxToObjIdx[bestMatchIdxs[i]] << endl;
				bestMatchIdxs[i] = matchIdxToObjIdx[bestMatchIdxs[i]];
			}
		}
	}



	void SetAlpha(double a){_alpha = a;}
	void SetSubclusterAlpha(double a){_emAlpha = a; }
	void SetThreshold(double t){ _thresh = t; }
	void SetPriorParameters(const map<string, Matrix>& params){_prior_params = params;} 		
	void SetReducePU(bool t){ _zoom_window = t; if(_zoom_window) cout << "Reducing PU with zoom window." << endl; }//draws rectangle around hard scattering to reduce # of rechits to cluster

	void _reset(){
		for(auto it = _obs.begin(); it != _obs.end(); it++){
			it->second = -1;
		}
		for(auto it = _vobs.begin(); it != _vobs.end(); it++){
			it->second.clear();
		}
		for(auto it = _vvobs.begin(); it != _vvobs.end(); it++){
			for(auto v : it->second){
				v.clear();
			}
			it->second.clear();
		}

	}

	private:
		TFile* _infile;
		string _oname;
		vector<TH1D*> _hists1D;
		vector<TH2D*> _hists2D;
		vector<TEllipse> _ellipses;
		vector<TMarker> _plot_particles;
		map<int, TEllipse> _jellipses;
		map<int, TMarker> _jcenters;
		map<int,map<int,TEllipse>> _subclellipses;
		map<int,map<int,TMarker>> _subclcenters;
		

		//value maps
		map<string, double> _obs;
		map<string, vector<double>> _vobs;
		map<string, vector<vector<double>>> _vvobs;
		vector<string> _types = {"BHC","BHCnoPU","genAK4", "genAK8", "genAK15","recoAK4", "recoAK8", "recoAK15"};
		vector<string> _genmatches = {"W", "Top", "b", "q", "gluon"};
		
	
		//ttree member variables
		TTree* _tree;

		vector<Jet> _phos; //photons for event
		vector<std::shared_ptr<BaseTree::node>> _trees;
		vector<Jet> _predJets, _genAK4jets, _recoAK4jets, _genparts, _genAK15jets, _genAK8jets, _recoAK8jets, _recoAK15jets;
		double _minTopPt, _minTopE, _minWPt;
		bool _smear;
		enum Strategy{
			//Delauney strategy - NlnN time - for 2pi cylinder
			NlnN = 0,
			//traditional strategy - N^2 time
			N2 = 1,
			//gmm only on reco'd jets
			gmmOnly = 2,
			//NlnN using rhs from reco AK4 jets
			NlnNonAK4 = 3,
		};
		//clustering strategy - N^2 or NlnN
		Strategy _strategy;
		//event selection type
		enum EvtSel{
			//default selection (generic hadronic)
			def = 0,	
			//boosted W selection
			singW = 1,
			//boosted top selection
			boostTop = 2,
			//QCD dijets selection
			QCDdijets = 3
		};
		EvtSel _sel;
		ReducedBaseSim* _base = nullptr;
		int _nEvts;
		std::unique_ptr<JetSimProducer> _prod;
		bool _data, _debug, _zoom_window;
		int _evti, _evtj;
		double _gev;
		double _c = 29.9792458; // speed of light in cm/ns
		double _radius; //radius of detector set by rhs in event (used for constructing jets)
		double _pvx, _pvy, _pvz;	
		double _alpha, _emAlpha, _thresh;
		//top decay info - 0 = fully had, 1 = semi lep, 2 = fully lep
		int _topDecayType;

		double _pt_thresh = 175.;


		//good gen objects available for matching 
		vector<Jet> _genW, _genb, _genTop, _genq, _genglu;
		int _evt2disp, _evt2disp_z;

		map<string, Matrix> _prior_params;

		double dR(double eta1, double phi1, double eta2, double phi2){
			//phi wraparound
			double dphi = (phi1-phi2);
			dphi = acos(cos(dphi));
			return sqrt((eta1-eta2)*(eta1-eta2) + dphi*dphi);
		}

		double CalcSize(const Matrix& cov, bool time = false){
			if(cov.nRows() != 3 || cov.nCols() != 3){
				cout << "Error: can't calculate size for matrix of size " << cov.nRows() << " x " << cov.nCols() << endl;
				return -1;
			}
			vector<double> eigvals;
			vector<Matrix> eigvecs;
			if(time){
				cov.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[2]);
			}
			else{
				Matrix cov2D(2,2);
				Get2DMat(cov,cov2D);	
				cov2D.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[1]);
			}
		}

		double CalcSubleadSize(const Matrix& cov, bool time = false){
			if(cov.nRows() != 3 || cov.nCols() != 3){
				cout << "Error: can't calculate size for matrix of size " << cov.nRows() << " x " << cov.nCols() << endl;
				return -1;
			}
			vector<double> eigvals;
			vector<Matrix> eigvecs;
			if(time){
				cov.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[1]);
			}
			else{
				Matrix cov2D(2,2);
				Get2DMat(cov,cov2D);	
				cov2D.eigenCalc(eigvals, eigvecs);
				//define jet size as length of major axis
				return sqrt(eigvals[0]);
			}
		}


		void Get2DMat(const Matrix& inmat, Matrix& outmat){
			if(!outmat.square()) return;
			if(outmat.nRows() != 2) return;
			outmat.reset();
			outmat.SetEntry(inmat.at(0,0),0,0);	
			outmat.SetEntry(inmat.at(0,1),0,1);	
			outmat.SetEntry(inmat.at(1,0),1,0);	
			outmat.SetEntry(inmat.at(1,1),1,1);
		}
		double Rotundity(const Matrix& cov){
			if(cov.nRows() != 3 || cov.nCols() != 3){
				cout << "Error: can't calculate rotundity for matrix of size " << cov.nRows() << " x " << cov.nCols() << endl;
				return -1;
			}
			Matrix cov2D(2,2);
			Get2DMat(cov,cov2D);	
			vector<Matrix> eigenvecs;
			vector<double> eigenvals;
			cov.eigenCalc(eigenvals, eigenvecs);
			int maxd = cov.nRows() - 1;
			double rot = 0;
			for(int i = 0; i < (int)eigenvals.size(); i++) rot += eigenvals[i];
			rot = eigenvals[maxd]/rot;
			//if(rot < 0.5 || rot > 1) cout << "rot: " << rot << endl;
			return rot;
		}


		TEllipse PlotEll(const Jet& jet){
			Matrix mu, cov;
			jet.GetClusterParams(mu, cov);
			//ellipse center
			double eta = mu.at(0,0); 
			double phi = mu.at(1,0);

			//get 2D matrix for jet size
			Matrix cov2D(2,2);
			Get2DMat(cov,cov2D);	
		
			double var_eta = cov2D.at(0,0);
			double var_phi = cov2D.at(1,1);
			double covar = cov2D.at(0,1);	
			//use analytic form of eigenvalues for 2x2 matrix - so correct r1, r2 can be passed to TEllipse
			double temp_a  = 0.5*(var_eta + var_phi);
        		double temp_b  = 0.5*sqrt((var_eta-var_phi)*(var_eta-var_phi) + 4*covar*covar);
        		double lambda_eta = temp_a + temp_b;
        		double lambda_phi = temp_a - temp_b;
			double x_r = sqrt(lambda_eta);
			double y_r = sqrt(lambda_phi); 
			
			double theta = 0;
			if(covar > 0) theta = atan((lambda_eta - var_eta)/covar);
			else theta = -atan((lambda_eta - var_eta)/-covar);
			theta = 180 * theta/(4*atan(1)); //put to degrees
			TEllipse el = TEllipse(eta, phi, x_r, y_r, 0, 360, theta);	
			return el;
		}

		
};
#endif
