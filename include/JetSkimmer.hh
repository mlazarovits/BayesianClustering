#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include "JetPoint.hh"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include "BayesCluster.hh"
#include "JetProducer.hh"
#include "PhotonProducer.hh"
#include <TFile.h>
#include <TGraph.h>
#include "TSystem.h"
#include "BaseTree.hh"

using node = BaseTree::node;
using procCat = BaseSkimmer::procCat;
class JetSkimmer : public BaseSkimmer{
	public:
		JetSkimmer(){
			_evti = 0;
			_evtj = 0;
			_gev = 1./10.;
		};
		virtual ~JetSkimmer(){ };

		//get rechits from file to cluster
		JetSkimmer(TFile* file) : BaseSkimmer(file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time		

			_prod = new JetProducer(file);
			_prod->SetIsoCut();
		

			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/jet_skims_"+_cms_label+".root";
			_gev = 1./10.;
				
			objE_clusterE->SetTitle("jetE_clusterE");
			objE_clusterE->SetName("jetE_clusterE");
			//true jet hists
			_hists1D.push_back(nTrueJets);
			_hists1D.push_back(TrueJet_pT); 
			_hists1D.push_back(TrueJet_nRhs); 
			_hists1D.push_back(TrueJet_EmE); 
			_hists1D.push_back(TrueJet_nConstituents);
			_hists1D.push_back(TrueJet_twoHardestpT);	
			_hists1D.push_back(nSubClusters_mm);
			
			_timeHists1D.push_back(PVtime);
			_timeHists1D.push_back(deltaT_jet);
			_timeHists1D.push_back(deltaT_pvGam);	
			//key figures of merit	
			_timeHists1D.push_back(diffDeltaT_recoGen);
			_timeHists1D.push_back(deltaT_pvGam_gen);	
			_timeHists1D.push_back(gamTime);
			_timeHists1D.push_back(geoEavg_sigmaDeltaTime_recoGen);
			_timeHists1D.push_back(geopTavg_sigmaDeltaTime_dijets);
			_timeHists1D.push_back(minpT_sigmaDeltaTime_dijets);
			_timeHists1D.push_back(Erh_sigmaDeltaTime_dijets);

			_timeHists2D.push_back(geoEavg_diffDeltaTime_recoGen);
			_timeHists2D.push_back(geopTavg_diffDeltaTime_dijets);	
			_timeHists2D.push_back(minpT_diffDeltaTime_dijets);	
			_timeHists2D.push_back(Erh_diffDeltaTime_dijets);	



			//_hists2D.push_back(erhs_trhs);		
			
			//predicted jets - from BHC
			//_hists1D.push_back(nClusters);
			//_hists1D.push_back(rhTime);
			//_hists1D.push_back(comptime);
			//_hists1D.push_back(PVtime_median_pred);
			//_hists1D.push_back(PVtime_eAvg_pred);
			//_hists1D.push_back(PVtime_mmAvg_pred);
			//_hists1D.push_back(PVdeltaT_jet_median_pred);
			//_hists1D.push_back(PVdeltaT_jet_eAvg_pred);
			//_hists1D.push_back(PVdeltaT_jet_mmAvg_pred);
			//_hists2D.push_back(e_nRhs);

			nrhs_comptime->SetName("nrhs_comptime");
			nrhs_comptime->SetTitle("nrhs_comptime");
			graphs.push_back(nrhs_comptime);


		};
		//ctor from rec hit collection - integrating into ntuplizer
		//needs to match order of timeRecoCats
		enum TimeStrategy{
			med = 0, 
			mmavg = 1, 
			eavg = 2
		};		
		//struct for different types of time reco (ie median, eAvg, mmAvg)
		struct timeRecoCat{
			//vector<TH1D*> hists1D;
			//vector<TH2D*> hists2D;
			
			string methodName;
			//histograms stored in each procCat as vector<vector<TH1D>>
			//procCats[p].hists1D[0] = nominal hists
			//procCats[p].hists1D[1] = leading hists (not included)
			//procCats[p].hists1D[2] = subleading hists (not included)
		
			vector<procCat> procCats;
			timeRecoCat(const vector<TH1D*>& in1dhists, const vector<TH2D*>& in2dhists, const TimeStrategy& ts, vector<procCat> pcs){
				//reset proc cat hists to these input hists
				for(int i = 0; i < pcs.size(); i++){
					procCats.push_back(pcs[i]);
					procCats[i].SetHists(in1dhists, in2dhists, false);
				}
				if(ts == med)
					methodName = "median";
				else if(ts == eavg)
					methodName = "eAvg";
				else if(ts == mmavg)
					methodName = "mmAvg";
				else
					cout << "Error: time strategy " << ts << " not supported." << endl;		
				string name;
				string addname;
				//for each histogram (variable or correlation)
				//	//make sure they have the right add-on name (ie leading, !lead, etc)
				//	TH1D* hist = (TH1D*)in1dhists[i]->Clone();
				//	name = hist->GetName();
				//	if(!methodName.empty()) name += "_"+methodName;
				//	hist->SetName(name.c_str());
				//	if(!methodName.empty()) hist->SetTitle(methodName.c_str());
				//	hists1D.push_back(hist);
					//do by process breakdown of sigma plots
					//if(name.find("sigma") != string::npos){
					for(int j = 0; j < procCats.size(); j++){
						for(int i = 0; i < (int)in1dhists.size(); i++){
							TH1D* hist = procCats[j].hists1D[0][i];
							name = hist->GetName();
							if(!procCats[j].plotName.empty()) name = name.substr(0,name.find("_"+procCats[j].plotName));
							if(!methodName.empty()) addname = "_"+methodName;
							if(!procCats[j].plotName.empty()) addname += "_"+procCats[j].plotName;
							hist->SetName((name+addname).c_str());
							if(!methodName.empty()) hist->SetTitle((addname.substr(1)).c_str());
							///hists1D.push_back(hist);
							//cout << "added hist " << hist->GetName() << endl;
						}
					//}
				}
				//for(auto hist : hists1D) cout << hist->GetName() << " " << hist->GetTitle() << endl;
//cout << "timeRecoCat n hists: " << procCats[0].hists1D[0].size() <<  " " << procCats[0].hists1D.size() << endl;
				//for each histogram
					for(int j = 0; j < procCats.size(); j++){
						for(int i = 0; i < (int)in2dhists.size(); i++){
							TH2D* hist = procCats[j].hists2D[0][i];
							name = hist->GetName();
							if(!procCats[j].plotName.empty()) name = name.substr(0,name.find("_"+procCats[j].plotName));
							if(!methodName.empty()) addname = "_"+methodName;
							if(!procCats[j].plotName.empty()) addname += "_"+procCats[j].plotName;
							hist->SetName((name+addname).c_str());
							if(!methodName.empty()) hist->SetTitle((addname.substr(1)).c_str());
							///hists2D.push_back(hist);
							//cout << "added hist " << hist->GetName() << " " << name << " " << addname << endl;
						}
					//}
				}
//cout << "2d timeRecoCat n hists: " << procCats[0].hists2D[0].size() <<  " " << procCats[0].hists2D.size() << endl;

			}
		
			void AddHist(TH1D* inhist){
				string name, addname;
				for(int p = 0; p < procCats.size(); p++){
					TH1D* hist = (TH1D*)inhist->Clone();
					name = hist->GetName();
					addname = "";
					if(!methodName.empty() && name.find(methodName) == string::npos) addname = "_"+methodName;
					if(!procCats[p].plotName.empty() && name.find(procCats[p].plotName) == string::npos) addname += "_"+procCats[p].plotName;
					hist->SetName((name+addname).c_str());
					if(!methodName.empty() && name.find(procCats[p].plotName) == string::npos) hist->SetTitle((addname.substr(1)).c_str());
					procCats[p].hists1D[0].push_back(hist);
				}
			}	
			void AddHist(TH2D* inhist){
				string name, addname;
				for(int p = 0; p < procCats.size(); p++){
					TH2D* hist = (TH2D*)inhist->Clone();
					name = hist->GetName();
					if(!methodName.empty() && name.find(methodName) == string::npos) addname = "_"+methodName;
					if(!procCats[p].plotName.empty() && name.find(procCats[p].plotName) == string::npos) addname += "_"+procCats[p].plotName;
					hist->SetName((name+addname).c_str());
					if(!methodName.empty() && name.find(procCats[p].plotName) == string::npos) hist->SetTitle((addname.substr(1)).c_str());
					procCats[p].hists2D[0].push_back(hist);
				}

			}	


		};

		void Skim();
		vector<TH1D*> _timeHists1D;	
		vector<TH2D*> _timeHists2D;	

		//true jet hists
		TH1D* nTrueJets = new TH1D("nTrueJets","nTrueJets",20,0,20);
		TH1D* rhTime = new TH1D("rhTime","rhTime",100,-30,30); 
		TH1D* TrueJet_pT = new TH1D("TrueJet_pT","TrueJet_pT",100,0,1000);
		TH1D* TrueJet_nRhs = new TH1D("TrueJet_nRhs","TrueJet_nRhs",25,0,100);
		TH1D* TrueJet_EmE = new TH1D("TrueJet_EmE","TrueJet_EmE",50,0,600);
		TH1D* TrueJet_nConstituents = new TH1D("TrueJet_nConstituents","TrueJet_nConstituents",20,0,50);
		TH1D* TrueJet_twoHardestpT =  new TH1D("TrueJet_twoHardestpT","TrueJet_twoHardestpT",100,0,1000);	
		TH1D* nSubClusters_mm = new TH1D("nSubClusters_mm","nSubClusters_mm",20,0,20);

		TH2D* e_nRhs = new TH2D("e_nRhs","e_nRhs",100,0,500,100,0,100);
		TH2D* erhs_trhs = new TH2D("erhs_trhs","erhs_trhs",100,0,4,100,-100,100);
		
		
		///////////////////// timeHists /////////////
		//0 - pv time
		TH1D* PVtime = new TH1D("PVtime", "PVtime",50,-10,10);	
		//1 - delta t between jets (pv time frame)
		TH1D* deltaT_jet = new TH1D("deltaT_jet", "deltaT_jet",50,-4,4);	
		//2 - reco delta t between pv and photon 
		TH1D* deltaT_pvGam = new TH1D("deltaT_pvGam_reco","deltaT_pvGam_reco",25,-4,4);	
		//3 - difference in deltaT_pvGam between gen and reco
		TH1D* diffDeltaT_recoGen = new TH1D("diffDeltaT_recoGen","diffDeltaT_recoGen",50,-3,3);
		//4 - gen deltaT bw photon and pv
		TH1D* deltaT_pvGam_gen = new TH1D("deltaT_pvGam_gen","deltaT_pvGam_gen",25,-10,10);	
		//5 - photon time
		TH1D* gamTime = new TH1D("gamTime_reco", "gamTime_reco",100,-10,10);	

		//these stay empty to be filled later (after hadding)	
		//6 - resolution of difference in reco - gen deltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_sigmaDeltaTime_recoGen = new TH1D("geoEavg_sigmaDeltaTime_recoGen","geoEavg_sigmaDeltaTime_recoGen",10,0,1000);
		//7 - resolution of difference between two jets for PV time as a function of their geo avg pT 	
		TH1D* geopTavg_sigmaDeltaTime_dijets = new TH1D("geopTavg_sigmaDeltaTime_dijets","geopTavg_sigmaDeltaTime_dijets",10,0,1000);
		//8 - resolution of difference between two jets for PV time as a function of the min pT 	
		TH1D* minpT_sigmaDeltaTime_dijets = new TH1D("minpT_sigmaDeltaTime_dijets","minpT_sigmaDeltaTime_dijets",10,0,1000);
		//9 - resolution of difference between two jets for PV time as a function of the sum ECAL energy	
		TH1D* Erh_sigmaDeltaTime_dijets = new TH1D("Erh_sigmaDeltaTime_dijets","Erh_sigmaDeltaTime_dijets",10,0,1000);
	
		//0 - 2D histogram for reco-gen resolution
		//may need to adjust binning/windows here
		TH2D* geoEavg_diffDeltaTime_recoGen = new TH2D("geoEavg_diffDeltaTime_recoGen","geoEavg_diffDeltaTime_recoGen;#sqrt{E^{pho}_{rh} #times E^{jets}_{rh}} (GeV);#Delta t^{PV,#gamma}_{reco, gen} (ns)",10,0,1000,10,-10,10);
		//1 - 2D histogram for dijets resolution - geometric avg of jet pT
		TH2D* geopTavg_diffDeltaTime_dijets = new TH2D("geopTavg_diffDeltaTime_dijets","geopTavg_diffDeltaTime_dijets;#sqrt{pT^{jet1} #times pT^{jet2}} (GeV); #Delta t^{PV}_{dijet}",10,0,1000,10,-10,10);	
		//2 - 2D histogram for dijets resolution - min E of jets
		TH2D* minpT_diffDeltaTime_dijets = new TH2D("minpT_diffDeltaTime_dijets","minpT_diffDeltaTime_dijets;min(pT^{jet1}, pT^{jet2}) (GeV); #Delta t^{PV}_{dijet}",10,0,1000,10,-10,10);	
		//3 - 2D histogram for dijets resolution - sum_rh E_rh of jets
		TH2D* Erh_diffDeltaTime_dijets = new TH2D("Erh_diffDeltaTime_dijets","Erh_diffDeltaTime_dijets;sum_{j} sum_{rh} E^{jet_{j}}_{rh} (GeV); #Delta t^{PV}_{dijet}",10,0,1000,10,-10,10);	

		//comparing predicted jets + true jets
		//TH2D* nSubClusters_nConstituents = new TH2D("nSubClusters_nConstituents", "nSubClusters_nConstituents",50,0,20,50,0,20);

		//predicted jet plots
		TH1D* nClusters = new TH1D("nClusters","nClusters",20,0,20);
		TH1D* PVtime_pred = new TH1D("PVtime_pred","PVtime_pred",100,-10,10);	
		//difference in tPV between two back-to-back jets
		//previous time definitions
		TH1D* deltaT_jet_pred = new TH1D("PVdeltaT_jet_pred","PVdeltaT_jet_pred",100,-10,10);	
		//comp time distribution
		TH1D* comptime = new TH1D("comptime","comptime",100,0,300);
		//comp time as a function of number of rechits per event
		TGraph* nrhs_comptime = new TGraph();


		vector<timeRecoCat> trCats;
		void MakeTimeRecoCatHists(){
			//don't want to separate lead/not lead histograms
			MakeProcCats(_oname, false);
			timeRecoCat trmed(_timeHists1D, _timeHists2D, med, _procCats);
			timeRecoCat treavg(_timeHists1D, _timeHists2D, eavg, _procCats);
			timeRecoCat trmmavg(_timeHists1D, _timeHists2D, mmavg, _procCats);
		
			trCats.push_back(trmed);
			trCats.push_back(trmmavg);	
			trCats.push_back(treavg);
		}

	

		void FillTrueJetHists(const vector<Jet>& jets){
			int njets = jets.size();	
			nTrueJets->Fill((double)njets);
		
			double eECAL = 0;
			int ijet = 0;
			vector<JetPoint> rhs;
			for(int j = 0; j < njets; j++){
				e_nRhs->Fill(jets[j].E(), jets[j].GetNRecHits());
				objE->Fill(jets[j].E());
				TrueJet_pT->Fill(jets[j].pt());
				TrueJet_nRhs->Fill(jets[j].GetNRecHits());
				ijet = jets[j].GetUserIdx();
				eECAL = _base->Jet_neEmEF->at(ijet) + _base->Jet_chEmEF->at(ijet);
				eECAL *= jets[j].E();
				TrueJet_EmE->Fill(eECAL);
				TrueJet_nConstituents->Fill(_base->Jet_nConstituents->at(ijet));			
			
				
				rhs.clear();
				rhs = jets[j].GetJetPoints();
				for(int r = 0; r < rhs.size(); r++){
					erhs_trhs->Fill(rhs[r].E(), rhs[r].t());
				}		
				if(njets < 2) continue;
				if(j == 0 || j == 1){	
					TrueJet_twoHardestpT->Fill(jets[j].pt());
					TrueJet_twoHardestpT->Fill(jets[j].pt());
				}
			}
		}

		//need to break down by procCat - see how this is done in PhotonSkimmer.cc
		void FillPVTimeHists(vector<Jet>& jets, int tr_idx, const Matrix& smear = Matrix(), double emAlpha = 0.5, double alpha = 0.1, double tres_c = 0.2, double tres_n = 0.3){
			int njets = jets.size();
			double jettime = -999;
			vector<double> jettimes;
			double gamtime = -999;
			double pvtime = -999;
			double deltaT_gampv = -999;
			double deltaT_gampv_gen = -999;
			double ptavg = 0;
			double geoEavg = 0;
			double Erh = 0;
			double Epho = 0;
			vector<JetPoint> rhs;
			TimeStrategy ts = TimeStrategy(tr_idx);
			//break down by process for this tr method
			int nProc = trCats[tr_idx].procCats.size();
			//cout << "method: " << tr_idx << " " << ts << " "<< trCats[tr_idx].methodName << endl;
			for(int p = 0; p < nProc; p++){
				//cout << " proc " << trCats[tr_idx].procCats[p].plotName << endl;
				for(int j = 0; j < njets; j++){
					jettime = CalcJetTime(ts, jets[j], smear, emAlpha, alpha, tres_c, tres_n);
					jets[j].SetJetTime(jettime);
					//cout << " jet #" << j << " time: " << jettime << endl;
					//fill jet time in pv frame - 0
					trCats[tr_idx].procCats[p].hists1D[0][0]->Fill(jettime);
					//add to pv time calculation
					jettimes.push_back(jettime);
					rhs = jets[j].GetJetPoints();
					for(int r = 0; r < rhs.size(); r++)
						Erh += rhs[r].E();
					rhs.clear();
				}
				//calculate jet time difference
				if(njets > 1){
					pair<Jet,Jet> hardjets;
					if(jets[0].pt() < jets[1].pt()){
						FindHardestJetPair(jets, hardjets);
					}
					//pt sorted
					else{
						hardjets = std::make_pair(jets[0],jets[1]);
					}
					double deltaT_jets = GetDeltaTime(hardjets);
					if(deltaT_jets != -999){
						ptavg = pow((hardjets.first.pt()*hardjets.second.pt()),0.5);	
						//fill deltaT_jets - 1
						trCats[tr_idx].procCats[p].hists1D[0][1]->Fill(deltaT_jets);
						//fill geopTavg vs deltaT jets - 1
						trCats[tr_idx].procCats[p].hists2D[0][1]->Fill(ptavg, deltaT_jets);
						//fill minpt vs deltaT jets - 2
						trCats[tr_idx].procCats[p].hists2D[0][2]->Fill(hardjets.second.pt(), deltaT_jets);
						//fill Erh vs deltaT jets - 3
						trCats[tr_idx].procCats[p].hists2D[0][3]->Fill(Erh, deltaT_jets);
					}
				}	
				//fill deltaT_pvGam - 2
				//only fill for two leading photons + weighted avg of jet time
				if(_phos.size() < 1) continue;
				//fill correct procCat
				vector<double> ids = trCats[tr_idx].procCats[p].ids;
				int phoidx = _phos[0].GetUserIdx();
				int phoid = _base->Photon_genLlpId->at(phoidx);
				vector<JetPoint> phorhs; 
				//make sure id is in current vector of ids (or ids does not contain -999)
				if(std::find(ids.begin(), ids.end(), phoid) != ids.end() || std::find(ids.begin(), ids.end(), -999) != ids.end()){
					//this assumes that the time for the jet was set previously with the respective method
					pvtime = CalcPVTime(ts, jets);
					gamtime = CalcJetTime(ts, _phos[0], smear, emAlpha, alpha, tres_c, tres_n);
					trCats[tr_idx].procCats[p].hists1D[0][5]->Fill(gamtime);
					deltaT_gampv = gamtime - pvtime;
					//get sum of pho rh energy
					phorhs = _phos[0].GetJetPoints();
					for(auto r : phorhs) Epho += r.E();
	
					trCats[tr_idx].procCats[p].hists1D[0][2]->Fill( deltaT_gampv );
					if(_data) continue; //no gen info with data
	
					//fill difference in deltaT_pvGam of reco and gen - 3
					deltaT_gampv_gen = CalcGenDeltaT(_phos[0]);
					trCats[tr_idx].procCats[p].hists1D[0][4]->Fill(deltaT_gampv_gen);
					//only for gen matches
					if(deltaT_gampv_gen != -999){
						trCats[tr_idx].procCats[p].hists1D[0][3]->Fill( deltaT_gampv - deltaT_gampv_gen);	
						//fill res (sigma from gaussian fit) for deltaT_recoGen as a function of ptAvg of jets that go into pv time calc - 4
						//sigma deltaT_recoGen as a function of geoEavg
						trCats[tr_idx].procCats[p].hists2D[0][0]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen);
					
					}	
	
				
				}

				//do same for subleading photon if it exists
				if(_phos.size() > 1){
					phoidx = _phos[1].GetUserIdx();
					phoid = _base->Photon_genLlpId->at(phoidx);
					//make sure id is in current vector of ids (or ids does not contain -999)
					if(std::find(ids.begin(), ids.end(), phoid) != ids.end() || std::find(ids.begin(), ids.end(), -999) != ids.end()){
						phorhs.clear();
						phorhs = _phos[1].GetJetPoints();
						for(auto r : phorhs) Epho += r.E();
						gamtime = CalcJetTime(ts, _phos[1], smear, emAlpha, alpha, tres_c, tres_n);
						trCats[tr_idx].procCats[p].hists1D[0][5]->Fill(gamtime);
						deltaT_gampv = gamtime - pvtime;
						trCats[tr_idx].procCats[p].hists1D[0][2]->Fill( deltaT_gampv );
						deltaT_gampv_gen = CalcGenDeltaT(_phos[1]);
						trCats[tr_idx].procCats[p].hists1D[0][4]->Fill(deltaT_gampv_gen);
						//only for gen matches
						if(deltaT_gampv_gen != -999){
							trCats[tr_idx].procCats[p].hists1D[0][3]->Fill( deltaT_gampv - deltaT_gampv_gen);
							trCats[tr_idx].procCats[p].hists2D[0][0]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen);
						}
					}
				}
			}	

		}

		double CalcAvgPt(const vector<Jet>& jets){
			double pt = 0;
			int njets = jets.size();
			for(int j = 0; j < njets; j++)
				pt += jets[j].pt();
			return pt/double(njets);
		}

		//should be deltaT = gam - pv
		double CalcGenDeltaT(const Jet& pho){
			//calc times differently for !sig and sig photons
			double dpho = -999;
			double c = 29.9792458; // speed of light in cm/ns
			int genidx, phoidx;
			phoidx = pho.GetUserIdx();
			//gen photon coordinates
			double px, py, pz, ptheta, peta, pphi, vx, vy, vz;
			//if no match
			if(_base->Photon_genLlpId->at(phoidx) == -1) return dpho;
			genidx = _base->Photon_genIdx->at(phoidx);
			peta = _base->Gen_eta->at(genidx);
			pphi = _base->Gen_phi->at(genidx);
			ptheta = 2*atan(exp(-1*peta));
			px = 120*sin(pphi);
			py = 120*cos(pphi);
			pz = 120/tan(ptheta);			
			
			double pvx = _base->PV_x;
			double pvy = _base->PV_y;
			double pvz = _base->PV_z;

			double beta;
			if(_base->Photon_genLlpId->at(phoidx) == 22){
				vx = _base->Photon_genSigMomVx->at(phoidx);
				vy = _base->Photon_genSigMomVy->at(phoidx);
				vz = _base->Photon_genSigMomVz->at(phoidx);
		
				double mompx, mompy, mompz, momE;			

				mompx = _base->Photon_genSigMomPx->at(phoidx);			
				mompy = _base->Photon_genSigMomPy->at(phoidx);			
				mompz = _base->Photon_genSigMomPz->at(phoidx);			
				momE = _base->Photon_genSigMomEnergy->at(phoidx);			

				beta = sqrt(mompx*mompx + mompy*mompy + mompz*mompz)/momE;
			
				//distance bw photon and production point (LLP)
				dpho = sqrt( (px - vx)*(px - vx) + (py - vy)*(py - vy) + (pz - vz)*(pz - vz) )/c;
		
				//distance bw LLP and PV
				dpho += sqrt( (vx - pvx)*(vx - pvx) + (vy - pvy)*(vy - pvy) + (vz - pvz)*(vz - pvz) )/(c*beta);	

			}
			//assume prompt production
			else{
				//distance bw photon and production point (PV)
				dpho = sqrt( (px - pvx)*(px - pvx) + (py - pvy)*(py - pvy) + (pz - pvz)*(pz - pvz) )/c;
			}
			
			return dpho;
		}



		void FillPredJetHists(const vector<node*>& trees){
			int njets = 0;
			vector<node*> cleaned_trees;
			for(int i = 0; i < trees.size(); i++){
				if(trees[i] == nullptr) continue;
				//check for mirrored point - would be double counted
				if(trees[i]->points->mean().at(1) > 2*acos(-1) || trees[i]->points->mean().at(1) < 0) continue;	
				FillModelHists(trees[i]->model);
				njets++;
				cleaned_trees.push_back(trees[i]);
			}
			nClusters->Fill(njets);
			FillPVHists_PredJets(cleaned_trees);
		}
	

		//this is for one jet
		//all hists referenced here are in hists1D
		void FillModelHists(BasePDFMixture* model){
			map<string, Matrix> params;
			vector<double> eigenvals, norms;
			vector<Matrix> eigenvecs;
			double theta, phi, r, id, npts, E_k;

			int nclusters = model->GetNClusters();
			
			nSubClusters->Fill(nclusters);
			model->GetNorms(norms);
		
			nClusters->Fill((double)nclusters);
			//k clusters = k jets in event -> subclusters are mixture model components
			for(int k = 0; k < nclusters; k++){
				E_k = norms[k]/_gev;

				params = model->GetPriorParameters(k);
				eta_center->Fill(params["mean"].at(0,0));
				phi_center->Fill(params["mean"].at(1,0));
				time_center->Fill(params["mean"].at(2,0));
		
				//calculate slopes from eigenvectors
				params["cov"].eigenCalc(eigenvals, eigenvecs);
				
				//total cluster energy
				clusterE->Fill(E_k);
			}
		}
		


		//find back to back jets
		void FillPVHists_PredJets(const vector<node*>& trees){
			int njets = (int)trees.size(); 
			double pi = acos(-1);

			double dtime, dphi, dr, phi1, t1, phi2, t2;
			//find pairs of jets to calculate resolution	
			//need to be back to back
			//time of subclusters is measured as center
			for(int i = 0; i < njets; i++){
				if(trees[i] == nullptr) continue;
				CalcMMAvgPhiTime(trees[i]->model, phi1, t1);
				//PVdeltaT_jet_mmAvg_pred->Fill(t1);	
				for(int j = i+1; j < njets; j++){
					if(trees[i] == nullptr) continue;
					CalcMMAvgPhiTime(trees[j]->model, phi2, t2);
					//dphi within [pi-0.1,pi+0.1]
					dphi = fabs(phi1 - phi2);
					if(dphi < pi-0.1 || dphi > pi+0.1) continue;
					//median time
					//energy-weighted average
					//mm average over subclusters
					//PVdeltaT_jet_mmAvg_pred->Fill(t1 - t2);	
				}
			}

		}
		void FindHardestJetPair(const vector<Jet>& injets, pair<Jet,Jet>& outjets){
			int njets = (int)injets.size(); 
			map<double,Jet> pt_jet;
			map<double, int> pt_idx;
			for(int i = 0; i < njets; i++){
				pt_jet[injets[i].pt()] = injets[i];
				pt_idx[injets[i].pt()] = i;
			}
			map<double,Jet>::reverse_iterator it = pt_jet.rbegin();	
			map<double,int>::reverse_iterator it_idx = pt_idx.rbegin();	
			outjets.first = it->second;
			outjets.first.SetUserIdx(it_idx->second);
			it++;
			it_idx++;
			outjets.second = it->second;
			outjets.second.SetUserIdx(it_idx->second);
			//if needing to find "true" pair
			//calculate dphi here
			//if dphi doesn't satisfy back-to-back pi requirement (see below) check next two jets (ie it->second and it++; it->second;)
			//continue until pair is found

		}
		double GetDeltaTime(pair<Jet,Jet>& jets){
			double pi = acos(-1);
			//not needed for GMSB
			if(_data){
				//dphi within [pi-0.1,pi+0.1]
				double phi1 = jets.first.phi_02pi();
				double phi2 = jets.second.phi_02pi();
				double dphi = fabs(phi1 - phi2);
				if(dphi > pi) dphi = 2*pi - dphi;
				if(dphi < pi-0.35 || dphi > pi+0.35) return -999;
			}
			double t1 = jets.first.time();
			double t2 = jets.second.time();
			return t1 - t2;
		}

	
		double CalcJetTime(const TimeStrategy& ts, Jet& jet, const Matrix& smear = Matrix(), double emAlpha = 0.5, double alpha = 0.1, double tres_c = 0.2, double tres_n = 0.3){
			//cout << "CalcJetTime method " << ts << endl;
			if(ts == med) return CalcMedianTime(jet);
			else if(ts == eavg) return CalcEAvgTime(jet);
			else if(ts == mmavg){
				GaussianMixture* gmm = _subcluster(jet, smear, emAlpha, alpha, tres_c, tres_n);
				nSubClusters_mm->Fill(gmm->GetNClusters());
				//cout << " nSubclusters: " << gmm->GetNClusters() << endl;
				return CalcMMAvgTime(gmm);
			}
			else return -999;

		}
		double CalcPVTime(const TimeStrategy& ts, vector<Jet>& jets){
			int njets = jets.size();
			double time = -999;
			if(ts == med){
				vector<double> times;
				for(int i = 0; i < njets; i++)
					times.push_back(jets[i].t());
				//even - return average of two median times
				if(njets % 2 == 0)
					time = (times[int(double(njets)/2.)] + times[int(double(njets)/2.)-1]/2.);
				//odd - return median
				else
					time = times[int(double(njets)/2.)];

			}
			else if(ts == eavg){
				double norm = 0;
				double t = 0;
				for(int j = 0; j < njets; j++){
					t += jets[j].E()*jets[j].t();
					norm += jets[j].E();
				}
				time = t/norm;
			}	
			else if(ts == mmavg){
				double norm = 0;
				double t = 0;
				for(int j = 0; j < njets; j++){
					t += jets[j].t();
				}
				time = t/double(njets);
			}	
			return time;
		}

		double CalcMedianTime(Jet& j){
			vector<double> times;
			double time = -999;
			vector<JetPoint> rhs = j.GetJetPoints();
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++)
				times.push_back(rhs[i].t());
			//even - return average of two median times
			if(nrhs % 2 == 0)
				time = (times[int(double(nrhs)/2.)] + times[int(double(nrhs)/2.)-1]/2.);
			//odd - return median
			else
				time = times[int(double(nrhs)/2.)];
			return time;
		}		


		double CalcEAvgTime(Jet& j){
			double t = 0;
			double norm = 0;
			double time = -999;
			vector<JetPoint> rhs = j.GetJetPoints();
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++){
				norm += rhs[i].E();
				t += rhs[i].E()*rhs[i].t();
			}
			time = t/norm;
			return time;
		}


		GaussianMixture* _subcluster(const Jet& jet, const Matrix& smear = Matrix(), double emAlpha = 0.5, double alpha = 0.1, double tres_c = 0.2, double tres_n = 0.3){
			vector<JetPoint> rhs = jet.GetJetPoints();
			vector<Jet> rhs_jet;
			for(int r = 0; r < rhs.size(); r++){
				rhs_jet.push_back( Jet(rhs[r]) );
				rhs_jet[r].SetVertex(jet.GetVertex());
			}
			BayesCluster* algo = new BayesCluster(rhs_jet);
			if(_smear && !smear.empty()) algo->SetDataSmear(smear);
			//set time resolution smearing
			if(_timesmear) algo->SetTimeResSmear(tres_c, tres_n*_gev);
			algo->SetThresh(1.);
			algo->SetAlpha(alpha);
			algo->SetSubclusterAlpha(emAlpha);
			algo->SetVerbosity(0);
			GaussianMixture* gmm = algo->SubCluster();
			return gmm;
		}



		double CalcMMAvgTime(BasePDFMixture* model){
			int kmax = model->GetNClusters();
			double t = 0;
			double norm = 0;
			map<string, Matrix> params;
			for(int k = 0; k < kmax; k++){
				params = model->GetPriorParameters(k);
				t += params["pi"].at(0,0)*params["mean"].at(2,0);
				norm += params["pi"].at(0,0);
				//cout << "cluster " << k << " has time " << params["mean"].at(2,0) << " and MM coeff " << params["pi"].at(0,0) << endl;
			}
			return t/norm; 
		}

		void CalcMMAvgPhiTime(BasePDFMixture* model, double& phi, double& t){
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double pi, ws;
			map<string, Matrix> params;
			for(int k = 0; k < kmax; k++){
				params = model->GetPriorParameters(k);
				phi += params["pi"].at(0,0)*params["mean"].at(1,0);
				t += params["pi"].at(0,0)*params["mean"].at(2,0);
				ws += params["pi"].at(0,0);
			}
			phi /= ws;
			t /= ws; 
		}
		
		void CalcEAvgPhiTime(BasePDFMixture* model, double& phi, double& t){
			int kmax = model->GetNClusters();
			phi = 0;
			t = 0;
			double pi, ws;
			map<string, Matrix> params;
			for(int i = 0; i < model->GetData()->GetNPoints(); i++){
				phi += model->GetData()->at(i).w()/_gev*model->GetData()->at(i).Value(1);
				t += model->GetData()->at(i).w()/_gev*model->GetData()->at(i).Value(2);
				ws += model->GetData()->at(i).w()/_gev;
			}
			phi /= ws;
			t /= ws; 
		}

		void WriteEmptyProfiles(TFile* ofile, timeRecoCat& tr){
			ofile->cd();
			string name, addname;
			int i, j, nbins, nprofs;
			string profname, histname, match;
			match = "_"+tr.methodName;
			//AddHist in timeRecoCat adds hist for each process - don't need to loop over it here
			for(int i = 0; i < tr.procCats[0].hists2D[0].size(); i++){
				histname = tr.procCats[0].hists2D[0][i]->GetName();
				//make sure taking info from 2D diff histogram
				if(histname.find("diffDeltaTime") == string::npos) continue;
				//check if data can be run
				if(histname.find("recoGen") != string::npos && _data) continue;
				//make sure profiles get written
				nprofs = tr.procCats[0].hists2D[0][i]->GetNbinsX();
				//histname = histname.substr(0,histname.find(match));		
				for(int k = 1; k < nprofs; k++){
					histname = tr.procCats[0].hists2D[0][i]->GetName();
					profname = "profile_"+histname;
					profname.insert(profname.find("_"+tr.methodName),"_bin"+std::to_string(k));
					nbins = tr.procCats[0].hists2D[0][i]->GetNbinsY();
					TH1D* prof = new TH1D(profname.c_str(), profname.c_str(), nbins, tr.procCats[0].hists2D[0][i]->GetYaxis()->GetBinLowEdge(1), tr.procCats[0].hists2D[0][i]->GetYaxis()->GetBinUpEdge(nbins));
					prof->GetXaxis()->SetTitle(tr.procCats[0].hists2D[0][i]->GetXaxis()->GetTitle());	
					tr.AddHist(prof);	
				}	
			}				
		}
		void WriteTimeRecoCatStack(TFile* ofile, const vector<timeRecoCat>& trs){
			ofile->cd();
			string name, dirname, histname;
			//write 1D hists
			//variables
			int nhists = trs[0].procCats[0].hists1D[0].size();
			for(int i = 0; i < nhists; i++){
				name = trs[0].procCats[0].hists1D[0][i]->GetName();
				dirname = name.substr(0,name.rfind("_"+trs[0].methodName));
			//cout << "i: " << i << " name " << name << " making dir " << dirname+"_stack" << endl;
				TDirectory* dir = ofile->mkdir((dirname+"_stack").c_str());
				for(int j = 0; j < trs.size(); j++){
					dir->cd();
					histname = trs[j].procCats[0].hists1D[0][i]->GetName();
					//write total method histogram outside process directory
					if(trs[j].procCats[0].hists1D[0][i] == nullptr) continue;
					if(trs[j].procCats[0].hists1D[0][i]->GetEntries() == 0 && (histname.find("sigma") == string::npos && histname.find("profile") == string::npos)){ continue; }
					//cout << "writing " << trs[j].procCats[0].hists1D[0][i]->GetName() << " " << trs[j].procCats[0].hists1D[0][i]->GetTitle() << " to " << dir->GetName() << endl;
					trs[j].procCats[0].hists1D[0][i]->Write();
					//make process breakdown directory
					TDirectory *dir2 = dir->mkdir((dirname+"_"+trs[j].methodName+"_procs").c_str());
					//cout << "  making dir " << dir2->GetName() << endl;
					dir2->cd();
					for(int p = 1; p < trs[j].procCats.size(); p++){
					//loop over processes
			//			cout << "    proc " << trs[j].procCats[p].plotName << " hist " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " entries " << trs[j].procCats[p].hists1D[0][i]->GetEntries() << endl;			
						histname = trs[j].procCats[p].hists1D[0][i]->GetName();
						if(trs[j].procCats[p].hists1D[0][i] == nullptr) continue;
						if(trs[j].procCats[p].hists1D[0][i]->GetEntries() == 0 && (histname.find("sigma") == string::npos && histname.find("profile") == string::npos)){ continue; }
						//cout << "  n hists " << trs[j].procCats[0].hists1D[0].size() << endl;
						//cout << "writing " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " to " << dir2->GetName() << endl;;
						trs[j].procCats[p].hists1D[0][i]->Write();
					

					} 
				}
			}
			//write 2D hists
			nhists = trs[0].procCats[0].hists2D[0].size();
			for(int i = 0; i < nhists; i++){
				name = trs[0].procCats[0].hists2D[0][i]->GetName();
				//only write sigma plots
				//if(name.find("sigma") == string::npos) continue;
				//first tr cat is always median - change to rfind("_")
				dirname = name.substr(0,name.rfind("_"));
			//cout << "i: " << i << " name " << name << " making dir " << dirname+"_stack" << endl;
				TDirectory* dir = ofile->mkdir((dirname+"_stack").c_str());
				dir->cd();
				for(int j = 0; j < trs.size(); j++){
					dir->cd();
					histname = trs[j].procCats[0].hists2D[0][i]->GetName();
					//write total method histogram outside process directory
					if(trs[j].procCats[0].hists2D[0][i] == nullptr) continue;
					if(trs[j].procCats[0].hists2D[0][i]->GetEntries() == 0 && (histname.find("sigma") == string::npos || histname.find("profile") == string::npos)){ continue; }
					//cout << "writing " << trs[j].procCats[0].hists1D[0][i]->GetName() << " " << trs[j].procCats[0].hists1D[0][i]->GetTitle() << " to " << dir->GetName() << endl;
					trs[j].procCats[0].hists2D[0][i]->Write();
					//write method as directory within directory
					TDirectory *dir2 = dir->mkdir((dirname+"_"+trs[j].methodName+"_procs").c_str());
					//cout << "  making dir " << dir2->GetName() << endl;
					dir2->cd();
					for(int p = 1; p < trs[j].procCats.size(); p++){
					//loop over processes
			//			cout << "    proc " << trs[j].procCats[p].plotName << " hist " << trs[j].procCats[p].hists2D[0][i]->GetName() << " " << trs[j].procCats[p].hists2D[0][i]->GetTitle() << " entries " << trs[j].procCats[p].hists2D[0][i]->GetEntries() << endl;			
						if(trs[j].procCats[p].hists2D[0][i] == nullptr) continue;
						if(trs[j].procCats[p].hists2D[0][i]->GetEntries() == 0){ continue; }//cout << "Histogram for proc " << trs[j].plotName << " not filled." << endl; continue; }
						//cout << "  n hists " << trs[j].procCats[0].hists1D[0].size() << endl;
						trs[j].procCats[p].hists2D[0][i]->Write();
					

					} 
				}
			}
			
		}


		void WriteHists(TFile* ofile){
			string name;

			ofile->cd();
			//for condor skims, histograms need to be written because these
			//can be hadded together (TCanvases can't)
			for(int i = 0; i < (int)_hists1D.size(); i++){
				//name = hists1D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRHist(hists1D[i], cv, name, name, "a.u.");
				//write cv to file			
			//	cv->SaveAs((fname+"/"+name+".pdf").c_str());
				//cv->Write();
				_hists1D[i]->Write();
			}
			for(int i = 0; i < (int)_hists2D.size(); i++){
				//name = hists2D[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDR2DHist(hists2D[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
			_hists2D[i]->Write();
			}
			//erhs_trhs->Write();		
			for(int i = 0; i < (int)graphs.size(); i++){
				//name = graphs[i]->GetName();
				//TCanvas* cv = new TCanvas(name.c_str(), "");
				//TDRGraph(graphs[i], cv, name, name, "a.u.");
				//write cv to file			
				//cv->Write();
				graphs[i]->Write();
			}
			for(int i = 0; i < trCats.size(); i++)
				WriteEmptyProfiles(ofile, trCats[i]);
			WriteTimeRecoCatStack(ofile, trCats);

			ofile->Close();


		}


		void SetStrategy(int i){
			if(i == 0) _strategy = NlnN;
			else if(i == 1) _strategy = N2;
			else if(i == 2) _strategy = MM;
			else return; 
		}


		private:
			enum Strategy{
				//Delauney strategy - NlnN time - for 2pi cylinder
				NlnN = 0,
				//traditional strategy - N^2 time
				N2 = 1,
				//mm only
				MM = 2
			};
		
		
			//clustering strategy - N^2 or NlnN
			Strategy _strategy;
			vector<TGraph*> graphs;
			vector<Jet> _phos; //photons for event

};
#endif
