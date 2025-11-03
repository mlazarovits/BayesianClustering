#ifndef JETSKIMMER_HH
#define JETSKIMMER_HH

#include <utility>
#include "JetPoint.hh"
#include "BaseSkimmer.hh"
#include "BasePDFMixture.hh"
#include "BayesCluster.hh"
#include "JetProducer.hh"
#include "PhotonProducer.hh"
#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include "TSystem.h"
#include "BaseTree.hh"

using procCat = BaseSkimmer::procCat;

class JetSkimmer : public BaseSkimmer{
	public:
		JetSkimmer(){
			_evti = 0;
			_evtj = 0;
			_gev = 1./10.;
			_swts.Init();
			_weight = 1.;
			_minRhE = 5;
			
			_cleansubcls = false;
			_smearMat = Matrix();
		
		};
		virtual ~JetSkimmer(){ };

		//get rechits from file to cluster
		JetSkimmer(TFile* file) : BaseSkimmer(file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time		
			SetupDetIDsEB( _detIDmap, _ietaiphiID );
			
			_prod = new JetProducer(file);
			_prod->SetIsoCut();
	
			//set producer to get jets with different kin reqs - can't use same file pointer ig?
			string fname = file->GetName();
			cout << "fname " << fname << endl;
			TFile* f2 = TFile::Open(fname.c_str());
			_scprod = new PhotonProducer(f2);	
			_scprod->SetMinPt(5);
			_scprod->SetMinRhE(0.5);

			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/jet_skims_"+_cms_label+".root";
			_gev = 1./10.;
			_swts.Init();
			_weight = 1;
			_minRhE = 5;
			_cleansubcls = false;
			_smearMat = Matrix();
			//set histogram weights
			//if(_data || fname.find("QCD") != string::npos){ _weight = 1.; }
			//if(_data){ _weight = 1.; }
			//if(_data || fname.find("QCD") != string::npos){ _weight = 1.; }
			//else{
			//	ifstream weights("info/EventWeights.txt", std::ios::in);
			//	string filein;
			//	string filename = file->GetName();
			//	double jet_weight, pho_weight;
			//	while( weights >> filein >> jet_weight >> pho_weight){
			//		if(filename.find(filein) == string::npos) continue;
			//		else{
			//			_weight = jet_weight;
			//			break;
			//		}
			//	}
			//} 


			
			InitHists();
		}
		JetSkimmer(string filelist) : BaseSkimmer(filelist){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time		
			SetupDetIDsEB( _detIDmap, _ietaiphiID );
		
			TChain* ch = MakeTChain(filelist);
                        if(ch == nullptr) return;	
			_prod = new JetProducer(ch);
			_prod->SetIsoCut();
	
			//set producer to get jets with different kin reqs - can't use same file pointer ig?
                        TChain* ch2 = MakeTChain(filelist);
			_scprod = new PhotonProducer(ch2);	
			_scprod->SetMinPt(5);
			_scprod->SetMinRhE(0.5);

			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/jet_skims_"+_cms_label+".root";
			_gev = 1./10.;
			_swts.Init();
			_weight = 1;
			_minRhE = 5;
			//set histogram weights
			//if(_data || fname.find("QCD") != string::npos){ _weight = 1.; }
			//if(_data){ _weight = 1.; }
			//if(_data || fname.find("QCD") != string::npos){ _weight = 1.; }
			//else{
			//	ifstream weights("info/EventWeights.txt", std::ios::in);
			//	string filein;
			//	string filename = file->GetName();
			//	double jet_weight, pho_weight;
			//	while( weights >> filein >> jet_weight >> pho_weight){
			//		if(filelist.find(filein) == string::npos) continue;
			//		else{
			//			_weight = jet_weight;
			//			break;
			//		}
			//	}
			//} 


			
			InitHists();
		}

		void SetObs(){ 
			_obj = "jet";
		};
		void WriteHeader(){ };
		void FillBranches(const Jet& bhc_obj){

		}

		void InitHists(){
			//true jet hists
			_hists1D.push_back(nTrueJets);
			_hists1D.push_back(TrueJet_pT); 
			_hists1D.push_back(TrueJet_nRhs); 
			_hists1D.push_back(TrueJet_EmE); 
			_hists1D.push_back(TrueJet_nConstituents);
			_hists1D.push_back(TrueJet_twoHardestpT);	
			_hists1D.push_back(TOFgam_rh_pv); 
			_hists1D.push_back(subclusterEfrac);
			_hists1D.push_back(subclDist_etaPhi);
			_hists1D.push_back(subclDist_time);
			_hists1D.push_back(etaSig);
			_hists1D.push_back(phiSig);
			_hists1D.push_back(timeSig);
			_hists1D.push_back(beta_k);
			_hists1D.push_back(W_ee_k);	
			_hists1D.push_back(W_pp_k);	
			_hists1D.push_back(W_tt_k);	
			_hists1D.push_back(t_physBkg);
			_hists1D.push_back(t_BH);
			_hists1D.push_back(t_spike);

		
			
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
			_timeHists1D.push_back(geoAvgEecal_sigmaDeltaTime_dijets);
			_timeHists1D.push_back(geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin1);
			_timeHists1D.push_back(geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin2);
			_timeHists1D.push_back(geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin3);
			_timeHists1D.push_back(geoEavg_meanDeltaTime_recoGen);
			_timeHists1D.push_back(geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin1);
			_timeHists1D.push_back(geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin2);
			_timeHists1D.push_back(geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin3);
			_timeHists1D.push_back(deltaT_pvGam_gen_Ebin1);	
			_timeHists1D.push_back(deltaT_pvGam_gen_Ebin2);	
			_timeHists1D.push_back(deltaT_pvGam_gen_Ebin3);	
			_timeHists1D.push_back(deltaT_pvGam_gen_Ebin4);	
			_timeHists1D.push_back(geoEavg_sigmaDeltaTime_gamPV);
			_timeHists1D.push_back(jetNrhs_Ebin1);
			_timeHists1D.push_back(jetNrhs_Ebin2);
			_timeHists1D.push_back(jetNrhs_Ebin3);
			_timeHists1D.push_back(jetNrhs_Ebin4);
			_timeHists1D.push_back(jetNrhs_Ebin5);
			_timeHists1D.push_back(jetNrhs_Ebin6);
			_timeHists1D.push_back(jetPt_Ebin1);
			_timeHists1D.push_back(jetPt_Ebin2);
			_timeHists1D.push_back(jetPt_Ebin3);
			_timeHists1D.push_back(jetPt_Ebin4);
			_timeHists1D.push_back(jetPt_Ebin5);
			_timeHists1D.push_back(jetPt_Ebin6);
			_timeHists1D.push_back(jetEta_Ebin1);
			_timeHists1D.push_back(jetEta_Ebin2);
			_timeHists1D.push_back(jetEta_Ebin3);
			_timeHists1D.push_back(jetEta_Ebin4);
			_timeHists1D.push_back(jetEta_Ebin5);
			_timeHists1D.push_back(jetEta_Ebin6);
			_timeHists1D.push_back(jetPhi_Ebin1);
			_timeHists1D.push_back(jetPhi_Ebin2);
			_timeHists1D.push_back(jetPhi_Ebin3);
			_timeHists1D.push_back(jetPhi_Ebin4);
			_timeHists1D.push_back(jetPhi_Ebin5);
			_timeHists1D.push_back(jetPhi_Ebin6);
			_timeHists1D.push_back(jetTime_Ebin1);
			_timeHists1D.push_back(jetTime_Ebin2);
			_timeHists1D.push_back(jetTime_Ebin3);
			_timeHists1D.push_back(jetTime_Ebin4);
			_timeHists1D.push_back(jetTime_Ebin5);
			_timeHists1D.push_back(jetTime_Ebin6);
			_timeHists1D.push_back(jetNSubclusters_Ebin1);
			_timeHists1D.push_back(jetNSubclusters_Ebin2);
			_timeHists1D.push_back(jetNSubclusters_Ebin3);
			_timeHists1D.push_back(jetNSubclusters_Ebin4);
			_timeHists1D.push_back(jetNSubclusters_Ebin5);
			_timeHists1D.push_back(jetNSubclusters_Ebin6);
			_timeHists1D.push_back(swCross_rhTimeNeg15toNeg9);
			_timeHists1D.push_back(swCross_rhTimeNeg6toNeg3);
			_timeHists1D.push_back(swCross_rhTimeNeg12);
			_timeHists1D.push_back(swCross_rhTimeNeg5);
			_timeHists1D.push_back(swCross_rhTime0);
			_timeHists1D.push_back(dRtrack_rhTimeNeg12);
			_timeHists1D.push_back(dRtrack_rhTimeNeg5);
			_timeHists1D.push_back(dRtrack_rhTime0);
			_timeHists1D.push_back(rhTime); 
			_timeHists1D.push_back(seedXtalEnergy_sigmaDeltaTime);		
			_timeHists1D.push_back(gamTime_maxErh); 
			_timeHists1D.push_back(seedXtalTime); 
			_timeHists1D.push_back(nSubclusters);
			_timeHists1D.push_back(nSubclusters_dijetSel);

			_timeHists2D.push_back(geoEavg_diffDeltaTime_recoGen);
			_timeHists2D.push_back(geopTavg_diffDeltaTime_dijets);	
			_timeHists2D.push_back(minpT_diffDeltaTime_dijets);	
			_timeHists2D.push_back(geoAvgEecal_diffDeltaTime_dijets);	
			_timeHists2D.push_back(geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin1);
			_timeHists2D.push_back(geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin2);
			_timeHists2D.push_back(geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin3);
			_timeHists2D.push_back(geoEavg_genDeltaTime_meanRecoGenDeltaT);
			_timeHists2D.push_back(genDeltaT_recoDeltaT_Ebin1);
			_timeHists2D.push_back(genDeltaT_recoDeltaT_Ebin2);
			_timeHists2D.push_back(genDeltaT_recoDeltaT_Ebin3);
			_timeHists2D.push_back(genDeltaT_recoDeltaT_Ebin4);
			_timeHists2D.push_back(recoGenDr_recoGenDeltaTRatio);
			_timeHists2D.push_back(genEnergy_recoGenDeltaTRatio);
			_timeHists2D.push_back(recoEnergy_genDeltaTpvGam);
			_timeHists2D.push_back(recoGenEnergyRatio_recoGenDeltaT);
			_timeHists2D.push_back(geoEavg_diffDeltaTime_gamPV);
			_timeHists2D.push_back(jetTime_Energy);
			_timeHists2D.push_back(rhTime_Energy);
			_timeHists2D.push_back(rhTime_eta);
			_timeHists2D.push_back(rhPhi_eta);
			_timeHists2D.push_back(rhEovP_dRtrack);
			_timeHists2D.push_back(rhTime_Energy_zoomed);
			_timeHists2D.push_back(swCross_rhTime);
			_timeHists2D.push_back(kWeird_rhTime);
			_timeHists2D.push_back(swCross_rhEnergy);	
			_timeHists2D.push_back(kWeird_rhEnergy);
			_timeHists2D.push_back(e4e1_kWeird_rhEge50rhTleNeg2);
			_timeHists2D.push_back(rhEta_rhPhi_rhTimeNeg15toNeg9);
			_timeHists2D.push_back(rhEta_rhPhi_rhTimeNeg6toNeg3);
			_timeHists2D.push_back(rhTime_rhEta_rhTimeNeg15toNeg9);
			_timeHists2D.push_back(rhTime_rhEta_rhTimeNeg6toNeg3);
			_timeHists2D.push_back(normNeighbors);	
			_timeHists2D.push_back(swCross_rhEnergy_timeNeg12);
			_timeHists2D.push_back(swCross_rhEnergy_timeNeg5);	
			_timeHists2D.push_back(swCross_rhEnergy_time0);	
			_timeHists2D.push_back(ENeighbors);	
			_timeHists2D.push_back(ENeighbors_timeNeg12);	
			_timeHists2D.push_back(ENeighbors_timeNeg5);	
			_timeHists2D.push_back(ENeighbors_time0);	
			_timeHists2D.push_back(rhEta_rhPhi_rhTimeNeg12);
			_timeHists2D.push_back(rhEta_rhPhi_rhTimeNeg5);
			_timeHists2D.push_back(rhEta_rhPhi_rhTime0);
			_timeHists2D.push_back(ENeighborsSkipCenter_timeNeg12);	
			_timeHists2D.push_back(ENeighborsSkipCenter_timeNeg5);	
			_timeHists2D.push_back(ENeighborsSkipCenter_time0);	
			_timeHists2D.push_back(ENeighborsJet);	
			_timeHists2D.push_back(seedXtalEnergy_diffDeltaTime); 
			_timeHists2D.push_back(AK4Jet_nRhs_nSubclusters);
			_timeHists2D.push_back(AK4Jet_nSuperclusters_nSubclusters);
			_timeHists2D.push_back(Photon_nRhs_nSubclusters);
			_timeHists2D.push_back(subclTime_Energy);
			_timeHists2D.push_back(subclTime_Energy_spikes);
			_timeHists2D.push_back(subclTime_Energy_BH);
			_timeHists2D.push_back(subclTime_Energy_physBkg);
			_timeHists2D.push_back(subclNNDiscr_Energy_spikes);
			_timeHists2D.push_back(subclNNDiscr_Energy_BH);
			_timeHists2D.push_back(subclNNDiscr_Energy_physBkg);
			_timeHists2D.push_back(subclNNDiscr_time_spikes);
			_timeHists2D.push_back(subclNNDiscr_time_BH);
			_timeHists2D.push_back(subclNNDiscr_time_physBkg);
			_timeHists2D.push_back(AK4Jet_nSuperclusters_nSubclusters_dijetSel);
			_timeHists2D.push_back(rhTime_Energy_highE);


		};
		//ctor from rec hit collection - integrating into ntuplizer
		//needs to match order of timeRecoCats
		enum TimeStrategy{
			med = 0, 
			mmavg = 1, 
			eavg = 2,
			emax = 3
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
				else if(ts == emax)
					methodName = "eMax";
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
		TH1D* TrueJet_pT = new TH1D("TrueJet_pT","TrueJet_pT",100,0,1000);
		TH1D* TrueJet_nRhs = new TH1D("TrueJet_nRhs","TrueJet_nRhs",25,0,100);
		TH1D* TrueJet_EmE = new TH1D("TrueJet_EmE","TrueJet_EmE",50,0,600);
		TH1D* TrueJet_nConstituents = new TH1D("TrueJet_nConstituents","TrueJet_nConstituents",20,0,100);
		TH1D* TrueJet_twoHardestpT =  new TH1D("TrueJet_twoHardestpT","TrueJet_twoHardestpT",100,0,1000);	
		TH1D* TOFgam_rh_pv = new TH1D("TOFgam_rh_pv","TOFgam_rh_pv",20,0,10); 
		TH1D* subclusterEfrac = new TH1D("subclusterEfrac","subclusterEfrac",50,0,1);
		TH1D* subclDist_etaPhi = new TH1D("subclDist_etaPhi","subclDist_etaPhi",50,0,0.2);
		TH1D* subclDist_time = new TH1D("subclDist_time","subclDist_time",50,-2,2);
		TH1D* etaSig = new TH1D("etaSig","etaSig",50,0,0.5);
		TH1D* phiSig = new TH1D("phiSig","phiSig",50,0,0.5);
		TH1D* timeSig = new TH1D("timeSig","timeSig",50,0,1);
		TH1D* beta_k = new TH1D("beta_k","beta_k",50,0,20);
		TH1D* W_ee_k = new TH1D("W_ee_k","W_ee_k",50,0,0.5);	
		TH1D* W_pp_k = new TH1D("W_pp_k","W_pp_k",50,0,0.5);	
		TH1D* W_tt_k = new TH1D("W_tt_k","W_tt_k",50,0,0.5);	
		TH1D* t_physBkg = new TH1D("t_physBkg","t_physBkg",50,-20,20);
		TH1D* t_BH = new TH1D("t_BH","t_BH",50,-20,20);
		TH1D* t_spike = new TH1D("t_spike","t_spike",50,-20,20);

		TH2D* e_nRhs = new TH2D("e_nRhs","e_nRhs",100,0,500,100,0,100);
		TH2D* erhs_trhs = new TH2D("erhs_trhs","erhs_trhs",100,0,4,100,-100,100);
		
		
		//bins for variable binning for resolution/mean plots
		//vector<double> xbins = {0, 200, 400, 600, 800, 1000, 1200}; 
		vector<double> xbins = {0, 200, 400, 600, 800, 1200}; 
		///////////////////// timeHists /////////////
		//0 - pv time
		TH1D* PVtime = new TH1D("jetTime_PV", "jetTime_PV",50,-10,10);	
		//1 - delta t between jets (pv time frame)
		TH1D* deltaT_jet = new TH1D("deltaT_jet", "deltaT_jet",50,-4,4);	
		//2 - reco delta t between pv and photon 
		TH1D* deltaT_pvGam = new TH1D("deltaT_gamPV_reco","deltaT_gamPV_reco",100,-1,1);
		//3 - difference in deltaT_pvGam between gen and reco
		TH1D* diffDeltaT_recoGen = new TH1D("diffDeltaT_recoGen","diffDeltaT_recoGen",100,-1,1);
		//4 - gen deltaT bw photon and pv
		TH1D* deltaT_pvGam_gen = new TH1D("deltaT_gamPV_gen","deltaT_gamPV_gen",25,0,12);	
		//5 - photon time (in pv frame for data, not for MC)
		TH1D* gamTime = new TH1D("gamTime_reco", "gamTime_reco",100,-1,1);	

		//these stay empty to be filled later (after hadding)	
		//6 - resolution of difference in reco - gen deltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_sigmaDeltaTime_recoGen = new TH1D("geoEavg_sigmaDeltaTime_recoGen","geoEavg_sigmaDeltaTime_recoGen",xbins.size()-1,&xbins[0]);
		//7 - resolution of difference between two jets for PV time as a function of their geo avg pT 	
		TH1D* geopTavg_sigmaDeltaTime_dijets = new TH1D("geopTavg_sigmaDeltaTime_dijets","geopTavg_sigmaDeltaTime_dijets",10,0,1000);
		//8 - resolution of difference between two jets for PV time as a function of the min pT 	
		TH1D* minpT_sigmaDeltaTime_dijets = new TH1D("minpT_sigmaDeltaTime_dijets","minpT_sigmaDeltaTime_dijets",10,0,1000);
		//9 - resolution of difference between two jets for PV time as a function of the sum ECAL energy	
		TH1D* geoAvgEecal_sigmaDeltaTime_dijets = new TH1D("geoAvgEecal_sigmaDeltaTime_dijets","geoAvgEecal_sigmaDeltaTime_dijets",xbins.size()-1,&xbins[0]);
		//10 - resolution of difference in reco - gen deltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin1 = new TH1D("geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin1","geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin1",xbins.size()-1,&xbins[0]);
		//11 - resolution of difference in reco - gen genDeltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin2 = new TH1D("geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin2","geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin2",xbins.size()-1,&xbins[0]);
		//12 - resolution of difference in reco - gen genDeltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin3 = new TH1D("geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin3","geoEavg_sigmaDeltaTime_recoGen_genDeltaTpvGambin3",xbins.size()-1,&xbins[0]);
		//13 - mean of difference in reco - gen genDeltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_meanDeltaTime_recoGen = new TH1D("geoEavg_meanDeltaTime_recoGen","geoEavg_meanDeltaTime_recoGen",xbins.size()-1,&xbins[0]);
		//14 - mean of difference in reco - gen genDeltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin1 = new TH1D("geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin1","geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin1",xbins.size()-1,&xbins[0]);
		//15 - mean of difference in reco - gen genDeltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin2 = new TH1D("geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin2","geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin2",xbins.size()-1,&xbins[0]);
		//16 - mean of difference in reco - gen genDeltaTs as a function of total E of rhs that go into PV time calculation
		TH1D* geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin3 = new TH1D("geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin3","geoEavg_meanDeltaTime_recoGen_genDeltaTpvGambin3",xbins.size()-1,&xbins[0]);
		//17 - gen deltaT bw photon and pv, reco E bin 1
		TH1D* deltaT_pvGam_gen_Ebin1 = new TH1D("deltaT_gamPV_gen_Ebin1","deltaT_gamPV_gen_Ebin1",25,0,12);	
		//18 - gen deltaT bw photon and pv, reco E bin 2
		TH1D* deltaT_pvGam_gen_Ebin2 = new TH1D("deltaT_gamPV_gen_Ebin2","deltaT_gamPV_gen_Ebin2",25,0,12);	
		//19 - gen deltaT bw photon and pv, reco E bin 3
		TH1D* deltaT_pvGam_gen_Ebin3 = new TH1D("deltaT_gamPV_gen_Ebin3","deltaT_gamPV_gen_Ebin3",25,0,12);	
		//20 - gen deltaT bw photon and pv, reco E bin 4
		TH1D* deltaT_pvGam_gen_Ebin4 = new TH1D("deltaT_gamPV_gen_Ebin4","deltaT_gamPV_gen_Ebin4",25,0,12);	
		//21 - resolution of difference in gam - pv
		TH1D* geoEavg_sigmaDeltaTime_gamPV = new TH1D("geoEavg_sigmaDeltaTime_gamPV","geoEavg_sigmaDeltaTime_gamPV",xbins.size()-1,&xbins[0]);

		//22 - nrhs per jet (E bin1)
		TH1D* jetNrhs_Ebin1 = new TH1D("jetNrhs_Ebin0to200","jetNrhs_Ebin0to200",50,0,100);
		//23 - nrhs per jet (E bin2)
		TH1D* jetNrhs_Ebin2 = new TH1D("jetNrhs_Ebin200to400","jetNrhs_Ebin200to400",50,0,100);
		//24 - nrhs per jet (E bin3)
		TH1D* jetNrhs_Ebin3 = new TH1D("jetNrhs_Ebin400to600","jetNrhs_Ebin400to600",50,0,100);
		//25 - nrhs per jet (E bin4)
		TH1D* jetNrhs_Ebin4 = new TH1D("jetNrhs_Ebin600to800","jetNrhs_Ebin600to800",50,0,100);
		//26 - nrhs per jet (E bin5)
		TH1D* jetNrhs_Ebin5 = new TH1D("jetNrhs_Ebin800to1000","jetNrhs_Ebin800to1000",50,0,100);
		//27 - nrhs per jet (E bin6)
		TH1D* jetNrhs_Ebin6 = new TH1D("jetNrhs_Ebin1000to1200","jetNrhs_Ebin1000to1200",50,0,100);
		//28 - jet pt (E bin1)
		TH1D* jetPt_Ebin1 = new TH1D("jetPt_Ebin0to200","jetPt_Ebin0to200",50,0,2000);
		//29 - jet pt (E bin2)
		TH1D* jetPt_Ebin2 = new TH1D("jetPt_Ebin200to400","jetPt_Ebin200to400",50,0,2000);
		//30 - jet pt (E bin3)
		TH1D* jetPt_Ebin3 = new TH1D("jetPt_Ebin400to600","jetPt_Ebin400to600",50,0,2000);
		//31 - jet pt (E bin4)
		TH1D* jetPt_Ebin4 = new TH1D("jetPt_Ebin600to800","jetPt_Ebin600to800",50,0,2000);
		//32 - jet pt (E bin5)
		TH1D* jetPt_Ebin5 = new TH1D("jetPt_Ebin800to1000","jetPt_Ebin800to1000",50,0,2000);
		//33 - jet pt (E bin6)
		TH1D* jetPt_Ebin6 = new TH1D("jetPt_Ebin1000to1200","jetPt_Ebin1000to1200",50,0,2000);
		//34 - jet eta (E bin1)
		TH1D* jetEta_Ebin1 = new TH1D("jetEta_Ebin0to200","jetEta_Ebin0to200",50,-1.5,1.5);
		//35 - jet eta (E bin2)
		TH1D* jetEta_Ebin2 = new TH1D("jetEta_Ebin200to400","jetEta_Ebin200to400",50,-1.5,1.5);
		//36 - jet eta (E bin3)
		TH1D* jetEta_Ebin3 = new TH1D("jetEta_Ebin400to600","jetEta_Ebin400to600",50,-1.5,1.5);
		//37 - jet eta (E bin4)
		TH1D* jetEta_Ebin4 = new TH1D("jetEta_Ebin600to800","jetEta_Ebin600to800",50,-1.5,1.5);
		//38 - jet eta (E bin5)
		TH1D* jetEta_Ebin5 = new TH1D("jetEta_Ebin800to1000","jetEta_Ebin800to1000",50,-1.5,1.5);
		//39 - jet eta (E bin6)
		TH1D* jetEta_Ebin6 = new TH1D("jetEta_Ebin1000to1200","jetEta_Ebin1000to1200",50,-1.5,1.5);
		//40 - jet phi (E bin1)
		TH1D* jetPhi_Ebin1 = new TH1D("jetPhi_Ebin0to200","jetPhi_Ebin0to200",50,-0.2,6.4);
		//41 - jet phi (E bin2)
		TH1D* jetPhi_Ebin2 = new TH1D("jetPhi_Ebin200to400","jetPhi_Ebin200to400",50,-0.2,6.4);
		//42 - jet phi (E bin3)
		TH1D* jetPhi_Ebin3 = new TH1D("jetPhi_Ebin400to600","jetPhi_Ebin400to600",50,-0.2,6.4);
		//43 - jet phi (E bin4)
		TH1D* jetPhi_Ebin4 = new TH1D("jetPhi_Ebin600to800","jetPhi_Ebin600to800",50,-0.2,6.4);
		//44 - jet phi (E bin5)
		TH1D* jetPhi_Ebin5 = new TH1D("jetPhi_Ebin800to1000","jetPhi_Ebin800to1000",50,-0.2,6.4);
		//45 - jet phi (E bin6)
		TH1D* jetPhi_Ebin6 = new TH1D("jetPhi_Ebin1000to1200","jetPhi_Ebin1000to1200",50,-0.2,6.4);
		//46 - jet time (E bin1)
		TH1D* jetTime_Ebin1 = new TH1D("jetTime_Ebin0to200","jetTime_Ebin0to200",50,-10,10);
		//47 - jet time (E bin2)
		TH1D* jetTime_Ebin2 = new TH1D("jetTime_Ebin200to400","jetTime_Ebin200to400",50,-10,10);
		//48 - jet time (E bin3)
		TH1D* jetTime_Ebin3 = new TH1D("jetTime_Ebin400to600","jetTime_Ebin400to600",50,-10,10);
		//49 - jet time (E bin4)
		TH1D* jetTime_Ebin4 = new TH1D("jetTime_Ebin600to800","jetTime_Ebin600to800",50,-10,10);
		//50 - jet time (E bin5)
		TH1D* jetTime_Ebin5 = new TH1D("jetTime_Ebin800to1000","jetTime_Ebin800to1000",50,-10,10);
		//51 - jet time (E bin6)
		TH1D* jetTime_Ebin6 = new TH1D("jetTime_Ebin1000to1200","jetTime_Ebin1000to1200",50,-10,10);
		//52 - n subclusters (E bin1)
		TH1D* jetNSubclusters_Ebin1 = new TH1D("jetNSubclusters_Ebin0to200","jetNSubclusters_Ebin0to200",50,-0,10);
		//53 - n subclusters (E bin2)
		TH1D* jetNSubclusters_Ebin2 = new TH1D("jetNSubclusters_Ebin200to400","jetNSubclusters_Ebin200to400",50,-0,10);
		//54 - n subclusters (E bin3)
		TH1D* jetNSubclusters_Ebin3 = new TH1D("jetNSubclusters_Ebin400to600","jetNSubclusters_Ebin400to600",50,-0,10);
		//55 - n subclusters (E bin4)
		TH1D* jetNSubclusters_Ebin4 = new TH1D("jetNSubclusters_Ebin600to800","jetNSubclusters_Ebin600to800",50,-0,10);
		//56 - n subclusters (E bin5)
		TH1D* jetNSubclusters_Ebin5 = new TH1D("jetNSubclusters_Ebin800to1000","jetNSubclusters_Ebin800to1000",50,-0,10);
		//57 - n subclusters (E bin6)
		TH1D* jetNSubclusters_Ebin6 = new TH1D("jetNSubclusters_Ebin1000to1200","jetNSubclusters_Ebin1000to1200",50,-0,10);
		//58 - swissCross, -15 < t < -9
		TH1D* swCross_rhTimeNeg15toNeg9 = new TH1D("swCross_rhTimeNeg15toNeg9","swCross_rhTimeNeg15toNeg9",50,-0.05,1);
		//59 - swissCross, -6 < t < -3
		TH1D* swCross_rhTimeNeg6toNeg3 = new TH1D("swCross_rhTimeNeg6toNeg3","swCross_rhTimeNeg6toNeg3",50,-0.05,1);
		//60 - swissCross, time ~ -12
		TH1D* swCross_rhTimeNeg12 = new TH1D("swCross_rhTimeNeg12","swCross_rhTimeNeg12",50,-0.05,1);
		//61 - swissCross, time ~ -5
		TH1D* swCross_rhTimeNeg5 = new TH1D("swCross_rhTimeNeg5","swCross_rhTimeNeg5",50,-0.05,1);
		//62 - swissCross, time ~ 0
		TH1D* swCross_rhTime0 = new TH1D("swCross_rhTime0","swCross_rhTime0",50,-0.05,1);
		//63 - dr bw rh and track, time ~ -12
		TH1D* dRtrack_rhTimeNeg12 = new TH1D("dRtrack_rhTimeNeg12","dRtrack_rhTimeNeg12",50,-0.01,0.1);
		//64 - dr bw rh and track, time ~ 0
		TH1D* dRtrack_rhTimeNeg5 = new TH1D("dRtrack_rhTimeNeg5","dRtrack_rhTimeNeg5",50,-0.01,0.1);
		//65 - dr bw rh and track, time ~ 0
		TH1D* dRtrack_rhTime0 = new TH1D("dRtrack_rhTime0","dRtrack_rhTime0",50,-0.01,0.1);
		//66 - rh time of jets in PV time
		TH1D* rhTime = new TH1D("rhTime","rhTime",100,-10,10); 
		//67 - profiled seed crystal energy vs time
		TH1D* seedXtalEnergy_sigmaDeltaTime = new TH1D("seedXtalEnergy_sigmaDeltaTime","seedXtalEnergy_sigmaDeltaTime",10,0,100);
		//68 - max e rh photon time
		TH1D* gamTime_maxErh = new TH1D("gamTime_maxErh","gamTime_maxErh",100,-1,1); 
		//69 - seed rh photon time
		TH1D* seedXtalTime = new TH1D("seedXtalTime","seedXtalTime",50,-2,2);	
		//70 - # subclusters from mm
		TH1D* nSubclusters = new TH1D("nSubclusters","nSubclusters",10,0,10);
		//71 - # subclusters from mm - dijet res jets
		TH1D* nSubclusters_dijetSel = new TH1D("nSubclusters_dijetSel","nSubclusters_dijetSel",10,0,10);
	
		//0 - 2D histogram for reco-gen resolution
		TH2D* geoEavg_diffDeltaTime_recoGen = new TH2D("geoEavg_diffDeltaTime_recoGen","geoEavg_diffDeltaTime_recoGen;#sqrt{E^{pho}_{rh} #times E^{jets}_{rh}} (GeV);#Delta t^{PV,#gamma}_{reco, gen} (ns)",xbins.size()-1,&xbins[0],100,-2,2);
		//1 - 2D histogram for dijets resolution - geometric avg of jet pT
		TH2D* geopTavg_diffDeltaTime_dijets = new TH2D("geopTavg_diffDeltaTime_dijets","geopTavg_diffDeltaTime_dijets;#sqrt{pT^{jet1} #times pT^{jet2}} (GeV); #Delta t^{PV}_{dijet}",10,0,1000,120,-4,4);
		//2 - 2D histogram for dijets resolution - min E of jets
		TH2D* minpT_diffDeltaTime_dijets = new TH2D("minpT_diffDeltaTime_dijets","minpT_diffDeltaTime_dijets;min(pT^{jet1}, pT^{jet2}) (GeV); #Delta t^{PV}_{dijet}",10,0,1000,120,-4,4);
		
		//3 - 2D histogram for dijets resolution - sum_rh E_rh of jets
		//variable binning
		TH2D* geoAvgEecal_diffDeltaTime_dijets = new TH2D("geoAvgEecal_diffDeltaTime_dijets","geoAvgEecal_diffDeltaTime_dijets;#sqrt{E^{jet 1}_{ECAL} #times E^{jet 2}_{ECAL}} (GeV); #Delta t^{PV}_{dijet}",xbins.size()-1,&xbins[0],120,-4,4);	
		//4 - 2D histogram for reco-gen resolution - genDeltaTpvGam ~ [3.5,4.5)
		TH2D* geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin1 = new TH2D("geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin1","geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin1;#sqrt{E^{pho}_{rh} #times E^{jets}_{rh}} (GeV);#Delta t^{PV,#gamma}_{reco, gen} (ns)",xbins.size()-1,&xbins[0],100,-2,2);
		//5 - 2D histogram for reco-gen resolution - genDeltaTpvGam ~ [4,8)
		TH2D* geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin2 = new TH2D("geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin2","geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin2;#sqrt{E^{pho}_{rh} #times E^{jets}_{rh}} (GeV);#Delta t^{PV,#gamma}_{reco, gen} (ns)",xbins.size()-1,&xbins[0],100,-2,2);
		//6 - 2D histogram for reco-gen resolution - genDeltaTpvGam ~ [8,12)
		TH2D* geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin3 = new TH2D("geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin3","geoEavg_diffDeltaTime_recoGen_genDeltaTpvGambin3;#sqrt{E^{pho}_{rh} #times E^{jets}_{rh}} (GeV);#Delta t^{PV,#gamma}_{reco, gen} (ns)",xbins.size()-1,&xbins[0],100,-2,2);
		//7 - mean of diff recoGen deltaT distribution as a function of geoEavg and gen deltaT
		TH2D* geoEavg_genDeltaTime_meanRecoGenDeltaT = new TH2D("geoEavg_genDeltaTime_meanRecoGenDeltaT","geoEavg_genDeltaTime_meanRecoGenDeltaT;geoEavg;genDeltaTime;meanRecoGenDeltaT",xbins.size()-1,&xbins[0],3,0,3);
		

		//below are for checking gen deltaT calculation
		//8 - gen delta time vs reco delta time for signal photons - 0 <= E < 100
		TH2D* genDeltaT_recoDeltaT_Ebin1 = new TH2D("genDeltaT_recoDeltaT_Ebin1","genDeltaT_recoDeltaT_Ebin1;genDeltaT_Ebin1;recoDeltaT;a.u.",25,0,15,25,0,15);
		//9 - gen delta time vs reco delta time for signal photons - 100 <= E < 400
		TH2D* genDeltaT_recoDeltaT_Ebin2 = new TH2D("genDeltaT_recoDeltaT_Ebin2","genDeltaT_recoDeltaT_Ebin2;genDeltaT_Ebin2;recoDeltaT;a.u.",25,0,15,25,0,15);
		//10 - gen delta time vs reco delta time for signal photons - 400 <= E < 700
		TH2D* genDeltaT_recoDeltaT_Ebin3 = new TH2D("genDeltaT_recoDeltaT_Ebin3","genDeltaT_recoDeltaT_Ebin3;genDeltaT_Ebin3;recoDeltaT;a.u.",25,0,15,25,0,15);
		//11 - gen delta time vs reco delta time for signal photons - 700 <= E
		TH2D* genDeltaT_recoDeltaT_Ebin4 = new TH2D("genDeltaT_recoDeltaT_Ebin4","genDeltaT_recoDeltaT_Ebin4;genDeltaT_Ebin4;recoDeltaT;a.u.",25,0,15,25,0,15);
		//12 - gen-matched dr vs ratio of reco to gen deltaT
		TH2D* recoGenDr_recoGenDeltaTRatio = new TH2D("recoGenDr_recoGenDeltaTRatio","recoGenDr_recoGenDeltaTRatio;recoGenDr;recoGenDeltaTRatio;a.u.",25,0,0.5,25,0,5);
		//13 - gen energy vs ratio of reco to gen deltaT
		TH2D* genEnergy_recoGenDeltaTRatio = new TH2D("genEnergy_recoGenDeltaTRatio","genEnergy_recoGenDeltaTRatio;genEnergy;recoGenDeltaTRatio;a.u.",25,0,1000,25,0,5);
		//14 - reco energy vs gen deltaT
		TH2D* recoEnergy_genDeltaTpvGam = new TH2D("recoEnergy_genDeltaTpvGam","recoEnergy_genDeltaTpvGam;recoEnergy;genDeltaTpvGam;a.u.",25,0,500,25,3.5,12);
		//15 - reco energy/gen energy vs reco deltaT - gen deltaT
		TH2D* recoGenEnergyRatio_recoGenDeltaT = new TH2D("recoGenEnergyRatio_recoGenDeltaT","recoGenEnergyRatio_recoGenDeltaT;recoGenEnergyRatio;recoGenDeltaT;a.u.",25,0,2,25,-10,10);
		
		//16 - 2D histogram for gamPV resolution 
		TH2D* geoEavg_diffDeltaTime_gamPV = new TH2D("geoEavg_diffDeltaTime_gamPV","geoEavg_diffDeltaTime_gamPV;#sqrt{E^{pho}_{rh} #times E^{jets}_{rh}} (GeV);#Delta t_{PV,#gamma} (ns)",xbins.size()-1,&xbins[0],100,-4,4);
		
		//17 - jet time vs jet E
		TH2D* jetTime_Energy = new TH2D("jetTime_Energy","jetTime_Energy;jetTime;Energy",50,-20,20,50,0,500);
		//18 - rh time vs rh E
		TH2D* rhTime_Energy = new TH2D("rhTime_Energy","rhTime_Energy;rhTime;Energy",50,-20,20,50,0,500);
		//19 - rh time v eta, time + E cuts
		TH2D* rhTime_eta = new TH2D("rhTime_eta","rhTime_eta;rhTime;eta",50,-20,20,50,-1.5,1.5);
		//20 - rh eta vs phi, time + E cuts
		TH2D* rhPhi_eta = new TH2D("rhPhi_eta","rhPhi_eta;rhPhi;eta",50,-3.4,3.4,50,-1.5,1.5);
		//21 - E/p vs drTrack, time + E cuts	
		TH2D* rhEovP_dRtrack = new TH2D("rhEovP_dRtrack","rhEovP_dRtrack;rhEovP;dRtrack",50,0,15,50,0,0.09);
		//22 - rh time vs rh E zoomed
		TH2D* rhTime_Energy_zoomed = new TH2D("rhTime_Energy_zoomed","rhTime_Energy_zoomed;rhTime;Energy_zoomed",50,-15,0,50,0,200);
		//23 - sw+ vs rh time
		TH2D* swCross_rhTime = new TH2D("swCross_rhTime","swCross_rhTime;swCross;rhTime",50,0,1.1,50,-20,20);
		//24 - kWeird vs rh time
		TH2D* kWeird_rhTime = new TH2D("kWeird_rhTime","kWeird_rhTime;kWeird;rhTime",50,0.01,0.08,50,-20,20);
		//25 - sw+ vs rh energy
		TH2D* swCross_rhEnergy = new TH2D("swCross_rhEnergy","swCross_rhEnergy;swissCross;rhEnergy",50,-0.05,1,50,0,500);	
		//26 - kWeird vs rh energy		
		TH2D* kWeird_rhEnergy = new TH2D("kWeird_rhEnergy","kWeird_rhEnergy;kWeird;rhEnergy",50,0.01,0.08,50,0,500);
		//27 - e4/e1 vs 0.02*log10(e1)+0.02 (kWeird == e4/e1 < 0.02*log10(e1)+0.02) for rh E > 50 && rh time < -2 ns
		TH2D* e4e1_kWeird_rhEge50rhTleNeg2 = new TH2D("e4e1_kWeird_rhEge50rhTleNeg2","e4e1_kWeird_rhEge50rhTleNeg2;e4e1;kWeird",50,-0.1,1.1,50,0.01,0.08);
		//28 - eta vs phi, -15 < t < -9 
		TH2D* rhEta_rhPhi_rhTimeNeg15toNeg9 = new TH2D("rhEta_rhPhi_rhTimeNeg15toNeg9","rhEta_rhPhi_rhTimeNeg15toNeg9;rhEta;rhPhi",50,-1.5,1.5,50,-3.4,3.4);
		//29 - eta vs phi, -6 < t < -3 
		TH2D* rhEta_rhPhi_rhTimeNeg6toNeg3 = new TH2D("rhEta_rhPhi_rhTimeNeg6toNeg3","rhEta_rhPhi_rhTimeNeg6toNeg3;rhEta;rhPhi",50,-1.5,1.5,50,-3.4,3.4);
		//30 - time vs eta, -15 < t < -9 
		TH2D* rhTime_rhEta_rhTimeNeg15toNeg9 = new TH2D("rhTime_rhEta_rhTimeNeg15toNeg9","rhTime_rhEta_rhTimeNeg15toNeg9;rhTime;rhEta",50,-20,20,50,-1.5,1.5);
		//31 - time vs eta, -6 < t < -3 
		TH2D* rhTime_rhEta_rhTimeNeg6toNeg3 = new TH2D("rhTime_rhEta_rhTimeNeg6toNeg3","rhTime_rhEta_rhTimeNeg6toNeg3;rhTime;rhEta",50,-20,20,50,-1.5,1.5);
		//32 - eta vs phi grid, overlaid neighbors normalized to center energy in 5x5 grid
		TH2D* normNeighbors = new TH2D("normNeighbors","normNeighbors;local ieta;local iphi",5,-2,3,5,-2,3);	
		//33 - sw+ vs rh energy, time ~ -12
		TH2D* swCross_rhEnergy_timeNeg12 = new TH2D("swCross_rhEnergy_timeNeg12","swCross_rhEnergy_timeNeg12;swCross;rhEnergy_timeNeg12",50,-0.05,1,50,0,500);	
		//34 - sw+ vs rh energy, time ~ -5
		TH2D* swCross_rhEnergy_timeNeg5 = new TH2D("swCross_rhEnergy_timeNeg5","swCross_rhEnergy_timeNeg5;swCross;rhEnergy_timeNeg5",50,-0.05,1,50,0,500);	
		//35 - sw+ vs rh energy, time ~ 0
		TH2D* swCross_rhEnergy_time0 = new TH2D("swCross_rhEnergy_time0","swCross_rhEnergy_time0;swCross;rhEnergy_time0",50,-0.05,1,50,0,500);	
		//36 - eta vs phi grid, overlaid neighbors energy in 5x5 grid
		TH2D* ENeighbors = new TH2D("ENeighbors","ENeighbors;local ieta;local iphi",5,-2,3,5,-2,3);	
		//37 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, time ~ -12
		TH2D* ENeighbors_timeNeg12 = new TH2D("ENeighbors_timeNeg12","ENeighbors_timeNeg12;local ieta;local iphi_timeNeg12",5,-2,3,5,-2,3);	
		//38 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, time ~ -5
		TH2D* ENeighbors_timeNeg5 = new TH2D("ENeighbors_timeNeg5","ENeighbors_timeNeg5;local ieta;local iphi_timeNeg5",5,-2,3,5,-2,3);	
		//39 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, time ~ -0
		TH2D* ENeighbors_time0 = new TH2D("ENeighbors_time0","ENeighbors_time0;local ieta;local iphi_time0",5,-2,3,5,-2,3);	
		//40 - eta vs phi, t ~ -12 
		TH2D* rhEta_rhPhi_rhTimeNeg12 = new TH2D("rhEta_rhPhi_rhTimeNeg12","rhEta_rhPhi_rhTimeNeg12;rhEta;rhPhi",50,-1.5,1.5,50,-3.4,3.4);
		//41 - eta vs phi, t ~ -5 
		TH2D* rhEta_rhPhi_rhTimeNeg5 = new TH2D("rhEta_rhPhi_rhTimeNeg5","rhEta_rhPhi_rhTimeNeg5;rhEta;rhPhi",50,-1.5,1.5,50,-3.4,3.4);
		//42 - eta vs phi, t ~ 0 
		TH2D* rhEta_rhPhi_rhTime0 = new TH2D("rhEta_rhPhi_rhTime0","rhEta_rhPhi_rhTime0;rhEta;rhPhi",50,-1.5,1.5,50,-3.4,3.4);
		//43 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, time ~ -12
		TH2D* ENeighborsSkipCenter_timeNeg12 = new TH2D("ENeighborsSkipCenter_timeNeg12","ENeighborsSkipCenter_timeNeg12;local ieta;local iphiSkipCenter_timeNeg12",5,-2,3,5,-2,3);	
		//44 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, time ~ -5
		TH2D* ENeighborsSkipCenter_timeNeg5 = new TH2D("ENeighborsSkipCenter_timeNeg5","ENeighborsSkipCenter_timeNeg5;local ieta;local iphiSkipCenter_timeNeg5",5,-2,3,5,-2,3);	
		//45 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, time ~ -0
		TH2D* ENeighborsSkipCenter_time0 = new TH2D("ENeighborsSkipCenter_time0","ENeighborsSkipCenter_time0;local ieta;local iphiSkipCenter_time0",5,-2,3,5,-2,3);	
		//46 - eta vs phi grid for whole jet in 11x11 grid
		TH2D* ENeighborsJet = new TH2D("ENeighborsJet","ENeighborsJet;local ieta;local iphi",61,-30,31,61,-30,31);	
		//47 - seed crystal energy vs time for profile
		TH2D* seedXtalEnergy_diffDeltaTime = new TH2D("seedXtalEnergy_diffDeltaTime","seedXtalEnergy_time;energy;time",10,0,100,50,-2,2);	
		//48 - # rhs vs # subclusters for jets
		TH2D* AK4Jet_nRhs_nSubclusters = new TH2D("AK4Jet_nRhs_nSubclusters","AK4Jet_nRhs_nSubclusters;AK4Jet_nRhs;nSubclusters;a.u.",100,15,115,15,0,15);
		//49 - # dr-matched superclusters vs # subclusters for jets
		TH2D* AK4Jet_nSuperclusters_nSubclusters = new TH2D("AK4Jet_nSuperclusters_nSubclusters","AK4Jet_nSuperclusters_nSubclusters;AK4Jet_nSuperclusters;nSubclusters;a.u.",15,0,15,15,0,15);
		//50 - # rhs vs # subclusters for photons
		TH2D* Photon_nRhs_nSubclusters = new TH2D("Photon_nRhs_nSubclusters","Photon_nRhs_nSubclusters;Photon_nRhs;nSubclusters;a.u.",80,0,80,5,0,5);
		//51 - subcl time vs subcl E
		TH2D* subclTime_Energy = new TH2D("subclTime_Energy","subclTime_Energy;subclTime;Energy",50,-20,20,50,0,500);
		//52 - subcl time vs subcl E - spikes
		TH2D* subclTime_Energy_spikes = new TH2D("subclTime_Energy_spikes","subclTime_Energy_spikes;subclTime;Energy_spikes",50,-20,20,50,0,500);
		//53 - subcl time vs subcl E - BH
		TH2D* subclTime_Energy_BH = new TH2D("subclTime_Energy_BH","subclTime_Energy_BH;subclTime;Energy_BH",50,-20,20,50,0,500);
		//54 - subcl time vs subcl E - physBkg
		TH2D* subclTime_Energy_physBkg = new TH2D("subclTime_Energy_physBkg","subclTime_Energy_physBkg;subclTime;Energy_physBkg",50,-20,20,50,0,500);
		//55 - subcl discriminator vs subcl E - spikes
		TH2D* subclNNDiscr_Energy_spikes = new TH2D("subclNNDiscr_Energy_spikes","subclNNDiscr_Energy_spikes;subclNNDiscr;Energy_spikes",50,-20,20,50,0,500);
		//56 - subcl discriminator vs subcl E - BH
		TH2D* subclNNDiscr_Energy_BH = new TH2D("subclNNDiscr_Energy_BH","subclNNDiscr_Energy_BH;subclNNDiscr;Energy_BH",50,-20,20,50,0,500);
		//57 - subcl discriminator vs subcl E - physBkg
		TH2D* subclNNDiscr_Energy_physBkg = new TH2D("subclNNDiscr_Energy_physBkg","subclNNDiscr_Energy_physBkg;subclNNDiscr;Energy_physBkg",50,-20,20,50,0,500);
		//58 - subcl discriminator vs subcl time - spikes
		TH2D* subclNNDiscr_time_spikes = new TH2D("subclNNDiscr_time_spikes","subclNNDiscr_time_spikes;subclNNDiscr;time_spikes",50,-20,20,50,0,500);
		//59 - subcl discriminator vs subcl time - BH
		TH2D* subclNNDiscr_time_BH = new TH2D("subclNNDiscr_time_BH","subclNNDiscr_time_BH;subclNNDiscr;time_BH",50,-20,20,50,0,500);
		//60 - subcl discriminator vs subcl time - physBkg
		TH2D* subclNNDiscr_time_physBkg = new TH2D("subclNNDiscr_time_physBkg","subclNNDiscr_time_physBkg;subclNNDiscr;time_physBkg",50,-20,20,50,0,500);
		//61 - # dr-matched superclusters vs # subclusters for dijet res jets
		TH2D* AK4Jet_nSuperclusters_nSubclusters_dijetSel = new TH2D("AK4Jet_nSuperclusters_nSubclusters_dijetSel","AK4Jet_nSuperclusters_nSubclusters_dijetSel;AK4Jet_nSuperclusters;nSubclusters_dijetSel;a.u.",15,0,15,15,0,15);
		//62 - rh time vs rh E - high E
		TH2D* rhTime_Energy_highE = new TH2D("rhTime_Energy_highE","rhTime_Energy_highE;rhTime;Energy_highE",50,-20,20,100,0,1000);


		vector<timeRecoCat> trCats;
		virtual void MakeTimeRecoCatHists(){
			//don't want to separate lead/not lead histograms
			MakeProcCats(_oname, false);
			timeRecoCat trmed(_timeHists1D, _timeHists2D, med, _procCats);
			timeRecoCat treavg(_timeHists1D, _timeHists2D, eavg, _procCats);
			timeRecoCat trmmavg(_timeHists1D, _timeHists2D, mmavg, _procCats);
			timeRecoCat tremax(_timeHists1D, _timeHists2D, emax, _procCats);
		
			trCats.push_back(trmed);
			trCats.push_back(trmmavg);	
			trCats.push_back(treavg);
			trCats.push_back(tremax);
		}


	
		void FillTruePhotonHists(vector<Jet>& phos){
			int nphos = phos.size();	
		
			vector<JetPoint> rhs;
			int rhidx;
			double pvx, pvy, pvz, tof;
			double rhx, rhy, rhz, maxE;
			JetPoint maxE_rh;
			for(int p = 0; p < nphos; p++){
				phos[p].GetJetPoints(rhs);
				pvx = phos[p].GetVertex().at(0);
				pvy = phos[p].GetVertex().at(1);
				pvz = phos[p].GetVertex().at(2);
				maxE = 0;

				int scidx = _base->Photon_scIndex->at(phos[p].GetUserIdx());
				int seedxtal_id = _base->SuperCluster_XtalSeedID->at(scidx);
				for(int r = 0; r < rhs.size(); r++){
					rhx = rhs[r].x();	
					rhy = rhs[r].y();	
					rhz = rhs[r].z();	
					tof = sqrt((rhx - pvx)*(rhx - pvx) + (rhy - pvy)*(rhy - pvy) + (rhz - pvz)*(rhz - pvz))/_c; 
					TOFgam_rh_pv->Fill(tof,_weight);
					
					if(rhs[r].E() > maxE){
						maxE = rhs[r].E();
						maxE_rh = rhs[r];
					}
					if(rhs[r].rhId() == seedxtal_id)
						trCats[0].procCats[0].hists1D[0][69]->Fill(rhs[r].t(), _weight);
						trCats[0].procCats[1].hists1D[0][69]->Fill(rhs[r].t(), _weight);
					
						trCats[0].procCats[1].hists2D[0][47]->Fill(rhs[r].E(), rhs[r].t(), _weight);		
						trCats[0].procCats[0].hists2D[0][47]->Fill(rhs[r].E(), rhs[r].t(), _weight);		
				}		
				trCats[0].procCats[0].hists1D[0][68]->Fill(maxE_rh.t(), _weight);
				trCats[0].procCats[1].hists1D[0][68]->Fill(maxE_rh.t(), _weight);
 
				//GaussianMixture* gmm = _subcluster(phos[p]);
				double nrhs = (double)phos[p].GetNRecHits();
				int n_k = phos[p].GetNConstituents(); //post PU cleaning
				trCats[0].procCats[0].hists2D[0][50]->Fill(nrhs, n_k, _weight);
				trCats[0].procCats[1].hists2D[0][50]->Fill(nrhs, n_k, _weight);
			}
		}


		//also subclusters jets
		void FillTrueJetHists(vector<Jet>& jets){
			int njets = jets.size();	
			nTrueJets->Fill((double)njets,_weight);
		
			double eECAL = 0;
			int ijet = 0;
			vector<JetPoint> rhs;
			for(int j = 0; j < njets; j++){
				e_nRhs->Fill(jets[j].E(), jets[j].GetNRecHits(), _weight);
				TrueJet_pT->Fill(jets[j].pt(), _weight);
				TrueJet_nRhs->Fill(jets[j].GetNRecHits(), _weight);
				ijet = jets[j].GetUserIdx();
				eECAL = _base->Jet_neEmEF->at(ijet) + _base->Jet_chEmEF->at(ijet);
				eECAL *= jets[j].E();
				TrueJet_EmE->Fill(eECAL, _weight);
				TrueJet_nConstituents->Fill(_base->Jet_nConstituents->at(ijet), _weight);
				
				jets[j].GetJetPoints(rhs);
				for(int r = 0; r < rhs.size(); r++){
					erhs_trhs->Fill(rhs[r].E(), rhs[r].t(), _weight);
				}		
				//GaussianMixture* gmm = _subcluster(jets[j]);
				//if(_cleansubcls) CleanSubclusters(gmm, jets[j]);
				double nrhs = (double)jets[j].GetNRecHits();
				int n_k = jets[j].GetNConstituents();
				trCats[0].procCats[0].hists2D[0][48]->Fill(nrhs, n_k, _weight);
				trCats[0].procCats[1].hists2D[0][48]->Fill(nrhs, n_k, _weight);
				int nscs = nSC_dRMatch(jets[j]);	
				trCats[0].procCats[0].hists2D[0][49]->Fill(nscs, n_k, _weight);
				trCats[0].procCats[1].hists2D[0][49]->Fill(nscs, n_k, _weight);
				
				trCats[TimeStrategy(1)].procCats[0].hists1D[0][70]->Fill(n_k, _weight);
				trCats[TimeStrategy(1)].procCats[1].hists1D[0][70]->Fill(n_k, _weight);
				
				//FillModelHists(gmm,jets[j]);
				FillModelHists(jets[j]);

				if(njets < 2) continue;
				if(j == 0 || j == 1){	
					TrueJet_twoHardestpT->Fill(jets[j].pt(), _weight);
					TrueJet_twoHardestpT->Fill(jets[j].pt(), _weight);
				}
			}
		}

		int nSC_dRMatch(const Jet& jet){ //put in E ratio req?
			double maxDr = 0.5;
			int nsc_tot = _SCs.size();
			int nsc = 0;
			double dr;
			for(int s = 0; s < nsc_tot; s++){
				dr = dR(jet.eta(), jet.phi(), _SCs[s].eta(), _SCs[s].phi());
				if(dr < maxDr)
					nsc++;
			} 
			return nsc;
		}
		
		//void FillModelHists(GaussianMixture* gmm, const Jet& jet){
		void FillModelHists(Jet jet){
			int n_k = jet.GetNConstituents();
			Matrix mu, cov;
			jet.GetClusterParams(mu, cov);

			//double Ek;
			//vector<double> norms;
			//gmm->GetNorms(norms);
			//map<string, Matrix> params;

			double subcl_dist_time = 0;
			double subcl_dist_etaphi = 0;
			double ec1, ec2, pc1, pc2, tc1, tc2, Ek;
			vector<pair<int, double>> subcl_predict;
			Jet subcl;
			for(int k = 0; k < n_k; k++){
				jet.GetConstituent(k, subcl);
				Ek = subcl.E();

				subclusterEfrac->Fill(Ek/jet.E());
				subcl.GetClusterParams(mu, cov);

				etaSig->Fill(cov.at(0,0));
				phiSig->Fill(cov.at(1,1));
				timeSig->Fill(cov.at(2,2));
	
				ec1 = mu.at(0,0);
				pc1 = mu.at(1,0);
				tc1 = mu.at(2,0);
				
				//posterior values of parameters from prior distributions
				//beta_k->Fill(params["scale"].at(0,0));
				//W_ee_k->Fill(params["scalemat"].at(0,0));
				//W_pp_k->Fill(params["scalemat"].at(1,1));
				//W_tt_k->Fill(params["scalemat"].at(2,2));
			
				//predicting det bkg MVA score for subclusters in jets
				//TODO - could change to grabbing rhs from SCs matched to jets 
				pair<int, double> class_discr = PredictSubcluster(subcl);
				subcl_predict.push_back(class_discr);
				int nclass = class_discr.first;
				if(nclass == 0)
					t_physBkg->Fill(tc1);
				if(nclass == 1)
					t_BH->Fill(tc1);
				if(nclass == 2)
					t_spike->Fill(tc1);

				//for(int kk = k+1; kk < n_k; kk++){
				//	params = gmm->GetLHPosteriorParameters(kk);
				//	ec2 = params["mean"].at(0,0);
				//	pc2 = params["mean"].at(1,0);
				//	tc2 = params["mean"].at(2,0);
				//	
				//	subcl_dist_time = (tc1 - tc2);
				//	subclDist_time->Fill(subcl_dist_time);
				//	
				//	double de = ec1 - ec2;
				//	double dp = pc1 - pc2;
				//	dp = acos(cos(dp)); //wraparound
				//	subcl_dist_etaphi = sqrt(de*de + dp*dp);
				//	subclDist_etaPhi->Fill(subcl_dist_etaphi);
				//}
			}
			//using subcl_predict predictions (class, score), sigclass == 0, minscore == 0.9, fully remove subclusters
			//TODO - verify signal class and minscore from ROC curve
			Jet detBkgCleanedJet;
			if(_cleansubcls) detBkgCleanedJet = jet.GenericClean(subcl_predict, 0, 0.9, true);
		}


		//need to break down by procCat - see how this is done in PhotonSkimmer.cc
		void FillPVTimeHists(vector<Jet> jets, vector<Jet> phos, int tr_idx){
			int njets = jets.size();
			double jettime = -999;
			vector<double> jettimes;
			double gamtime = -999;
			double pvtime = -999;
			double gamtime_pv = -999;
			double deltaT_gampv = -999;
			double deltaT_gampv_gen = -999;
			double ptavg, geoEavg, Epho, Erh, Ejets;
			int phoidx, genidx, phoid;
			double dphi_phoJets = -999;
			double pi = acos(-1);
			TimeStrategy ts = TimeStrategy(tr_idx);
			//break down by process for this tr method
			int nProc = trCats[tr_idx].procCats.size();
			double dtime = 1;
			for(int p = 0; p < nProc; p++){
				//cout << " proc " << trCats[tr_idx].procCats[p].plotName << endl;
				Erh = 0;
				ptavg = 0;
				geoEavg = 0;
			//cout << "\nmethod: " << tr_idx << " " << ts << " "<< trCats[tr_idx].methodName << " proc " << p << " with " << njets << " jets" << endl;
				for(int j = 0; j < njets; j++){
					//if(ts == mmavg) cout << "\njet #" << j << endl;
					//jets from RunCluster already have PU components remove and energy-weighted time calculated hit-by-hit with PU weights applied 
					jettime = CalcJetTime(ts, jets[j]);
					jets[j].SetJetTime(jettime);
					//fill jet map
					vector<JetPoint> rhs;
					jets[j].GetJetPoints(rhs);
					if(tr_idx == 0){
						vector<double> neighborEs_jet; 
						std::vector<std::pair<int, int>> icoords_jet;
						GetNeighborE(rhs, -1, icoords_jet, neighborEs_jet,false,61);
						for(int e = 0; e < (int)neighborEs_jet.size(); e++){
							trCats[tr_idx].procCats[p].hists2D[0][46]->Fill(icoords_jet[e].first, icoords_jet[e].second, neighborEs_jet[e]);
						}
					}
					//cout << " jet #" << j << " time: " << jettime << endl;
					//fill jet time in pv frame - 0
					trCats[tr_idx].procCats[p].hists1D[0][0]->Fill(jettime, _weight);
					for(int r = 0; r < rhs.size(); r++)
						Erh += rhs[r].E();
					rhs.clear();
					//add to pv time calculation
					jettimes.push_back(jettime);
				}
				//calculate jet time difference
				double Erh1, Erh2;
				if(njets > 1){
					pair<Jet,Jet> hardjets;
					int pair = FindJetPair(jets, hardjets);
					if(pair == -999) continue; //jets did not pass selection
					double deltaT_jets = GetDeltaTime(hardjets);

					ptavg = pow((hardjets.first.pt()*hardjets.second.pt()),0.5);
					Erh1 = 0;
					Erh2 = 0;	
					double bestp, bestdr, kweird;
					vector<JetPoint> rhs; hardjets.first.GetJetPoints(rhs);
					vector<double> neighborEs, neighborEs_noCenter, neighborEs_time;
					std::vector<std::pair<int, int>> icoords, icoords_noCenter, icoords_time;
				

					//fill subcl hists for dijet resolution
					int nk = hardjets.first.GetNConstituents();
					int nscs = nSC_dRMatch(hardjets.first);	
					trCats[tr_idx].procCats[p].hists2D[0][61]->Fill(nscs, nk, _weight);
					trCats[tr_idx].procCats[p].hists1D[0][71]->Fill(nk, _weight);
					Jet constit;
					for(int k = 0; k < nk; k++){
						hardjets.first.GetConstituent(k, constit);
						trCats[tr_idx].procCats[p].hists2D[0][51]->Fill(constit.t(), constit.E(), _weight);
					}
					nk = hardjets.second.GetNConstituents();
					nscs = nSC_dRMatch(hardjets.second);	
					trCats[tr_idx].procCats[p].hists2D[0][61]->Fill(nscs, nk, _weight);
					trCats[tr_idx].procCats[p].hists1D[0][71]->Fill(nk, _weight);
					for(int k = 0; k < nk; k++){
						hardjets.second.GetConstituent(k, constit);
						trCats[tr_idx].procCats[p].hists2D[0][51]->Fill(constit.t(), constit.E(), _weight);
					}

					for(int r = 0; r < rhs.size(); r++){
						if(rhs[r].E() < _minRhE) continue;
						//if(rhs[r].t() > -1) continue;
						Erh1 += rhs[r].E();
						//only need for one time reco method (doesnt depend on this)
						if(tr_idx == 0){
							trCats[tr_idx].procCats[p].hists2D[0][18]->Fill(rhs[r].t(), rhs[r].E(), _weight);
							trCats[tr_idx].procCats[p].hists2D[0][62]->Fill(rhs[r].t(), rhs[r].E(), _weight);
							trCats[tr_idx].procCats[p].hists2D[0][22]->Fill(rhs[r].t(), rhs[r].E(), _weight);
							TrackMatched(rhs[r], bestdr, bestp);
							if(rhs[r].E() > 50 && rhs[r].t() < -2){
								trCats[tr_idx].procCats[p].hists2D[0][19]->Fill(rhs[r].t(), rhs[r].eta(), _weight);
								trCats[tr_idx].procCats[p].hists2D[0][20]->Fill(rhs[r].phi(), rhs[r].eta(), _weight);
								if(bestdr != 999) trCats[tr_idx].procCats[p].hists2D[0][21]->Fill(rhs[r].E()/bestp, bestdr,_weight);
								trCats[tr_idx].procCats[p].hists2D[0][27]->Fill(1 - _base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), kweird,_weight);
							}
							if(rhs[r].t() > -15 && rhs[r].t() < -9){
								trCats[tr_idx].procCats[p].hists1D[0][58]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][28]->Fill(rhs[r].eta(), rhs[r].phi(),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][30]->Fill(rhs[r].time(), rhs[r].eta(),_weight);
							}
							if(rhs[r].t() > -6 && rhs[r].t() < -3){
								trCats[tr_idx].procCats[p].hists1D[0][59]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][29]->Fill(rhs[r].eta(), rhs[r].phi(),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][31]->Fill(rhs[r].time(), rhs[r].eta(),_weight);
							}
							kweird = 0.02*log10(rhs[r].E()) + 0.02;
							trCats[tr_idx].procCats[p].hists2D[0][23]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].t(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][24]->Fill(kweird, rhs[r].t(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][25]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][26]->Fill(kweird, rhs[r].E(),_weight);
							//2D map of neighbors
							GetNeighborE(rhs, r, icoords, neighborEs,false);
							GetNeighborE(rhs, r, icoords_noCenter, neighborEs_noCenter,true);
							for(int e = 0; e < (int)neighborEs.size(); e++){
								trCats[tr_idx].procCats[p].hists2D[0][32]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]/rhs[r].E()*_weight);
								trCats[tr_idx].procCats[p].hists2D[0][36]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
							}
							//time slice hists
							if(rhs[r].t() > -12-dtime && rhs[r].t() < -12+dtime){
								trCats[tr_idx].procCats[p].hists1D[0][60]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][33]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
								for(int e = 0; e < (int)neighborEs.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][37]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
								}
								for(int e = 0; e < (int)neighborEs_noCenter.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][43]->Fill(icoords_noCenter[e].first, icoords_noCenter[e].second, neighborEs_noCenter[e]*_weight);
								}
								trCats[tr_idx].procCats[p].hists1D[0][63]->Fill(bestdr,_weight);
								trCats[tr_idx].procCats[p].hists2D[0][40]->Fill(rhs[r].eta(),rhs[r].phi(),_weight);
							}
							if(rhs[r].t() > -5-dtime && rhs[r].t() < -5+dtime){
								trCats[tr_idx].procCats[p].hists1D[0][61]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][34]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
								for(int e = 0; e < (int)neighborEs.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][38]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
								}
								for(int e = 0; e < (int)neighborEs_noCenter.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][44]->Fill(icoords_noCenter[e].first, icoords_noCenter[e].second, neighborEs_noCenter[e]*_weight);
								}
								trCats[tr_idx].procCats[p].hists1D[0][64]->Fill(bestdr,_weight);
								trCats[tr_idx].procCats[p].hists2D[0][41]->Fill(rhs[r].eta(),rhs[r].phi(),_weight);
							}
							if(rhs[r].t() > -dtime && rhs[r].t() < dtime){
								trCats[tr_idx].procCats[p].hists1D[0][62]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][35]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
								for(int e = 0; e < (int)neighborEs.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][39]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
								}
								for(int e = 0; e < (int)neighborEs_noCenter.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][45]->Fill(icoords_noCenter[e].first, icoords_noCenter[e].second, neighborEs_noCenter[e]*_weight);
								}
								trCats[tr_idx].procCats[p].hists1D[0][65]->Fill(bestdr,_weight);
								trCats[tr_idx].procCats[p].hists2D[0][42]->Fill(rhs[r].eta(),rhs[r].phi(),_weight);
							}
						
						}
					}
					hardjets.second.GetJetPoints(rhs);
					for(int r = 0; r < rhs.size(); r++){
						if(rhs[r].E() < _minRhE) continue;
						//if(rhs[r].t() > -1) continue;
						Erh2 += rhs[r].E();
						//only need for one time reco method (doesnt depend on this)
						if(tr_idx == 0){
							trCats[tr_idx].procCats[p].hists2D[0][18]->Fill(rhs[r].t(), rhs[r].E(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][62]->Fill(rhs[r].t(), rhs[r].E(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][22]->Fill(rhs[r].t(), rhs[r].E(),_weight);
							if(rhs[r].E() > 50 && rhs[r].t() < -2){
								TrackMatched(rhs[r], bestdr, bestp);
								trCats[tr_idx].procCats[p].hists2D[0][19]->Fill(rhs[r].t(), rhs[r].eta(),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][20]->Fill(rhs[r].phi(), rhs[r].eta(),_weight);
								if(bestdr != 999) trCats[tr_idx].procCats[p].hists2D[0][21]->Fill(rhs[r].E()/bestp, bestdr,_weight);
								trCats[tr_idx].procCats[p].hists2D[0][27]->Fill(1 - _base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), kweird,_weight);
							}
							if(rhs[r].t() > -15 && rhs[r].t() < -9){
								trCats[tr_idx].procCats[p].hists1D[0][58]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][28]->Fill(rhs[r].eta(), rhs[r].phi(),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][30]->Fill(rhs[r].time(), rhs[r].eta(),_weight);
							}
							if(rhs[r].t() > -6 && rhs[r].t() < -3){
								trCats[tr_idx].procCats[p].hists1D[0][59]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][29]->Fill(rhs[r].eta(), rhs[r].phi(),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][31]->Fill(rhs[r].time(), rhs[r].eta(),_weight);
							}
							kweird = 0.02*log10(rhs[r].E()) + 0.02;
							trCats[tr_idx].procCats[p].hists2D[0][23]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].t(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][24]->Fill(kweird, rhs[r].t(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][25]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
							trCats[tr_idx].procCats[p].hists2D[0][26]->Fill(kweird, rhs[r].E(),_weight);
							//2D map of neighbors/
							GetNeighborE(rhs, r, icoords, neighborEs,false);
							GetNeighborE(rhs, r, icoords_noCenter, neighborEs_noCenter,true);
							for(int e = 0; e < (int)neighborEs.size(); e++){
								trCats[tr_idx].procCats[p].hists2D[0][32]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]/rhs[r].E()*_weight);
								trCats[tr_idx].procCats[p].hists2D[0][36]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
							}
							//time slice hists
							if(rhs[r].t() > -12-dtime && rhs[r].t() < -12+dtime){
								trCats[tr_idx].procCats[p].hists1D[0][60]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][33]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
								for(int e = 0; e < (int)neighborEs.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][37]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
								}
								for(int e = 0; e < (int)neighborEs_noCenter.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][43]->Fill(icoords_noCenter[e].first, icoords_noCenter[e].second, neighborEs_noCenter[e]*_weight);
								}
								trCats[tr_idx].procCats[p].hists1D[0][63]->Fill(bestdr,_weight);
                                                                trCats[tr_idx].procCats[p].hists2D[0][40]->Fill(rhs[r].eta(),rhs[r].phi(),_weight);
							}
							if(rhs[r].t() > -5-dtime && rhs[r].t() < -5+dtime){
								trCats[tr_idx].procCats[p].hists1D[0][61]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][34]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
								for(int e = 0; e < (int)neighborEs.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][38]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
								}
								for(int e = 0; e < (int)neighborEs_noCenter.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][44]->Fill(icoords_noCenter[e].first, icoords_noCenter[e].second, neighborEs_noCenter[e]*_weight);
								}
								trCats[tr_idx].procCats[p].hists1D[0][64]->Fill(bestdr,_weight);
                                                                trCats[tr_idx].procCats[p].hists2D[0][41]->Fill(rhs[r].eta(),rhs[r].phi(),_weight);
							}
							if(rhs[r].t() > -dtime && rhs[r].t() < dtime){
								trCats[tr_idx].procCats[p].hists1D[0][62]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()),_weight);
								trCats[tr_idx].procCats[p].hists2D[0][35]->Fill(_base->ECALRecHit_swCross->at(rhs[r].GetUserIdx()), rhs[r].E(),_weight);
								for(int e = 0; e < (int)neighborEs.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][39]->Fill(icoords[e].first, icoords[e].second, neighborEs[e]*_weight);
								}
								for(int e = 0; e < (int)neighborEs_noCenter.size(); e++){
									trCats[tr_idx].procCats[p].hists2D[0][45]->Fill(icoords_noCenter[e].first, icoords_noCenter[e].second, neighborEs_noCenter[e]*_weight);
								}

								trCats[tr_idx].procCats[p].hists1D[0][65]->Fill(bestdr,_weight);
                                                                trCats[tr_idx].procCats[p].hists2D[0][42]->Fill(rhs[r].eta(),rhs[r].phi(),_weight);
							}
						}
					}
	//cout << "deltaT_jets " << deltaT_jets << endl;
					//fill deltaT_jets - 1
					trCats[tr_idx].procCats[p].hists1D[0][1]->Fill(deltaT_jets, _weight);
					//fill geopTavg vs deltaT jets - 1
					trCats[tr_idx].procCats[p].hists2D[0][1]->Fill(ptavg, deltaT_jets, _weight);
					//fill minpt vs deltaT jets - 2
					trCats[tr_idx].procCats[p].hists2D[0][2]->Fill(hardjets.second.pt(), deltaT_jets, _weight);
					//fill geoAvgEecal vs deltaT jets - 3
					trCats[tr_idx].procCats[p].hists2D[0][3]->Fill(sqrt(Erh1*Erh2), deltaT_jets, _weight);
					//trCats[tr_idx].procCats[p].hists2D[0][3]->Fill(Erh, deltaT_jets);
					//fill jet property hists
					if(xbins[0] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[1]){
						//cout << " bin 1 - E " << sqrt(Erh1*Erh2) << " nrhs " << hardjets.first.GetNRecHits() << " pt " << hardjets.first.pt() << " eta " << hardjets.first.eta() << " phi " << hardjets.first.phi() << " time " << hardjets.first.time() << " n subclusters " << hardjets.first.GetNConstituents() << endl;
						trCats[tr_idx].procCats[p].hists1D[0][22]->Fill(hardjets.first.GetNRecHits());
						trCats[tr_idx].procCats[p].hists1D[0][22]->Fill(hardjets.second.GetNRecHits());
						
						trCats[tr_idx].procCats[p].hists1D[0][28]->Fill(hardjets.first.pt());
						trCats[tr_idx].procCats[p].hists1D[0][28]->Fill(hardjets.second.pt());
						
						trCats[tr_idx].procCats[p].hists1D[0][34]->Fill(hardjets.first.eta());
						trCats[tr_idx].procCats[p].hists1D[0][34]->Fill(hardjets.second.eta());
						
						trCats[tr_idx].procCats[p].hists1D[0][40]->Fill(hardjets.first.phi());
						trCats[tr_idx].procCats[p].hists1D[0][40]->Fill(hardjets.second.phi());
						
						trCats[tr_idx].procCats[p].hists1D[0][46]->Fill(hardjets.first.time());
						trCats[tr_idx].procCats[p].hists1D[0][46]->Fill(hardjets.second.time());
					
						if(ts == mmavg){	
							trCats[tr_idx].procCats[p].hists1D[0][52]->Fill(hardjets.first.GetNConstituents());
							trCats[tr_idx].procCats[p].hists1D[0][52]->Fill(hardjets.second.GetNConstituents());
						}

					}
					else if(xbins[1] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[2]){
						//cout << " bin 2 - E " << sqrt(Erh1*Erh2) << " nrhs " << hardjets.first.GetNRecHits() << " pt " << hardjets.first.pt() << " eta " << hardjets.first.eta() << " phi " << hardjets.first.phi() << " time " << hardjets.first.time() << " n subclusters " << hardjets.first.GetNConstituents() << endl;
						trCats[tr_idx].procCats[p].hists1D[0][23]->Fill(hardjets.first.GetNRecHits());
						trCats[tr_idx].procCats[p].hists1D[0][23]->Fill(hardjets.second.GetNRecHits());
						
						trCats[tr_idx].procCats[p].hists1D[0][29]->Fill(hardjets.first.pt());
						trCats[tr_idx].procCats[p].hists1D[0][29]->Fill(hardjets.second.pt());
						
						trCats[tr_idx].procCats[p].hists1D[0][35]->Fill(hardjets.first.eta());
						trCats[tr_idx].procCats[p].hists1D[0][35]->Fill(hardjets.second.eta());
						
						trCats[tr_idx].procCats[p].hists1D[0][41]->Fill(hardjets.first.phi());
						trCats[tr_idx].procCats[p].hists1D[0][41]->Fill(hardjets.second.phi());
						
						trCats[tr_idx].procCats[p].hists1D[0][47]->Fill(hardjets.first.time());
						trCats[tr_idx].procCats[p].hists1D[0][47]->Fill(hardjets.second.time());
					
						if(ts == mmavg){	
							trCats[tr_idx].procCats[p].hists1D[0][53]->Fill(hardjets.first.GetNConstituents());
							trCats[tr_idx].procCats[p].hists1D[0][53]->Fill(hardjets.second.GetNConstituents());
						}
					
					}
					else if(xbins[2] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[3]){
						//cout << " bin 3 - E " << sqrt(Erh1*Erh2) << " nrhs " << hardjets.first.GetNRecHits() << " pt " << hardjets.first.pt() << " eta " << hardjets.first.eta() << " phi " << hardjets.first.phi() << " time " << hardjets.first.time() << " n subclusters " << hardjets.first.GetNConstituents() << endl;
						trCats[tr_idx].procCats[p].hists1D[0][24]->Fill(hardjets.first.GetNRecHits());
						trCats[tr_idx].procCats[p].hists1D[0][24]->Fill(hardjets.second.GetNRecHits());

						trCats[tr_idx].procCats[p].hists1D[0][30]->Fill(hardjets.first.pt());
						trCats[tr_idx].procCats[p].hists1D[0][30]->Fill(hardjets.second.pt());
						
						trCats[tr_idx].procCats[p].hists1D[0][36]->Fill(hardjets.first.eta());
						trCats[tr_idx].procCats[p].hists1D[0][36]->Fill(hardjets.second.eta());
						
						trCats[tr_idx].procCats[p].hists1D[0][42]->Fill(hardjets.first.phi());
						trCats[tr_idx].procCats[p].hists1D[0][42]->Fill(hardjets.second.phi());
						
						trCats[tr_idx].procCats[p].hists1D[0][48]->Fill(hardjets.first.time());
						trCats[tr_idx].procCats[p].hists1D[0][48]->Fill(hardjets.second.time());
					
						if(ts == mmavg){	
							trCats[tr_idx].procCats[p].hists1D[0][54]->Fill(hardjets.first.GetNConstituents());
							trCats[tr_idx].procCats[p].hists1D[0][54]->Fill(hardjets.second.GetNConstituents());
						}
					}
					else if(xbins[3] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[4]){
						//cout << " bin 4 - E " << sqrt(Erh1*Erh2) << " nrhs " << hardjets.first.GetNRecHits() << " pt " << hardjets.first.pt() << " eta " << hardjets.first.eta() << " phi " << hardjets.first.phi() << " time " << hardjets.first.time() << " n subclusters " << hardjets.first.GetNConstituents() << endl;
						trCats[tr_idx].procCats[p].hists1D[0][25]->Fill(hardjets.first.GetNRecHits());
						trCats[tr_idx].procCats[p].hists1D[0][25]->Fill(hardjets.second.GetNRecHits());
						
						trCats[tr_idx].procCats[p].hists1D[0][31]->Fill(hardjets.first.pt());
						trCats[tr_idx].procCats[p].hists1D[0][31]->Fill(hardjets.second.pt());
						
						trCats[tr_idx].procCats[p].hists1D[0][37]->Fill(hardjets.first.eta());
						trCats[tr_idx].procCats[p].hists1D[0][37]->Fill(hardjets.second.eta());
						
						trCats[tr_idx].procCats[p].hists1D[0][43]->Fill(hardjets.first.phi());
						trCats[tr_idx].procCats[p].hists1D[0][43]->Fill(hardjets.second.phi());
						
						trCats[tr_idx].procCats[p].hists1D[0][49]->Fill(hardjets.first.time());
						trCats[tr_idx].procCats[p].hists1D[0][49]->Fill(hardjets.second.time());
					
						if(ts == mmavg){	
							trCats[tr_idx].procCats[p].hists1D[0][55]->Fill(hardjets.first.GetNConstituents());
							trCats[tr_idx].procCats[p].hists1D[0][55]->Fill(hardjets.second.GetNConstituents());
						}

					}
					else if(xbins[4] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[5]){
						//cout << " bin 5 - E " << sqrt(Erh1*Erh2) << " nrhs " << hardjets.first.GetNRecHits() << " pt " << hardjets.first.pt() << " eta " << hardjets.first.eta() << " phi " << hardjets.first.phi() << " time " << hardjets.first.time() << " n subclusters " << hardjets.first.GetNConstituents() << endl;
						trCats[tr_idx].procCats[p].hists1D[0][26]->Fill(hardjets.first.GetNRecHits());
						trCats[tr_idx].procCats[p].hists1D[0][26]->Fill(hardjets.second.GetNRecHits());
						
						trCats[tr_idx].procCats[p].hists1D[0][32]->Fill(hardjets.first.pt());
						trCats[tr_idx].procCats[p].hists1D[0][32]->Fill(hardjets.second.pt());
						
						trCats[tr_idx].procCats[p].hists1D[0][38]->Fill(hardjets.first.eta());
						trCats[tr_idx].procCats[p].hists1D[0][38]->Fill(hardjets.second.eta());
						
						trCats[tr_idx].procCats[p].hists1D[0][44]->Fill(hardjets.first.phi());
						trCats[tr_idx].procCats[p].hists1D[0][44]->Fill(hardjets.second.phi());
						
						trCats[tr_idx].procCats[p].hists1D[0][50]->Fill(hardjets.first.time());
						trCats[tr_idx].procCats[p].hists1D[0][50]->Fill(hardjets.second.time());
					
						if(ts == mmavg){	
							trCats[tr_idx].procCats[p].hists1D[0][56]->Fill(hardjets.first.GetNConstituents());
							trCats[tr_idx].procCats[p].hists1D[0][56]->Fill(hardjets.second.GetNConstituents());
						}

					}
					/*
					else if(xbins[5] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[6]){
						//cout << " bin 6 - E " << sqrt(Erh1*Erh2) << " nrhs " << hardjets.first.GetNRecHits() << " pt " << hardjets.first.pt() << " eta " << hardjets.first.eta() << " phi " << hardjets.first.phi() << " time " << hardjets.first.time() << " n subclusters " << hardjets.first.GetNConstituents() << endl;
						trCats[tr_idx].procCats[p].hists1D[0][27]->Fill(hardjets.first.GetNRecHits());
						trCats[tr_idx].procCats[p].hists1D[0][27]->Fill(hardjets.second.GetNRecHits());
						
						trCats[tr_idx].procCats[p].hists1D[0][33]->Fill(hardjets.first.pt());
						trCats[tr_idx].procCats[p].hists1D[0][33]->Fill(hardjets.second.pt());
						
						trCats[tr_idx].procCats[p].hists1D[0][39]->Fill(hardjets.first.eta());
						trCats[tr_idx].procCats[p].hists1D[0][39]->Fill(hardjets.second.eta());
						
						trCats[tr_idx].procCats[p].hists1D[0][45]->Fill(hardjets.first.phi());
						trCats[tr_idx].procCats[p].hists1D[0][45]->Fill(hardjets.second.phi());
						
						trCats[tr_idx].procCats[p].hists1D[0][51]->Fill(hardjets.first.time());
						trCats[tr_idx].procCats[p].hists1D[0][51]->Fill(hardjets.second.time());
					
						if(ts == mmavg){	
							trCats[tr_idx].procCats[p].hists1D[0][57]->Fill(hardjets.first.GetNConstituents());
							trCats[tr_idx].procCats[p].hists1D[0][57]->Fill(hardjets.second.GetNConstituents());
						}

					}
					*/
					else{ }
					trCats[tr_idx].procCats[p].hists2D[0][17]->Fill(hardjets.first.time(), Erh1);
					trCats[tr_idx].procCats[p].hists2D[0][17]->Fill(hardjets.second.time(), Erh2);
					
				}	
				//this assumes that the time for the jet was set previously with the respective method
				pvtime = CalcPVTime(ts, jets);
				//cout << "pv time " << pvtime << endl;
				Ejets = 0;
				vector<Jet> jets_highpt;
				vector<JetPoint> jrhs;
				for(auto j : jets){
					if(j.pt() > 50) jets_highpt.push_back(j);
					j.GetJetPoints(jrhs);
					for(auto r : jrhs){
						Ejets += r.E();
						trCats[tr_idx].procCats[p].hists1D[0][66]->Fill(r.t());
					}

				}
				double pvtime_highpt = CalcPVTime(ts, jets_highpt);
				//only fill for two leading photons + weighted avg of jet time
				if(phos.size() < 1) continue;
				vector<JetPoint> phorhs; 
				if(_data){
					//want difference in PV frame (no shift to detector)
					gamtime = CalcJetTime(ts, phos[0], false);
					deltaT_gampv = gamtime - pvtime;
					cout << "data - gam time - pv time " << deltaT_gampv << endl;
					//should only be one process in data
					trCats[tr_idx].procCats[p].hists1D[0][5]->Fill(gamtime, _weight);
					trCats[tr_idx].procCats[p].hists1D[0][2]->Fill(deltaT_gampv, _weight);
					//gampv resolution
					if(phos[0].pt() > 70){
						dphi_phoJets = deltaPhi(phos[0], jets);
						if(dphi_phoJets > pi-0.35 && dphi_phoJets < pi+0.35){
							Epho = 0;
							phos[0].GetJetPoints(phorhs);
							for(auto r : phorhs) Epho += r.E();
						//cout << "x " << sqrt(Epho*Ejets) << " y " << deltaT_gampv << endl;
							trCats[tr_idx].procCats[p].hists2D[0][16]->Fill(sqrt(Epho*Ejets) , deltaT_gampv, _weight);
						}
					}

					//do same for subleading photon if it exists
					if(phos.size() > 1){
						//want difference in PV frame (no shift to detector)
						gamtime = CalcJetTime(ts, phos[1], false);
						deltaT_gampv = gamtime - pvtime;
					cout << "data - gam time - pv time " << deltaT_gampv << endl;
					//cout << "gampv time " << deltaT_gampv << endl;
						trCats[tr_idx].procCats[p].hists1D[0][5]->Fill(gamtime, _weight);
						trCats[tr_idx].procCats[p].hists1D[0][2]->Fill(deltaT_gampv, _weight);
						//gampv resolution
						if(phos[1].pt() > 70){
							dphi_phoJets = deltaPhi(phos[1], jets);
							if(dphi_phoJets > pi-0.35 && dphi_phoJets < pi+0.35){
								Epho = 0;
								phos[1].GetJetPoints(phorhs);
								for(auto r : phorhs) Epho += r.E();
								trCats[tr_idx].procCats[p].hists2D[0][16]->Fill(sqrt(Epho*Ejets), deltaT_gampv, _weight);
							}
						}
					}			

				}
				else{
					//fill correct procCat
					vector<double> ids = trCats[tr_idx].procCats[p].ids;
					phoidx = phos[0].GetUserIdx();
					//cout << "leading phoidx " << phoidx << endl;
					genidx = _base->Photon_genIdx->at(phoidx);
					if(genidx == -1) phoid = -1; //gen match for MC
					else phoid = _base->Gen_susId->at(genidx);
					
					//cout << "leading phoid " << phoid << endl;
					//cout << (std::find(ids.begin(), ids.end(), phoid) != ids.end()) << " null id " <<  (std::find(ids.begin(), ids.end(), -999) != ids.end()) << endl;
					//make sure id is in current vector of ids (or ids does not contain -999)
					if(std::find(ids.begin(), ids.end(), phoid) != ids.end() || std::find(ids.begin(), ids.end(), -999) != ids.end()){
						//get sum of pho rh energy
						phos[0].GetJetPoints(phorhs);
						Epho = 0;
						for(auto r : phorhs) Epho += r.E();
						//cout << "LEAD CALC GAMTIME - tridx: " << tr_idx << " p " << p << " E " << Epho << endl;
						//want the photon time at the PV so we can compare it to the jet time 
						gamtime_pv = CalcJetTime(ts, phos[0], false);
						deltaT_gampv = gamtime_pv - pvtime;
						if(phos[0].pt() > 70) trCats[tr_idx].procCats[p].hists2D[0][16]->Fill(sqrt(Epho*Ejets) , deltaT_gampv, _weight);
						//deltaT of pho time - pv time (with both times in pv frame)
						trCats[tr_idx].procCats[p].hists1D[0][2]->Fill(deltaT_gampv, _weight);
						//cout << "ts " << ts << " gampv time " << deltaT_gampv << endl;
						//want the photon time at the detector so we can compare it to the gen version
						gamtime = CalcJetTime(ts, phos[0], true);
						//cout << "LEAD CALC GAMTIME END" << endl;
						trCats[tr_idx].procCats[p].hists1D[0][5]->Fill(gamtime_pv, _weight);
						deltaT_gampv = gamtime - pvtime;
					
						//fill difference in deltaT_pvGam of reco and gen - 3
						//cout << "calc gen delta t" << endl;
						deltaT_gampv_gen = CalcGenDeltaT(phos[0]);
					//cout << "LEAD tr idx: " << tr_idx << " pho id " << phoid << " gen deltaT: " << deltaT_gampv_gen << " reco deltaT: " << deltaT_gampv << " gamtime: " << gamtime << " pvtime: " << pvtime << " Epho: " << Epho << endl;
						//only for gen matches
						if(deltaT_gampv_gen != -999){
							trCats[tr_idx].procCats[p].hists1D[0][4]->Fill(deltaT_gampv_gen, _weight);
							//will be centered on zero - for resolution comparison with/without pv time need minimum cut on jet pt for pv time
							if(!(xbins[0] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[1])){
								deltaT_gampv = gamtime - pvtime;	
								trCats[tr_idx].procCats[p].hists1D[0][3]->Fill(deltaT_gampv - deltaT_gampv_gen, _weight);
							}
							//fill res (sigma from gaussian fit) for deltaT_recoGen as a function of ptAvg of jets that go into pv time calc - 4
							//sigma deltaT_recoGen as a function of geoEavg
							trCats[tr_idx].procCats[p].hists2D[0][0]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
							//if photon time (gamtime - pvtime) is in different delayed time windows, fill different 2D hists
							if(deltaT_gampv_gen >= 3.5 && deltaT_gampv_gen < 4.5)
								trCats[tr_idx].procCats[p].hists2D[0][4]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
							if(deltaT_gampv_gen >= 4.5 && deltaT_gampv_gen < 8)
								trCats[tr_idx].procCats[p].hists2D[0][5]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
							if(deltaT_gampv_gen >= 8 && deltaT_gampv_gen < 12)
								trCats[tr_idx].procCats[p].hists2D[0][6]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
							//fill gen deltaT vs reco deltaT in energy bins
							if(Epho >= 0 && Epho < 100){
								trCats[tr_idx].procCats[p].hists2D[0][8]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
								trCats[tr_idx].procCats[p].hists1D[0][17]->Fill(deltaT_gampv, _weight); 
							}
							if(Epho >= 100 && Epho < 400){
								trCats[tr_idx].procCats[p].hists2D[0][9]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
								trCats[tr_idx].procCats[p].hists1D[0][18]->Fill(deltaT_gampv, _weight); 
							}
							if(Epho >= 400 && Epho < 700){
								trCats[tr_idx].procCats[p].hists2D[0][10]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
								trCats[tr_idx].procCats[p].hists1D[0][19]->Fill(deltaT_gampv, _weight); 
							}
							if(Epho >= 700){
								trCats[tr_idx].procCats[p].hists2D[0][11]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
								trCats[tr_idx].procCats[p].hists1D[0][20]->Fill(deltaT_gampv, _weight); 
							}
							
							trCats[tr_idx].procCats[p].hists2D[0][12]->Fill(CalcGenDr(phos[0]), deltaT_gampv/deltaT_gampv_gen, _weight);
							trCats[tr_idx].procCats[p].hists2D[0][13]->Fill(GenEnergy(phos[0]), deltaT_gampv/deltaT_gampv_gen, _weight);
							trCats[tr_idx].procCats[p].hists2D[0][14]->Fill(Epho, deltaT_gampv_gen, _weight);
							trCats[tr_idx].procCats[p].hists2D[0][15]->Fill(Epho/GenEnergy(phos[0]), deltaT_gampv - deltaT_gampv_gen, _weight);

						}	
	
					
					}

					//do same for subleading photon if it exists
					if(phos.size() > 1){
						phoidx = phos[1].GetUserIdx();
						genidx = _base->Photon_genIdx->at(phoidx);
						if(genidx == -1) phoid = -1;
						else phoid = _base->Gen_susId->at(genidx);
						phorhs.clear();
					      //make sure id is in current vector of ids (or ids does not contain -999)
					      //cout << "!leading phoid " << phoid << " " << (std::find(ids.begin(), ids.end(), phoid) != ids.end()) << " null id " <<  (std::find(ids.begin(), ids.end(), -999) != ids.end()) << endl;
						if(std::find(ids.begin(), ids.end(), phoid) != ids.end() || std::find(ids.begin(), ids.end(), -999) != ids.end()){
							phos[1].GetJetPoints(phorhs);
							Epho = 0;
							for(auto r : phorhs) Epho += r.E();
							//cout << "SUBLEAD CALC GAMTIME - tridx: " << tr_idx << " p " << p << " Epho " << Epho << endl;
							//cout << "SUBLEAD GAMTIME" << endl;	
							gamtime_pv = CalcJetTime(ts, phos[1], false);
							deltaT_gampv = gamtime_pv - pvtime;
							if(phos[1].pt() > 70) trCats[tr_idx].procCats[p].hists2D[0][16]->Fill(sqrt(Epho*Ejets) , deltaT_gampv, _weight);
							trCats[tr_idx].procCats[p].hists1D[0][2]->Fill(deltaT_gampv, _weight);
							//cout << "gampv time " << deltaT_gampv << endl;
							//cout << "SUBLEAD GAMTIME END" << endl	
							gamtime = CalcJetTime(ts, phos[1], true);
							trCats[tr_idx].procCats[p].hists1D[0][5]->Fill(gamtime_pv, _weight);
							deltaT_gampv = gamtime - pvtime;
							
							deltaT_gampv_gen = CalcGenDeltaT(phos[1]);
					//cout << "SUBLEAD tr idx: " << tr_idx << " pho id " << phoid << " gen deltaT: " << deltaT_gampv_gen << " reco deltaT: " << deltaT_gampv << " gamtime: " << gamtime << " pvtime: " << pvtime << endl;
							trCats[tr_idx].procCats[p].hists1D[0][4]->Fill(deltaT_gampv_gen, _weight);
							//only for gen matches
							if(deltaT_gampv_gen != -999){
								if(!(xbins[0] <= sqrt(Erh1*Erh2) && sqrt(Erh1*Erh2) < xbins[1])){
									deltaT_gampv = gamtime - pvtime;	
									trCats[tr_idx].procCats[p].hists1D[0][3]->Fill(deltaT_gampv - deltaT_gampv_gen, _weight);
								}
								trCats[tr_idx].procCats[p].hists2D[0][0]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
								if(deltaT_gampv_gen >= 3.5 && deltaT_gampv_gen < 4.5)
									trCats[tr_idx].procCats[p].hists2D[0][4]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
								if(deltaT_gampv_gen >= 4.5 && deltaT_gampv_gen < 8)
									trCats[tr_idx].procCats[p].hists2D[0][5]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
								if(deltaT_gampv_gen >= 8 && deltaT_gampv_gen < 12)
									trCats[tr_idx].procCats[p].hists2D[0][6]->Fill(sqrt(Epho*Erh), deltaT_gampv - deltaT_gampv_gen, _weight);
								//fill gen deltaT vs reco deltaT in energy bins
								if(Epho >= 0 && Epho < 100){
									trCats[tr_idx].procCats[p].hists2D[0][8]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
									trCats[tr_idx].procCats[p].hists1D[0][17]->Fill(deltaT_gampv, _weight); 
								}
								if(Epho >= 100 && Epho < 400){
									trCats[tr_idx].procCats[p].hists2D[0][9]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
									trCats[tr_idx].procCats[p].hists1D[0][18]->Fill(deltaT_gampv, _weight); 
								}
								if(Epho >= 400 && Epho < 700){
									trCats[tr_idx].procCats[p].hists2D[0][10]->Fill(deltaT_gampv_gen, deltaT_gampv, _weight); 
									trCats[tr_idx].procCats[p].hists1D[0][19]->Fill(deltaT_gampv, _weight); 
								}
								if(Epho >= 700){
									trCats[tr_idx].procCats[p].hists2D[0][11]->Fill(deltaT_gampv_gen, deltaT_gampv); 
									trCats[tr_idx].procCats[p].hists1D[0][20]->Fill(deltaT_gampv, _weight); 
								}
								trCats[tr_idx].procCats[p].hists2D[0][12]->Fill(CalcGenDr(phos[1]), deltaT_gampv/deltaT_gampv_gen, _weight);
								trCats[tr_idx].procCats[p].hists2D[0][13]->Fill(GenEnergy(phos[1]), deltaT_gampv/deltaT_gampv_gen, _weight);
								trCats[tr_idx].procCats[p].hists2D[0][14]->Fill(Epho, deltaT_gampv_gen, _weight);
								trCats[tr_idx].procCats[p].hists2D[0][15]->Fill(Epho/GenEnergy(phos[1]), deltaT_gampv - deltaT_gampv_gen, _weight);
	
							}
						}
					}
				
				}
			//cout << "\n" << endl;
			}
			//cout << "\n" << endl;
		}

		double CalcAvgPt(const vector<Jet>& jets){
			double pt = 0;
			int njets = jets.size();
			for(int j = 0; j < njets; j++)
				pt += jets[j].pt();
			return pt/double(njets);
		}

		double GenEnergy(const Jet& pho){
			//calc times differently for !sig and sig photons
			int genidx, phoidx, phoid;
			//gen photon coordinates
			double geneta, genphi, phoeta, phophi;
			//if no match
			phoidx = pho.GetUserIdx();
			genidx = _base->Photon_genIdx->at(phoidx);
			if(genidx == -1) phoid = -1;
			else phoid = _base->Gen_susId->at(genidx);
			if(phoid == -1) return -999;
			return _base->Gen_energy->at(genidx);

		}


		double CalcGenDr(const Jet& pho){
			//calc times differently for !sig and sig photons
			int genidx, phoidx, phoid;
			//gen photon coordinates
			double geneta, genphi, phoeta, phophi;
			//if no match
			phoidx = pho.GetUserIdx();
			genidx = _base->Photon_genIdx->at(phoidx);
			if(genidx == -1) phoid = -1;
			else phoid = _base->Gen_susId->at(genidx);
			if(phoid == -1) return -999;
			//geneta = _base->Gen_eta->at(genidx);
			//genphi = _base->Gen_phi->at(genidx);

			//need to correct for displacement from 0 point
			double gvx = _base->Gen_vx->at(genidx);
                        double gvy = _base->Gen_vy->at(genidx);
                        double gvz = _base->Gen_vz->at(genidx);

                        //double rx = 129*cos(_base->Photon_phi->at(phoidx));
                        //double ry = 129*sin(_base->Photon_phi->at(phoidx));
                        //double rtheta = 2*atan2(1,exp(_base->Photon_eta->at(phoidx)));
                        //double rz = 129/tan(rtheta);

			int scidx = _base->Photon_scIndex->at(phoidx);
			double rx = _base->SuperCluster_x_calo->at(scidx);
			double ry = _base->SuperCluster_y_calo->at(scidx);
			double rz = _base->SuperCluster_z_calo->at(scidx);

			//phoeta = pho.eta();
			//phophi = pho.phi();
			//double deta = geneta - phoeta;
			//double dphi = genphi - phophi;

			double dx = rx - gvx;
                        double dy = ry - gvy;
                        double dz = rz - gvz;

			double reta = asinh(dz/sqrt(dx*dx + dy*dy));
                        double rphi = atan2(dy,dx);
		
			double deta = _base->Gen_eta->at(genidx) - reta;
			double dphi = _base->Gen_phi->at(genidx) - rphi;
		
			
			if(dphi > acos(-1)) dphi -= 2*acos(-1);
			if(dphi < -acos(-1)) dphi += 2*acos(-1);
			return sqrt( deta*deta + dphi*dphi );

		}

		//should be deltaT = gam - pv
		double CalcGenDeltaT(const Jet& pho){
		//cout << "CalcGenDeltaT - start" << endl;
			//calc times differently for !sig and sig photons
			double dpho = -999;
			int genidx, phoidx, phoid;
			//gen photon coordinates
			double genx, geny, genz, gentheta, geneta, genphi, vx, vy, vz;
			//if no match
			phoidx = pho.GetUserIdx();
			genidx = _base->Photon_genIdx->at(phoidx);
			//cout << "phoidx " << phoidx << " genidx " << genidx << endl;
			if(genidx == -1) phoid = -1;
			else phoid = _base->Gen_susId->at(genidx);
			if(phoid == -1) return dpho;

			//photon production vertex coordinates
			vx = _base->Gen_vx->at(genidx);
			vy = _base->Gen_vy->at(genidx);
			vz = _base->Gen_vz->at(genidx);		
			//momentum components
			double px = _base->Gen_px->at(genidx);
			double py = _base->Gen_py->at(genidx);
			double pz = _base->Gen_pz->at(genidx);		

			//energy
			double e = _base->Gen_energy->at(genidx);

			//find path length from photon production vertex to detector
			double R = 129; //radius of ECAL
			double vt = sqrt(vx*vx + vy*vy); //(l) distance to prod vertex in transverse plane from (0, 0)
			double pt = _base->Gen_pt->at(genidx); //momentum in transverse plane (for betaT)
			double dot = px*vx + py*vy; //pT \dot vT
			double L = sqrt(R*R - vt*vt + (dot/pt)*(dot/pt)) - dot/pt; //flight path from prod vertex to detector
			
			double betaT = pt/e; //beta in transverse plane (would be 1 for photons with z component of momentum)
			double tof = L/(_c*betaT); //time of flight to detector along calculated path length

			double genx_ECAL = vx + (px/e)*_c*tof;
			double geny_ECAL = vy + (py/e)*_c*tof;
			double genz_ECAL = vz + (pz/e)*_c*tof; 

			double ftheta = atan2(sqrt(genx_ECAL*genx_ECAL + geny_ECAL*geny_ECAL),genz_ECAL);
			double feta = -log(tan(ftheta/2.));
			//calculate TOF for that path (transverse plane only)
			//use TOF to propagate x, y, z from gen vertex position (vi) to detector (geni_ECAL)


			double pvx = _base->PV_x;
			double pvy = _base->PV_y;
			double pvz = _base->PV_z;

			double beta;
			if(phoid == 22){
				int momidx = _base->Photon_genSigXMomId->at(phoidx);
				//cout << "signal - phoid " << phoid << " phoidx " << phoidx << " genidx " << genidx << " momidx " << momidx << endl;
				//TODO: remove when ntuples are fixed (10/9/24)
				if(momidx < 0) return dpho;
				//check gen pdgids
				//want production vertex of photon (where LLP -> photon)
				//not production vertex of mother (close to PV)
				vx = _base->Gen_vx->at(genidx);
				vy = _base->Gen_vy->at(genidx);
				vz = _base->Gen_vz->at(genidx);
		
				//production vertex of mother particle (should be close to PV)
				double momvx, momvy, momvz;
				momvx = _base->Gen_vx->at(momidx);
				momvy = _base->Gen_vy->at(momidx);
				momvz = _base->Gen_vz->at(momidx);

				double mompx, mompy, mompz, momE;			

				mompx = _base->Gen_px->at(momidx);			
				mompy = _base->Gen_py->at(momidx);			
				mompz = _base->Gen_pz->at(momidx);			
				momE = _base->Gen_energy->at(momidx);			

				double p = sqrt(mompx*mompx + mompy*mompy + mompz*mompz);
				double m = _base->Gen_mass->at(momidx);
				double gam1 = sqrt(1 + (p/(m))*(p/(m)));
				double beta1 = sqrt(1 - 1/(gam1*gam1));

				//beta/c = p/E = v
				beta = sqrt(mompx*mompx + mompy*mompy + mompz*mompz)/momE;
				//cout << "beta: " << beta << " vel: " << beta*_c << endl;
				//check gen photon energy	
				//cout << "photon energy: " << _base->Photon_energy->at(phoidx) << endl;
				//distance bw photon and production point (where LLP decays to photon)
				dpho = sqrt( (genx_ECAL - vx)*(genx_ECAL - vx) + (geny_ECAL - vy)*(geny_ECAL - vy) + (genz_ECAL - vz)*(genz_ECAL - vz) )/_c;
		
				//distance bw LLP decay point and PV
				//LLP is produced close to PV (should take into account?)
				//LLPdecay - LLPprod?
				dpho += sqrt( (vx - momvx)*(vx - momvx) + (vy - momvy)*(vy - momvy) + (vz - momvz)*(vz - momvz) )/(_c*beta);	

			}
			//assume prompt production
			else{
				//cout << "prompt" << endl;
				//distance bw photon and production point (PV)
				dpho = sqrt( (genx_ECAL - pvx)*(genx_ECAL - pvx) + (geny_ECAL - pvy)*(geny_ECAL - pvy) + (genz_ECAL - pvz)*(genz_ECAL - pvz) )/_c;
			}
			
		//cout << "CalcGenDeltaT - end" << endl;
			return dpho;
		}

		//dphi bw photon and jet system
		double deltaPhi(const Jet& pho, const vector<Jet>& jets){
			double r = 129;
			double jx = 0;
			double jy = 0;
			//vector sum of jets
			for(auto j : jets){
				jx += r*cos(j.phi());
				jy += r*sin(j.phi());
			}
			double jphi = atan2(jx,jy);
			double dphi = fabs(pho.phi() - jphi);
			double pi = acos(-1);
			if(dphi > pi) return 2*pi - dphi;
			else return dphi;	
			
		}
		
		int FindJetPair(vector<Jet>& injets, pair<Jet,Jet>& outjets){
			//sort
			sort(injets.begin(), injets.end(), ptsort);
			Jet jet1 = injets[0];
			Jet jet2 = injets[1];	
		
			//if needing to find "true" pair
			//calculate dphi here
			double pi = acos(-1);
			double ptasym = 0.4;
			//have higher min pt for dijet selection
			double minpt = 100;
			if(jet1.pt() < minpt || jet2.pt() < minpt) return -999;


				
			//not needed for GMSB
			if(_data){
				//dphi within [pi-0.1,pi+0.1]
				double phi1 = jet1.phi_02pi();
				double phi2 = jet2.phi_02pi();
				double dphi = fabs(phi1 - phi2);
				if(dphi > pi) dphi = 2*pi - dphi;
				//cout << "phi1 = " << phi1 << " phi2 = " << phi2 << " dphi = " << dphi << endl;
				//cout << "pt2 = " << jet2.pt() << " pt1 = " << jet1.pt() << " pt2/pt1 = " << jet2.pt()/jet1.pt() << endl;
				if(dphi < pi-0.35 || dphi > pi+0.35){
					//cout << "failed dphi cut" << endl;
					return -999;
				}
				//pt asymmetry cut
				if(jet2.pt() / jet1.pt() < 1 - ptasym){
					//cout << "failed pt cut with thresh " << 1 - ptasym << endl;
					return -999;
				}
			}
			else{ //only do for bkg MC
				if(_oname.find("SMS") == string::npos){
					//make MC selection as similar to data as possible
					//dphi within [pi-0.1,pi+0.1]
					double phi1 = jet1.phi_02pi();
					double phi2 = jet2.phi_02pi();
					double dphi = fabs(phi1 - phi2);
					if(dphi > pi) dphi = 2*pi - dphi;
					//cout << "phi1 = " << phi1 << " phi2 = " << phi2 << " dphi = " << dphi << endl;
					//cout << "pt2 = " << jet2.pt() << " pt1 = " << jet1.pt() << " pt2/pt1 = " << jet2.pt()/jet1.pt() << endl;
					if(dphi < pi-0.35 || dphi > pi+0.35){
						//cout << "failed dphi cut" << endl;
						return -999;
					}
					//use only ECAL energy for asymmetry cut
					double e1 = 0;
					double e2 = 0;
					vector<JetPoint> rhs;
					jet1.GetJetPoints(rhs);
					for(auto r : rhs) e1 += r.e();
					jet2.GetJetPoints(rhs);
					for(auto r : rhs) e2 += r.e();
					if(e2 < e1){
						if(e2 / e1 < 1 - ptasym) return -999; 
					}
					else{
						if(e1 / e2 < 1 - ptasym) return -999; 
					}
				}
			}

			outjets.first = jet1;
			outjets.second = jet2;
			return 0;
		}
		double GetDeltaTime(pair<Jet,Jet>& jets){
			double t1 = jets.first.time();
			double t2 = jets.second.time();
			return t1 - t2;
		}

		//if pho, add TOF back to time (shift to detector)	
		double CalcJetTime(const TimeStrategy& ts, Jet& jet, bool pho = false){
			//cout << "CalcJetTime method " << ts << endl;
			double time = -999;
			vector<JetPoint> rhs; jet.GetJetPoints(rhs);
			if(ts == med) time = CalcMedianTime(jet);
			else if(ts == eavg) time = CalcEAvgTime(jet);
			else if(ts == mmavg){
				time = CalcMMAvgTime(jet, pho);
				//cout << "mmavg time " << time << " with " << jet.GetNConstituents() << " subcls" << endl;
			}
			else if(ts == emax) time = CalcMaxTime(jet);
			else cout << "Error: invalid time reconstruction method specified for calculating jet time" << endl;
			//if photon, shift time to detector by centroid as defined above
			if(pho){
				BayesPoint center;
				if(ts == med) center = CalcMedianCenter(jet);
				else if(ts == eavg) center = CalcEAvgCenter(jet);
				else if(ts == mmavg){
					center = CalcEAvgCenter(jet); //assume a BHC jet has been passed with appropriate weights on rh energies
				}
				else if(ts == emax) center = CalcMaxCenter(jet);
				double rtheta = atan2(129,center.at(2));
				if(center.at(2) < 0) rtheta = atan2(129,-center.at(2))+acos(-1)/2.;
				double reta = -log(tan(rtheta/2));
				double rphi = atan2(center.at(1),center.at(0));
				if(rphi < 0) rphi += 2*acos(-1);
				//jet.Print();
				//cout << "PV frame time: " << time << endl;
				//else return time;
				BayesPoint pv = jet.GetVertex();
				double dx = center.at(0) - pv.at(0);
				double dy = center.at(1) - pv.at(1);
				double dz = center.at(2) - pv.at(2);
				double t_shift = sqrt(dx*dx + dy*dy + dz*dz)/_c;
				//cout << "center for " << ts << " eta " << reta << " phi " << rphi << " x " << center.at(0) << " y " << center.at(1) << " z " << center.at(2) << " time " << time << " shift " << t_shift << endl;
				//shift from PV to detector face
				time += t_shift;
			}
			//cout << "return time " << time << " for jet with energy " << jet.E() << endl;	
			return time;

		}
		double CalcPVTime(const TimeStrategy& ts, vector<Jet>& jets){
			int njets = jets.size();
			double time = -999;
			if(njets < 1) return time;
			if(ts == med){
				vector<double> times;
				for(int i = 0; i < njets; i++)
					times.push_back(jets[i].t());
				//even - return average of two median times
				if(njets % 2 == 0)
					time = (times[int(double(njets)/2.)] + times[int(double(njets)/2.)-1])/2.;
				//odd - return median
				else
					time = times[int(double(njets)/2.)];

			}
			//ECAL energy weighted avg
			else if(ts == eavg || ts == emax){
				double norm = 0;
				double t = 0;
				double e = 0;
				vector<JetPoint> rhs;
				for(int j = 0; j < njets; j++){
					jets[j].GetJetPoints(rhs);
					for(int r = 0; r < rhs.size(); r++)
						e += rhs[r].E();
					t += e*jets[j].t();
					norm += e;
					e = 0;
					rhs.clear();
				}
				//cout << "pv time: " << t/norm << endl;
				time = t/norm;
			}	
			else if(ts == mmavg){
				double t = 0;
				for(int j = 0; j < njets; j++){
					//cout << "jet time (pv frame) " << jets[j].t() << endl;
					t += jets[j].t();
				}
				time = t/double(njets);
			}	
			return time;
		}

		double CalcMaxTime(Jet& j){
			map<double, double> eneTime;
			vector<JetPoint> rhs;
			j.GetJetPoints(rhs);
			for(int i = 0; i < rhs.size(); i++)
				eneTime[rhs[i].E()] = rhs[i].t();
			return eneTime.rbegin()->second;
		}

		BayesPoint CalcMaxCenter(Jet& j){
			map<double, BayesPoint> eneLoc;
			vector<JetPoint> rhs; j.GetJetPoints(rhs);
			for(int i = 0; i < rhs.size(); i++){
				BayesPoint pt({rhs[i].x(), rhs[i].y(), rhs[i].z()});
				//cout << "x " << rhs[i].x() << " y " << rhs[i].y() << " z " << rhs[i].z() << " eta " << rhs[i].eta() << " phi " << rhs[i].phi() << " phi02pi " << rhs[i].phi_02pi() << " rphi " << atan2(rhs[i].y(), rhs[i].x()) << endl;
				eneLoc[rhs[i].E()] = pt;
			}
			return eneLoc.rbegin()->second;
		}

		double CalcMedianTime(Jet& j){
			vector<double> times;
			double time = -999;
			vector<JetPoint> rhs; j.GetJetPoints(rhs);
			int nrhs = rhs.size();
			for(int i = 0; i < nrhs; i++)
				times.push_back(rhs[i].t());
			//make sure times are sorted
			sort(times.begin(), times.end());
			//even - return average of two median times
			if(nrhs % 2 == 0)
				time = (times[int(double(nrhs)/2.)] + times[int(double(nrhs)/2.)-1])/2.;
			//odd - return median
			else
				time = times[int(double(nrhs)/2.)];
			return time;
		}		

		BayesPoint CalcMedianCenter(Jet& j){
			map<double,BayesPoint> timeLoc;
			BayesPoint center(3);
			double time = -999;
			vector<double> times;
			vector<JetPoint> rhs; j.GetJetPoints(rhs);
			int nrhs = rhs.size();
			BayesPoint pt(3);
			//cout << "MEDIAN CENTER" << endl;
			//cout << nrhs << " rhs" << endl;
			for(int i = 0; i < nrhs; i++){
				pt.SetValue(rhs[i].x(),0);
				pt.SetValue(rhs[i].y(),1);
				pt.SetValue(rhs[i].z(),2);
				//cout << "i " << i << " x " << rhs[i].x() << " y " << rhs[i].y() << " z " << rhs[i].z() << " eta " << rhs[i].eta() << " phi " << rhs[i].phi() << " phi02pi " << rhs[i].phi_02pi() << " rphi " << atan2(rhs[i].y(), rhs[i].x()) << endl;
				timeLoc[rhs[i].t()] = pt;
				//cout << "i " << i << " x " << rhs[i].x() << " y " << rhs[i].y() << " z " << rhs[i].z() << " t " << rhs[i].t() << endl;
			}
			//cout << "# times " << timeLoc.size() << endl;
			for(auto it = timeLoc.begin(); it != timeLoc.end(); it++){
				//cout << "time " << it->first << " size " << times.size() << endl;
				times.push_back(it->first);
			}
			//even - return average of two median times
			if(nrhs % 2 == 0){
				BayesPoint center(3);
				//cout << "idx1 " << int(double(nrhs)/2.) << " idx2 " << int(double(nrhs)/2.)-1 << endl;
				double t1 = times[int(double(nrhs)/2.)];
				double t2 = times[int(double(nrhs)/2.)-1];
				//cout << "t1 " << t1 << " t2 " << t2 << endl;
				center.SetValue((timeLoc[t1].at(0) + timeLoc[t2].at(0))/2.,0);
				center.SetValue((timeLoc[t1].at(1) + timeLoc[t2].at(1))/2.,1);
				center.SetValue((timeLoc[t1].at(2) + timeLoc[t2].at(2))/2.,2);
				//cout << "even" << endl;
				//timeLoc[t1].Print(); timeLoc[t2].Print();				

				return center;	
			}
			//odd - return median
			else{
				time = times[int(double(nrhs)/2.)];
				//cout << "odd" << endl; timeLoc[time].Print();
				return timeLoc[time];
			}
		}		



		BayesPoint CalcEAvgCenter(Jet& j){
			double norm = 0;
			double x = 0;
			double y = 0;
			double z = 0;
			BayesPoint center(3);
			vector<JetPoint> rhs; j.GetJetPoints(rhs);
			int nrhs = rhs.size();
			double minpt = 10;
			for(int i = 0; i < nrhs; i++){
				if(rhs[i].E() < minpt) continue;
				norm += rhs[i].E();
				x += rhs[i].E()*rhs[i].x();
				y += rhs[i].E()*rhs[i].y();
				z += rhs[i].E()*rhs[i].z();
			}
			//cout << "CalcEAvgCenter - eta " << j.eta() << " phi " << j.phi() << endl;
			center.SetValue(x/norm,0);
			center.SetValue(y/norm,1);
			center.SetValue(z/norm,2);
			return center;
		}
		
		double CalcEAvgTime(Jet& j){
			double t = 0;
			double norm = 0;
			vector<JetPoint> rhs; j.GetJetPoints(rhs);
			int nrhs = rhs.size();
			double minpt = 10;
			for(int i = 0; i < nrhs; i++){
				if(rhs[i].E() < minpt) continue;
				norm += rhs[i].E();
				t += rhs[i].E()*rhs[i].t();
			//cout << "time " << rhs[i].t() << " energy " << rhs[i].E() << endl;
			}
			//cout << "return " << t/norm << endl;
			return t/norm;
		}



		GaussianMixture* _subcluster(Jet& jet){
			vector<JetPoint> rhs; jet.GetJetPoints(rhs);
			vector<Jet> rhs_jet;
			for(int r = 0; r < rhs.size(); r++){
				rhs_jet.push_back( Jet(rhs[r], jet.GetVertex()) );
			}
			BayesCluster* algo = new BayesCluster(rhs_jet);
			if(_smear && !_smearMat.empty()) algo->SetDataSmear(_smearMat);
			//set time resolution smearing
			//if(_timesmear) algo->SetTimeResSmear(_tres_c, _tres_n*_gev);
			algo->SetMeasErrParams(_cell, _tresCte, _tresStoch*_gev, _tresNoise*_gev); 
			algo->SetThresh(_thresh);
			//algo->SetAlpha(_alpha); //isn't used should probably remove this line
			algo->SetSubclusterAlpha(_emAlpha);
			algo->SetPriorParameters(_prior_params);
			algo->SetVerbosity(0);
			GaussianMixture* gmm = algo->SubCluster();
			//set constituents
			vector<double> norms;
			gmm->GetNorms(norms);
			for(int k = 0; k < gmm->GetNClusters(); k++){
				map<string, Matrix> params;
				gmm->GetLHPosteriorParameters(k, params);
				double E_k = norms[k]/_gev;
				Matrix mu = params["mean"];
				Matrix cov = params["cov"];
				double pi = params["pi"].at(0,0);
				//cout << "cluster k " << k << " pi " << pi << " center" << endl;
				//mu.Print();
				Jet subcl(mu,cov,E_k,pi,jet.GetVertex());
				jet.AddConstituent(subcl);
			}
			return gmm;
		}

		//returns {class, max discr val}
		//pair<int,double> PredictSubcluster(BasePDFMixture* model, int k, const Jet& jet){
		pair<int,double> PredictSubcluster(const Jet& jet){
			//nclass options (index of predicted class)
			//0 == 1 == phys bkg
			//1 == 2 == BH
			//2 == 3 == spike
			vector<JetPoint> rhs; jet.GetJetPoints(rhs);
			vector<double> ovalues; //discriminator output value, pass-by-ref
			map<string, double> obs; //k maps for k subclusters
			//features = vector of strings of features to use - set with SetNNFeatures
			//obs = map<string, double> observations for each feature per subcluster
			//makeDNNinputs too (if using later) 
			//TODO - make sure bhc jet itself is PU cleaned (ie weights applied to rhs)
			//can apply CNN weights as an additional weight to rhs in SCs matched to jet
			MakeCNNInputGrid(rhs, obs);
			CNNPredict(obs,ovalues);
			//TODO - update with discriminator score cut
			auto max_el = max_element(ovalues.begin(), ovalues.end());
			double ovalue = *max_el;
			int nclass = std::distance(ovalues.begin(), max_el);
			return make_pair(nclass, ovalue);
		}	



		//removes/downweights subclusters that are tagged (probID) as detector bkg by MVA
		//saves new prob-weight in jet
		//if remove, probID = 0 for detector bkg
		//NN json needs to have been set with SetNNModel(string jsonname)
		void CleanSubclusters(BasePDFMixture* model, Jet& jet){
			int kmax = model->GetNClusters();
			pair<int, double> class_discr;
			int nclass; double score;
			//cout << kmax << " # subclusters " << endl;
			Jet subcl;
			for(int k = 0; k < kmax; k++){
				//cout << "k " << k << " of " << kmax;
				jet.GetConstituent(k, subcl);
				class_discr = PredictSubcluster(subcl);
				nclass = class_discr.first;
				score = class_discr.second;

				double e = subcl.E();
				double t = subcl.t();

				//downweight - pushback prob == physbkg value
				//remove - if nclass != phys bkg, probID = 0
				//remove for now bc need to custom hack frugally-deep to return output of all neurons in the last layer
				//if(nclass == 0 && score < 0.8){ //need a minimial phys bkg score to stay in jet
				//if(nclass != 0 && score > 0.99){ //need a maximal bh/spike score to get cleaned
				//	jet.GetConstituent(k).SetCoefficient(0.); //pi*physbkg_val
				//	//Jet subcl = jet.GetConstituent(k);
				//	//subcl.SetCoefficient(0.); //pi*physbkg_val
				//}
				//fill dijet subcl hists
				if(nclass == 0){ //phys bkg
					trCats[0].procCats[0].hists2D[0][54]->Fill(t, e, _weight);
					trCats[0].procCats[0].hists2D[0][57]->Fill(class_discr.second, e, _weight);
					trCats[0].procCats[0].hists2D[0][60]->Fill(class_discr.second, t, _weight);
					
					trCats[0].procCats[1].hists2D[0][54]->Fill(t, e, _weight);
					trCats[0].procCats[1].hists2D[0][57]->Fill(class_discr.second, e, _weight);
					trCats[0].procCats[1].hists2D[0][60]->Fill(class_discr.second, t, _weight);
				}	
				else if(nclass == 1 && score > 0.99){ //bh
					trCats[0].procCats[0].hists2D[0][53]->Fill(t, e, _weight);
					trCats[0].procCats[0].hists2D[0][56]->Fill(class_discr.second, e, _weight);
					trCats[0].procCats[0].hists2D[0][59]->Fill(class_discr.second, t, _weight);
					
					trCats[0].procCats[1].hists2D[0][53]->Fill(t, e, _weight);
					trCats[0].procCats[1].hists2D[0][56]->Fill(class_discr.second, e, _weight);
					trCats[0].procCats[1].hists2D[0][59]->Fill(class_discr.second, t, _weight);
				}	
				else if(nclass == 2 && score > 0.99){ //spikes
					trCats[0].procCats[0].hists2D[0][52]->Fill(t, e, _weight);
					trCats[0].procCats[0].hists2D[0][55]->Fill(class_discr.second, e, _weight);
					trCats[0].procCats[0].hists2D[0][58]->Fill(class_discr.second, t, _weight);
					
					trCats[0].procCats[1].hists2D[0][52]->Fill(t, e, _weight);
					trCats[0].procCats[1].hists2D[0][55]->Fill(class_discr.second, e, _weight);
					trCats[0].procCats[1].hists2D[0][58]->Fill(class_discr.second, t, _weight);
				}	
			}
		}



		double CalcMMAvgTime(Jet& jet, bool pho){
			int kmax = jet.GetNConstituents();
			//cout << "CalcMMAvgTime - jet has " << kmax << " subcls" << endl;
			double t = 0;
			double norm = 0;
			double pik, tk;
			Matrix muk, covk;
			
			return jet.t();
			/*	
			for(int k = 0; k < kmax; k++){
				subcl = jet.GetConstituent(k);
				pik = subcl.GetCoefficient(); 
				subcl.GetClusterParams(muk,covk);
				tk = pik*muk.at(2,0);
				t += tk;
				norm += pik;
				//cout << "k: " << k << " pi: " << params["pi"].at(0,0) << " mean " << params["mean"].at(2,0) << endl;
				if(pho){
					//cout << "pho mm time: " << t/norm << endl;
					return t/norm;
				}
				//cout << "cluster " << k << " has time " << params["mean"].at(2,0) << " and MM coeff " << params["pi"].at(0,0) << endl;
			}
			//cout << "jet mm time: " << t/norm << endl;
			return t/norm; 
			*/
		}
		
		BayesPoint CalcMMAvgCenter(Jet& jet, bool pho, vector<double> probIDs = {}){
			int kmax = jet.GetNConstituents();
			double norm = 0;
			double x = 0;
			double y = 0;
			double z = 0;
			double phi, eta, theta;
			BayesPoint center(3);
			map<string, Matrix> params;
			double xk, yk, zk, pik;
			Matrix muk, covk;
			Jet subcl;	
			for(int k = 0; k < kmax; k++){
				jet.GetConstituent(k, subcl);
				pik = subcl.GetCoefficient();
				subcl.GetClusterParams(muk,covk);
				eta = muk.at(0,0);
				phi = muk.at(1,0);
				theta = 2*atan2(1,exp(eta));
				//detector y (in y,z plane) or r (in x,y plane) is 1.29 m
				xk = pik*129*cos(phi);
				yk = pik*129*sin(phi);
				//if(theta > acos(-1)/2.) z = params["pi"].at(0,0)*-129/tan(theta - acos(-1)/2.); 
				//else z += params["pi"].at(0,0)*129/tan(theta);	
				zk = pik*129/tan(theta);
				//clean out det bkg
				if(probIDs.size() == kmax){
					xk *= probIDs[k];
					yk *= probIDs[k];
					zk *= probIDs[k];
				}
				x += xk;
				y += yk;
				z += zk;
				//for checking calculations
				double rtheta = atan2(129,z);
				//if(z < 0) rtheta = atan2(129.,z)+acos(-1)/2.;
				double reta = -log(tan(rtheta/2));
				double rphi = atan2(y,x);
				norm += pik;
				if(pho){
		//cout << "MMCenter - pi: " << params["pi"].at(0,0) << " phi: " << phi << " eta: " << eta << " x: " << x/norm << " y: " << y/norm << " z: " << z/norm << " reta: " << reta << " rphi: " << rphi << " theta: " << theta << " rtheta: " << rtheta << endl;	
					center.SetValue(xk/norm,0);
					center.SetValue(yk/norm,1);
					center.SetValue(zk/norm,2);
					return center;
				}
			}
			center.SetValue(x/norm,0);
			center.SetValue(y/norm,1);
			center.SetValue(z/norm,2);
			return center;
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
				for(int k = 1; k < nprofs+1; k++){
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
					if(trs[j].procCats[0].hists1D[0][i] == nullptr) continue;
					histname = trs[j].procCats[0].hists1D[0][i]->GetName();
					//write total method histogram outside process directory
					if(trs[j].procCats[0].hists1D[0][i]->GetEntries() == 0 && ((histname.find("sigma") == string::npos && histname.find("mean") == string::npos) && histname.find("profile") == string::npos)){ continue; }
					//cout << "writing " << trs[j].procCats[0].hists1D[0][i]->GetName() << " " << trs[j].procCats[0].hists1D[0][i]->GetTitle() << " to " << dir->GetName() << endl;
					trs[j].procCats[0].hists1D[0][i]->Write();
					//make process breakdown directory
					TDirectory *dir2 = dir->mkdir((dirname+"_"+trs[j].methodName+"_procStack").c_str());
					//cout << "  making dir " << dir2->GetName() << endl;
					dir2->cd();
					for(int p = 1; p < trs[j].procCats.size(); p++){
					//loop over processes
			//			cout << "    proc " << trs[j].procCats[p].plotName << " hist " << trs[j].procCats[p].hists1D[0][i]->GetName() << " " << trs[j].procCats[p].hists1D[0][i]->GetTitle() << " entries " << trs[j].procCats[p].hists1D[0][i]->GetEntries() << endl;			
						histname = trs[j].procCats[p].hists1D[0][i]->GetName();
						if(trs[j].procCats[p].hists1D[0][i] == nullptr) continue;
						if(trs[j].procCats[p].hists1D[0][i]->GetEntries() == 0 && ((histname.find("sigma") == string::npos && histname.find("mean") == string::npos) && histname.find("profile") == string::npos)){ continue; }
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
					if(trs[j].procCats[0].hists2D[0][i]->GetEntries() == 0 && (histname.find("sigma") == string::npos && histname.find("profile") == string::npos && histname.find("meanRecoGenDeltaT") == string::npos)){ continue; }
					//check if data can be run
					if(histname.find("recoGen") != string::npos && _data) continue;
					//cout << "writing " << trs[j].procCats[0].hists1D[0][i]->GetName() << " " << trs[j].procCats[0].hists1D[0][i]->GetTitle() << " to " << dir->GetName() << endl;
					trs[j].procCats[0].hists2D[0][i]->Write();
					//write method as directory within directory
					TDirectory *dir2 = dir->mkdir((dirname+"_"+trs[j].methodName+"_procStack").c_str());
					//cout << "  making dir " << dir2->GetName() << endl;
					dir2->cd();
					for(int p = 1; p < trs[j].procCats.size(); p++){
					//loop over processes
						//if(dirname.find("meanRecoGenDeltaT") != string::npos) cout << "    proc " << trs[j].procCats[p].plotName << " hist " << trs[j].procCats[p].hists2D[0][i]->GetName() << " " << trs[j].procCats[p].hists2D[0][i]->GetTitle() << " entries " << trs[j].procCats[p].hists2D[0][i]->GetEntries() << endl;			
						if(trs[j].procCats[p].hists2D[0][i] == nullptr) continue;
						if(trs[j].procCats[p].hists2D[0][i]->GetEntries() == 0 && dirname.find("meanRecoGenDeltaT") == string::npos){ continue; }//cout << "Histogram for proc " << trs[j].plotName << " not filled." << endl; continue; }
						//check if data can be run
						histname = trs[j].procCats[p].hists2D[0][i]->GetName();
						if(histname.find("recoGen") != string::npos && _data) continue;
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
			for(int i = 0; i < trCats.size(); i++)
				WriteEmptyProfiles(ofile, trCats[i]);
			WriteTimeRecoCatStack(ofile, trCats);

			ofile->Close();


		}

		//void TrackMatched(BasePDFMixture* model, int subcl, double& bestdr, double& bestp){
		void TrackMatched(JetPoint& rh, double& bestdr, double& bestp){
			//do track matching
			double bestTrackDr = 999;
			//double maxTrackDr;
			double dr, teta, tphi, de;
			unsigned int detid;

			double dphi = -999;
			//double bestde_dr;
			int trackidx;
			double ec = rh.eta();
			double pc = rh.phi();

			//loop through tracks to get best match to this subcluster (tracks are matched to superclusters, if not idx < 0)
			//Track_scIndexs[i][j] is for track i that matched to supercluster j (tracks can match to multiple SCs)
			int nTracks = _base->Track_scIndexs->size();
			for(int t = 0; t < nTracks; t++){
				if(_base->Track_scIndexs->at(t).at(0) < 0) continue; //not matched to any SC
				
				int nSCs = _base->Track_scIndexs->at(t).size();
				for(int sc = 0; sc < nSCs; sc++){
					int sc_idx = _base->Track_scIndexs->at(t).at(sc);
					//get eta, phi of supercluster sc that track is matched to
					double sc_eta = _base->SuperCluster_eta->at(sc_idx);
					double sc_phi = _base->SuperCluster_phi->at(sc_idx);
					double sc_phi_02pi = sc_phi;
					if(sc_phi_02pi < 0) sc_phi_02pi += 2*acos(-1);
					else if(sc_phi_02pi > 2*acos(-1)) sc_phi_02pi -= 2*acos(-1); 
					else sc_phi_02pi = sc_phi;

					dphi = fabs(pc - sc_phi_02pi);
                                	dphi = acos(cos(dphi));

					dr = sqrt((sc_eta - ec)*(sc_eta - ec) + dphi*dphi);

	
					//E = p for photons
					if(dr < bestTrackDr){
						bestTrackDr = dr;
						bestp = _base->Track_p->at(t);
					}
					
				}
				bestdr = bestTrackDr;
			}

		}


		
		void SetMinRhE_PV(double m){ _minRhE = m; }	
		void SetCleanSubclusters(bool s){ _cleansubcls = s; }	
		void SetMistClean(bool m){ _prod->SetMistClean(m); }
		private:
			vector<Jet> _SCs; //superclusters for event - for # subcluster matching
			double _minRhE;
			Matrix _smearMat;
			PhotonProducer* _scprod;
			bool _cleansubcls;
};
#endif
