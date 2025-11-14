#ifndef SUPERCLUSTERSKIMMER_HH
#define SUPERCLUSTERSKIMMER_HH

#include "JetPoint.hh"
#include <TFile.h>
#include "BaseSkimmer.hh"
#include "PhotonProducer.hh"
#include "TSystem.h"
#include <math.h>
#include <fstream>

using procCat = BaseSkimmer::procCat;
class SuperClusterSkimmer : public BaseSkimmer{
	public:
		SuperClusterSkimmer(){
			SetObs();
			InitMapTree();
			_evti = 0;
			_evtj = 0;
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_gev = 1/30.;
			_applyFrac = false;
		};
		virtual ~SuperClusterSkimmer(){ };

		//get rechits from file to cluster
		//SuperClusterSkimmer(TFile* file) : BaseSkimmer(file){
		SuperClusterSkimmer(string filename) : BaseSkimmer(filename){
			SetObs();
			InitMapTree();
			
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time
			//this is just the type of producer, there is a GetSuperCluster fcn in the base producer class
				
			_prod = new PhotonProducer(_ch);
			_prod->SetTimeCalibrationTool(_timecalib);	
			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/photon_skims_"+_cms_label+".root";
			
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_gev = 1/30.;
			_applyFrac = false;

			SetupDetIDs();
		}
		/*
		SuperClusterSkimmer(string filelist) : BaseSkimmer(filelist){
			SetObs();
			InitMapTree();
			//this is just the type of producer, there is a GetSuperCluster fcn in the base producer class
			TChain* ch = MakeTChain(filelist);
			if(ch == nullptr) return;
			_prod = new PhotonProducer(ch);
			_prod->SetTimeCalibrationTool(_timecalib);	
			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/sc_skims_"+_cms_label+".root";
			
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 0.1;
			_emAlpha = 0.5;
			_gev = 1/30.;
			_applyFrac = false;

			
			SetupDetIDs();
		}
		*/

		void SuperClusterAddBranches(){
			_obj = "SC";
			_obsnames.push_back("EovP_track");
			_obsnames.push_back("dR_track");
			_obsnames.push_back("trueLabel"); //CR designation
			_obsnames.push_back("predLabel"); //CNN prediction
			_obsnames.push_back("predScore_physBkg"); //CNN prediction
			_obsnames.push_back("predScore_BH"); //CNN prediction
			_obsnames.push_back("predScore_spike"); //CNN prediction
			_obsnames.push_back("swCrossPrime");
			_obsnames.push_back("swCrossCMS");

			vAddBranch("nRHs_grid","# rhs in CNN input grid");
			vAddBranch("nRHs_ogSC","# of rhs in original SC");
			vAddBranch("nRHs_bhcSC","# of rhs in BHC SC");
			vAddBranch("seedTimeSignificance","time significance from seed rh in original SC");
			vAddBranch("weightedTimeSignificance","time significance with inverse variance weights");
			vAddBranch("seedTime","time from seed rh in original SC");

			//rh obs
			vvAddBranch("rh_iEta","ieta coord of rh in 7x7 CNN grid");
			vvAddBranch("rh_iPhi","iphi coord of rh in 7x7 CNN grid");
			vvAddBranch("rh_Energy","energy of rh in 7x7 CNN grid");
			vvAddBranch("rh_Weight","weight of rh in 7x7 CNN grid");
		}

		//302 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, phys bkg 
		TH2D* ENeighbors_physBkg = new TH2D("ENeighbors_physBkg","ENeighbors_physBkg;local ieta;local iphi_physBkg",5,-2,3,5,-2,3);	
		//303 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, BH 
		TH2D* ENeighbors_BH = new TH2D("ENeighbors_BH","ENeighbors_BH;local ieta;local iphi_BH",5,-2,3,5,-2,3);	
		//304 - eta vs phi grid, overlaid neighbors energy in 5x5 grid, spike 
		TH2D* ENeighbors_spikes = new TH2D("ENeighbors_spikes","ENeighbors_spikes;local ieta;local iphi_spikes",5,-2,3,5,-2,3);	
		//305 - eta vs phi grid no center, overlaid neighbors energy in 5x5 grid, phys bkg 
		TH2D* ENeighborsSkipCenter_physBkg = new TH2D("ENeighborsSkipCenter_physBkg","ENeighborsSkipCenter_physBkg;local ieta;local iphi_physBkg",5,-2,3,5,-2,3);	
		//306 - eta vs phi grid no center, overlaid neighbors energy in 5x5 grid, BH 
		TH2D* ENeighborsSkipCenter_BH = new TH2D("ENeighborsSkipCenter_BH","ENeighborsSkipCenter_BH;local ieta;local iphi_BH",5,-2,3,5,-2,3);	
		//307 - eta vs phi grid no center, overlaid neighbors energy in 5x5 grid, spike 
		TH2D* ENeighborsSkipCenter_spikes = new TH2D("ENeighborsSkipCenter_spikes","ENeighborsSkipCenter_spikes;local ieta;local iphi_spikes",5,-2,3,5,-2,3);	
		//308 - eta vs phi grid
		TH2D* ENeighbors = new TH2D("ENeighbors","ENeighbors;local ieta;local iphi",61,-30,31,61,-30,31);

		enum weightScheme{
			noWeight = 0,
			Eweight = 1,
			logEweight = 2
		};
		

	

			
		void Skim();

		void AddSample(TFile* file);

		void SetIsoCuts(){ _isocuts = true; }
		bool _isocuts;
		//set skip for outstream
		void SetSkip(int i){ _oskip = i; }
		int _oskip;
		double _timeoffset;
		bool _BHcluster;
		int _clusterSize;
		void ApplyFractions(bool a){ _applyFrac = a; if(_applyFrac) cout << "Applying RH fractions" << endl; }
		bool _applyFrac;


		void WritePlotCat1D(TFile* ofile, const procCat& pc){
			ofile->cd();
			string name;
			vector<vector<TH1D*>> hists1D = pc.hists1D;
			//write 1D hists
			for(int i = 0; i < (int)hists1D.size(); i++){
				for(int j = 0; j < hists1D[i].size(); j++){
					if(pc.hists1D[i][j] == nullptr) continue;
					name = pc.hists1D[i][j]->GetName();
					if(pc.hists1D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram: " << name << " not filled." << endl; continue; }
					pc.hists1D[i][j]->Write();
				}
			}
		}
		void WritePlotCat2D(TFile* ofile, const procCat& pc){
			ofile->cd();
			string name;
			vector<vector<TH2D*>> hists2D = pc.hists2D;
			
			//write 2D hists
			for(int i = 0; i < (int)hists2D.size(); i++){
				for(int j = 0; j < hists2D[i].size(); j++){
					if(pc.hists2D[i][j] == nullptr) continue;
					//write total hist to file
					name = pc.hists2D[i][j]->GetName();
					//if ends in "_" remove
					if(strcmp(&name[name.size() - 1],"_") == 0) name.pop_back();
					name += "2D";
					pc.hists2D[i][j]->SetName(name.c_str());		
					pc.hists2D[i][j]->SetTitle(pc.plotName.c_str());
					if(pc.hists2D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram: " << name << " not filled." << endl; continue; }
					pc.hists2D[i][j]->Write();
				}
			}

		}



		void WritePlotCatStack(TFile* ofile, const vector<procCat>& pcs){
			ofile->cd();
			string name;
			//number of histogram categories (ie leading, !leading, etc)
			int nhistCats = pcs[0].hists1D.size();
			//write 1D hists
			//variables
			for(int j = 0; j < (int)_hists1D.size(); j++){
				//lead, not lead, etc.
				name = _hists1D[j]->GetName();
				TDirectory* dir = ofile->mkdir((name+"_stack").c_str());
				dir->cd(); 
				for(int i = 0; i < nhistCats; i++){
					if(pcs[0].hists1D[i][j] == nullptr) continue;
					//needs to be reset for each category
					//name = _hists1D[j]->GetName();
					if(!pcs[0].histcatnames[i].empty()) name += "_"+pcs[0].histcatnames[i];
					//cout << "	category: " << pcs[0].histcatnames[i] << endl; 
					//proc
					for(int k = 0; k < pcs.size(); k++){
						if(pcs[k].hists1D[i][j] == nullptr) continue;
						if(pcs[k].hists1D[i][j]->GetEntries() == 0){ continue; }//cout << "Histogram for proc " << pcs[k].plotName << " not filled." << endl; continue; }
						//cout << "		adding proc " << pcs[k].plotName << " to plot with hist " << pcs[k].hists1D[i][j]->GetName() << endl;
						pcs[k].hists1D[i][j]->SetTitle(pcs[k].plotName.c_str());
						pcs[k].hists1D[i][j]->Write();
					}
				}
			}
			
		}


		void WriteHists(TFile* ofile){
			//normalize histograms	
			for(int i = 0; i < (int)_procCats.size(); i++){
				for(int j = 0; j < _procCats[i].hists1D.size(); j++){
					//relative fraction histograms
					//nSubClusters
					//_procCats[i].hists1D[j][0]->Scale(1./_procCats[i].hists1D[j][0]->Integral());
					//ellipsoid center coordinates
					_procCats[i].hists1D[j][1]->Scale(1./_procCats[i].hists1D[j][1]->Integral());
					_procCats[i].hists1D[j][2]->Scale(1./_procCats[i].hists1D[j][2]->Integral());
					_procCats[i].hists1D[j][3]->Scale(1./_procCats[i].hists1D[j][3]->Integral());
					//theta + azimuthal angles
					_procCats[i].hists1D[j][7]->Scale(1./_procCats[i].hists1D[j][7]->Integral());
					_procCats[i].hists1D[j][8]->Scale(1./_procCats[i].hists1D[j][8]->Integral());
					//rotundity
					_procCats[i].hists1D[j][10]->Scale(1./_procCats[i].hists1D[j][10]->Integral());
					_procCats[i].hists1D[j][11]->Scale(1./_procCats[i].hists1D[j][11]->Integral());
					//velocity
					_procCats[i].hists1D[j][14]->Scale(1./_procCats[i].hists1D[j][14]->Integral());
				}
			}

			WritePlotCat1D(ofile, _procCats[0]);
			for(int i = 0; i < (int)_procCats.size(); i++)
				WritePlotCat2D(ofile, _procCats[i]);
			vector<procCat> id_cats(_procCats.begin()+1, _procCats.end());
			
			if(_procCats.size() > 1) WritePlotCatStack(ofile, id_cats);

			ofile->Close();

		}


		void TrackMatched(Matrix mu, double& bestdr, double& bestp){
			//do track matching
			bestdr = 999;
			bestp = 999;
			double bestTrackDr = 999;
			//double maxTrackDr;
			double dr, teta, tphi, de;
			unsigned int detid;

			double dphi = -999;
			//double bestde_dr;
			int trackidx;
			double ec = mu.at(0,0);
			double pc = mu.at(1,0);

			double track_phi = -999;
			double track_eta = -999;

			//loop through tracks to get best match to this subcluster (tracks are matched to superclusters, if not idx < 0)
			//Track_scIndexs[i][j] is for track i that matched to supercluster j (tracks can match to multiple SCs)
			int nTracks = _base->Track_scIndexs->size();
			for(int t = 0; t < nTracks; t++){
				if(_base->Track_scIndexs->at(t).at(0) < 0) continue; //not matched to any SC
				
				int nSCs = _base->Track_scIndexs->at(t).size();
//cout << "track #" << t << " matched to " << nSCs << " superclusters" << endl;
				for(int sc = 0; sc < nSCs; sc++){
					int sc_idx = _base->Track_scIndexs->at(t).at(sc);
					//get eta, phi of supercluster sc that track is matched to
					double sc_eta = _base->SuperCluster_eta->at(sc_idx);
					double sc_phi = _base->SuperCluster_phi->at(sc_idx);
					double sc_phi_02pi = sc_phi;
					if(sc_phi_02pi < 0) sc_phi_02pi += 2*acos(-1);
					else if(sc_phi_02pi > 2*acos(-1)) sc_phi_02pi -= 2*acos(-1); 
					else sc_phi_02pi = sc_phi;
//cout << "   track #" << t << " matched to supercluster " << sc << " with eta " << sc_eta << " and phi " << sc_phi << endl;

					dphi = fabs(pc - sc_phi_02pi);
                                	dphi = acos(cos(dphi));

					dr = sqrt((sc_eta - ec)*(sc_eta - ec) + dphi*dphi);

	
					//E = p for photons
					if(dr < bestTrackDr){
						bestTrackDr = dr;
						track_phi = sc_phi;
						track_eta = sc_eta;
						bestp = _base->Track_p->at(t);
					}
					
				}
				//cout << "best dr " << bestdr << " " << bestTrackDr << " track eta " << track_eta << " SC eta " << ec << " track phi " << track_phi << " SC phi " << pc << endl;
				bestdr = bestTrackDr;
			}

		}

		void SetObs(){
			SuperClusterAddBranches();
			//sample
			_inputs.push_back("sample");
			//event
			_inputs.push_back("event");
			//event weight
			_inputs.push_back("event_weight");
			//supercl
			_inputs.push_back("object");
		}
		void WriteHeader(){
			//CNN inputs
			for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++)
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++)
					_inputs.push_back("CNNgrid_cell"+to_string(i)+"_"+to_string(j));
			//label
			_inputs.push_back("label");
			for(auto s : _inputs){
				if(s != "label") _csvfile << s << ","; 
				else _csvfile << s << endl;
			}

		}

		void FillBranches(const Jet& bhc_obj){
			double E_tot = bhc_obj.E();
			vFillBranch(E_tot, "Energy");		


			double ec = bhc_obj.eta();
			double pc = bhc_obj.phi();
			double tc = bhc_obj.t();
			if(isnan(pc)) cout << "pc is nan" << endl;
			if(isinf(pc)) cout << "pc is inf" << endl;
			if(pc < 0 || pc > 2*acos(-1)) cout << "pc out of bounds " << pc << endl;
			vFillBranch(ec, "EtaCenter");		
			vFillBranch(pc, "PhiCenter");		
			vFillBranch(tc, "TimeCenter");		
			
			Matrix cov = bhc_obj.GetCovariance();
			vFillBranch(cov.at(0,0), "EtaVar");		
			vFillBranch(cov.at(1,1), "PhiVar");		
			vFillBranch(cov.at(2,2), "TimeVar");
			vFillBranch(cov.at(0,1), "EtaPhiCov");
			vFillBranch(cov.at(0,2), "EtaTimeCov");
			vFillBranch(cov.at(1,2), "PhiTimeCov");


			
			vector<double> spikeObs;
			PointCollection* points = GetPointsFromJet(bhc_obj);
			vFillBranch((double)points->GetNPoints(),"nRHs_bhcSC");

			SpikeObs(points, spikeObs);
			double swCP = spikeObs[0];
			vFillBranch(swCP, "swCrossPrime");
			vector<Jet> rh_jets;
			bhc_obj.GetJets(rh_jets);
			double sw = swissCross(rh_jets);
			vFillBranch(sw, "swCrossCMS");	
		
			//time signifiance	
			double timesig = CalcTimeSignificance(points);
			vFillBranch(timesig, "timeSignificance");			
	
		}



	double sqrtcov(double c){
		if(c > 0)
			return sqrt(c);
		else
			return -sqrt(-c);

	}


	double CalcCov(const Matrix& cov, int i, int j, bool norm = true){
		if(i >= cov.nRows()) return -999;
		if(j >= cov.nCols()) return -999;
		if(cov.empty()) return -999;
		double denom = 1;
		if(norm){
			denom = sqrt(cov.at(i,i)*cov.at(j,j));
		}
		//return sqrtcov(cov.at(i,j)/denom);
		return cov.at(i,j)/denom;
	}


	void SpikeObs(PointCollection* pc, vector<double>& obs){
		obs.clear();
		int npts = pc->GetNPoints();
		double wmax = 0;
		BayesPoint xmax;
		for(int i = 0; i < npts; i++){
			if(pc->at(i).w() > wmax){
				wmax = pc->at(i).w();
				xmax = pc->at(i);
			}
		}
		//find 4 closest neighbors to xmax
		//map<double, Point> distToMax;
		double dist, deta, dphi;
		//distance threshold to include rh in sw+'
		double distThresh = 0.02;
		int nNeighbors = 0;
		double sumNeighbors = 0;
		for(int i = 0; i < npts; i++){
			//skip xmax
			if(xmax == pc->at(i)) continue;
			deta = pc->at(i).at(0) - xmax.at(0);
			dphi = pc->at(i).at(1) - xmax.at(1);
			if(dphi > 4*atan(1)) dphi -= 8*atan(1);
			dist = sqrt(deta*deta + dphi*dphi);
			if(dist < distThresh){
				sumNeighbors += pc->at(i).w();
				nNeighbors++;
			}
			//distToMax[dist] = points->at(i);
		}
		//swCP
		obs.push_back(sumNeighbors/wmax);

	}

	double swissCross(const vector<Jet>& jets){
		//find seed crystal (highest weight, E = w/_gev)
		PointCollection sc(1); //save cmsswIds with associated weights for rec hit id
		JetPoint rh;
		for(int j = 0; j < jets.size(); j++){
			//should only have 1 rh per jet
			if(jets[j].GetNRecHits() > 1) continue;
			jets[j].GetJetPointAt(0,rh);
			BayesPoint pt(1);
			pt.SetValue(int(rh.rhId()),0);
			pt.SetWeight(rh.GetWeight());
			sc += pt;	
		}
		sc.Sort();
		double e1, e4;
		e1 = sc.at(sc.GetNPoints()-1).w()/_gev;
		int e1id = sc.at(sc.GetNPoints()-1).at(0);
		int e1_iphi = _detIDmap[e1id].i1;
		int e1_ieta = _detIDmap[e1id].i2;
		
		//find up, left, right, down crystals (e4)
		//eta on x-axis - shift in eta = left/right
		//phi on y-axis - shift in phi = up/down
		int e4upid = offsetBy(e1_ieta, e1_iphi, 0, 1);
		int e4downid = offsetBy(e1_ieta, e1_iphi, 0, -1);
		int e4rightid = offsetBy(e1_ieta, e1_iphi, 1, 0);
		int e4leftid = offsetBy(e1_ieta, e1_iphi, -1, 0);
		e4 = 0;
		for(int i = 0; i < sc.GetNPoints(); i++){
			if(sc.at(i).at(0) == e4upid && e4upid != -1)
				e4 += sc.at(i).w()/_gev;
			if(sc.at(i).at(0) == e4downid && e4downid != -1)
				e4 += sc.at(i).w()/_gev;
			if(sc.at(i).at(0) == e4leftid && e4leftid != -1)
				e4 += sc.at(i).w()/_gev;
			if(sc.at(i).at(0) == e4rightid && e4rightid != -1)
				e4 += sc.at(i).w()/_gev;
		}
		return 1 - e4/e1;

	}

	//from CMSSW: https://cmssdt.cern.ch/lxr/source/DataFormats/EcalDetId/src/EBDetId.cc
	UInt_t offsetBy(int ieta, int iphi, int nrStepsEta, int nrStepsPhi) const {
		int newEta = ieta + nrStepsEta;
		if (newEta * ieta <= 0) {
		  if (ieta < 0) {
		    newEta++;
		  } else if (ieta > 0) {
		    newEta--;
		  }
		}
		int newPhi = iphi + nrStepsPhi;
		while (newPhi > 360)
		  newPhi -= 360;
		while (newPhi <= 0)
		  newPhi += 360;
		pair<int, int> newpair = make_pair(newEta, newPhi);	
		if (validDetId(newEta, newPhi)) {
		  UInt_t id = _ietaiphiID.at(newpair);
		  return id;
		} else {
		  return -1;
		}
	}

	//from CMSSW: https://cmssdt.cern.ch/lxr/source/DataFormats/EcalDetId/interface/EBDetId.h
	/// check if a valid index combination
	static bool validDetId(int i, int j) {
		int max_ieta = 85;
		int max_iphi = 360;
		int min_ieta = 1;
		int min_iphi = 1;
		return i != 0 && (std::abs(i) <= max_ieta) && (j >= min_iphi) && (j <= max_iphi);
	}






	//recreate R9 variable
	double sc3x3E(vector<JetPoint>& rhs){
		//find seed xtal
		double maxE = 0;
		JetPoint seed;
		for(auto rh : rhs){
			if(rh.E() > maxE){ seed = rh; maxE = seed.E(); } 
		}
		//find 3x3 grid centered on xmax
		int detid = seed.rhId();
		int ieta = _detIDmap[detid].i2;
		int iphi = _detIDmap[detid].i1;
		double gridE = 0;
			
		int id, ie, ip, ide, idp;
		for(auto rh : rhs){
			id = rh.rhId();
			ie = _detIDmap[detid].i2;
			ip = _detIDmap[detid].i1;
			ide = ieta - ie;
			idp = iphi - ip;
			if(fabs(ide) < 2 && fabs(idp) < 2){
				gridE += rh.E();
			}
		}
		return gridE;	
	}



	//return BH cluster size, rhs to use in time discriminant
	int MakeBHFilterCluster(Jet& sc, vector<JetPoint>& bh_rhs){
		bh_rhs.clear();
		//get rhs
		vector<JetPoint> rhs; sc.GetJetPoints(rhs);
		//find seed xtal
		double maxE = 0;
		pair<int, int> iSeed; //ieta, iphi of seed
		unsigned int seedId;
		for(int i = 0; i < rhs.size(); i++){
			if(rhs[i].E() > maxE){
				maxE = rhs[i].E();
				iSeed = make_pair(_detIDmap[rhs[i].rhId()].i2, _detIDmap[rhs[i].rhId()].i1);
				seedId = rhs[i].rhId();
			}
		}
		//seed must be at least 5 GeV
		if(maxE < 5) return -1;
		//build SC band
		double deta = acos(-1)/180;
		double dphi = acos(-1)/180;
		double maxEta = 0.4;
		double maxPhi = 0.16;

		int maxiEta = round(maxEta/deta);
		int maxiPhi = round(maxPhi/dphi); 

		unsigned int rhid;
		int ieta, iphi;
		int nRhs = 0;
		for(int i = 0; i < rhs.size(); i++){
			rhid = rhs[i].rhId();
			//skip seed
			if(rhid == seedId) continue;
			//check if rh is in grid
			ieta = _detIDmap[rhid].i2;
			iphi = _detIDmap[rhid].i1;
			if(fabs(ieta - iSeed.first) > maxiEta/2) continue;
			if(fabs(iphi - iSeed.second) > maxiPhi/2) continue;
			//if in grid, check that meets energy threshold
			if(rhs[i].E() <= 1) continue;
			//make sure no rhs with E too high is found on outer phi edges
			if(fabs(iphi - iSeed.second) > 1){
				if(rhs[i].E() > 2) continue;
			}
			//grab rhs in central and neighbor bands for time discriminant
			if(fabs(iphi - iSeed.second < 2)) bh_rhs.push_back(rhs[i]);
			//get cluster size from central band
			if(fabs(iphi - iSeed.second == 0) && rhs[i].E() > 2) nRhs++;
		}
		return nRhs;
	}	
	double BHTimeDiscriminant(vector<JetPoint>& rhs){
		double td = 0;
		double fsep;
		int ieta;
		for(auto rh : rhs){
			ieta = _detIDmap[rh.rhId()].i2;
			fsep = -0.5 * (sqrt( 130*130 + 9*(double)ieta*(double)ieta) - 3*(double)ieta)/30;
			td += (rh.t() - fsep)*log10(rh.E());
		}
		return td;
	}

	bool BHclusterIso(vector<JetPoint>& rhs){
		unsigned int rhid, seedId;
		int ieta, iphi;
		//find seed xtal
		double maxE = 0;
		pair<int, int> iSeed; //ieta, iphi of seed
		for(int i = 0; i < rhs.size(); i++){
			if(rhs[i].E() > maxE){
				maxE = rhs[i].E();
				iSeed = make_pair(_detIDmap[rhs[i].rhId()].i2, _detIDmap[rhs[i].rhId()].i1);
				seedId = rhs[i].rhId();
			}
		}
		for(int i = 0; i < rhs.size(); i++){
			rhid = rhs[i].rhId();
			//check if rh is in grid
			ieta = _detIDmap[rhid].i2;
			iphi = _detIDmap[rhid].i1;
			//if any rhs in neighbor bands, return false (failed isolation) 
			if(fabs(iphi - iSeed.second == 1)) return false;
		}
		return true;
	}
	int GetTrainingLabel(int nobj, const Jet& bhc_sc, const Jet& og_sc){
		//labels
		//unmatched = -1

		//sig vs bkg - PhotonSkimmer
		//signal = 0
		//iso bkg = 4

		//iso vs !iso - PhotonSkimmer
		//iso sig = 5
		//!iso bkg = 6

		//phys bkg vs det bkg
		//phys bkg = 1 
		//BH = 2
		//spike = 3
	
		double ec, pc, tc;
		double E = bhc_sc.E(); 
		Matrix mu, cov;
		bhc_sc.GetClusterParams(mu, cov);


		ec = mu.at(0,0);
		pc = mu.at(1,0);
		//tc = mu.at(2,0);
		//use seed time from original SC
		unsigned int seed_id = _base->SuperCluster_XtalSeedID->at(og_sc.GetUserIdx());
		//time significance to enhance purity of det bkg CRs
		//using seed of original SC for time and time significance
		PointCollection* points = GetPointsFromJet(og_sc);
		double tsig_seed = CalcTimeSignificance(points, seed_id);
		double tsig = CalcTimeSignificance(points);
		tc = -999;
		for(int n = 0; n < points->GetNPoints(); n++){
			if(points->at(n).GetUserIdx() == seed_id){
				tc = points->at(n).at(2);
				break;
			}
		}
		//if tc is -999 it means that the seed rh was removed in processing (due to either time req or eta req)
		vFillBranch(tc,"seedTime");
		vFillBranch(tsig_seed,"seedTimeSignificance");
		vFillBranch(tsig,"weightedTimeSignificance");	
	
		//do track matching for spikes
		og_sc.GetClusterParams(mu, cov);
		double bestdr, bestde_dr;
		TrackMatched(mu,bestdr, bestde_dr);
		if(bestde_dr != 999) bestde_dr = og_sc.e()/bestde_dr;
		vFillBranch(bestdr,"dR_track");
		vFillBranch(bestde_dr,"EovP_track");
		
		bool trksum, ecalrhsum, htowoverem;
		int iso;
		//for BH definition
		//phi center is either at ~0, ~pi, ~2pi (within ~10 xtals)
		bool pcFilter = (pc < 0.1 || (acos(-1) - 0.1 < pc && pc < acos(-1) + 0.1) || 2*acos(-1) - 0.1 < pc );
		bool pcFilter_wide = (pc < 0.3 || (acos(-1) - 0.3 < pc && pc < acos(-1) + 0.3) || 2*acos(-1) - 0.3 < pc );
		bool spikeTime = (tc <= -10 && tc != -999);
		bool bhTime = (-7 < tc && tc <= -2);
		bool physBkgTime = (-0.5 < tc && tc <= 0.5); 
		bool detBkgTimeSig = (tsig < -3);
		bool physBkgTimeSig = (-1 < tsig && tsig < 1);
		bool spikeTrackMatch = (bestdr <= 0.02);
		bool notSpikeTrackVeto = (bestdr > 0.03);
		bool evtfilters = _base->Flag_BadChargedCandidateFilter && _base->Flag_BadPFMuonDzFilter && _base->Flag_BadPFMuonFilter && _base->Flag_EcalDeadCellTriggerPrimitiveFilter && _base->Flag_HBHENoiseFilter && _base->Flag_HBHENoiseIsoFilter && _base->Flag_ecalBadCalibFilter && _base->Flag_goodVertices && _base->Flag_hfNoisyHitsFilter;
		bool bh_filter = _base->Flag_globalSuperTightHalo2016Filter;

		int label = -1;
		bool passGJetsObjSel = GJetsCR_ObjSel(og_sc,false);
cout << "time " << tc << " phi " << pc << " time sig " << tsig << " dr track " << bestdr << " gjets evt sel " << _passGJetsEvtSel << " gjets obj sel " << passGJetsObjSel << endl;	
		//in data - could be spikes or BH
		if(_data){
			//early times, phi left/right for BH
			//if subcl is spike
//cout << "time center " << tc << " phi center " << pc << " dr to track " << bestdr << " pcfilter BH " << pcFilter << endl;
			int phoidx = _base->SuperCluster_PhotonIndx->at(nobj);
			if(phoidx == -1){
				//not not matched to an e/gamma candidate so cant be det bkg
cout << "no photon match" << endl; 
				label = -1;
				iso = -1;
			}
			else{ //matched to e/gamma candidate
                		trksum = _base->Photon_trkSumPtSolidConeDR04->at(phoidx) < 6.0;
                		ecalrhsum = _base->Photon_ecalRHSumEtConeDR04->at(phoidx) < 10.0;
                		htowoverem = _base->Photon_hadTowOverEM->at(phoidx) < 0.02;
                		iso = trksum && ecalrhsum && htowoverem;
cout << "pass iso? " << iso << " Photon_trkSumPtSolidConeDR04 " << _base->Photon_trkSumPtSolidConeDR04->at(phoidx) << " Photon_ecalRHSumEtConeDR04 " << _base->Photon_ecalRHSumEtConeDR04->at(phoidx) << " Photon_hadTowOverEM " << _base->Photon_hadTowOverEM->at(phoidx) << endl;
                		if(!iso) label = -1; //not isolated photon - won't make it into analysis anyway
				//do track seed matching for photon (bh, phys bkg) and e (spike) candidates
cout << "pass pixel seed " << _base->Photon_pixelSeed->at(phoidx) << endl;
				//spike
				//if(spikeTrackMatch && spikeTime && !pcFilter_wide && detBkgTimeSig && evtfilters && iso){
				if(spikeTrackMatch && spikeTime && !pcFilter_wide && detBkgTimeSig && iso){
					label = 3;
				}
				//for physics bkg + BH match to photons + apply isolation criteria
				//if subcl is BH - need to match to photon and apply isolation
				//if(bhTime && pcFilter && iso && notSpikeTrackVeto && detBkgTimeSig && evtfilters){	
				if(bhTime && pcFilter && iso && notSpikeTrackVeto && detBkgTimeSig){	
					label = 2;
				}
				//if subcl is not BH or spike (ie prompt, 'physics' bkg) - need to match to photon and apply isolation
				if(physBkgTime && iso && physBkgTimeSig && !pcFilter && notSpikeTrackVeto && _passGJetsEvtSel && passGJetsObjSel && evtfilters && bh_filter){
					label = 1;
				}

			}
			vFillBranch((double)iso,"isoPresel");
			vFillBranch((double)passGJetsObjSel,"PassGJetsCR_Obj");
			double pixseed = phoidx == -1 ? -1 : (double)_base->Photon_pixelSeed->at(phoidx); 
			vFillBranch(pixseed,"passPixelSeed");
		}
		else // not data
			label = -1;

		return label;
	}
	int _nBH_hist = 0;
	int _nSpike_hist = 0;

};		
#endif
