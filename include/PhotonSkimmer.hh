#ifndef PHOTONSKIMMER_HH
#define PHOTONSKIMMER_HH

#include "JetPoint.hh"
#include <TFile.h>
#include "BaseSkimmer.hh"
#include "PhotonProducer.hh"
#include "JetProducer.hh"
#include "TSystem.h"
#include <math.h>
#include <fstream>

using procCat = BaseSkimmer::procCat;
class PhotonSkimmer : public BaseSkimmer{
	public:
		PhotonSkimmer(){
			SetObs();
			InitMapTree();
			_evti = 0;
			_evtj = 0;
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 1e-300;
			_emAlpha = 1e-5;
			_gev = 1/30.;
			_applyFrac = false;


			_weight = 1;
			_cell = 0;
			_tresCte = 0;
			_tresNoise = 0;
			_tresStoch = 0;

		};
		virtual ~PhotonSkimmer(){ };

		//get rechits from file to cluster
		//PhotonSkimmer(TFile* file) : BaseSkimmer(file){
		PhotonSkimmer(string filename) : BaseSkimmer(filename){
			SetObs();
			InitMapTree();
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time
			_prod = make_unique<PhotonProducer>(filename);
			_fname = filename;

			
			_base = _prod->GetBase();
			_nEvts = _base->fChain->GetEntries();
			_evti = 0;
			_evtj = _nEvts;
			_oname = "plots/photon_skims_"+_cms_label+".root";
			
			_isocuts = false;
			_oskip = 10;
			_thresh = 1.;
			_alpha = 1e-300;
			_emAlpha = 1e-5;
			_gev = 1/30.;
			_applyFrac = false;
			
			SetupDetIDs();
		}

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
		void SetThresh(double t){ _thresh = t; }
		void SetBHCAlpha(double a){ _alpha = a; }
		void SetEMAlpha(double a){ _emAlpha = a; }
		double _thresh, _alpha, _emAlpha; 
		void ApplyFractions(bool a){ _applyFrac = a; if(_applyFrac) cout << "Applying RH fractions" << endl; }
		bool _applyFrac;


		void PhotonAddBranches(){
			_obj = "Photon";
			_obsnames.push_back("Pt");
			//sw+
			_obsnames.push_back("swCP");
			//subcl major length
			_obsnames.push_back("majorLength");
			//subcl minor length
			_obsnames.push_back("minorLength");
			_obsnames.push_back("rot2D");
			_obsnames.push_back("phiE2D");
			//subcl max pt / total E
			_obsnames.push_back("maxOvtotE");
			//R9
			_obsnames.push_back("R9");
			//Sietaieta
			_obsnames.push_back("Sietaieta");
			//Siphiiphi
			_obsnames.push_back("Siphiiphi");
			//Smajor
			_obsnames.push_back("Smajor");
			//Sminor
			_obsnames.push_back("Sminor");
			//hcalTowerSumEtConeDR04
			_obsnames.push_back("hcalTowerSumEtConeDR04");
			//trkSumPtSolidConeDR04
			_obsnames.push_back("trkSumPtSolidConeDR04");
			//trkSumPtHollowConeDR04
			_obsnames.push_back("trkSumPtHollowConeDR04");
                	//ecalRHSumEtConeDR04
                	_obsnames.push_back("ecalRHSumEtConeDR04");
                	//hadTowOverEM
                	_obsnames.push_back("hadTowOverEM");
			
			vAddBranch("predScore_isoBkg","score for iso bkg designation from DNN"); //DNN prediction
			vAddBranch("predScore_nonIsoBkg","score for non iso bkg designation from DNN"); //DNN prediction
			vAddBranch("trueLabel","true DNN label of photon");
			
		}
		//obsnames will save as branches
		//inputs will save to training CSV


		void SetObs(){
			PhotonAddBranches();
			//sample
			_inputs.push_back("sample");
			//event
			_inputs.push_back("event");
			//event weight
			_inputs.push_back("event_weight");
			//supercl
			_inputs.push_back("object");
			for(int o = 0; o < _obsnames.size(); o++){
				_inputs.push_back(_obsnames[o]);
			}
			_inputs.push_back("2017_presel");
			//label
			_inputs.push_back("label");
		}

		void WriteHeader(){
			for(auto s : _inputs){
				if(s != "label") _csvfile << s << ","; 
				else _csvfile << s << endl;
			}
		}


		void InitObs(map<string, double>& obs){
			for(int i = 0; i < _inputs.size(); i++)
				obs[_inputs[i]] = -999;
		}


		void FillBranches(const Jet& bhc_obj, string tag = ""){ }

		void FillBranches(map<string, double> obs){
			for(auto it = obs.begin(); it != obs.end(); it++){
				if(it->second == -999) continue; //don't fill branches that haven't had their value updated
				vFillBranch(it->second, it->first);
			}

		}

		void FillJetObs(const Jet& bhc_obj, map<string, double>& obs){
			double E_tot = bhc_obj.E();
			obs.at("Energy") = E_tot;
			
			double ec = bhc_obj.eta();
			double pc = bhc_obj.phi();
			double tc = bhc_obj.t();
			if(isnan(pc)) cout << "pc is nan" << endl;
			if(isinf(pc)) cout << "pc is inf" << endl;
			if(pc < 0 || pc > 2*acos(-1)) cout << "pc out of bounds " << pc << endl;
			obs.at("EtaCenter") = ec;		
			obs.at("PhiCenter") = pc;		
			obs.at("TimeCenter") = tc;		
			
			Matrix cov = bhc_obj.GetCovariance();
			obs.at("EtaVar") = cov.at(0,0);		
			obs.at("PhiVar") = cov.at(1,1);		
			obs.at("TimeVar") = cov.at(2,2);		
			obs.at("EtaPhiCov") = cov.at(0,1);		
			obs.at("EtaTimeCov") = cov.at(0,2);		
			obs.at("PhiTimeCov") = cov.at(1,2);


			//rotundity - 2D
			//take upper 2x2 submatrix from covariance
			Matrix space_mat(2,2);
			Get2DMat(cov,space_mat);
			vector<Matrix> eigenvecs_space;
			vector<double> eigenvals_space;
			space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
			double majLength_2D = sqrt(eigenvals_space[1]);
			if(eigenvals_space[1] < 0) cout << "negative eigenvalue " << eigenvals_space[1] << endl;
			double minLength_2D; 
			if(eigenvals_space[0] < 0) minLength_2D = -sqrt(-eigenvals_space[0]);
			else minLength_2D = sqrt(eigenvals_space[0]);	
			double phi2D = PhiEll(space_mat);
			double rot2D = Rotundity(space_mat);
			
			obs.at("majorLength") = majLength_2D;
			obs.at("minorLength") = minLength_2D;
			obs.at("rot2D") = rot2D;
			obs.at("phiE2D") = phi2D;

			std::unique_ptr<PointCollection> points = GetPointsFromJet(bhc_obj);
			points->Sort();
			double maxE = points->at(points->GetNPoints() - 1).w();
			obs.at("maxOvtotE") = maxE/E_tot;
			
			vector<double> spikeObs;
			SpikeObs(points.get(), spikeObs);
			double swCP = spikeObs[0];
			obs.at("swCP") = swCP;

			//time significance not set

		}




	void MakeCovMat(PointCollection* pc, Matrix& outcov, const weightScheme& ws){
		if(!outcov.square()) return;
		if(outcov.nRows() != pc->Dim()) return;
		
		//set weights to logE
		//og w_i = _gev*E_i
		int npts = pc->GetNPoints();
		int maxd = pc->Dim();
		double E_tot = 0.;
		//CMSSW value for super clusters
		double w0 = 4.7;
		//zero suppression involves a hitsAndFractions transfer factor (noZS in cmssw)
		//not user here (is 1)
		for(int i = 0; i < npts; i++){
			E_tot += pc->at(i).w()/_gev;
		}
		PointCollection pcnew;
		BayesPoint pt;
		double denom = 0;
		BayesPoint mean = BayesPoint(maxd);
		double meta = 0;
		double mphi = 0;
		double mtime = 0;
		for(int i = 0; i < npts; i ++){
			pt = pc->at(i);
			if(ws == 0) pt.SetWeight( 1.0 );
			//already e-weighted
			if(ws == 2) pt.SetWeight( log( w0 + (pc->at(i).w()/_gev)/E_tot ) );
			pcnew.AddPoint(pt);
		}
		//center at 0 - [-pi, pi]
		pcnew.Center();
		//calculate weighted mean - with phi wraparound
		mean.SetValue(pcnew.Centroid(0),0);
		mean.SetValue(pcnew.Centroid(1),1);
		mean.SetValue(pcnew.Centroid(2),2);
		double ent;
		double ent_pt, dd1, dd2;
		double pi = acos(-1);
		for(int d1 = 0; d1 < maxd; d1++){
			for(int d2 = d1; d2 < maxd; d2++){ 
				ent = 0;
				for(int i = 0; i < npts; i++){
					dd1 = pcnew.at(i).at(d1) - mean.at(d1);
					//phi wraparound
					if(d1 == 1){
						//if(dd1 > pi) dd1 = 2*pi - dd1;
						dd1 = acos(cos(dd1));
						if(dd1 > pi) cout << "dd1: " << dd1 << endl;
						
					}
					dd2 = pcnew.at(i).at(d2) - mean.at(d2);
					if(d2 == 1){
						//if(dd2 > pi) dd2 = 2*pi - dd2;
						dd2 = acos(cos(dd2));
						if(dd2 > pi) cout << "dd2: " << dd2 << endl;
					}
					ent_pt = pcnew.at(i).w() * (dd1)*(dd2) / pcnew.Sumw();
					//if(d1 == 1) ent_pt *= sqrt(phiCorrectionFactor);
					//if(d2 == 1) ent_pt *= sqrt(phiCorrectionFactor);
					ent += ent_pt;
				}
				outcov.SetEntry(ent,d1,d2);
				if(d1 != d2) outcov.SetEntry(ent,d2,d1);
			}
		}
		//eta time sign convention
		//if(mean.at(0) < 0){
		//	//time sign does NOT match eta sign
		//	//flip sign of eta-time entry
		//	outcov.SetEntry(-outcov.at(0,2),0,2);	
		//	outcov.SetEntry(-outcov.at(2,0),2,0);	
		//}



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


	double Rotundity(Matrix& inmat){
		vector<Matrix> eigenvecs;
		vector<double> eigenvals;
		inmat.eigenCalc(eigenvals, eigenvecs);
		int maxd = inmat.nRows() - 1;
		double rot = 0;
		for(int i = 0; i < (int)eigenvals.size(); i++) rot += eigenvals[i];
		rot = eigenvals[maxd]/rot;
		if(isnan(rot))
			cout << "rot is nan - eigenvals[0] " <<  eigenvals[0] << " eigenvals[1] " << eigenvals[1] << endl;
		//if(rot < 0.5 || rot > 1) cout << "rot: " << rot << endl;
		return rot;
	}

	double PhiEll(Matrix& inmat){
		vector<Matrix> eigenvecs;
		vector<double> eigenvals;
		inmat.eigenCalc(eigenvals, eigenvecs);
		int maxd = inmat.nRows() - 1;
		double v_x = eigenvecs[maxd].at(0,0);	
		double v_y = eigenvecs[maxd].at(0,1);	
		//azimuthal angle with lead eigenvector (from 2D spatial submatrix)
		double phi = atan2( v_y , v_x );
		return phi;
	}

	void SpikeObs(PointCollection* pc, vector<double>& obs){
		obs.clear();
		pc->Sort();
		int npts = pc->GetNPoints();
		double wmax = pc->at(npts-1).w();
		BayesPoint xmax = pc->at(npts-1);
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
		//find seed crystal (highest weight, E = w*_gev)
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

	//used for sig/bkg MVA, iso/!iso MVA
	int GetTrainingLabel(int phoidx, const Jet& bhc_pho, const Jet& og_pho){
		//labels
		//unmatched = -1

		//signal = 0
		//iso bkg = 4
		//iso bkg (gen level) = 5
		//!iso bkg = 6

		//phys bkg vs det bkg - these are done in SuperClusterSkimmer
		//phys bkg = 1 
		//BH = 2
		//spike = 3
		
		//bool trksum, ecalrhsum, htowoverem, iso;	
		bool evtfilters = _base->Flag_BadChargedCandidateFilter && _base->Flag_BadPFMuonDzFilter && _base->Flag_BadPFMuonFilter && _base->Flag_EcalDeadCellTriggerPrimitiveFilter && _base->Flag_HBHENoiseFilter && _base->Flag_HBHENoiseIsoFilter && _base->Flag_ecalBadCalibFilter && _base->Flag_goodVertices && _base->Flag_hfNoisyHitsFilter;
		bool bh_filter = _base->Flag_globalSuperTightHalo2016Filter;
		int label = -1;
		bool isoBkgSel = GJetsCR_ObjSel(og_pho, true);
		bool nonIsoBkgSel = DijetsCR_ObjSel(og_pho);
		bool passMinPt = (og_pho.pt() > _minPhoPt_CRsel);

		vFillBranch((double)isoBkgSel,"PassGJetsCR_Obj");
		//MC
		if(!_data){
                	//trksum = _base->Photon_trkSumPtSolidConeDR04->at(phoidx) < 6.0;
                	//ecalrhsum = _base->Photon_ecalRHSumEtConeDR04->at(phoidx) < 10.0;
                	//htowoverem = _base->Photon_hadTowOverEM->at(phoidx) < 0.02;
                	//iso = trksum && ecalrhsum && htowoverem;
			if(_base->Photon_genIdx->at(phoidx) != -1){
				int genidx = _base->Photon_genIdx->at(phoidx);
                		//needs to be isolated
				//non isolated bkg
				if(_oname.find("SMS") != string::npos){
					//photon from C2
					if(_base->Gen_susId->at(genidx) == 22)
						label = 0;
					//photon from hard subprocess - isolated sig 
					else if(_base->Gen_status->at(genidx) == 23)
						label = 5;
					else if(_base->Gen_motherIdx->at(genidx) != -1 && _base->Gen_status->at(_base->Gen_motherIdx->at(genidx)) == 23)
						label = 5;
					else
						label = 1; //removal of GMSB !sig photons is done in data processing for NN

				}
				//isolated bkg
				if(isoBkgSel && _passGJetsEvtSel && evtfilters && bh_filter)
					label = 4;
				else if(!isoBkgSel && nonIsoBkgSel && _passDijetsEvtSel && evtfilters && bh_filter)
					label = 6;
				else
					label = -1;
			}
			else //no gen match
				label = -1;
		}
		//else in data - could be spikes or BH
		else{
			//iso bkg = 4
			
			//obj selection for iso bkg
			//non isolated bkg
			//isolated bkg
			//isolated bkg
			if(isoBkgSel && _passGJetsEvtSel && evtfilters && bh_filter)
				label = 4;
			else if(!isoBkgSel && _passDijetsEvtSel && evtfilters && bh_filter && passMinPt)
				label = 6;
			else label = -1; 
		}
		return label;
	}



	private:
		string _fname; 

};
#endif
