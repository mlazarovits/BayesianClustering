#ifndef BaseProducer_HH
#define BaseProducer_HH

#include "ReducedBase.hh"
#include "RandomSample.hh"
#include "Jet.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TCanvas.h"
#include <fstream>

class BaseProducer{
	public:
		BaseProducer(){ 
			_gev = 1;
			_isocut = false;
			_minpt = 0;
			_mineme = 0;
			_minnrhs = 0;
			_minrhE = 0.5;
			_maxrhE = -999;
			_minobjeta = 1.5;
			_year = 2018;
			_data = false;
			_calibmap = nullptr;
			_applyFrac = false;
			_spikes = false;
			_timesmear = false;
			//if(gSystem->AccessPathName("info/KUCMS_GJets_v14_met50_rhE5_Cali.root")){
			//	cout << "Calibration map file " << "info/KUCMS_GJets_v14_met50_rhE5_Cali.root" << " does not exist." << endl;
			//	return;
			//}
			//TFile* calibfile = TFile::Open("info/KUCMS_GJets_v14_met50_rhE5_Cali.root");
			//SetTimeCalibrationMap(calibfile);
			SetupDetIDsEB();
		};
		BaseProducer(TFile* file){
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			if(gSystem->AccessPathName(file->GetName())){ cout << "Error: file " << file->GetName() << " doesn't exist." << endl; return; }
			TTree* tree = (TTree*)file->Get("tree/llpgtree");
			_base = new ReducedBase(tree);
			_nEvts = _base->fChain->GetEntries();
			//default to 1 GeV = 1 entry -> gev = 1
			_gev = 1;
			_isocut = false;
			_minpt = 30;
			_mineme = 20;
			_minnrhs = 15;
			_minrhE = 0.5;
			_maxrhE = -999;
			_minobjeta = 1.5;
			_applyFrac = false;
			_spikes = false;
			_timesmear = false;
			
			//set year
			string name = file->GetName();
			if(name.find("2017") != string::npos) _year = 2017;
			else if(name.find("2018") != string::npos) _year = 2018;
			else if(name.find("2022") != string::npos) _year = 2022;
			
			//set if data
			if(name.find("SIM") == string::npos) _data = true;
			else _data = false;
			_calibmap = nullptr;
			//SetTimeCalibrationMap(calibfile);
			SetupDetIDsEB();
			
			if(name.find("_v20_")) useFilters = true;


		}
		virtual ~BaseProducer(){ 
			delete _base;
			delete _calibmap;	
		};


		bool useFilters = false;
		//returns vector of rec hits (as Jets) for each event (vector of vectors)
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt) = 0;
		virtual void GetRecHits(vector<Jet>& rhs, int evt){};
		virtual void GetSimRecHits(vector<Jet>& rhs, int evt){};
		virtual void GetRecHits(vector<JetPoint>& rhs, int evt, int obj) = 0;
		virtual void GetRecHits(vector<Jet>& rhs, int evt, int obj){};
		virtual void GetPrimaryVertex(BayesPoint& vtx, int evt) = 0;
		virtual void GetGenJets(vector<Jet>& genjets, int evt){}; 

		//amp = ampeff = amp/sigma
		double SmearRecHitTime(double amp, double time){
			double n = 14;
			double s = 0;
			double c = 0.1083;

			//sigma^2
			double sigma = (n/amp)*(n/amp) + (s*s)/(amp) + 2*c*c;
			sigma = sqrt(sigma/2); //calculated for 2 rhs, need for 1
			RandomSample rs;
			rs.SetRange(time-5*sigma,time+5*sigma);
			double newtime = rs.SampleGaussian(time,sigma,1).at(0);
			return newtime;
		}
		void SetTimeSmear(bool t){_timesmear = t;}
		bool _timesmear;
		void GetTrueJets(vector<Jet>& jets, int evt, double gev = -1);
		void GetTruePhotons(vector<Jet>& phos, int evt, double gev = -1);
		int GetTrueSuperClusters(vector<Jet>& phos, int evt, double gev = -1);

		bool _isocut;
		void SetIsoCut(){ _isocut = true; }		

		enum beamHaloFilter{
			notApplied = 0,
			applied = 1,
			invApplied = 2
		};
		beamHaloFilter _bh;
		void SetBeamHaloFilter(int bh){ 
			_bh = beamHaloFilter(bh); 
			if(_bh == notApplied) cout << "Not applying beam halo filter." << endl;
			else if(_bh == applied) cout << "Applying beam halo filter." << endl;
			else cout << "Applying inverse beam halo filter." << endl;
		} 

		ReducedBase* GetBase(){ return _base; }

		void PrintPreselection(){
			cout << "Default energy transfer factor: " << _gev << endl;
			cout << "Minimum pt: " << _minpt << endl;
			cout << "Minimum ECAL energy: " << _mineme << endl;
			cout << "Minimum object eta: " << _minobjeta << endl;
			cout << "Minimum rh (barrel only) energy: " << _minrhE << endl;
        		cout << "Minimum # of in-time rhs: " << _minnrhs << endl;
			cout << "Rechit time smear? " << _timesmear << endl;
		}


		ReducedBase* _base = nullptr;
		int _nEvts;

		//energy weight transfer factor
		//w = g*E -> g = N/GeV
		void SetTransferFactor(double g){ _gev = g; }
		double _gev;


		void SetMinPt(double p){ _minpt = p; }
		double _minpt;
		void SetMinNrhs(double p){ _minnrhs = p; }
		double _minnrhs;
		void SetMinEmE(double p){ _mineme = p; }
		double _mineme;
		void SetMinRhE(double r){ _minrhE = r; }
		double _minrhE;
		void SetMaxRhE(double r){ _maxrhE = r; }
		double _maxrhE; 
		void SetMinObjEta(double e){ _minobjeta = e; }
		double _minobjeta;
		void ApplyFractions(bool a){ _applyFrac = a; }
		bool _applyFrac;
		void RejectSpikes(bool s){ _spikes = s; }
		bool _spikes;

		double deltaR2(double e1, double p1, double e2, double p2){
			double de = e1 - e2;
			double dp = fabs(p1 - p2);
			double pi = acos(-1);
			if(dp > pi)
				dp -= 2*pi;
			return de*de + dp*dp;
		}


		double hypo(double x, double y, double z){
			return sqrt(x*x + y*y + z*z);
		}
		double _c = 29.9792458;	
		bool _data;
		int _year;

		TH2D* _calibmap;
		void SetTimeCalibrationMap(TFile* f){
			if(!f){ cout << "File for calibration map not set." << endl; return; }
			_calibmap = (TH2D*)f->Get("AveXtalRatioRecTimeEBMap");
		};
		
		double GetTimeCalibrationFactor(int ieta, int iphi){
			if(!_calibmap){ cout << "Calibration map not set." << endl; return -999; }
			if(iphi < 1 || iphi > 360){cout << "Invalid iphi: " << iphi << endl; return -999; }
			if(ieta < -85 || ieta > 85){cout << "Invalid ieta: " << ieta << endl; return -999; }
			return _calibmap->GetBinContent(ieta+86, iphi);
		};


		double GetTimeCalibrationFactor(unsigned int rhid){
			//transform from (rh) -> (ieta, iphi)
			unsigned int ieta = _detIDMap[rhid].i2;
			unsigned int iphi = _detIDMap[rhid].i1;
			return GetTimeCalibrationFactor(ieta, iphi);
		}
		struct DetIDStruct {
                        DetIDStruct() {}
                        DetIDStruct(const int ni1, const int ni2, const Int_t nTT, const Int_t & necal) : i1(ni1), i2(ni2), TT(nTT), ecal(necal){}
                        int i1;
                        int i2;
                        Int_t TT; 
                        Int_t ecal; 
                };


		map<UInt_t, DetIDStruct> _detIDMap;

		//this function and the corresponding DetIDStruct (above) are courtesy of Jack King 
		//https://github.com/jking79/LLPgammaAnalyzer/blob/master/macros/KUCMS_Skimmer/KUCMSHelperFunctions.hh	
		void SetupDetIDsEB(){
		    const std::string detIDConfigEB("info/fullinfo_detids_EB.txt");
		    std::ifstream infile( detIDConfigEB, std::ios::in);
		
		    UInt_t cmsswId, dbID;
		    int hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
		    std::string pos;
		
		    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM){
		        _detIDMap[cmsswId] = {iphi,ieta,TT25,0};
		    }
		}

};
#endif
