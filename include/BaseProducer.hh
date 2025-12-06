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
#include "KUCMS_TimeCalibration.hh"
#include <fstream>

class BaseProducer{
	public:
		BaseProducer() :
			_gev(1),
			_isocut(false),
			_minpt(0),
			_mineme(0),
			_minnrhs(0),
			_minrhE(0.5),
			_maxrhE(-999),
			_minobjeta(3.),
			_year(2018),
			_data(false),
			_calib(true),
			_applyFrac(false),
			_spikes(false),
			_timesmear(false),
			_spatial_corr(true),
			_mctype(-1),
			_mistcuts(false),
			_nEvts(-1),
			//pointers
			_timecalibTool(make_unique<KUCMS_TimeCalibration>()),
			_base(nullptr),
			_ch(nullptr)
		{  };
		//BaseProducer(TTree* tree) :
		BaseProducer(string file) :
			//default to 1 GeV = 1 entry -> gev = 1
			_gev(1),
			_isocut(false),
			_minpt(30),
			_mineme(30),
			_minnrhs(15),
			_minrhE(0.5),
			_maxrhE(-999),
			_minobjeta(3.),
			_applyFrac(false),
			_spikes(false),
			_timesmear(false),
			_spatial_corr(true),
			_mctype(-1),
			_mistcuts(false),
			_calib(true),
			//pointers
			_timecalibTool(make_unique<KUCMS_TimeCalibration>())
		{
			MakeBase(file);	
			_nEvts = _base->fChain->GetEntries();
			
			//set year
			string name = file;
			if(name.find("2016") != string::npos) _year = 2016;
			else if(name.find("2017") != string::npos) _year = 2017;
			else if(name.find("2018") != string::npos) _year = 2018;
			else if(name.find("2022") != string::npos) _year = 2022;
			
			//cout << "name " << name << endl;
			if(name.find("SIM") == string::npos) _data = true;
			else _data = false;
			
			//set if data
			cout << "name " << name << " year " << _year << " data " << _data << endl;
			
			if(_year == 2016)
				_timecalibTag = "r2_ul18"; //use 2017 time reco tag
			else if(_year == 2017)
				_timecalibTag = "r2_ul17";
			else if(_year == 2018)
				_timecalibTag = "r2_ul18";
			else if(_year == 2022)
				_timecalibTag = "r2_ul18"; //use 2018 time reco tag
			else
				cout << "Time calibration tag not set for year " << _year << endl;
			if(_data){
				_mctype = 1;
			}
			else{
				_mctype = 0;
			}
			if(!_data){
				if(_year == 2018)
					_timecalibTag += "_mc";
				else
					cout << "Time calibration tag not set for " << _year << " MC." << endl;
			}	
			cout << "Using time calibration + smearing tag " << _timecalibTag << endl;

		}

		BaseProducer(TFile* file) : 
			//default to 1 GeV = 1 entry -> gev = 1
			_gev(1),
			_isocut(false),
			_minpt(30),
			_mineme(30),
			_minnrhs(15),
			_minrhE(0.5),
			_maxrhE(-999),
			_minobjeta(1.5),
			_applyFrac(false),
			_spikes(false),
			_timesmear(false),
			_spatial_corr(true),
			_mctype(-1),
			_mistcuts(false),
			_calib(true),
			//pointers
			_timecalibTool(make_unique<KUCMS_TimeCalibration>())
		{	
			//jack does rh_adjusted_time = rh_time - (d_rh - d_pv)/c = rh_time - d_rh/c + d_pv/c
			//tof = (d_rh-d_pv)/c
			//in ntuplizer, stored as rh time

			//grab rec hit values
			//x, y, z, time (adjusted), energy, phi, eta
			if(gSystem->AccessPathName(file->GetName())){ cout << "Error: file " << file->GetName() << " doesn't exist." << endl; return; }
			TTree* tree = (TTree*)file->Get("tree/llpgtree");
			_base = new ReducedBase(tree);
			_nEvts = _base->fChain->GetEntries();
			//set year
			string name = file->GetName();
			if(name.find("2017") != string::npos) _year = 2017;
			else if(name.find("2018") != string::npos) _year = 2018;
			else if(name.find("2022") != string::npos) _year = 2022;
			
			//set if data
			if(name.find("SIM") == string::npos) _data = true;
			else _data = false;

			if(_year == 2018)
				_timecalibTag = "r2_ul18";
			else if(_year == 2017)
				_timecalibTag = "r2_ul17";
			else
				cout << "Time calibration tag not set for year " << _year << endl;
			if(_data){
				_mctype = 1;
			}
			else{
				_mctype = 0;
			}
			if(!_data){
				if(_year == 2018)
					_timecalibTag += "_mc";
				else
					cout << "Time calibration tag not set for " << _year << " MC." << endl;
			}	
			cout << "Using time calibration + smearing tag " << _timecalibTag << endl;
		}
		virtual ~BaseProducer(){ 
			delete _base;
		};
		void MakeBase(string f){
			if(gSystem->AccessPathName(f.c_str())){ 
				cout << "Error: file " << f << " doesn't exist." << endl; 
			}
			_ch = make_unique<TChain>("tree/llpgtree");
			_ch->SetTitle(f.c_str());
			_ch->Add(f.c_str()); //skip non-recoverable files
			_base = new ReducedBase(_ch.get());
		}

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
			_timecalibTool->setSmearTag(_timecalibTag);
			return _timecalibTool->getSmearedTime(time, amp);
		}
		void SetTimeSmear(bool t){_timesmear = t;}
		bool _timesmear;
		void SetMistClean(bool m){ _mistcuts = m; }
		bool _mistcuts; //true = cut out mist, false = dont cut out mist
		void GetTrueJets(vector<Jet>& jets, int evt, double gev = -1);
		void GetTruePhotons(vector<Jet>& phos, int evt, double gev = -1);
		int GetTrueSuperClusters(vector<Jet>& phos, int evt, double gev = -1);

		bool _isocut;
		void SetIsoCut(){ _isocut = true; }		

		bool _spatial_corr;
		void SetSpatialCorr(bool b){ _spatial_corr = b; }

		ReducedBase* GetBase(){ return _base; }

		unique_ptr<KUCMS_TimeCalibration> _timecalibTool;
		unique_ptr<TChain> _ch;
		ReducedBase* _base; //needs to be raw because ROOT
		string _timecalibTag;
		int _mctype;

		void PrintPreselection(){
			cout << "Default energy transfer factor: " << _gev << endl;
			cout << "Minimum pt: " << _minpt << endl;
			cout << "Minimum ECAL energy: " << _mineme << endl;
			cout << "Minimum rh (barrel only) energy: " << _minrhE << endl;
			if(_maxrhE != -999) cout << "Maximum rh (barrel only) energy: " << _maxrhE << endl;
        		cout << "Minimum # of in-time rhs: " << _minnrhs << endl;
			cout << "Rechit time smear? " << _timesmear << endl;
			cout << "Spatially correcting times to PV location? " << _spatial_corr << endl;
		}

		map<unsigned int, double> _rhIdToRes;
		void GetRhIdResMap(map<unsigned int, double>& omap){
			omap.clear();
			omap = _rhIdToRes;
		}

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

		bool _calib;
		void SetTimeCalibration(bool c){
			_calib = c;
		};
		
		double GetTimeCalibrationFactor(unsigned int rhid, int run){
			return _timecalibTool->getCalibration(rhid, run, _timecalibTag);
		}
		struct DetIDStruct {
                        DetIDStruct() {}
                        DetIDStruct(const int ni1, const int ni2, const Int_t nTT, const Int_t & necal) : i1(ni1), i2(ni2), TT(nTT), ecal(necal){}
                        int i1;
                        int i2;
                        Int_t TT; 
                        Int_t ecal; 
                };



};
#endif
