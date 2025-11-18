#ifndef CLUSTERANALYZER_HH
#define CLUSTERANALYZER_HH
#include "BayesCluster.hh"
#include "Jet.hh"
#include "TSystem.h"
#include <utility>
#include <string>
#include <math.h>
#include "fdeep/fdeep.hpp"

using std::make_pair;
using std::to_string;

constexpr static double _SOL = 29.9792458; //in cm/ns
//functions for calculating PV time, cluster observables, etc
struct ClusterObj{
	Jet _jet;
	BayesPoint _PV;
	int _idx;
	ClusterObj(){ 
		_idx = -1;
	};
	ClusterObj(Jet jet){ 
		_jet = jet; _PV = _jet.GetVertex(); _jet.SortConstituents();
		CalculateObjTimes(); 
	}

	void SetUserIndex(int i){ _idx = i; }
	int GetUserIndex(){ return _idx; }

        std::map<UInt_t,pair<int,int>> _detIDmap;
        void SetupDetIDs( std::map<UInt_t,pair<int,int>> DetIDMap){
		_detIDmap = DetIDMap;
	}

	vector<bool> _PUscores;
	//will fully remove PU clusters - ie set all w_nk = 0 for PU subcluster k
	void CleanOutPU(){
		Jet cleanedJet = _jet.CleanOutPU(_PUscores, false);
		_jet = cleanedJet;
		//recalculate observables
		CalculateObjTimes(); 
	}
	void CleanOutDetBkg(double minscore, int sigclass = 0, bool remove = false){
		//_jet.NeuralNetClean(_detBkgScores, sigclass, minscore, remove);
		//recalculate observables
		//CalculateObjTimes(); 
	} 


	float GetEtaCenter(){ return (float)_jet.eta(); }
	float GetPhiCenter(){ return (float)_jet.phi_std(); }
	float GetEtaVar(){ return (float)_jet.GetCovariance().at(0,0); }
	float GetPhiVar(){ return (float)_jet.GetCovariance().at(1,1); }
	float GetTimeVar(){ return (float)_jet.GetCovariance().at(2,2); }
	float GetEtaPhiCov(){ return (float)_jet.GetCovariance().at(0,1); }
	float GetEtaTimeCov(){ return (float)_jet.GetCovariance().at(0,2); }
	float GetPhiTimeCov(){ return (float)_jet.GetCovariance().at(1,2); }
	void GetMajMinLengths(double& majlen, double& minlen){
		majlen = -1;
		minlen = -1;

		Matrix space_mat(2,2);
		Matrix cov = _jet.GetCovariance();
		vector<double> eigvals;
		vector<Matrix> eigvecs; 
		Get2DMat(cov,space_mat);
		space_mat.eigenCalc(eigvals, eigvecs);
		majlen = sqrt(eigvals[1]);
		if(eigvals[0] < 0) minlen = -sqrt(-eigvals[0]);
		else minlen = sqrt(eigvals[0]);	
	}
	double GetEnergy(){ return _jet.E(); }	
	double GetPt(){ return _jet.pt(); }

	int GetNSubclusters(){ return _jet.GetNConstituents(); }
	void GetSubclusters(vector<Jet>& subcls){
		_jet.GetConstituents(subcls);
	}
	float GetSubclusterEtaCenter(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.eta(); 
	}
	float GetSubclusterPhiCenter(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.phi(); 
	}
	float GetSubclusterEtaVar(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.GetCovariance().at(0,0); 
	}
	float GetSubclusterPhiVar(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.GetCovariance().at(1,1); 
	}
	float GetSubclusterTimeVar(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.GetCovariance().at(2,2); 
	}
	float GetSubclusterEtaPhiCov(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.GetCovariance().at(0,1); 
	}
	float GetSubclusterEtaTimeCov(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.GetCovariance().at(0,2); 
	}
	float GetSubclusterPhiTimeCov(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.GetCovariance().at(1,2); 
	}
	float GetSubclusterEnergy(int k){
		Jet subcl;
		_jet.GetConstituent(k, subcl);
		return (float)subcl.E(); 
	}
	//per-object information
	double GetObjTime_PV(){ return _pvTime; }
	double GetObjTime_Det(){ return _detTime; }
	void GetRecHitWeights(map<unsigned int, float>& ws){
		ws.clear();
		Jet subcl;
		JetPoint rh, rrh;
		for(int r = 0; r < _jet.GetNRecHits(); r++){
			_jet.GetJetPointAt(r,rh);
			unsigned int id = rh.rhId();
			ws[id] = rh.GetWeight(); //weight = resp_rk*isPU_k
			//cout << "rh # " << r << " has weight " << rh.GetWeight() << " and id " << id << " for " << _jet.GetNConstituents() << " subclusters" << endl;
		}
	}

	//per-subcluster information - vector of length k for k subclusters
	//0 - PU, 1 - not PU
	void GetPUScores(vector<bool>& scores){
		scores.clear(); scores = _PUscores;
	}
	void GetDetBkgScores(vector<vector<float>>& detBkgScores){ detBkgScores.clear(); detBkgScores = _detBkgScores; }
	void GetPhotonIDScores(vector<float>& photonIDScores){ photonIDScores.clear(); photonIDScores = _photonIDScores; }
	void SetCNNModel(string model){ _detbkgmodel = fdeep::load_model(model,true,fdeep::dev_null_logger); }
	void SetDNNModel(string model){ _photonidmodel = fdeep::load_model(model,true,fdeep::dev_null_logger); }
	void SetCNNModel(fdeep::model model){ _detbkgmodel = model; }
	void SetDNNModel(fdeep::model model){ _photonidmodel = model; }
	

	void SetIsolationInfo(map<string, double> obs){
		_isomap.clear();
		_isomap = obs;
	}


	void CalculatePUScores(){
		_jet.CleanOutPU(_PUscores, true);
	}
	void CalculatePhotonIDScores(double pt, const map<string, double>& isomap){
		_photonIDScores.clear();
		map<string, double> nn_input_map;
		MakeDNNInputs(nn_input_map, pt, isomap);
		vector<float> input_sample;
		//put map into vector with same order as input features (which is the same order that the model sees)	
		for(int f = 0; f < _nnfeatures.size(); f++){
			input_sample.push_back((float)nn_input_map[_nnfeatures[f]]);
		}

		int size = input_sample.size();
		fdeep::tensor input_tensor = fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(size)), input_sample);
		
		//predict_class returns predicted class number and value of max output neuron
		fdeep::tensors result = _photonidmodel.predict({input_tensor});
		for(int i = 0; i < result.size(); i++){
			vector<float> reti = result[i].to_vector();
			for(int j = 0; j < reti.size(); j++)
				_photonIDScores.push_back(reti[j]);
		}
	}



	void CalculateDetBkgScores(bool pho){
		_detBkgScores.clear();
		fdeep::tensor_shape tensor_shape(_ngrid, _ngrid, 1);
		fdeep::tensor input_tensor(tensor_shape, 0.0f);
		//make grid for each subclusters in jet
		JetPoint center;
		GetCenterXtal(center);
		map<string, double> inputs;
		if(pho){
			MakeCNNGrid(-1, center, inputs);
			//transform grid to input_sample
			for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++){
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
					double val = inputs["CNNgrid_cell"+to_string(i)+"_"+to_string(j)];
					input_tensor.set(fdeep::tensor_pos(i+(_ngrid-1)/2, j+(_ngrid-1)/2, 0), val);
				}							
			}
			//predict_class returns predicted class number and value of max output neuron
			fdeep::tensors result = _detbkgmodel.predict({input_tensor});
                        for(int i = 0; i < result.size(); i++){
                                vector<float> reti = result[i].to_vector();
                                _detBkgScores.push_back(reti);
                        }
		}
		else{
			for(int k = 0; k < _jet.GetNConstituents(); k++){
				MakeCNNGrid(k, center, inputs);
				//transform grid to input_sample
				for(int i = -(_ngrid-1)/2; i < (_ngrid-1)/2+1; i++){
					for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
						double val = inputs["CNNgrid_cell"+to_string(i)+"_"+to_string(j)];
						input_tensor.set(fdeep::tensor_pos(i+(_ngrid-1)/2, j+(_ngrid-1)/2, 0), val);
					}							
				}
				//predict_class returns predicted class number and value of max output neuron
				fdeep::tensors result = _detbkgmodel.predict({input_tensor});
                        	for(int i = 0; i < result.size(); i++){
                        	        vector<float> reti = result[i].to_vector();
                                	_detBkgScores.push_back(reti);
                        	}
			}
		}
	}
	
	private:
		map<string, double> _isomap;
		double _pvTime = -999; //time at PV
		double _detTime = -999;
		void CalculateObjTimes(){ 
			//obj time = sum_i w_i * t_i / w_i
			double t = 0;
			double norm = 0;
			PointCollection locs;
			Jet subcl;
			JetPoint jrh;
			for(int r = 0; r < _jet.GetNPoints(); r++){
				_jet.GetJetPointAt(r,jrh);
				double w = jrh.GetWeight()*jrh.E(); //weight = r_nk*isPU_k
				BayesPoint loc({jrh.x(), jrh.y(), jrh.z(), jrh.t()}); //time is at detector face
				loc.SetWeight(w);
				locs += loc;	
			}
			_detTime = locs.Centroid(3);
			//include geo factor to PV for PV time
			BayesPoint center({locs.Centroid(0), locs.Centroid(1), locs.Centroid(2)});
			//
			double dx = center.at(0) - _PV.at(0);
			double dy = center.at(1) - _PV.at(1);
			double dz = center.at(2) - _PV.at(2);
			double d_rh_pv = sqrt(dx*dx + dy*dy + dz*dz)/_SOL;
			_pvTime = _detTime - d_rh_pv;
		}

		vector<vector<float>> _detBkgScores; //is a vector of vectors s.t. _detBkgScores[subcl_idx] = {score_physbkg, score_bh, score_spike} 
		vector<float> _photonIDScores; //is a vector of vectors s.t. _photonIDScores[subcl_idx] = {score_isobkg, score_nonisobkg} 
		fdeep::model _detbkgmodel = fdeep::load_model("config/json/small3CNN_EMultr_2017and2018.json",true,fdeep::dev_null_logger);
		fdeep::model _photonidmodel = fdeep::load_model("config/json/med16DNN_MCtrained_photonID.json",true,fdeep::dev_null_logger);
		int _ngrid = 7;
		vector<string> _nnfeatures = {"EtaSig", "PhiSig", "EtaPhiCov", "majorLength", "minorLength", "hcalTowerSumEtConeDR04OvPt", "trkSumPtSolidConeDR04OvPt", "trkSumPtHollowConeDR04OvPt", "hadTowOverEMOvPt", "ecalRHSumEtConeDR04OvPt"};
		
		void MakeDNNInputs(map<string, double>& obs, double phopt, const map<string, double>& isomap){
			obs.clear();
			//get shape variables
			Matrix cov = _jet.GetCovariance();
			obs["EtaSig"] = sqrt(cov.at(0,0));
			obs["PhiSig"] = sqrt(cov.at(1,1));
			obs["EtaPhiCov"] = cov.at(0,1);
			double majlen, minlen;
			GetSpaceMajMinLens(minlen, majlen);
			obs["majorLength"] = majlen;
			obs["minorLength"] = minlen;

			//get photon isolation
			obs["hcalTowerSumEtConeDR04OvPt"] = isomap.at("hcalTowerSumEtConeDR04") / phopt;
			obs["trkSumPtSolidConeDR04OvPt"] = isomap.at("trkSumPtSolidConeDR04") / phopt;
			obs["trkSumPtHollowConeDR04OvPt"] = isomap.at("trkSumPtHollowConeDR04") / phopt;
			obs["hadTowOverEMOvPt"] = isomap.at("hadTowOverEM") / phopt;
			obs["ecalRHSumEtConeDR04OvPt"] = isomap.at("hadTowOverEM") / phopt;
			
		}


		void GetCenterXtal(JetPoint& center){
			//get center of pts in ieta, iphi -> max E point
			vector<JetPoint> rhs;
			_jet.GetJetPoints(rhs);
			double maxE = 0;
			for(int r = 0; r < rhs.size(); r++){
				if(rhs[r].E() > maxE){
					maxE = rhs[r].E();
					center = rhs[r];
				}
			}
		}
		void MakeCNNGrid(int k, JetPoint& center, map<string,double>& mapobs){
			if(_detIDmap.size() == 0){
				cout << "ERROR: detIDmap not setup for this ClusterObj. Please use SetupDetIDs( std::map<UInt_t,pair<int,int>> DetIDMap) to set the detID map for the CNN grid creation. Not running CNN detector bkg classification." << endl;
				return;
			}
			mapobs.clear();
			map<pair<int,int>, vector<double>> grid;
			//make sure ngrid is odd to include center crystal
			if(_ngrid % 2 == 0)
				_ngrid++;

			int ngrid_boundary = (_ngrid-1)/2;
			//set default channel values to 0
			//{E, t, r, Er}
			for(int i = -ngrid_boundary; i < ngrid_boundary+1; i++)
				for(int j = -ngrid_boundary; j < ngrid_boundary+1; j++)
					grid[make_pair(i,j)] = {0., 0., 0.};


			//get ngrid x ngrid around center point 
			int ieta, iphi;
			int rh_ieta = _detIDmap.at(center.rhId()).first;
			int rh_iphi = _detIDmap.at(center.rhId()).second;
			int deta, dphi;
			Jet subcl;
			if(k == -1) //do over whole object
				subcl = _jet;
			else
				_jet.GetConstituent(k, subcl);
			vector<JetPoint> rhs;
			subcl.GetJetPoints(rhs);
			for(int j = 0; j < rhs.size(); j++){
				ieta = _detIDmap.at(rhs[j].rhId()).first;
				iphi = _detIDmap.at(rhs[j].rhId()).second;
				//do eta flip
				if(rh_ieta < 0)
					deta = -(ieta - rh_ieta);
				else
					deta = ieta - rh_ieta;
				dphi = iphi - rh_iphi; 
				//needs wraparound
				if(dphi > 180)
					dphi = 360 - dphi;
				else if(dphi < -180)
					dphi = -(360 + dphi);
				if(fabs(deta) <= ngrid_boundary && fabs(dphi) <= ngrid_boundary){
					//posterior is weighted s.t. sum_k post_nk = w_n = E*_gev, need to just have unweighted probs since E is already here
					//r_nk = post_nk/w_n s.t. sum_k (post_nk/w_n) = w_n/w_n = 1
					grid[make_pair(deta, dphi)] = {rhs[j].E()}; //subcluster jets already have rh energies "projected in" to their subclusters

				}
			}
			pair<int, int> icoords_grid;
			for(int i = -(_ngrid-1)/2.; i < (_ngrid-1)/2+1; i++){
				for(int j = -(_ngrid-1)/2; j < (_ngrid-1)/2+1; j++){
					icoords_grid = make_pair(i,j);
					mapobs["CNNgrid_cell"+to_string(i)+"_"+to_string(j)] = grid[icoords_grid][0];
				}
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
		void GetSpaceMajMinLens(double& minlen, double& majlen){
			Matrix cov = _jet.GetCovariance();
			Matrix space_mat(2,2);
			Get2DMat(cov,space_mat);
			vector<Matrix> eigenvecs_space;
			vector<double> eigenvals_space;
			space_mat.eigenCalc(eigenvals_space, eigenvecs_space);
			majlen = sqrt(eigenvals_space[1]);
			if(eigenvals_space[1] < 0) cout << "negative eigenvalue " << eigenvals_space[1] << endl;
			if(eigenvals_space[0] < 0) minlen = -sqrt(-eigenvals_space[0]);
			else minlen = sqrt(eigenvals_space[0]);	

		}



};
//wrapper class for clustering
class ClusterAnalyzer{
	public:
		ClusterAnalyzer();
		virtual ~ClusterAnalyzer(){
			delete _algo;
		};

		//needs PV info for geometric corrections and correct momentum calculations of clustered elements
		void SetPV(double pvx, double pvy, double pvz){
			vector<double> pv = {pvx, pvy, pvz};
			_PV = BayesPoint(pv);
		}
		void SetDetectorCenter(double x, double y, double z){ _detCenter = BayesPoint({x,y,z});}
		//sets transfer factor for energy weighting in clustering
		void SetTransferFactor(double g){ _gev = g; }
		//adding data to clustering algorithm
		void AddRecHit(double rhx, double rhy, double rhz, double rhE, double rht, int rhId, bool invalidTime = false);
		int GetNRecHits(){ return _rhs.size(); }
		void ClearRecHitList();

		void SetVerbosity(int v){ _verb = v; }

		int RunClustering(ClusterObj& retobj, bool pho = true);
		int NoClusterRhs(ClusterObj& retobj, bool pho = true);

		void SetDetIDs(std::map<UInt_t,pair<int,int>> detidmap){
			_detIDmap.clear();
			_detIDmap = detidmap;
		}
		void SetCNNModel(string model){ 
			if(_verb > -1) cout << "Using model " << model << " for instrumental background." << endl; 
			_detbkgmodel = fdeep::load_model(model); }
		void SetDNNModel(string model){
			if(_verb > -1) cout << "Using model " << model << " for photon ID." << endl; 
			_photonidmodel = fdeep::load_model(model); }

	private:
		vector<Jet> _rhs;
		double _gev; 
		double _radius; //detector radius in m
		BayesPoint _PV;
		BayesPoint _detCenter;
		BayesCluster* _algo;
		void _treesToObjs(vector<node*>& trees, vector<ClusterObj>& objs);
		std::map<UInt_t,pair<int,int>> _detIDmap;
		fdeep::model _detbkgmodel = fdeep::load_model("config/json/small3CNN_EMultr_2017and2018.json",true,fdeep::dev_null_logger);
		fdeep::model _photonidmodel = fdeep::load_model("config/json/med16DNN_MCtrained_photonID.json",true,fdeep::dev_null_logger);
		int _verb;

};
static bool Esort(ClusterObj j1, ClusterObj j2){ return (j1._jet.E() > j2._jet.E()); }
static bool ptsort(ClusterObj j1, ClusterObj j2){ return (j1._jet.pt() > j2._jet.pt()); }
#endif
