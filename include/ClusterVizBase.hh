#ifndef CLUSTERVIZBASE_HH
#define CLUSTERVIZBASE_HH

#include <TGraph.h>
#include <TCanvas.h>
#include "BasePDFMixture.hh"
//#include "GaussianMixture.hh"
#include "BaseCluster.hh"
#include "PointCollection.hh"
#include <string>
#include <vector>

using std::vector;
using std::string;

class ClusterVizBase{
	public:
		ClusterVizBase(){ 
			_model = nullptr;
			_post = Matrix();
			_points = new PointCollection();	
			_n = 0; //number of points
			_k = 0; //number of clusters
			_fname = "";	
			_cvs = {}; 
			_verb = 0;
		};
		ClusterVizBase(BaseCluster* algo, string fname = "test"){
			_model = algo->GetModel();
			_fname = fname;
			_points = algo->GetData();
			_n = _points->GetNPoints();
			_k = algo->GetNClusters();
			_post = _model->GetPosterior();
			_verb = 0;
		};
		virtual ~ClusterVizBase(){ 
		//	cvs.clear();
		//	delete _model;
		};
	
		void UpdatePosterior(){
			_post = _model->GetPosterior();
		
		}	
		virtual void AddPlot(string plotName = "test") = 0;
		virtual void Write() = 0;
		virtual void SeeData() = 0;	

		void SetVerbosity(int v){ _verb = v; }
		//for calculating energies s.t. E_n = w_n*k
		void SetTransfFactor(double k){ _transf = k; }
	
	//protected:
		BasePDFMixture* _model;
		Matrix _post;
		PointCollection* _points;	
		int _n; //number of points
		int _k; //number of clusters
		string _fname;	
		vector<TCanvas*> _cvs;
		int _verb; 
		double _transf;




};

#endif
