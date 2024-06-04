#ifndef VIZCLUSTERVIZ3D_HH
#define	VIZCLUSTERVIZ3D_HH

#include "ClusterVizBase.hh"
#include "VarEMCluster.hh"
#include <TCanvas.h>
#include <string.h>

using std::string;

class VarClusterViz3D : public ClusterVizBase{
	public:
		VarClusterViz3D(){
			_model = nullptr;
			_post = Matrix();
			_points = new PointCollection();	
			_n = 0; //number of points
			_k = 0; //number of clusters
			_fname = "";	
			_cvs = {}; 
		}
		VarClusterViz3D(const VarClusterViz3D& viz);
		VarClusterViz3D(VarEMCluster* algo, string fname = "test");
		virtual ~VarClusterViz3D(){ };

		void Write(){ };
		void AddPlot(string filename = "test"){ };
	
		void WriteJson(string filename = "test");
		void UpdatePosterior(){  
			if(_n == 0){
				return;
			}
			_post = _model->GetPosterior();
			_k = _model->GetNClusters();
		}	
		void SeeData();

		BasePDFMixture* GetModel(){ return _model; }
		string GetFileName(){ return _fname; }

	private:
		BayesPoint _scale;
		BayesPoint _shift;
};

#endif
