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
			m_model = nullptr;
			m_post = Matrix();
			m_points = new PointCollection();	
			m_n = 0; //number of points
			m_k = 0; //number of clusters
			m_fname = "";	
			m_cvs = {}; 
		};
		ClusterVizBase(BaseCluster* algo, string fname = "test"){
			m_model = algo->GetModel();
			m_fname = fname;
			m_points = algo->GetData();
			m_n = m_points->GetNPoints();
			m_k = algo->GetNClusters();
			m_post = m_model->GetPosterior();
			m_algo = algo;
		};
		virtual ~ClusterVizBase(){ 
		//	cvs.clear();
		//	delete m_model;
		};
	
		void UpdatePosterior(){
			m_post = m_model->GetPosterior();
		
		}	
		virtual void AddPlot(string plotName = "test") = 0;
		virtual void Write() = 0;
		virtual void SeeData() = 0;	

	
	//protected:
		BaseCluster* m_algo;
		Int_t m_palette[100];
		BasePDFMixture* m_model;
		Matrix m_post;
		PointCollection* m_points;	
		int m_n; //number of points
		int m_k; //number of clusters
		string m_fname;	
		vector<TCanvas*> m_cvs; 




};

#endif
