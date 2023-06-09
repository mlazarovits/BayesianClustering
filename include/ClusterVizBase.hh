#ifndef CLUSTERVIZBASE_HH
#define CLUSTERVIZBASE_HH

#include <TGraph.h>
#include <TCanvas.h>
#include "GaussianMixture.hh"
#include "PointCollection.hh"
#include <string>
#include <vector>

using std::vector;
using std::string;

class ClusterVizBase{
	public:
		ClusterVizBase(){ 
		};
		ClusterVizBase(GaussianMixture* model, string fname = "test"){
			m_model = model;
			m_fname = fname;
			m_points = m_model->GetData();
			m_n = m_points.GetNPoints();
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
		
		//virtual void DrawGausParams() = 0;

	
	protected:
		GaussianMixture* m_model;
		Matrix m_post;
		PointCollection m_points;	
		int m_n; //number of points
		string m_fname;	
		vector<TCanvas*> m_cvs; 




};

#endif
