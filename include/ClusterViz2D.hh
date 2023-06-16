#ifndef CLUSTERVIZ2D_HH
#define	CLUSTERVIZ2D_HH

#include "ClusterVizBase.hh"
#include "BasePDFMixture.hh"
#include <TGraph2D.h>
#include <string.h>

using std::string;

class ClusterViz2D{
	public:
		ClusterViz2D(){ };
		ClusterViz2D(BasePDFMixture* model, string fname = "test");
		virtual ~ClusterViz2D(){ };
		
		void AddPlot(string plotName = "test");
		void UpdatePosterior(){
			m_post = m_model->GetPosterior();
		}
		void Write();
		void SetPalette(int k);

	private:
		BasePDFMixture* m_model;
		Matrix m_post;
		PointCollection m_points;	
		int m_n; //number of points
		int m_k; //number of clusters
		string m_fname;	
		vector<TCanvas*> m_cvs; 
		Int_t m_palette[100];
};

#endif
