#ifndef VIZCLUSTERVIZ2D_HH
#define	VIZCLUSTERVIZ2D_HH

#include "VarGaussianMixture.hh"
#include <TCanvas.h>
#include <string.h>

using std::string;

class VarClusterViz2D{
	public:
		VarClusterViz2D(){ };
		VarClusterViz2D(VarGaussianMixture* model, string fname = "test");
		//VarClusterViz2D(VarGaussianMixture* model, string fname = "test");
		virtual ~VarClusterViz2D(){ };
		
		void AddPlot(string plotName = "test");
		void UpdatePosterior(){
			m_post = m_model->GetPosterior();
		}	
		void Write();
		void SetPalette(int k);

	private:
		VarGaussianMixture* m_model;
		Matrix m_post;
		PointCollection m_points;	
		int m_n; //number of points
		int m_k; //number of clusters
		string m_fname;	
		vector<TCanvas*> m_cvs; 
		Int_t m_palette[100];
};

#endif
