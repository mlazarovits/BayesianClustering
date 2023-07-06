#ifndef VIZCLUSTERVIZ3D_HH
#define	VIZCLUSTERVIZ3D_HH

#include "VarGaussianMixture.hh"
#include <TCanvas.h>
#include <string.h>

using std::string;

class VarClusterViz3D{
	public:
		VarClusterViz3D(){
			m_model = nullptr;
			m_post = Matrix();
			m_points = PointCollection();	
			m_n = 0; //number of points
			m_k = 0; //number of clusters
			m_fname = "";	
			m_cvs = {}; 
		}
		VarClusterViz3D(const VarClusterViz3D& viz);
		VarClusterViz3D(VarGaussianMixture* model, string fname = "test");
		//VarClusterViz3D(VarGaussianMixture* model, string fname = "test");
		virtual ~VarClusterViz3D(){ };
	
		//add plot at specified t (z) value	
		//time window?
		void AddPlot(double t, string plotName = "test");
		void AddAnimationPlots(string dirName = "test");
		void UpdatePosterior(){  
			if(m_n == 0){
				return;
			}
			m_post = m_model->GetPosterior();
		}	
		void Write();
		void SetPalette(int k);

		VarGaussianMixture* GetModel(){ return m_model; }
		string GetFileName(){ return m_fname; }

	private:
		VarGaussianMixture* m_model;
		Matrix m_post;
		PointCollection m_points;	
		int m_n; //number of points
		int m_k; //number of clusters
		double m_deltaT; //time slices for individual plots
		string m_fname;	
		vector<TCanvas*> m_cvs; 
		Int_t m_palette[100];
};

#endif
