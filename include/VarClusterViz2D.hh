#ifndef VIZCLUSTERVIZ2D_HH
#define	VIZCLUSTERVIZ2D_HH

#include "BasePDFMixture.hh"
#include "VarEMCluster.hh"
#include <TCanvas.h>
#include <string.h>

using std::string;

class VarClusterViz2D{
	public:
		VarClusterViz2D(){
			m_model = nullptr;
			m_post = Matrix();
			m_points = new PointCollection();	
			m_n = 0; //number of points
			m_k = 0; //number of clusters
			m_fname = "";	
			m_cvs = {}; 
		}
		VarClusterViz2D(const VarClusterViz2D& viz);
		VarClusterViz2D(VarEMCluster* algo, string fname = "test");
		//VarClusterViz2D(VarGaussianMixture* model, string fname = "test");
		virtual ~VarClusterViz2D(){ };
		
		void AddPlot(string plotName = "test");
		void UpdatePosterior(){  
			if(m_n == 0){
				return;
			}
			m_post = m_algo->GetPosterior();
		}	
		void Write();
		void SetPalette(int k);

		BasePDFMixture* GetModel(){ return m_model; }
		string GetFileName(){ return m_fname; }

	private:
		BaseCluster* m_algo;
		BasePDFMixture* m_model;
		Matrix m_post;
		PointCollection* m_points;	
		int m_n; //number of points
		int m_k; //number of clusters
		string m_fname;	
		vector<TCanvas*> m_cvs; 
		Int_t m_palette[100];
};

#endif
