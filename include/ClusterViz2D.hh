#ifndef CLUSTERVIZ2D_HH
#define	CLUSTERVIZ2D_HH

#include "ClusterVizBase.hh"
#include "GaussianMixture.hh"
#include <TGraph2D.h>
#include <string.h>

using std::string;

class ClusterViz2D : public ClusterVizBase{
	public:
		ClusterViz2D(){ };
		ClusterViz2D(GaussianMixture* model, string fname = "test");
		virtual ~ClusterViz2D(){ };
		
		void AddPlot(string plotName = "test");
		void Write();
		void SetPalette(int k);

	private:
		string m_fname;
		Int_t m_palette[100];
		vector<TCanvas*> m_cvs; 
};

#endif
