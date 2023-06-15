#ifndef CLUSTERVIZ2D_HH
#define	CLUSTERVIZ2D_HH

#include "ClusterVizBase.hh"
#include "BasePDFMixture.hh"
#include <TGraph2D.h>
#include <string.h>

using std::string;

class ClusterViz2D : public ClusterVizBase{
	public:
		ClusterViz2D(){ };
		ClusterViz2D(BasePDFMixture* model, string fname = "test");
		//ClusterViz2D(VarGaussianMixture* model, string fname = "test");
		virtual ~ClusterViz2D(){ };
		
		void AddPlot(string plotName = "test");
		void Write();
		void SetPalette(int k);

	private:
		Int_t m_palette[100];
};

#endif
