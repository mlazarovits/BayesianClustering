#ifndef CLUSTERVIZ2D_HH
#define	CLUSTERVIZ2D_HH

#include "ClusterVizBase.hh"
#include "BasePDFMixture.hh"
#include "EMCluster.hh"
#include <TGraph2D.h>
#include <string.h>

using std::string;

class ClusterViz2D : public ClusterVizBase{
	public:
		ClusterViz2D(){ };
		ClusterViz2D(EMCluster* algo, string fname = "test");
		virtual ~ClusterViz2D(){ };
		
		void AddPlot(string plotName = "test");
		void SeeData();
		void Write();
		void SetPalette(int k);

		Int_t _palette[100];
};

#endif
