# BayesianClustering
Repository for generic EM/hierarchical clustering algorithm (to be used for jet clustering at CMS).


### Formatting
Plot formatting is done in `ClusterViz2D` but also done in `rootlogon.C` to correctly render colors when loading through the ROOT interpreter. Make sure that these formats are synced and are using the correct palette for the appropriate number of clusters in the GMM.

To do:
- `Matrix.cc`
	- make sure all dims match for mult, add
	- if dims aren't set, set them; otherwise make sure they match
- `GaussianMixture.cc`
	- make sure there is data set before Initialize()
