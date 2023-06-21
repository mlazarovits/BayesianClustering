# BayesianClustering
Repository for generic EM/hierarchical clustering algorithm (to be used for jet clustering at CMS).


### Model Initialization
- means of Gaussians are randomly initialized while the covariance matrices are initialized to the identity matrix
- if you are having trouble with the variational GMM, play around with the initial values of the parameters
	- it helps if $\beta_k$ and $\nu_k$ are the same order of magnitude 

### Formatting
Plot formatting is done in `ClusterViz2D` but also done in `rootlogon.C` to correctly render colors when loading through the ROOT interpreter. Make sure that these formats are synced and are using the correct palette for the appropriate number of clusters in the GMM.

To do:
- `GaussianMixture.cc`
	- make sure there is data set before Initialize()
- set jet application in it's own folder or repository
