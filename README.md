# BayesianClustering
Repository for generic EM/hierarchical clustering algorithm (to be used for jet and photon clustering at CMS).


### Model Initialization
- means of Gaussians are initialized via K-means while the covariance matrices are initialized to the identity matrix
- if you are having trouble with the variational GMM, play around with the initial values of the parameters
	- it helps if $\beta_k$ and $\nu_k$ are the same order of magnitude 
- default initial values of the Gaussian priors are set to:
	- $\beta_0 = 0.001$
	- $\vec{m}_0 = \vec{0}$
	- $\nu_0 = (D - 1) + 0.001$
	- $\W_0 = \mathbb{I}$

### Formatting
There are muliple visualization classes:
- `ClusterViz2D` produces a root file output for the vanilla EM algorithm
- `VarClusterViz2D` produces a root file output for the variational EM algorithm
- `VarClusterViz3D` produces a json output for the variational EM algorithm
	- this is used in the standalone variational EM algorithm and the full algorithm
	- run either `macros/Viz3DVarEM.py` (var EM) or `macros/Viz3DFullCluster.py` (full algo) over the json output
	- for var EM, one json (plot) per iteration is produced into a directory
	- for the full algorithm, one json is produced, with individual plots corresponding descending levels in the binary agglomeration tree


### (Jet) Clustering (Physics-specific)
- `JetProducer` either grabs rechits from specified ROOT file (NANOAOD/Ntuple) or simulates rechits to cluster
- `Jet` is 1 cluster unit (it can be anything, doesn't have to be a physics jet)
	- it can be composite (made up of multiple Jets) or singular (ie one rechit)
	- these are the units that are clustered in `JetCluster` 
- `Cluster` calls the specified model (either Gaussian Mixture or Bayesian Gaussian Mixture) and runs the corresponding EM algorithm to find the cluster parameters
	- TODO: eventually will introduct Bayesian Hierarchical clustering to combine jets into cluster hierarchy and use varGMM EM to re-estimate clusters at each step


