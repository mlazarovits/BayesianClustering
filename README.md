# BayesianClustering
Repository for generic EM/hierarchical clustering algorithm (to be used for jet and photon clustering at CMS).

### Compiling
- to run on LPC: `make lpc`
- to run on local machine:
	- make sure paths in Makefile are updated accordingly
	- `make`

### Dependencies
- ROOT for producers and skimmers
- Python for macros for visualization
	- at least v3.x
- [eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [boost](https://www.boost.org/doc/libs/1_82_0/libs/math/doc/html/special.html)
- [nlohmman-json](https://github.com/nlohmann/json) to make jsons for python macros
- [cgal](https://www.cgal.org/) for Voronoi diagrams
	- on LPC: make sure the path to the CGAL lib is added to `$LD_LIBRARY_PATH`
		- this can be done within the user's `~/.bash_profile` by adding the following line
		```
		export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib
		```
- [fastjet](https://fastjet.fr/) for detector simulation
	- at least v3.4.1 for compatibility with C++17
	- at least v3.4.2 for compatibility with C++20 

### Model Initialization
- means of Gaussians are initialized via K-means while the covariance matrices are initialized to the identity matrix
	- this is also the case for the variational EM algorithm, except once the prior parameters are set, the parameters are updated to seed the algorithm (initial M0-step, then alternate between E-M)
- if you are having trouble with the variational GMM, play around with the initial values of the parameters
	- it helps if $\beta_k$ and $\nu_k$ are the same order of magnitude 
- default initial values of the Gaussian prior (Normal-Wishart) are set to the least informative priors (Bishop and [Breheny](https://myweb.uiowa.edu/pbreheny/uk/teaching/701/notes/3-28.pdf)):
	- $\beta_0 = 0.001$
	- $\vec{m}_0 = \vec{0}$
	- $\nu_0 = D$  for $D$ dimensional data
	- $W_0 = \nu_0^{-1}\mathbb{I}$
- any data smearing needs to be set before the data is set
	- this is especially true for the BHC model because the leaf $r_k$ values are calculated when the data is set
- the algorithm will remove subclusters from a mixture model if the number of effective points ($\alpha + N_k$) falls below a specified threshhold `algo->SetThreshold(t)`
	- the mininum number of subclusters in a mixture model though is 1
	- if the algorithm removes all subclusters, it will set the number of subclusters to 1

### Formatting
There are muliple visualization classes:
- `ClusterViz2D` produces a root file output for the vanilla EM algorithm
- `VarClusterViz2D` produces a root file output for the variational EM algorithm
- `VarClusterViz3D` produces a json output for the variational EM algorithm
	- this is used in the standalone variational EM algorithm and the full algorithm
	- run either `macros/Viz3DVarEM.py` (var EM) or `macros/Viz3DFullCluster.py` (full algo) over the json output
	- for var EM, one json (plot) per iteration is produced into a directory
	- for the full algorithm, one json is produced, with individual plots corresponding descending levels in the binary agglomeration tree
- `FullViz3D.cc` provides visualization capabilities for the full algorithm
	- this produces a single json file (and corresponding directory for pdfs to be deposited)
	- to produce pdfs, run `python macros/Viz3DFullCluster.py -j [json_file]`
		- this will produce pdfs and open the highest level clustering in a browser

### (Jet) Clustering (Physics-specific)
- `JetProducer` either grabs rechits from specified ROOT file (NANOAOD/Ntuple) or simulates rechits to cluster
- `Jet` is 1 cluster unit (it can be anything, doesn't have to be a physics jet)
	- it can be composite (made up of multiple Jets) or singular (ie one rechit)
	- these are the units that are clustered in `JetCluster` 
- `Cluster` calls the full algorithm and runs the corresponding EM algorithm to find the cluster parameters
	- it returns a vector of `node`s, where each `node` is a jet
	- a jet's (node's) clustering history can be traced back individually through its parents
- `FindSubjets` only calls the variational EM algorithm with Gaussian mixture model
	- can be done in $\eta - \phi$ space or $X - Y - Z$ space, both with a time dimension
- `Jet` and `JetPoint` (the massless counterpart of `Jet`, a version of a RecHit) both store the phi ($\phi$) coordinate on a $[-\pi, \pi]$ range
	- `phi_02pi()` for both classes returns phi on a range from $[0, 2\pi]$ but still keeps the original coordinate value
	- for `Jet`s and `JetPoint`s that are instantiated from a momentum four-vector with no external $\phi$ or $\eta$ provided via the `SetPhi()` or `SetEta()` functions, these coordinates will be calculated as they are in [FastJet](https://fastjet.fr/repo/doxygen-3.4.2/PseudoJet_8cc_source.html) for `Jet`s and ROOT for `JetPoint`s ([eta](https://root.cern.ch/doc/master/eta_8h_source.html#l00048) and [phi](https://root.cern.ch/doc/master/GenVector_2Cartesian3D_8h_source.html#l00117))
	- `Jet`s are instantiated with a momentum four-vector $(p_x, p_y, p_z, E)$ (in that order) while a `JetPoint` is instantiated from a spatial fourvector $(x, y, z, t)$
		- this is because a `JetPoint` (or RecHit) has no mass and no associated vertex (only spatial and energy information) and therefore a momentum four-vector cannot be specified (unless you assume a zero-mass scenario and provide the coordinates for the vertex)
- for the Bayesian Hierarchical Clustering (BHC) a distance constraint of $\pi/2$ is defaulted to in `BayesCluster`
	- the distance constraint sets the merge posterior ($r_k$) of a potential merge if the constituents are too far apart
	- see implementation in `MergeTree::DistanceConstraint()`

### Plotting
- to add a new histogram to the skims, need to declare it in BaseSkimmer and push it back to the vector in ctor
	- fill histogram in the derived object skimmer
- plots are saved as histograms (not PDFs or TCanvases) so they can be hadded together
	- if running the skimmer on condor (see below)
- to format histograms and save them as TCanvases, run `root -l -b -q 'macros/HistFormat.C("[skim.root]")'` from command line
- to add a sample to an overlaid (stack) plot (ie when looking at data as a proxy for background) you can run that sample (like JetHT) separately then hadd the total root files to the ones with signal and other backgrounds
	- make sure to add sample to MakeIDHists() in `PhotonSkimmer.hh`

#### Condor
- the skimmer can be run on condor (on the LPC) with the following steps:
	- `python2 generateSubmission.py --inputFile [file]` generates the submission script for condor
		- needs to be run in python2 because on the LPC in CMSSW_10_X_X PyROOT is not available in python3
- run the python scripts and submit scripts from the condor folder

### References and Acknowledgements
- [Bayesian Hierarchical Clustering](https://www2.stat.duke.edu/~kheller/bhcnew.pdf)
- The seeding for the Bayesian Hierarchical Clustering is based on the Voronoi diagram implementation of finding geometric nearest neighbors as done in [FastJet](https://arxiv.org/abs/1111.6097)
