# BayesianClustering
Repository for generic EM/hierarchical clustering algorithm (to be used for jet and photon clustering at CMS).

### Compiling Executables
- to run on LPC: `make lpc`
	- because of the fact that this command makes (or updates) directories with condor tarballs, this command cannot be run with multithreading (ie `-j` flag not allowed)
- to run on local machine:
	- make sure paths in Makefile are updated accordingly
	- `make`
- for compiling in a larger package, like [KUCMSNtupleizer](https://github.com/jking79/KUCMSNtupleizer/tree/master) with some classes already defined, like the KUCMSTimeCalibration class, make sure to use the most updated versions of these classes from [here](https://github.com/jking79/KUCMSNtupleizer/tree/master/KUCMSSkimmer/KUCMSTimeCaliFiles)


### Dependencies
Once all of the below packages are downloaded, the corresponding paths in the Makefile must be updated to point to the relevant packages
- g++ for compilation of executables and linking of shared objects
	- this framework was developed with g++11, moving to a more recent g++ version may have some compilation/linking issues
- ROOT for producers and skimmers
- Python for macros for visualization
	- at least v3.x
- [KUCMSTimeCalibration code](https://github.com/jking79/KUCMSNtupleizer.git) (until skimmers are split out...)
All external packages below are header only: 
- [eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [boost](https://www.boost.org/doc/libs/1_82_0/libs/math/doc/html/special.html)
- [nlohmman-json](https://github.com/nlohmann/json) to make jsons for python macros
- [cgal](https://www.cgal.org/) for Voronoi diagrams
	- on LPC: make sure the path to the CGAL lib is added to `$LD_LIBRARY_PATH`
		- this can be done within the user's `~/.bash_profile` by adding the following line
		```
		export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/cgal/installation
		```
- [fastjet](https://fastjet.fr/) for detector simulation
	- at least v3.4.1 for compatibility with C++17
	- at least v3.4.2 for compatibility with C++20 
- [plotly](https://plotly.com/graphing-libraries/) for producing 3D variational clustering plots and gifs
- [frgually-deep](https://github.com/Dobiasd/frugally-deep/blob/master/INSTALL.md) for importing neural networks trained in Tensorflow to C++
- [Functional-Plus](https://github.com/Dobiasd/FunctionalPlus) required for frugally-deep 

### Branches
- `main` branch is latest stable release
- always work on `dev` branch (current ntuple version is v24) 
- `dev_v15` is compatible with v15, v16 of KUCMS ntuples
- `dev_v14` is compatible with v14 of KUCMS ntuples

### Use as an external package
- make sure `lib/libBayesianClustering.so` exists
	- if not, it can be made with `make lib` or `make lpclib` if on LPC
- you MUST include the following flags in the compilation of the main executable that `libBayesianClustering` is being linked against
	- `-O3 -march=native -ffast-math -funroll-loops -flto`
	- this is to ensure that `Eigen::MatrixXd` objects are aligned properly
- be sure to add the path to the repository to `$LD_LIBRARY_PATH` so the library can be dynamically found
	```
	export $LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/BayesianClustering/lib
	```
- In Makefile of overall framework where you want to include `BayesianClustering`:
	- Add the following flags to `CXXFLAGS` (or equivalent)
	```
	CXXFLAGS += -I/path/to/BayesianClustering/include -frounding-math
	```
	- Add the following flags to `GLIBS` (or equivalent)
	```
	GLIBS += -L/path/to/BayesianClustering/lib -lBayesianClustering
	```
- If there are errors related to CGAL (ie undefined references)
	- Add the following flags to `GLIBS` (or equivalent)
	```
	GLIBS += -L/path/to/CGALvX.X-install/X.X/lib -lCGAL
	```
- If there are errors related to boost (ie undefined references)
	- Add the following flags to `CXXFLAGS` (or equivalent)
	```
	CXXFLAGS += -I/path/to/boost_vX.X.X-install/X.X.X/include
	```
- If there are errors related to how CGAL needs boost (ie `warning: libboost_thread.so.1.63.0, needed by /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib//libCGAL.so, not found`)
	- Add the following flags to `GLIBS` (or equivalent) *before* the libraries for CGAL
	```
	GLIBS   += -Wl,-rpath-link,/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0/lib/ -lboost_thread
	```
- Then, in your executable be sure to include the header `#include BayesianClustering/BayesianClustering.hh`

### Interfacing with `BayesianClustering`
A user can interface with the clustering algorithm through the API class `ClusterAnalyzer.` After declaring an instance of the class and setting any advance values (like `SetDetectorCenter(x, y, z)` or `SetTransferFactor(t)`), the user can loop over the points for clustering. Within this for loop, the user can call `AddRecHit(x, y, z, E, t, rhId)` to add a point to the list for clustering. After all the points have been added and the clustering is run, be sure to call `ClearRecHitList()` to clear the list of data points and start again, if desired.

When all the relevant data points have been added, the user can call `RunClustering()` to run the clustering algorithm and produce a `ClusterObj`, returned by this method. This object will have methods to calculate and return various observables such as the object time at the detector face (`GetObjTime_Det()`), the object time at the primary vertex (`GetObjTime_PV()`), the object's time significance given a resolution map in the form `map<int rhId, double res>` (`GetObjTimeSig()`). Each of these methods must be called after the corresponding calculation method. So, for example, the `GetObjTime()` methods need to be called after `CalculateObjTime()` and the `GetObjTimeSig()` needs to be called after `CalculateObjTimeSig(map<rhId, res>)`. As default, the clustering algorithm run in the `RunClustering()` method is the NlnN version of the full hierarchical clustering algorithm. The priors and other hyperparameters for this algorithm are hard-coded into the `RunClustering()` method. 


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
- it is assumed that the prior parameters are given in the final scaled + shifted data frame of reference

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
- to create profiles and derived histograms from 2D histograms, run `root -l -b -q 'macros/Profile2D.C("[skim.root]")'
	- this will edit the passed histograms to fill the respective, empty histograms that depend on profiles
- to format histograms and save them as TCanvases, run `root -l -b -q 'macros/HistFormatJets.C("[skim.root]")'` or the photon/sim version (`macros/HistFormat.C("[skim.root")'`) from command line
- to quickly make PDFs for a subset of histograms, run `root -l -b -q 'macros/MakePDFs.C("[input_skim_formatted.root]","[output_dir]","[hist_name_match_string]")'`
- to add a sample to an overlaid (stack) plot (ie when looking at data as a proxy for background) you can run that sample (like JetHT) separately then hadd the total root files to the ones with signal and other backgrounds
	- make sure to add sample to MakeIDHists() in `PhotonSkimmer.hh`
	- because of this, for JetSkimming, the 1D profiles and time sigma plots are *not* filled (only the 2D diffDeltaTime plots and the other 1D variables are filled) so the user can hadd the 2D histograms together
		- the user should hadd the skims together first then run `macros/Profile2D.C("file")` then `macros/HistFormatJets.C("file")` for the correct plots


### Condor
- the skimmer can be run on condor (on the LPC) with the following steps:
	- `python2 generateSubmission.py --inputFile [file]` generates the submission script for condor
		- needs to be run in python2 because on the LPC in CMSSW_10_X_X PyROOT is not available in python3
- run the python scripts and submit scripts from the condor folder
- to combine CSV files for MVA inputs, run `source combineCSVs.h` from the `condor/` folder, passing the directory to the `out/*.csv` path as a positional argument



### Including classes as dictionaries in ROOT
Sometimes, you may need to add more classes to the ROOT compiler. For example, ROOT cannot fill TBranches with `vector<vector<type>>`. In this case, you will have to add the required class as a dictionary. First, you need to create a `LinkDef.h` file defining the classes. An example of this for the `vector<vector<int>>` class is in the base directory of this repository. Then, you will have to create the dictionary by running the following command
```
 rootcling -v4 -f vecDict.cxx  -rmf libvecDict.rootmap -rml libvecDict.so LinkDef.h
```
Once the dictionary is made, you then need to compile it as a shared library to include in the compilation of the framework. This is done with
```
g++ -shared -o libvecDict.so vecDict.cxx `root-config --cflags --libs`
```
After this step, I usually move the `libvecDict.so` and `libmydict_rdict.pcm` files to the `lib` directory to keep things clean. Then, you will want to add the following commands to all the constructors of the classes that require the new dictionary-based classes
```
gSystem->Load("lib/libvecDict.so");
```
I also include this line in the `rootlogon.C` script so it can be read by the interpreter.



### Troubleshooting
- if experiencing undefined behavior, crashes on copy assignment operators, destructors, etc, make SURE the package you are linking BayesianClustering against is compiled with these flags:
	`-O3 -march=native -ffast-math -funroll-loops -flto`
Undefined behavior happens because if the vector<Matrix> or Matrix object was allocated without proper alignment, the destructor may try to free memory that is misaligned, which could cause a crash in `free()` in a compiler-generated destructor call

### References and Acknowledgements
- [Bayesian Hierarchical Clustering](https://www2.stat.duke.edu/~kheller/bhcnew.pdf)
- The seeding for the Bayesian Hierarchical Clustering is based on the Voronoi diagram implementation of finding geometric nearest neighbors as done in [FastJet](https://arxiv.org/abs/1111.6097)
