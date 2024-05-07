#include ROOT cflags and libraries
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTGLIBS   = $(shell root-config --glibs)
#include pythia cflags and libraries
local: PYTHIACFLAGS = $(shell /Users/margaretlazarovits/pythia8307/bin/pythia8-config --cflags)
local: PYTHIAGLIBS  = $(shell /Users/margaretlazarovits/pythia8307/bin/pythia8-config --libs) 
lpc:   PYTHIACFLAGS = -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pythia8/309-c48e277ae4ccb32ab17f6e0b5f0c5d07/include
lpc:   PYTHIAGLIBS  = -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pythia8/309-c48e277ae4ccb32ab17f6e0b5f0c5d07/lib -Wl,-rpath,/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pythia8/309-c48e277ae4ccb32ab17f6e0b5f0c5d07/lib -lpythia8 -ldl

#specify compiler
CXX         = g++
CXX	   += -frounding-math

#set c(xx)flags 
CXXFLAGS    = $(ROOTCFLAGS)
CXXFLAGS   += $(PYTHIACFLAGS) 
 
#add libraries
GLIBS       = $(ROOTGLIBS)
GLIBS      += $(PYTHIAGLIBS)

#add eigen include path
local: CXXFLAGS    += -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/
#add digamma include path
local: CXXFLAGS    += -I/opt/homebrew/Cellar/boost/1.82.0_1/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-f76596f4b83666ac3468f34a5f342677/include/ 
#add jsoncpp flags
local: CXXFLAGS    += -I/opt/homebrew/Cellar/nlohmann-json/3.11.2/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/json/3.10.2-a6d86565b09ec3d0e02bf7b52c31bbfc/include/ 
#add CGAL flags and libraries
local: CXXFLAGS    += -I/opt/homebrew/Cellar/cgal/5.6/include/
#include necessary CGAL libraries BEFORE the include file so the compile knows about them
lpc:   GLIBS        += -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gmp-static/6.2.1-f4591b847fcbe5753bfc5d2b02f57089/lib/ -lgmp
lpc:   CXXFLAGS     += -I$(CMSSW_BASE)/src/CGAL-5.6.1/include
#lpc:   CXXFLAGS     += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/CGAL-5.6.1/include/
#boost for CGAL
#example from KUEWkinoAnalysis makefile #lpc:   GLIBS += -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_5/external/slc7_amd64_gcc700/lib/ -lvdt -lboost_program_options -lboost_filesystem -lboost_regex -lboost_system
#need to dynamically link the boost shared library so CGAL can use it
#use -Wl comma separated to pass to linker
lpc:   GLIBS       += -Wl,-rpath-link,/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-f76596f4b83666ac3468f34a5f342677/lib/ -lboost_thread
#include FastJet cxxflags and libraries
local: CXXFLAGS  += $(shell ~/fastjet-install/bin/fastjet-config --cxxflags)
local: GLIBS     += $(shell ~/fastjet-install/bin/fastjet-config --libs)
lpc: CXXFLAGS  += $(shell /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fastjet/3.4.1-b5a7b930eb5755ed2b2b87b323687b41/bin/fastjet-config  --cxxflags)
lpc: GLIBS     += $(shell /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fastjet/3.4.1-b5a7b930eb5755ed2b2b87b323687b41/bin/fastjet-config --libs)

#add to make lib ig?
lpc: CXXFLAGS += -fPIC
#make list of libraries for shared library
lib: SGLIBS          = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
lib: SGLIBS         += -lRooFit -lRooFitCore
lib: SGLIBS         += $(shell ~/fastjet-install/bin/fastjet-config --libs)
lib: SGLIBS         += $(PYTHIAGLIBS)

lpclib: SGLIBS        = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
lpclib: SGLIBS       += -lRooFit -lRooFitCore
#fastjet
lpclib: SGLIBS       += $(shell /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0/bin/fastjet-config --libs)
#pythia
lpclib: SGLIBS       += $(PYTHIAGLIBS)
lpclib: SGLIBS        += -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gmp-static/6.2.1-f4591b847fcbe5753bfc5d2b02f57089/lib/ -lgmp
#lpclib: SGLIBS       += -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib/
#lpclib: SGLIBS       += -lCGAL -lCGAL_Core

#specify local paths
INCLUDEDIR  = ./include/
SRCDIR      = ./src/
#make sure compiler knows where local include files are
CXX	   += -I$(INCLUDEDIR) -I.
OUTOBJ	    = ./obj/

#specify local source, include, and object files
CC_FILES        = $(wildcard src/*.cc)
HH_FILES        = $(wildcard include/*.hh)
OBJ_FILES       = $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))

SOBJ_FILES = $(filter-out ./obj/BasicDetectorSim.o, $(OBJ_FILES))

#specify what to make
#all: GMM.x varGMM.x jetAlgo.x photonAlgo.x FullClusterSingle.x FullClusterSkim.x detectorSim.x 
all: FullClusterSingle.x FullClusterSkim.x detectorSimNtuples.x detectorSimSkimmer.x 
local: all
lpc:   all configtar lpclib simconfigtar
lib: lib/libBayesCluster.so
lpclib: lib/libBayesCluster.so

#executables
GMM.x: $(SRCDIR)GMM.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o GMM.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch GMM.x

varGMM.x: $(SRCDIR)varGMM.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o varGMM.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch varGMM.x

jetAlgo.x: $(SRCDIR)jetAlgo.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o jetAlgo.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch jetAlgo.x

photonAlgo.x: $(SRCDIR)photonAlgo.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o photonAlgo.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch photonAlgo.x

FullClusterSingle.x: $(SRCDIR)FullClusterSingle.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o FullClusterSingle.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch FullClusterSingle.x

FullClusterSkim.x: $(SRCDIR)FullClusterSkim.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o FullClusterSkim.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch FullClusterSkim.x

detectorSimNtuples.x: $(SRCDIR)detectorSimNtuples.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o detectorSimNtuples.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch detectorSimNtuples.x

detectorSimSkimmer.x: $(SRCDIR)detectorSimSkimmer.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o detectorSimSkimmer.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch detectorSimSkimmer.x

configtar:
	cp FullCluster*.x config/
	tar -czf condor/config.tgz config/

simconfigtar:
	cp detectorSim*.x configSim/
	tar -czf condorSim/configSim.tgz configSim/
	tar -czf condorSimNtuples/configSim.tgz configSim/

#make shared library
lib/libBayesCluster.so: $(SOBJ_FILES)
	mkdir -p lib
	$(CXX) -shared -o lib/libBayesianClustering.so $(SOBJ_FILES) $(SGLIBS) -fPIC
	touch lib/libBayesianClustering.so

#where to put object files
$(OUTOBJ)%.o: src/%.cc include/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

#clean options
clean:
	rm -f $(OUTOBJ)*.o
	rm -f *.x
	rm -f AutoDict*
	rm -f *.d
	rm -f lib/libBayesCluster.so


