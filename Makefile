#include ROOT cflags and libraries
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTGLIBS   = $(shell root-config --glibs)
#include pythia cflags and libraries
local: PYTHIACFLAGS = $(shell /Users/margaretlazarovits/pythia8307/bin/pythia8-config --cflags)
local: PYTHIAGLIBS  = $(shell /Users/margaretlazarovits/pythia8307/bin/pythia8-config --libs) 
lpc:   PYTHIACFLAGS = -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/240/include
lpc:   PYTHIAGLIBS  = -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/240/lib -Wl,-rpath,/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pythia8/240/lib -lpythia8 -ldl

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
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/eigen/1ae2849542a7892089f81f2ee460b510cdb0a16d/include/eigen3/
#add digamma include path
local: CXXFLAGS    += -I/opt/homebrew/Cellar/boost/1.82.0_1/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0/include/ 
#add jsoncpp flags
local: CXXFLAGS    += -I/opt/homebrew/Cellar/nlohmann-json/3.11.2/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/json/3.7.3/include/ 
#add CGAL flags and libraries
local: CXXFLAGS    += -I/opt/homebrew/Cellar/cgal/5.6/include/
#CGAL is at v4.2 on LPC - not a header-only version ig
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/include/
lpc:   GLIBS       += -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib/
lpc:   GLIBS       += -lCGAL -lCGAL_Core
#boost for CGAL
#example from KUEWkinoAnalysis makefile #lpc:   GLIBS += -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_5/external/slc7_amd64_gcc700/lib/ -lvdt -lboost_program_options -lboost_filesystem -lboost_regex -lboost_system
#need to dynamically link the boost shared library so CGAL can use it
#use -Wl comma separated to pass to linker
lpc:   GLIBS       += -Wl,-rpath-link,/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0/lib/ -lboost_thread
#include FastJet cxxflags and libraries
local: CXXFLAGS  += $(shell ~/fastjet-install/bin/fastjet-config --cxxflags)
local: GLIBS     += $(shell ~/fastjet-install/bin/fastjet-config --libs)
lpc: CXXFLAGS  += $(shell /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0/bin/fastjet-config --cxxflags)
lpc: GLIBS     += $(shell /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/fastjet/3.3.0/bin/fastjet-config --libs)

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
lpclib: SGLIBS       += -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib/
lpclib: SGLIBS       += -lCGAL -lCGAL_Core

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


