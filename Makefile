#include ROOT cflags and libraries
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTGLIBS   = $(shell root-config --glibs)

#set c(xx)flags and libraries
CXXFLAGS    = $(ROOTCFLAGS)

#add eigen include path
#CXXFLAGS    += -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/
##add digamma include path
#CXXFLAGS    += -I/opt/homebrew/Cellar/boost/1.81.0_1/include/
#GLIBS       = $(ROOTGLIBS)
##add jsoncpp flags
#CXXFLAGS    += $(shell pkg-config --cflags /opt/homebrew/Cellar/jsoncpp/1.9.5/lib/pkgconfig/jsoncpp.pc)
#GLIBS       += $(shell pkg-config --libs /opt/homebrew/Cellar/jsoncpp/1.9.5/lib/pkgconfig/jsoncpp.pc)
##add nlohmann-json flags
#CXXFLAGS   += -I/opt/homebrew/Cellar/nlohmann-json/3.11.2/include/
#specify compiler
CXX         = g++

#add libraries
GLIBS       = $(ROOTGLIBS)

#add eigen include path
local: CXXFLAGS    += -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/eigen/1ae2849542a7892089f81f2ee460b510cdb0a16d/include/eigen3/
#add digamma include path
local: CXXFLAGS    += -I/opt/homebrew/Cellar/boost/1.81.0_1/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.67.0/include/ 
#add jsoncpp flags
local: CXXFLAGS    += -I/opt/homebrew/Cellar/nlohmann-json/3.11.2/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/json/3.7.3/include/ 
#add CGAL flags
local: CXXFLAGS    += -I/opt/homebrew/Cellar/cgal/5.6/include/
lpc:   CXXFLAGS    += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/include/

#specify local paths
INCLUDEDIR  = ./include/
SRCDIR      = ./src/
#make sure compiler knows where local include files are
CXX	   += -I$(INCLUDEDIR) -I.
OUTOBJ	    = ./obj/

#specify local source, include, and object files
CC_FILES    = $(wildcard src/*.cc)
HH_FILES    = $(wildcard include/*.hh)
OBJ_FILES   = $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))

local: all
lpc:   all
#specify what to make
all: GMM.x varGMM.x jetAlgo.x fullAlgo.x photonAlgo.x

#executables
GMM.x: $(SRCDIR)GMM.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o GMM.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch GMM.x

varGMM.x: $(SRCDIR)varGMM.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o varGMM.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch varGMM.x

fullAlgo.x: $(SRCDIR)fullAlgo.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o fullAlgo.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch fullAlgo.x

jetAlgo.x: $(SRCDIR)jetAlgo.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o jetAlgo.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch jetAlgo.x

photonAlgo.x: $(SRCDIR)photonAlgo.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o photonAlgo.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch photonAlgo.x
#where to put object files
$(OUTOBJ)%.o: src/%.cc include/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

#clean options
clean:
	rm -f $(OUTOBJ)*.o
	rm -f *.x
	rm -f AutoDict*
	rm -f *.d



