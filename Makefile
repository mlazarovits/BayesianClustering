#include ROOT cflags and libraries
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTGLIBS   = $(shell root-config --glibs)

#set c(xx)flags and libraries
CXXFLAGS    = $(ROOTCFLAGS)

#add eigen include path
CXXFLAGS    += -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/
#add digamma include path
CXXFLAGS    += -I/opt/homebrew/Cellar/boost/1.81.0_1/include
GLIBS       = $(ROOTGLIBS)
#specify compiler
CXX         = g++

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


#specify what to make
all: GMM.x

#executables
GMM.x: $(SRCDIR)GMM.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o GMM.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch GMM.x

#where to put object files
$(OUTOBJ)%.o: src/%.cc include/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

#clean options
clean:
	rm -f $(OUTOBJ)*.o
	rm -f *.x
	rm -f AutoDict*
	rm -f *.d



