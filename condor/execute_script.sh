#!/bin/bash
#wget --quiet --no-check-certificate http://stash.osgconnect.net/+zflowers/backup/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2

tar -xzf ./config.tgz

export SCRAM_ARCH=slc7_amd64_gcc700
#sets CMSSW environment variables
source /cvmfs/cms.cern.ch/cmsset_default.sh
#sets up CMSSW environment from a sandbox
source ./config/cmssw_setup_connect.sh

cmssw_setup sandbox-CMSSW_10_6_5-6403d6f.tar.bz2

#don't need restframes library
#source ./config/setup_RestFrames_connect.sh 
 
#don't need lhapdf library - but may need to do this for cgal, etc.
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/lib/
./config/FullAlgo.x "$@"
