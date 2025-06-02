#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################

import sys
import os
import argparse
import shutil
import submissionHelper as SH
from fractions import Fraction


# Create workspace and condor submit files.
def generateSubmission(args):
    # Ensure that the directory includes "/" at the end.
    odir = args.directory
    if odir[-1] != "/":
    	odir += "/"
    
    print("Directory for condor submission: {0}".format(odir))
    print("------------------------------------------------------------")
    # Create output directory for condor results if it does not exist.
    SH.makeDir(odir)
    
    sampleOptions = ['ttbar','QCD']
    sampleName = ""
    if args.ttbar:
    	sampleName = "ttbar"
    elif args.QCD:
    	sampleName = "QCD"
    else:
                print("Sample must be provided, current options are",sampleOptions)
                exit()
    
    
    dirname = odir+sampleName
    ofilename = "condorSimNtuples_"+sampleName
    if args.output is not None:
                ofilename = ofilename+"_"+args.output
                dirname = dirname+"_"+args.output
    #put algo config in file name
    print("Preparing sample directory: {0}".format(dirname))
    ##### Create a workspace (remove existing directory) #####
    if os.path.exists(dirname):
    	print("Removing existing directory: {0}".format(dirname))
    	shutil.rmtree(dirname)
    print("ofilename",ofilename) 
    # Create directories for work area.
    SH.createWorkArea(dirname)
    
    # grab relevant flags
    eventnums = SH.eventsSplit(int(args.nevts), args.split)
    flags = '-v '+str(args.verbosity)+' --nevts '+str(args.nevts)+' --spikeProb '+str(args.spikeProb) + ' --eThresh ' +str(args.eThresh)
    if(args.ttbar):
        flags += ' --ttbar'
    if(args.QCD):
        flags += ' --QCD'
    if(args.sigDelayed):
        flags += ' --sigDelayed'
    if(args.sigBoosted):
        flags += ' --sigBoosted'
    if(args.pileup):
        flags += ' --pileup'

    flags += ' --energyCte '+str(args.energyCte)

    ##### Create condor submission script in src directory #####
    condorSubmitFile = dirname + "/src/submit.sh"
    subf = open(condorSubmitFile, "w")
    print("outputfile name "+ofilename)
    SH.writeSubmissionBase(subf, dirname, ofilename)
    SH.writeQueueList(subf, ofilename, eventnums, flags)
    #subf.close()
    if eventnums == 0 or eventnums is None:
    	return
    
    print("------------------------------------------------------------")
    print("Submission ready, to run use:")
    print("condor_submit "+condorSubmitFile)

def main():
    # options
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", "-d", default="Output", help="working directory for condor submission")
    #Ntuple file to run over
    parser.add_argument('--output','-o',help='output label')
    parser.add_argument('--nevts',help='number of events to simulate (default = 100)',default=100)
    parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
    parser.add_argument('--verbosity','-v',help="verbosity",default=0)
    parser.add_argument('--ttbar',help="run ttbar process",action='store_true')
    parser.add_argument('--QCD',help="run QCD process",action='store_true')
    parser.add_argument('--sigDelayed',help="run sigDelayed process",action='store_true')
    parser.add_argument('--sigBoosted',help="run sigBoosted process",action='store_true')
    parser.add_argument('--pileup','-pu',help="run pileup process",action='store_true')
    parser.add_argument('--spikeProb',help='set probability of spike occuring (default = 0, off)',default = 0)	
    parser.add_argument('--energyCte',help='set energy smearing constant (default = 0.26)',default = 0.26)
    parser.add_argument('--eThresh',help='set energy threshold for rechit reco (default = 0.5)',default = 0.5)
    args = parser.parse_args()
    
    generateSubmission(args)

if __name__ == "__main__":
    main()
