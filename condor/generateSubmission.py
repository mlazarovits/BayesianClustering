#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################

import sys
import os
import shutil
import submissionHelper as SH
import argparse

odir = 'Output/'
#Load input arguments ()
parser = argparse.ArgumentParser()
#Ntuple file to run over
parser.add_argument('--inputFile','-i',help='Ntuple file to create skims from')
#which object to analyze (jets or photons currently supported)
objName = sys.argv[4]
if(objName != zip("jet","photon")):
	print("Error: only jet and photon supported currently. Input:",objName)
	exit()
parser.add_argument('--object','-o',help='which object to skim (currently only jets or photons supported)',choices=["jets","photons"])
strategy = sys.argv[5]
parser.add_argument('--strategy','-st',help='if skimming jets, which strategy to use for BHC (NlnN = 0 default, N2 = 1)',default=0,type=int,choices=["1","0"])
parser.add_argument('--verbosity','-v',help='verbosity',default=0,type=int)
parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
#add algorithm parameters - alpha, emAlpha, verbosity, thresh
parser.add_argument('--alpha','-a',help="alpha for BHC",default=0.1)
parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=0.5)
parser.add_argument('--verbosity','-v',help="verbosity",default=0)
parser.add_argument('--thresh','-t',help='threshold for GMM clusters')
args = parser.parse_args()

inputFile = args.inputFile
objName = args.object
strategyName = args.strategy
print("Skimming for",objName,"in file",args.inputFile)


#organize output by sample, object (ie jets or photons), and strategy (for jets only - NlnN or N2)
#find .root and then find the / before that (if it exists) - everything in between is the file name
sampleName = inputFile[ inputFile.rfind("/")+1 : inputFile.find(".root") ]
sampleNameShort = sampleName[ : sampleName.find("_") ]

dirname = ""
ofilename = ""
if(objName == "jets"):
	dirname = odir+"/"+sampleNameShort+"/"+sampleName+"/"+objName+"/"+strategyName
	ofilename = sampleNameShort+"_"+objName+"_"+strategyName
else:
	dirname = odir+"/"+sampleNameShort+"/"+sampleName+"/"+objName
	ofilename = sampleNameShort+"_"+objName+"_"+strategyName
#####################Create a workspace (remove existing directory)############
if os.path.exists(dirname):
	print("removing existing directory",dirname)
	shutil.rmtree(dirname)
os.mkdirs(dirname)
os.mkdirs(dirname+"/src")
os.mkdirs(dirname+"/log")
os.mkdirs(dirname+"/out")


####################Create submission script in src directory##################
subf = open(dirname+"/src/submit.sh","w")
SH.writeSubmissionBase( subf, dirname )
if(args.split != 0):
	eventnums = SH.eventsSplit(inFile, args.split)
flags = '--alpha '+args.alpha+' --EMalpha '+args.EMalpha+' -v '+args.verbosity+' -t '+args.thresh+" -s "+args.strategy 
SH.writeQueueList(subf, inputFile, ofilename, eventnums, flags)

print("submission ready, to run use:")
#need to be in directory with the execution script to run
print("pushd ../ && condor_submit MakeNtuple_LPC/Output/"+dataSetName+"_"+yearTag+"/src/submit.sh")


















































