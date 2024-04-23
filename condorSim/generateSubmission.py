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

	if args.inputSample == "ttbar":
	        inputFile = "simNtuples_ttbar.root"
	if args.inputSample == "ttbar_QCD":
	        inputFile = "simNtuples_ttbarQCD.root"
	else:
                print("Sample "+args.inputSample+" not found")
                exit()
	inputFile = "root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/"+inputFile
	
	#organize output by sample, object (ie jets or photons), and strategy (for jets only - NlnN or N2)
	#find .root and then find the / before that (if it exists) - everything in between is the file name
	sampleName = inputFile[ inputFile.rfind("/")+1 : inputFile.find(".root") ]
	sampleNameShort = sampleName[  sampleName.find("_")+1 : ]

	dirname = odir+sampleNameShort
	ofilename = "condorSim_"+sampleNameShort
	if args.output is not None:
                ofilename = ofilename+"_"+args.output
                dirname = dirname+"_"+args.output
	#put algo config in file name
	kname = "%.3f" % args.alpha
	kname = kname.replace(".","p")
	paramsname = "_bhcAlpha"+kname
	kname = "%.3f" % args.EMalpha
	kname = kname.replace(".","p")
	paramsname += "_emAlpha"+kname
	thresh = str(args.thresh).replace(".","p")
	paramsname += "_thresh"+thresh
	k = float(sum(Fraction(s) for s in args.gev.split()))
	kname = "%.3f" % k
	kname = kname.replace(".","p")
	paramsname += "_NperGeV"+kname
	ofilename += paramsname
	dirname += paramsname
	print("Preparing sample directory: {0}".format(dirname))
	##### Create a workspace (remove existing directory) #####
	if os.path.exists(dirname):
		print("Removing existing directory: {0}".format(dirname))
		shutil.rmtree(dirname)

	# Create directories for work area.
	SH.createWorkArea(dirname)

	# grab relevant flags
	eventnums = SH.eventsSplit(inputFile, args.split)
	flags = '--alpha '+str(args.alpha)+' --EMalpha '+str(args.EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" --gev "+str(args.gev)+' --minpt '+str(args.minpt)+' --minNrhs '+str(args.minnrhs)+' --minemE '+str(args.minemE)+' --minRhE '+str(args.minRhE)+' -s '+str(args.strategy)
	if(args.noSmear):
		flags += ' --noSmear'
	if(args.timeSmear):
		flags += ' --timeSmear'
	if(args.applyFrac):
		flags += ' --applyFrac'

	##### Create condor submission script in src directory #####
	condorSubmitFile = dirname + "/src/submit.sh"
	subf = open(condorSubmitFile, "w")
	print("outputfile name "+ofilename)
	SH.writeSubmissionBase(subf, dirname, ofilename, inputFile)
	SH.writeQueueList(subf, inputFile, ofilename, eventnums, flags)
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
	parser.add_argument('-inputSample','-i',help='Ntuple sample to create skims from',required=True,choices=['ttbar','ttbar_QCD'])
	parser.add_argument('--output','-o',help='output label')
	parser.add_argument('--strategy','-st',help='which strategy to use for BHC (NlnN = 0 default, N2 = 1)',default=0,type=int,choices=[1,0])
	parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
	parser.add_argument('--verbosity','-v',help="verbosity",default=0)
	#add algorithm parameters - alpha, emAlpha, verbosity, thresh
	parser.add_argument('--alpha','-a',help="alpha for BHC",default=0.1)
	parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=0.5)
	parser.add_argument('--thresh','-t',help='threshold for GMM clusters',default=1.)
	parser.add_argument('--gev',help='energy transfer factor',default=1/10)
	parser.add_argument('--minpt',help='min object pt',default=10.)
	parser.add_argument('--minnrhs',help='min object nrhs',default=2)
	parser.add_argument('--minemE',help='min object ECAL energy',default=0)
	parser.add_argument('--minRhE',help='min rechit ECAL energy',default=0.5)
	parser.add_argument('--noSmear',help="turn off spatial smearing",default=False,action='store_true')
	parser.add_argument('--timeSmear',help="turn on time smearing",default=False,action='store_true')
	parser.add_argument('--applyFrac',help="apply fractions from hitsAndFractions list to rh energies for photons",default=False,action='store_true')
	args = parser.parse_args()

	generateSubmission(args)

if __name__ == "__main__":
    main()
