#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################

import sys
import os
import argparse
import shutil
import submissionHelper as SH


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
	
	inputFile = args.inputFile
	objName = args.object
	strategyName = ""
	if args.strategy == 0:
	        strategyName = "NlnN"
	elif args.strategy == 1:
	        strategyName = "N2"
	elif args.strategy == 2:
		strategyName = "GMMonly"
	else:
	        print("Invalid strategy",args.strategy,"specified")
	        exit()
	
	print "Skimming for",objName,"in file",args.inputFile,"with strategy",strategyName 
	
	
	#organize output by sample, object (ie jets or photons), and strategy (for jets only - NlnN or N2)
	#find .root and then find the / before that (if it exists) - everything in between is the file name
	sampleName = inputFile[ inputFile.rfind("/")+1 : inputFile.find(".root") ]
	sampleNameShort = sampleName[ : sampleName.find("_") ]
	
	dirname = ""
	ofilename = ""
	if(objName == "jets"):
	        dirname = odir+sampleNameShort+"/"+sampleName+"/"+objName+"/"+strategyName
	        ofilename = sampleNameShort+"_"+objName+"_"+strategyName
	#strategy doesn't apply to photons (GMM only)
	else:
	        dirname = odir+"/"+sampleNameShort+"/"+sampleName+"/"+objName
	        ofilename = dirname+"/out/"+sampleNameShort+"_"+objName

	
	print("Preparing sample directory: {0}".format(dirname))
	##### Create a workspace (remove existing directory) #####
	if os.path.exists(dirname):
		print("Removing existing directory: {0}".format(dirname))
		shutil.rmtree(dirname)

	# Create directories for work area.
	SH.createWorkArea(dirname)

	# grab relevant flags
	eventnums = SH.eventsSplit(inputFile, args.split)
	flags = '--alpha '+str(args.alpha)+' --EMalpha '+str(args.EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" -s "+str(args.strategy)

	##### Create condor submission script in src directory #####
	condorSubmitFile = dirname + "/src/submit.sh"
	subf = open(condorSubmitFile, "w")
	oFileName = "condorSkim"
	SH.writeSubmissionBase(subf, dirname, oFileName, inputFile)
	SH.writeQueueList(subf, inputFile, oFileName, eventnums, flags)
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
	parser.add_argument('--inputFile','-i',help='Ntuple file to create skims from (absolute path)',required=True)
	#which object to analyze (jets or photons currently supported)
	parser.add_argument('--object','-o',help='which object to skim (currently only jets or photons supported)',choices=["jets","photons"],required=True)
	parser.add_argument('--strategy','-st',help='if skimming jets, which strategy to use for BHC (NlnN = 0 default, N2 = 1, GMM only = 2)',default=0,type=int,choices=[2,1,0])
	parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
	parser.add_argument('--verbosity','-v',help="verbosity",default=0)
	#add algorithm parameters - alpha, emAlpha, verbosity, thresh
	parser.add_argument('--alpha','-a',help="alpha for BHC",default=0.1)
	parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=0.5)
	parser.add_argument('--thresh','-t',help='threshold for GMM clusters',default=1.)
	args = parser.parse_args()

	generateSubmission(args)

if __name__ == "__main__":
    main()
