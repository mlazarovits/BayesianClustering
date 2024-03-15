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

	if args.inputSample == "GMSB_L500_ctau1000":	
		inputFile = "GMSB_AOD_v14_GMSB_L-500TeV_Ctau-1000cm_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GMSB_L350_ctau200":	
		inputFile = "GMSB_AOD_v14_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GMSB_L150_ctau200":	
		inputFile = "GMSB_AOD_v14_GMSB_L-150TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "JetHT":
		inputFile = "JetHT_Met150_AOD_v14_JetHT_AOD_Run2018CRun2018D-15Feb2022_UL2018-v1.root"
	elif args.inputSample == "GJets_HT400To600":
		inputFile = "GJets_AOD_v14_GJets_HT-400To600_AODSIM_RunIISummer20UL18RECO-106X_upgrade2018.root"
	elif args.inputSample == "GJets_HT600ToInf":
		inputFile = "GJets_AOD_v14_GJets_HT-600ToInf_AODSIM_RunIISummer20UL18RECO-106X_upgrade2018.root" 
	elif args.inputSample == "GMSB_L100_ctau0p1":
		inputFile = "GMSB_AOD_v14_GMSB_L-100TeV_Ctau-0_1cm_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GJets_HT400To600_v15":
		inputFile = "GJets_R17_v15_GJets_HT-400To600_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GJets_HT600ToInf_v15":
		inputFile = "GJets_R17_v15_GJets_HT-600ToInf_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GJets_HT400To600_2017_v16":
		inputFile = "GJets_R17_v16_GJets_HT-400To600_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GJets_HT400To600_2018_v16":
		inputFile = "GJets_R18_v16_GJets_HT-400To600_AODSIM_RunIISummer20UL18RECO-106X_upgrade2018.root"
	elif args.inputSample == "GMSB_L300ctau600_v15":
		inputFile = "GMSB_R17_v15_GMSB_L-300TeV_Ctau-600cm_AODSIM_RunIIFall17DRPremix.root"
	elif args.inputSample == "GMSB_L150ctau0p1_v15":
		inputFile = "GMSB_R17_v15_GMSB_L-150TeV_Ctau-0_1cm_AODSIM_RunIIFall17DRPremix.root"

	else:
		print("Sample "+args.inputSample+" not found")
	#to use xrootd path cannot be relative
	#find any ../ and remove it and the dir before it
	inputFile = "root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/"+inputFile
	
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
	
	printstring = "Skimming for "+objName+" in file "+inputFile
	if(objName == "jets"):
		printstring += " with strategy "+strategyName 
	print printstring

	#organize output by sample, object (ie jets or photons), and strategy (for jets only - NlnN or N2)
	#find .root and then find the / before that (if it exists) - everything in between is the file name
	sampleName = inputFile[ inputFile.rfind("/")+1 : inputFile.find(".root") ]
	sampleNameShort = sampleName[ : sampleName.find("_") ]

	dirname = odir+sampleNameShort+"/"+sampleName+"/"+objName
	ofilename = args.inputSample+"_"+objName
	if(objName == "jets"):
		dirname += "/"+strategyName
		ofilename += "_"+strategyName
	#strategy doesn't apply to photons (GMM only)

	if args.output is not None and objName not in args.output:
		ofilename = ofilename+"_"+args.output 
		dirname = dirname+"_"+args.output
	ofilename = "condor_"+ofilename

	print("Preparing sample directory: {0}".format(dirname))
	##### Create a workspace (remove existing directory) #####
	if os.path.exists(dirname):
		print("Removing existing directory: {0}".format(dirname))
		shutil.rmtree(dirname)

	# Create directories for work area.
	SH.createWorkArea(dirname)

	# grab relevant flags
	eventnums = SH.eventsSplit(inputFile, args.split)
	flags = '--alpha '+str(args.alpha)+' --EMalpha '+str(args.EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" -s "+str(args.strategy)+" --gev "+str(args.gev)+' --minpt '+str(args.minpt)+' --minNrhs '+str(args.minnrhs)+' --minemE '+str(args.minemE)+' --minRhE '+str(args.minRhE)
	if(args.noSmear):
		flags += ' --noSmear'
	if(args.timeSmear):
		flags += ' --timeSmear'
	if(args.applyFrac):
		flags += ' --applyFrac'
	if(objName == "jets"):
		flags += " --object 0"
	if(objName == "photons"):
		flags += " --object 1"

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
	parser.add_argument('--inputSample','-i',help='Ntuple sample to create skims from',required=True,choices=['GMSB_L500_ctau1000','GMSB_L350_ctau200','GMSB_L150_ctau200','GMSB_L100_ctau0p1','JetHT','GJets_HT400To600','GJets_HT600ToInf','GJets_HT400To600_v15','GJets_HT600ToInf_v15','GMSB_L300ctau600_v15','GMSB_L150ctau0p1_v15','GJets_HT400To600_2017_v16','GJets_HT400To600_2018_v16'])
	parser.add_argument('--output','-o',help='output label')
	parser.add_argument('--year',help='year of sample',default=2017)
	#which object to analyze (jets or photons currently supported)
	parser.add_argument('--object',help='which object to skim (currently only jets or photons supported)',choices=["jets","photons"],required=True)
	parser.add_argument('--strategy','-st',help='if skimming jets, which strategy to use for BHC (NlnN = 0 default, N2 = 1, GMM only = 2)',default=0,type=int,choices=[2,1,0])
	parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
	parser.add_argument('--verbosity','-v',help="verbosity",default=0)
	#add algorithm parameters - alpha, emAlpha, verbosity, thresh
	parser.add_argument('--alpha','-a',help="alpha for BHC",default=0.1)
	parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=0.5)
	parser.add_argument('--thresh','-t',help='threshold for GMM clusters',default=1.)
	parser.add_argument('--gev',help='energy transfer factor',default=1./30.)
	parser.add_argument('--minpt',help='min object pt',default=30.)
	parser.add_argument('--minnrhs',help='min object nrhs',default=15)
	parser.add_argument('--minemE',help='min object ECAL energy',default=20)
	parser.add_argument('--minRhE',help='min rechit ECAL energy',default=0.5)
	parser.add_argument('--noSmear',help="turn off spatial smearing",default=False,action='store_true')
	parser.add_argument('--timeSmear',help="turn on time smearing",default=False,action='store_true')
	parser.add_argument('--applyFrac',help="apply fractions from hitsAndFractions list to rh energies for photons",default=False,action='store_true')
	args = parser.parse_args()

	generateSubmission(args)

if __name__ == "__main__":
    main()
