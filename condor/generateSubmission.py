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
    beta0 = 1e-5
    m0 = [0, 0, 0]
    W0diag = [0.33333333,0.33333333,33.333333e-5]
    nu0 = 3
    EMalpha = 1e-5
    # Ensure that the directory includes "/" at the end.
    odir = args.directory
    if odir[-1] != "/":
    	odir += "/"
    
    print("Directory for condor submission: {0}".format(odir))
    print("------------------------------------------------------------")
    # Create output directory for condor results if it does not exist.
    SH.makeDir(odir)

    ver = "v31"
    sel = args.selection#"MET100"
    yr = str(args.year)[-2:]
    inputFilterList = ""

    data = False
    if args.era is not None:
        data = True

    #to use xrootd path cannot be relative
    #find any ../ and remove it and the dir before it
    inputPathList = "/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/BayesianClustering/filelists/"
    #loop through all files in filelist dir
    condor_subs = []
    for filelist in os.listdir(inputPathList):
        if args.inputSample+"_" not in filelist:
            continue
        if sel not in filelist:
            continue
        if ver not in filelist:
            continue
        if data and "Run"+args.year+args.era not in filelist:
            continue
        if "SMS" in args.inputSample:
            if args.mGl is not None and args.mGl not in filelist:
                continue
            if args.mN2 is not None and args.mGl not in filelist:
                continue
            if args.mN1 is not None and args.mGl not in filelist:
                continue
        else:
            if args.HT is not None and args.HT not in filelist and not data:
                continue
            if args.year is not None and args.year not in filelist:
                continue
        if data:
            if args.era not in filelist:
                continue
        print("filelist",filelist)
        inputFileList = inputPathList+filelist

        objName = args.object
        
        printstring = "Skimming for "+objName+" in file "+inputFileList
        #if(objName == "jets"):
        #	printstring += " with strategy "+strategyName 
        print(printstring)
        
        #organize output by sample, object (ie jets or photons), and strategy (for jets only - NlnN or N2)
        #find .root and then find the / before that (if it exists) - everything in between is the file name
        match = "kucmsntuple_"
        sampleName = inputFileList[ inputFileList.find(match)+len(match) : inputFileList.find("_list") ]
        sampleNameShort = sampleName[ : sampleName.find("_")  ]

        dirname = odir+sampleNameShort+"/"+sampleName+"/"+objName
        ofilename = args.inputSample+"_"+objName
        
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
        
        if(args.gev == -999):
            if(args.object == "jets"):
                args.gev = 1./10.
            else:
                args.gev = 1./30.
        

        # grab relevant flags
        filearr = []
        if args.nfilesfirst != -1 or args.nfileslast != -1:
            filearr = SH.filesSplit(inputFileList, args.nfilesfirst, args.nfileslast, args.maxnevts)
        else:
            filearr = SH.filesSplit(inputFileList, args.maxnevts, args.maxnfiles)
        flags = '--EMalpha '+str(EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" --gev "+str(args.gev)+' --minpt '+str(args.minpt)+' --minNrhs '+str(args.minnrhs)+' --minemE '+str(args.minemE)+' --minRhE '+str(args.minRhE)+' --BHFilter '+str(args.beamHaloFilter)+' --beta0 '+str(beta0)+' --m0 '
        for m in m0:
            flags += str(m)+' '
        flags += '--W0diag '
        for w in W0diag:
            flags += str(w)+' '
        flags += '--nu0 '+str(nu0)

        if(args.smear):
        	flags += ' --smear'
        if(args.timeSmear and "AODSIM" in inputFileList):
        	flags += ' --timeSmear'
        if(args.timeSmear and "AODSIM" not in inputFileList):
            print("---timeSmear only for MC")
        if(args.applyFrac):
	        flags += ' --applyFrac'
        if(args.noCalibrate):
        	flags += ' --noCalibrate'
        if(args.rejectSpikes):
            flags += ' --rejectSpikes'
        if(args.noIso):
            flags += ' --noIso'
            flags += ' --maxmet_cr '+str(args.maxmet_cr)+' --minphopt_cr '+str(args.minphopt_cr)+' --minht_cr '+str(args.minht_cr)+' --minjetpt_cr '+str(args.minjetpt_cr)
        if(args.maxRhE != -999):
            flags += ' --maxRhE '+str(args.maxRhE)
        if(args.noSpatCorr):
            flags += ' --noSpatCorr'
        if(args.cleanSubclusters and "AODSIM" not in inputFileList):
            flags += ' --cleanSubclusters'
        if(args.cleanSubclusters and "AODSIM" in inputFileList):
            print("--cleanSubclusters only for data") 
        if(args.cleanMist and "AODSIM" not in inputFileList):
            flags += ' --cleanMist'
        if(args.cleanMist and "AODSIM" in inputFileList):
            print("--cleanMist only for data") 
        if(args.noLumiMask):
            flags += ' --noLumiMask'
        

        if(objName == "jets"):
        	flags += " --object 0"
        elif(objName == "superclusters"):
        	flags += " --object 1"
        elif(objName == "photons"):
        	flags += " --object 2"
        else:
            print("Object",objName,"not recognized")
            exit()

        ##### Create condor submission script in src directory #####
        condorSubmitFile = dirname + "/src/submit.sh"
        subf = open(condorSubmitFile, "w")
        print("outputfile name "+ofilename)
        SH.writeSubmissionBase(subf, dirname, ofilename, args.max_mat, args.max_idle, args.request_memory)
        #need to remove local lpc path for actual args
        inputFileList = inputFileList[inputFileList.rfind("/",0,inputFileList.rfind("/"))+1:]
        print("inputfilelist",inputFileList)
        SH.writeQueueList(subf, ofilename, filearr, flags)
        condor_subs.append(condorSubmitFile)
        print()
    if len(condor_subs) < 1:
        print("No file lists found for input sample",args.inputSample,"and selection",sel)
        exit()
    elif len(condor_subs) < 2:
        print("------------------------------------------------------------")
        print("Submission ready, to run use:")
        print("condor_submit "+condor_subs[0]+"\n")
    else:
        #write bash script
        mult_bash_name = odir+sampleNameShort+"/"
        if "SMS" in args.inputSample:
            mult_bash_name += args.inputSample+"_"
        mult_bash_name += objName+"_"+args.output+"_MultiSub.sh"
        mult_bash = open(mult_bash_name, "w")
        SH.writeMultiSubScript(mult_bash, condor_subs)
        print("------------------------------------------------------------")
        print("Submission ready, to run use:")
        print("source "+mult_bash_name+"\n")



def main():
    # options
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", "-d", default="Output", help="working directory for condor submission")
    parser.add_argument("--max_mat",help='max_materialization condor option (default: off)',default=-1)
    parser.add_argument("--max_idle",help='max_idle condor option (default: off)',default=-1)
    parser.add_argument("--request_memory",help='memory to request from condor scheduler in bits (default = 2048)',default=-1)
    parser.add_argument('--maxnevts',help="maximum number of events to run over",default=-999,type=int)
    parser.add_argument('--maxnfiles',help="maximum number of files total",default=-999,type=int)
    parser.add_argument('--nfilesfirst',help='index of first file to process',default=-1,type=int)
    parser.add_argument('--nfileslast',help='index of last file to process',default=-1,type=int)
    #Ntuple file to run over
    parser.add_argument('-inputSample','-i',help='Ntuple sample to create skims from',required=True,choices=['DisplacedJet','JetHT','EGamma','DoubleEG','GJets','QCD','MET','SMS-GlGl','SMS-SqSq','SMS-GlGlZ'])
    parser.add_argument('--mGl',help='gluino mass for signal',default=None)
    parser.add_argument('--mN2',help='neutralino2 mass for signal',default=None)
    parser.add_argument('--mN1',help='neutralino1 mass for signal',default=None)
    parser.add_argument('--HT',help='ht bin',default=None)
    parser.add_argument('--era',help='era (run) for PDs',default=None)
    parser.add_argument('--selection',help='ntuple preselection or gluino mass',required=True)#choices=['MET100','AL1IsoPho','AL1IsoPhoMET100','MRL_MET100','MRL_None'],required=True)
    parser.add_argument('--output','-o',help='output label')
    parser.add_argument('--year',help='year of sample',default=None,choices=['2016','2017','2018','2022'])
    #which object to analyze (jets or photons currently supported)
    parser.add_argument('--object',help='which object to skim',choices=["jets","superclusters","photons"],required=True)
    #parser.add_argument('--strategy','-st',help='if skimming jets, which strategy to use for BHC (NlnN = 0 default, N2 = 1, GMM only = 2)',default=0,type=int,choices=[2,1,0])
    parser.add_argument('--verbosity','-v',help="verbosity",default=-1)
    
    #add algorithm parameters - emAlpha, priors, verbosity, thresh
    
    parser.add_argument('--thresh','-t',help='threshold for GMM clusters',default=1.)
    parser.add_argument('--gev',help='energy transfer factor (default = 1/10 for jets, 1/30 for photons + superclusters',default=1/10)
    parser.add_argument('--minpt',help='min object pt',default=30.)
    parser.add_argument('--minnrhs',help='min object nrhs',default=15)
    parser.add_argument('--minemE',help='min object ECAL energy',default=30)
    parser.add_argument('--minRhE',help='min rechit ECAL energy',default=0.5)
    parser.add_argument('--maxRhE',help='max rechit ECAL energy (-999 = off)',default=-999)
    parser.add_argument('--smear',help="turn on spatial smearing (turns off measurement error)",default=False,action='store_true')
    parser.add_argument('--timeSmear',help="turn on time smearing",default=False,action='store_true')
    parser.add_argument('--applyFrac',help="apply fractions from hitsAndFractions list to rh energies for photons",default=False,action='store_true')
    parser.add_argument('--noCalibrate',help="turn off channel-by-channel time calibration",default=False,action='store_true')
    parser.add_argument('--rejectSpikes',help="reject spikes based on swiss cross cut (default = false, off)",default=False,action='store_true')
    parser.add_argument('--beamHaloFilter',help="apply BH filter (0: off - default, 1: on, 2: inverse)",default=0)
    parser.add_argument('--noIso',help='turn off isolation for photons in preselection',default=False,action='store_true')
    parser.add_argument('--noLumiMask',help='do not apply lumi mask',default=False,action='store_true')
    parser.add_argument('--isoBkg',help='turn on event selection for isolated background',default=False,action='store_true')
    parser.add_argument('--maxmet_cr',help='max MET for isolated background event selection',default=100)
    parser.add_argument('--minphopt_cr',help='min photon pt for isolated background event selection',default=50)
    parser.add_argument('--minht_cr',help='min HT for isolated background event selection',default=50)
    parser.add_argument('--minjetpt_cr',help='min jet pt for isolated background event selection',default=50)
    parser.add_argument('--noSpatCorr',help='turn off spatial corrections for rechit times to put in PV frame (jets only)',default=False,action='store_true')
    parser.add_argument('--cleanSubclusters',help='clean subclusters from jet time (data, jets only)',default=False,action='store_true')
    parser.add_argument('--cleanMist',help='extreme \'mist\' cleaning in rh time + energy (data, jets only)',default=False,action='store_true')
    args = parser.parse_args()

    generateSubmission(args)

if __name__ == "__main__":
    main()
