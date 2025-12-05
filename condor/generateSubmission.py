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

    ver = "v31"
    sel = args.selection#"MET100"
    yr = str(args.year)[-2:]
    reco_date = {}
    reco_date["2016_MET"] = "-21Feb2020_UL2016_HIPM-v1" 
    reco_date["2017"] = "_17Nov2017"
    reco_date["2017_MET"] = "-09Aug2019_UL2017_rsb-v1"
    reco_date["2017_DEG"] = "-09Aug2019_UL2017-v1"
    reco_date["2018_DEG"] = "-15Feb2022"
    reco_date["2018_MET"] = "-15Feb2022_UL2018-v1"
    reco_date["2022_MET"] = "-27Jun2023-v2"
    reco_date["2018_MC"] = "RunIISummer20UL18"
    reco_date["2017_MC"] = "RunIIFall17DRPremix"
    reco_date["2018"] = ""
    inputFilterList = ""
    if "GJets" in args.inputSample:
        inputFileList = "kucmsntuple_GJets_R"+yr+"_"+sel+"_"+ver+"_GJets_HT-"+args.HT+"_TuneCP5_AODSIM_"+reco_date[args.year+"_MC"]+"RECO_list.txt"
    elif "MET" in args.inputSample:
        inputFileList = "kucmsntuple_MET_R"+yr+"_"+sel+"_"+ver+"_MET_AOD_Run20"+yr+args.era+reco_date[args.year+"_MET"]+"_list.txt"
    elif "JetHT" in args.inputSample:
    	inputFileList = "kucmsntuple_JetHT_R"+yr+"_"+sel+"_"+ver+"_JetHT_AOD_Run20"+yr+args.era+reco_date[args.year+"_MET"]+"_list.txt"
    elif "QCD" in args.inputSample:
    	inputFileList = "kucmsntuple_QCD_R"+yr+"_"+sel+"_"+ver+"_QCD_HT"+args.HT+"_TuneCP5_AODSIM_"+reco_date[args.year+"_MC"]+"RECO_list.txt"
    elif "EGamma" in args.inputSample:
        if "AL1SelEle" in sel:
            inputFileList = "kucmsntuple_EGamma_R"+yr+"_"+sel+"_"+ver+"_EGamma_AOD_Run20"+yr+args.era+reco_date[args.year+"_DEG"]+"_UL2018-v1_list.txt"
        else:
            inputFileList = "kucmsntuple_EGamma_R"+yr+"_"+sel+"_"+ver+"_EGamma_AOD_Run20"+yr+args.era+"_list.txt"
    elif "DoubleEG" in args.inputSample:
        inputFileList = "kucmsntuple_DoubleEG_R"+yr+"_"+sel+"_"+ver+"_DoubleEG_AOD_Run20"+yr+args.era+reco_date[args.year+"_DEG"]+"_list.txt"
    elif "gogo" in args.inputSample:
        inputFileList = "kucmsntuple_SMS_Sig_"+sel+"_"+ver+"_SMS-GlGl_AODSIM_mGl-"+args.mGl+"_mN2-"+args.mN2+"_mN1-"+args.mN1+"_list.txt"
    elif "sqsq" in args.inputSample:
        inputFileList = "kucmsntuple_SMS_Sig_"+sel+"_"+ver+"_SMS-SqSq_AODSIM_mGl-"+args.mGl+"_mN2-"+args.mN2+"_mN1-"+args.mN1+"_list.txt"
    else:
    	print("Sample "+args.inputSample+" not found")
    	exit()
    #to use xrootd path cannot be relative
    #find any ../ and remove it and the dir before it
    inputPathList = "/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/BayesianClustering/filelists/"
    #inputPathList = "filelists/"
    inputFileList = inputPathList+inputFileList

    objName = args.object
    #strategyName = "GMMonly" #only option for CMS jets/photons
    
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
   
    #prior parameters
    priorname = ""
    #add emAlpha to output name
    emAlphastr = str(args.EMalpha)
    emAlphastr = emAlphastr.replace(".","p")

    betastr = str(args.beta0)
    betastr = betastr.replace(".","p")
    priorname = priorname+"_beta0-"+betastr

    #m0
    #check dims
    if len(args.m0) != 3:
        print("Error: m0 must be a vector of length 3")
        exit()
    mstr = ""
    for i, m in enumerate(args.m0):
        m = float(m)
        m = round(m,3)
        m = str(m)
        m = m.replace(".","p")
        if i == len(args.m0)-1:
            mstr += m
        else:
            mstr += m+"-"
    priorname = priorname+"_m0-"+mstr

    #W0
    if len(args.W0diag) != 3:
        print("Error: W0diag must be a vector of length 3")
        exit()
    Wstr = ""
    for i, w in enumerate(args.W0diag):
        w = float(w)
        w = round(w,3)
        w = str(w)
        w = w.replace(".","p")
        if i == len(args.W0diag)-1:
            Wstr += w
        else:
            Wstr += w+"-"
    priorname = priorname+"_W0diag-"+Wstr
    
    #nu0
    nustr = str(args.nu0)
    nustr = nustr.replace(".","p")
    priorname = priorname+"_nu0-"+nustr


    #add npergev to output name
    gevstr = str(args.gev)
    gevstr = gevstr.replace(".","p")
    priorname += "_NperGeV-"+gevstr
    
    #add emalpha to output name
    emalphastr = str(args.EMalpha)
    emalphastr = emalphastr.replace(".","p")
    priorname += "_emAlpha-"+emalphastr

    #if(objName == "jets"):
    #	dirname += "/"+strategyName
    #	ofilename += "_"+strategyName
    #strategy doesn't apply to photons (GMM only)
    
    if args.output is not None and objName not in args.output:
    	ofilename = ofilename+"_"+args.output 
    	dirname = dirname+"_"+args.output
    #dirname += priorname
    #ofilename += priorname
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
    filearr = SH.filesSplit(inputFileList, args.maxnevts, args.maxnfiles)
    flags = '--EMalpha '+str(args.EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" --gev "+str(args.gev)+' --minpt '+str(args.minpt)+' --minNrhs '+str(args.minnrhs)+' --minemE '+str(args.minemE)+' --minRhE '+str(args.minRhE)+' --BHFilter '+str(args.beamHaloFilter)+' --beta0 '+str(args.beta0)+' --m0 '
    for m in args.m0:
        flags += str(m)+' '
    flags += '--W0diag '
    for w in args.W0diag:
        flags += str(w)+' '
    flags += '--nu0 '+str(args.nu0)

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
    #SH.writeQueueList(subf, inputFileList, ofilename, eventnums, flags)
    SH.writeQueueList(subf, ofilename, filearr, flags)
    
    print("------------------------------------------------------------")
    print("Submission ready, to run use:")
    print("condor_submit "+condorSubmitFile)

def main():
    # options
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", "-d", default="Output", help="working directory for condor submission")
    parser.add_argument("--max_mat",help='max_materialization condor option (default: off)',default=-1)
    parser.add_argument("--max_idle",help='max_idle condor option (default: off)',default=-1)
    parser.add_argument("--request_memory",help='memory to request from condor scheduler in bits (default = 2048)',default=-1)
    parser.add_argument('--maxnevts',help="maximum number of events to run over",default=-999,type=int)
    parser.add_argument('--maxnfiles',help="maximum number of files total",default=-999,type=int)
    #Ntuple file to run over
    parser.add_argument('-inputSample','-i',help='Ntuple sample to create skims from',required=True,choices=['JetHT','EGamma','DoubleEG','GJets','QCD','MET','gogo','sqsq'])
    parser.add_argument('--mGl',help='gluino mass for signal',default='2000')
    parser.add_argument('--mN2',help='neutralino2 mass for signal',default='1950')
    parser.add_argument('--mN1',help='neutralino1 mass for signal',default='1')
    parser.add_argument('--HT',help='ht bin',default='100to200')
    parser.add_argument('--era',help='era (run) for PDs',default='RunF')
    parser.add_argument('--selection',help='ntuple preselection or gluino mass',required=True)#choices=['MET100','AL1IsoPho','AL1IsoPhoMET100','MRL_MET100','MRL_None'],required=True)
    parser.add_argument('--output','-o',help='output label')
    parser.add_argument('--year',help='year of sample',default='2017',choices=['2016','2017','2018','2022'])
    #which object to analyze (jets or photons currently supported)
    parser.add_argument('--object',help='which object to skim',choices=["jets","superclusters","photons"],required=True)
    #parser.add_argument('--strategy','-st',help='if skimming jets, which strategy to use for BHC (NlnN = 0 default, N2 = 1, GMM only = 2)',default=0,type=int,choices=[2,1,0])
    parser.add_argument('--verbosity','-v',help="verbosity",default=-1)
    
    #add algorithm parameters - emAlpha, priors, verbosity, thresh
    parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=1e-5)
    parser.add_argument('--beta0',help="beta0 prior",default=1e-5)
    parser.add_argument('--m0',help="m0 prior (must be n-dim entries)",default=[0,0,0],nargs="+")
    parser.add_argument('--W0diag',help="diagonal elements of W0 prior (must be n-dim entries)",default=[0.33333333,0.33333333,33.333333e-5],nargs="+")
    parser.add_argument('--nu0',help="nu0 prior",default=3)
    
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
