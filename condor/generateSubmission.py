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

    if "GJets_HT400to600" in args.inputSample:
    	inputFile = "GJets_R17_MET100_v24_GJets_HT-400To600_AODSIM_RunIIFall17DRPremix.root"
    elif "GJets_HT100to200" in args.inputSample:
    	inputFile = "GJets_R17_MET100_v24_GJets_HT-100To200_AODSIM_RunIIFall17DRPremix.root"
    elif "GJets_HT200to400" in args.inputSample:
    	inputFile = "GJets_R17_MET100_v24_GJets_HT-200To400_AODSIM_RunIIFall17DRPremix.root"
    elif "GJets_HT40to100" in args.inputSample:
    	inputFile = "GJets_R17_MET100_v24_GJets_HT-40To100_AODSIM_RunIIFall17DRPremix.root"
    elif "GJets_HT600toInf" in args.inputSample:
    	inputFile = "GJets_R17_MET100_v24_GJets_HT-600ToInf_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-250_Ctau-10" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-250TeV_Ctau-10cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-300_Ctau-400" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-300TeV_Ctau-400cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-350_Ctau-0p1" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v22_GMSB_L-350TeV_Ctau-0_1cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-350_Ctau-200" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v22_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-400_Ctau-800" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-400TeV_Ctau-800cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-300_Ctau-1000" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-300TeV_Ctau-1000cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-300_Ctau-600" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-300TeV_Ctau-600cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-350_Ctau-10" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-350TeV_Ctau-10cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-350_Ctau-800" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v22_GMSB_L-350TeV_Ctau-800cm_AODSIM_RunIIFall17DRPremix.root"
    elif "GMSB_L-400_Ctau-200" in args.inputSample:
        inputFile = "GMSB_R17_MET100_v21_GMSB_L-400TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix.root"
    elif "MET_RunE" in args.inputSample:
    	inputFile = "MET_R17_MET100_v21_MET_AOD_Run2017E_17Nov2017.root"
    elif "JetHT_RunF_2017" in args.inputSample:
    	inputFile = "JetHT_R17_MET100_v21_JetHT_AOD_Run2017F_17Nov2017.root"
    elif "JetHT_RunC_2018" in args.inputSample:
    	inputFile = "JetHT_R18_MET100_v22_JetHT_AOD_Run2018C_15Feb2022_UL2018.root"
    elif "QCD_HT500to700" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT500to700_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT1000to1500" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT1000to1500_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT100to200" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT100to200_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT1500to2000" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT1500to2000_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT2000toInf" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT2000toInf_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT300to500" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT300to500_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT500to700" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT500to700_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT700to1000" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT700to1000_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT200to300" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT200to300_AODSIM_RunIIFall17DRPremix.root"
    elif "QCD_HT50to100" in args.inputSample:
    	inputFile = "QCD_R17_MET100_v24_QCD_HT50to100_AODSIM_RunIIFall17DRPremix.root"
    elif "EGamma_RunF" in args.inputSample:
    	inputFile = "DEG_R17_MET100_v21_DoubleEG_AOD_Run2017F_09Aug2019_UL2017.root"
    elif "SMS-GlGl" in args.inputSample:
        inputFile = "SMS_GlGl_v23_justin_mc_noFilter_AODSIM_private.root"
    else:
    	print("Sample "+args.inputSample+" not found")
    	exit()
    #to use xrootd path cannot be relative
    #find any ../ and remove it and the dir before it
    if "PhoSlim" in args.inputSample:
        inputPath = "root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/"
        inputFile = inputFile.replace("MET100_v21","AL1IsoPho_v22")
    else:
        inputPath = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSNtuples/"
    inputFile = inputPath+inputFile

    objName = args.object
    #strategyName = "GMMonly" #only option for CMS jets/photons
    
    printstring = "Skimming for "+objName+" in file "+inputFile
    #if(objName == "jets"):
    #	printstring += " with strategy "+strategyName 
    print(printstring)
    
    #organize output by sample, object (ie jets or photons), and strategy (for jets only - NlnN or N2)
    #find .root and then find the / before that (if it exists) - everything in between is the file name
    sampleName = inputFile[ inputFile.rfind("/")+1 : inputFile.find(".root") ]
    sampleNameShort = sampleName[ : sampleName.find("_") ]
    
    dirname = odir+sampleNameShort+"/"+sampleName+"/"+objName
    ofilename = args.inputSample+"_"+objName
    
    #prior parameters
    priorname = ""
    #add emAlpha to output name
    emAlphastr = str(args.EMalpha)
    emAlphastr = emAlphastr.replace(".","p")
    ofilename = ofilename+"_emAlpha"+emAlphastr

    #beta0
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
        print('m',m)
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
        print('w',w)
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

    #if(objName == "jets"):
    #	dirname += "/"+strategyName
    #	ofilename += "_"+strategyName
    #strategy doesn't apply to photons (GMM only)
    
    if args.output is not None and objName not in args.output:
    	ofilename = ofilename+"_"+args.output 
    	dirname = dirname+"_"+args.output
    dirname += priorname
    ofilename += priorname
    ofilename = "condor_"+ofilename
    print("ofilename",ofilename)
    print("dirname",dirname)
    exit()
    
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
    eventnums = SH.eventsSplit(inputFile, args.split)
    flags = '--EMalpha '+str(args.EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" --gev "+str(args.gev)+' --minpt '+str(args.minpt)+' --minNrhs '+str(args.minnrhs)+' --minemE '+str(args.minemE)+' --minRhE '+str(args.minRhE)+' --BHFilter '+str(args.beamHaloFilter)
    if(args.noSmear):
    	flags += ' --noSmear'
    if(args.timeSmear):
    	flags += ' --timeSmear'
    if(args.applyFrac):
	    flags += ' --applyFrac'
    if(args.noCalibrate):
    	flags += ' --noCalibrate'
    if(args.rejectSpikes):
        flags += ' --rejectSpikes'
    if(args.noIso):
        flags += ' --noIso'
    if(args.maxRhE != -999):
        flags += ' --maxRhE '+str(args.maxRhE)

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
    parser.add_argument('-inputSample','-i',help='Ntuple sample to create skims from',required=True,choices=['GJets_HT400to600','GJets_HT100to200','GJets_HT200to400','GJets_HT40to100','GJets_HT600toInf','GMSB_L-350_Ctau-200','GMSB_L-250_Ctau-10','GMSB_L-300_Ctau-400','GMSB_L-350_Ctau-0p1','GMSB_L-400_Ctau-800','GMSB_L-300_Ctau-1000','GMSB_L-300_Ctau-600','GMSB_L-350_Ctau-10','GMSB_L-350_Ctau-800','GMSB_L-400_Ctau-200','MET_RunE','JetHT_RunF_2017','JetHT_RunC_2018','EGamma_RunF','QCD_HT200to1500','QCD_HT100to200','QCD_HT1500to2000','QCD_HT2000toInf','QCD_HT200to300','QCD_HT50to100','QCD_HT700to1000','QCD_HT300to500','QCD_HT500to700','QCD_HT200to300','QCD_HT50to100','QCD_HT1000to1500','GJets_HT400to600_PhoSlim','GJets_HT100to200_PhoSlim','GJets_HT200to400_PhoSlim','GJets_HT40to100_PhoSlim','GJets_HT600toInf_PhoSlim','MET_RunE_PhoSlim','JetHT_RunF_PhoSlim','EGamma_RunF_PhoSlim','QCD_HT200to1500_PhoSlim','QCD_HT100to200_PhoSlim','QCD_HT1500to2000_PhoSlim','QCD_HT2000toInf_PhoSlim','QCD_HT200to300_PhoSlim','QCD_HT50to100_PhoSlim','QCD_HT700to1000_PhoSlim','QCD_HT300to500_PhoSlim','QCD_HT500to700_PhoSlim','QCD_HT200to300_PhoSlim','QCD_HT50to100_PhoSlim','QCD_HT1000to1500_PhoSlim','SMS-GlGl'])
    parser.add_argument('--output','-o',help='output label')
    parser.add_argument('--year',help='year of sample',default=2017)
    #which object to analyze (jets or photons currently supported)
    parser.add_argument('--object',help='which object to skim',choices=["jets","superclusters","photons"],required=True)
    #parser.add_argument('--strategy','-st',help='if skimming jets, which strategy to use for BHC (NlnN = 0 default, N2 = 1, GMM only = 2)',default=0,type=int,choices=[2,1,0])
    parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
    parser.add_argument('--verbosity','-v',help="verbosity",default=0)
    
    #add algorithm parameters - emAlpha, priors, verbosity, thresh
    parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=0.5)
    parser.add_argument('--beta0',help="beta0 prior",default=1e-3)
    parser.add_argument('--m0',help="m0 prior (must be n-dim entries)",default=[0,0,0],nargs="+")
    parser.add_argument('--W0diag',help="diagonal elements of W0 prior (must be n-dim entries)",default=[0.33333333,0.33333333,0.33333333],nargs="+")
    parser.add_argument('--nu0',help="nu0 prior",default=3)
    
    parser.add_argument('--thresh','-t',help='threshold for GMM clusters',default=1.)
    parser.add_argument('--gev',help='energy transfer factor (default = 1/10 for jets, 1/30 for photons + superclusters',default=-999)
    parser.add_argument('--minpt',help='min object pt',default=30.)
    parser.add_argument('--minnrhs',help='min object nrhs',default=15)
    parser.add_argument('--minemE',help='min object ECAL energy',default=30)
    parser.add_argument('--minRhE',help='min rechit ECAL energy',default=0.5)
    parser.add_argument('--maxRhE',help='max rechit ECAL energy (-999 = off)',default=-999)
    parser.add_argument('--noSmear',help="turn off spatial smearing",default=False,action='store_true')
    parser.add_argument('--timeSmear',help="turn on time smearing",default=False,action='store_true')
    parser.add_argument('--applyFrac',help="apply fractions from hitsAndFractions list to rh energies for photons",default=False,action='store_true')
    parser.add_argument('--noCalibrate',help="turn off channel-by-channel time calibration",default=False,action='store_true')
    parser.add_argument('--rejectSpikes',help="reject spikes based on swiss cross cut (default = false, off)",default=False,action='store_true')
    parser.add_argument('--beamHaloFilter',help="apply BH filter (0: off - default, 1: on, 2: inverse)",default=0)
    parser.add_argument('--noIso',help='turn off isolation for photons in preselection',default=False,action='store_true')
    args = parser.parse_args()

    generateSubmission(args)

if __name__ == "__main__":
    main()
