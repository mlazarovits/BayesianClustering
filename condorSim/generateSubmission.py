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
    
    #make sure ntuple names are updated for latest version otherwise skimmer might crash
    if args.inputSample == "ttbar":
            inputFile = "condorSimNtuples_ttbar_defaultv9p8"
    elif args.inputSample == "QCD":
            inputFile = "condorSimNtuples_QCD_defaultv9p1"
    #elif args.inputSample == "QCD_noSpatialSmear":
    #        inputFile = "condorSimNtuples_QCD_defaultv4_noSpatialSmear.root"
    #elif args.inputSample == "QCD_noSpatialSmear_highEnergySmear":
    #        inputFile = "condorSimNtuples_QCD_defaultv4_noSpatialSmear_highEnergySmear.root"
    #elif args.inputSample == "ttbar_v3_noEnergySmear":
    #        inputFile = "condorSimNtuples_ttbar_v3_noEnergySmear.root"
    #elif args.inputSample == "ttbar_v3_noEnergySmear_clusterFromRecoParticles":
    #        inputFile = "condorSimNtuples_ttbar_v3_noEnergySmear_clusterFromRecoParticles.root"
    #elif args.inputSample == "ttbar_noSpatialSmear":
    #        inputFile = "condorSimNtuples_ttbar_defaultv4_noSpatialSmear.root"
    #elif args.inputSample == "ttbar_noSpatialSmear_highEnergySmear":
    #        inputFile = "condorSimNtuples_ttbar_defaultv4_noSpatialSmear_highEnergySmear.root"
    else:
                print("Sample "+args.inputSample+" not found")
                exit()
    if(args.zeroSup != '0.5'):
        zeroSupStr = zeroSupStr.replace(".","p")
        inputFile += "_eThresh-"+zeroSupStr

    inputFile += ".root"
    inputFile = "root://cmseos.fnal.gov//store/user/mlazarov/SimNtuples/"+inputFile

    print("inputFile",inputFile)
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
    kname = "%.3f" % float(args.alpha)
    kname = kname.replace(".","p")
    paramsname = "_bhcAlpha"+kname
    kname = "%.3f" % float(args.EMalpha)
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


    priorname = ""
    #beta0
    betastr = str(args.beta0)
    betastr = betastr.replace(".","p")
    priorname = priorname+"beta0-"+betastr

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

    ofilename += "_"+priorname+"_"+args.strategy
    dirname += "_"+priorname+"_"+args.strategy
    print("Preparing sample directory: {0}".format(dirname))
    ##### Create a workspace (remove existing directory) #####
    if os.path.exists(dirname):
    	print("Removing existing directory: {0}".format(dirname))
    	shutil.rmtree(dirname)
    
    # Create directories for work area.
    SH.createWorkArea(dirname)
    
    # grab relevant flags
    eventnums = SH.eventsSplit(inputFile, args.split)
    if eventnums is None:
        return
    flags = '--alpha '+str(args.alpha)+' --EMalpha '+str(args.EMalpha)+' -v '+str(args.verbosity)+' -t '+str(args.thresh)+" --gev "+str(k)+' --minpt '+str(args.minpt)+' --minE '+str(args.minE)+' --minNrhs '+str(args.minnrhs)+' --minRhE '+str(args.minRhE)+' --minNconsts '+str(args.minNconsts)+" --minTopPt "+str(args.minTopPt)

    flags += ' --beta0 '+str(args.beta0)+' --m0 '
    for m in args.m0:
        flags += str(m)+' '
    flags += '--W0diag '
    for w in args.W0diag:
        flags += str(w)+' '
    flags += '--nu0 '+str(args.nu0)
    
    if(args.smear):
    	flags += ' --smear'
    if(args.timeSmear):
    	flags += ' --timeSmear'
    if(args.applyFrac):
    	flags += ' --applyFrac'
    if(args.checkMerges):
        flags += ' --checkMerges'
    flags += ' --tResCte '+str(args.tResCte)
    flags += ' --tResStoch '+str(args.tResStoch)
    flags += ' --tResNoise '+str(args.tResNoise) 

    strategyMap = {}
    strategyMap["NlnN"] = 0
    strategyMap["N2"] = 1
    strategyMap["GMMonly"] = 2
    strategyMap["NlnNonAK4"] = 3

    flags += ' --strategy '+str(strategyMap[args.strategy])

    ##### Create condor submission script in src directory #####
    condorSubmitFile = dirname + "/src/submit.sh"
    subf = open(condorSubmitFile, "w")
    print("outputfile name "+ofilename)
    SH.writeSubmissionBase(subf, dirname, ofilename)
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
    parser.add_argument('--inputSample','-i',help='Ntuple sample to create skims from',required=True,choices=['ttbar','QCD'])
    parser.add_argument('--zeroSup',help='min rechit energy at ntuple level for reco (zero suppression)',default='0.5')
    parser.add_argument('--output','-o',help='output label')
    parser.add_argument('--strategy','-st',help='which strategy to use for BHC (default = NlnN)',default='NlnN',choices=['NlnN','N2','GMMonly','NlnNonAK4'])
    parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
    parser.add_argument('--verbosity','-v',help="verbosity",default=0)
    #add algorithm parameters - alpha, emAlpha, verbosity, thresh
    parser.add_argument('--alpha','-a',help="alpha for BHC",default=0.1)
    parser.add_argument('--EMalpha','-EMa',help="alpha for GMM (EM algo)",default=1.)
    parser.add_argument('--beta0',help="beta0 prior",default=1e-3)
    parser.add_argument('--m0',help="m0 prior (must be n-dim entries)",default=[0,0,0],nargs="+")
    parser.add_argument('--W0diag',help="diagonal elements of W0 prior (must be n-dim entries)",default=[0.33333333,0.33333333,33.333333e-5],nargs="+")
    parser.add_argument('--nu0',help="nu0 prior",default=3)
    parser.add_argument('--thresh','-t',help='threshold for GMM clusters',default=1.)
    parser.add_argument('--gev',help='energy transfer factor',default='1/10')
    parser.add_argument('--tResCte',help='set time smearing constant parameter in ns',default=0.1727)
    parser.add_argument('--tResStoch',help='set time smearing stochastic parameter in ns',default=0.5109)
    parser.add_argument('--tResNoise',help='set time smearing noise parameter in ns',default=2.106)
    parser.add_argument('--minnrhs',help='min object nrhs',default=2)
    parser.add_argument('--minpt',help='min gen jet pt',default=0.)
    parser.add_argument('--minE',help='min gen jet energy',default=0)
    parser.add_argument('--minTopPt',help='min gen top pt',default=0.)
    parser.add_argument('--minRhE',help='min rechit energy',default=0.5)
    parser.add_argument('--smear',help="turn on spatial smearing",default=False,action='store_true')
    parser.add_argument('--checkMerges',help="turn on subcluster merging for BHC jets",default=False,action='store_true')
    parser.add_argument('--timeSmear',help="turn on time smearing",default=False,action='store_true')
    parser.add_argument('--applyFrac',help="apply fractions from hitsAndFractions list to rh energies for photons",default=False,action='store_true')
    parser.add_argument('--minNconsts',help='set minimum number of constituents for gen jets',default=5)
    args = parser.parse_args()
    
    generateSubmission(args)

if __name__ == "__main__":
    main()
