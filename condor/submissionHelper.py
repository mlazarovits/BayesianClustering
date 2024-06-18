#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################

import os
import re
import ROOT
import numpy as np

# Create directory if it does not exist.
def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

# Create directories for work area.
def createWorkArea(sampleDir):
	os.makedirs(sampleDir)
	os.mkdir(sampleDir+"/src")
	os.mkdir(sampleDir+"/log")
	os.mkdir(sampleDir+"/out")

def getDataSetName(pathToList):
        tmp = pathToList.split('/')
        tmp = tmp[-1].split('.')
        #print(tmp[0])
        return tmp[0]

# Write the header of the condor submit file
def writeSubmissionBase(subf, dirname, ofilename, infile):
        subf.write("universe = vanilla\n")
        subf.write("executable = execute_script.sh\n")
        subf.write("output = ./"+dirname+"/log/job.$(Process).out\n")
        subf.write("error = ./"+dirname+"/log/job.$(Process).err\n")
        subf.write("log = ./"+dirname+"/log/job.log\n")
        #include tarball with CMSSW environment
        #subf.write("transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2, config.tgz, \n")
        subf.write("transfer_input_files = /uscms/home/mlazarov/nobackup/sandboxes/sandbox-CMSSW_13_0_13.tar.bz2, config.tgz, \n")
        subf.write("MY.wantOS=\"el9\"\n")
        subf.write("should_transfer_files = YES\n")
        subf.write("when_to_transfer_output = ON_EXIT\n")
        if "jet" in dirname:
            subf.write("request_memory=4096\n")
        outnames = []
        outnames.append(ofilename+".$(Process).root")
        #if photons, write csv file for MVA
        if "photon" in dirname:
            outnames.append(ofilename+".$(Process).csv")
        outname = ""
        for o in outnames:
            outname += o+", "
        #remove last comma and space
        outname = outname[:-2]
        #subf.write("transfer_output_files = "+outname+"\n")
        subf.write("transfer_output_files = "+outname+"\n")
        # need to supply absolute path for remap
        #absCWD = os.path.abspath(".") # these cwd give the wrong abs path, there is something special in the environment
        #absCWD = os.getcwd()
        absCWD = os.popen('pwd').readline().rstrip()
        #print("abs path is "+ absCWD)
        remap = ""
        for o in outnames:
            remap += o+"="+absCWD+"/"+dirname+"/out/"+o+";"
            #print("remap is "+ remap)
	        #print("outname is "+outname)
        #remove last semicolon 
        remap = remap[:-1]
        print("remap",remap)
        subf.write("transfer_output_remaps = \""+remap+"\"\n")	

#splits by event number
def eventsSplit(infile, nChunk):
    if nChunk == 0:
        nChunk += 1
    print("Splitting each file into "+str(nChunk)+" jobs ")
    #should split by event number in file
    rfile = ROOT.TFile.Open(infile)
    tree = rfile.Get("tree/llpgtree")
    nevts = tree.GetEntries()
    evts = range(nevts+1)
    #return array of pairs of evtFirst and evtLast to pass as args into the exe to run
    #make sure to count first event in every chunk after first
    arr = [[min(i)-1, max(i)] for i in np.array_split(evts,nChunk)]
    #set first entry to 0
    arr[0][0] = 0
    return arr



# Write each job to the condor submit file.
def writeQueueList( subf, inFile, ofilename, evts, flags ):
    if evts == 0 or evts is None:
        print("No events found")
        return
    #.root is set in exe
    outFileArg = ofilename+".$(Process)"
    
    #infile should only be file name (no path)
    #inFile = inFile[inFile.rfind("/")+1:]
    #inFile = "root://cmsxrootd.fnal.gov/"+inFile
    
    jobCtr=0
    for e in evts:
            inFileArg = " -i "+inFile
            Args = "Arguments ="+inFileArg+" "+flags+" --evtFirst "+str(e[0])+" --evtLast "+str(e[1])+" -o "+outFileArg+"\n"
            subf.write("\n\n\n")
            subf.write("###### job"+str(jobCtr)+ "######\n")
            subf.write(Args)
            subf.write("Queue\n")
            jobCtr=jobCtr+1
