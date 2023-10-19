#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################
import ROOT
import os
import re
import numpy as np
def writeSubmissionBase( subf, dirname, infile ):
	subf.write("universe = vanilla\n")
	subf.write("executable = execute_script.sh\n")
	subf.write("output = ./"+dirname+"/log/job.$(Process).out\n")
	subf.write("error = ./"+dirname+"/log/job.$(Process).err\n")
	subf.write("log = ./"+dirname+"/log/job.log\n")
	#include tarball with CMSSW environment
	subf.write("transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2, config.tgz, "+infile+"\n")
	subf.write("should_transfer_files = YES\n")
	subf.write("when_to_transfer_output = ON_EXIT\n")
	outname = dirname+"/out/skim.$(Process).root"
	subf.write("transfer_output_files = "+outname+"\n")
	# need to supply absolute path for remap
	#absCWD = os.path.abspath(".") # these cwd give the wrong abs path, there is something special in the environment
	#absCWD = os.getcwd()
	absCWD = os.popen('pwd').readline().rstrip() 
	#print("abs path is "+ absCWD)
	remap= absCWD+"/"+outname
	#print("remap is "+ remap)
	subf.write("transfer_output_remaps = \""+outname+"="+remap+"\"\n")

#splits by event number
def eventsSplit(infile, nChunk):
	if nChunk == 0:
		nChunk += 1
	print("Splitting each file into "+str(nChunk)+" jobs ")
	#should split by event number in file
	rfile = ROOT.TFile(infile)
	nevts = rfile.Get("tree/llpgtree").GetEntries()
	evts = range(nevts+1)
	#return array of pairs of evtFirst and evtLast to pass as args into the exe to run
	arr = [[min(i), max(i)] for i in np.array_split(evts,nChunk)] 
	return arr

def writeQueueList( subf, inFile, ofilename, evts, flags ):
	configArgs = " --skim"
	outFileArg = " -o "+ofilename+".$(Process).root"

	#infile should only be file name (no path)
	inFile = inFile[inFile.rfind("/")+1:]
	
	jobCtr=0
	for e in evts:
		inFileArg = " -i "+inFile
		Args = "Arguments ="+inFileArg+" "+configArgs+" "+flags+" --evtFirst "+str(e[0])+" --evtLast "+str(e[1])+"\n"
		subf.write("\n\n\n")
		subf.write("###### job"+str(jobCtr)+ "######\n")
		subf.write(Args)
		subf.write("Queue\n")	
		jobCtr=jobCtr+1

