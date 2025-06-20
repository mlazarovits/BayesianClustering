
import argparse

################
parser = argparse.ArgumentParser()
parser.add_argument("--directory","-d",help="Directory with src/submit.sh to resubmit",required=True)
parser.add_argument("--jobs","-j",help="job numbers of jobs to resubmit",required=True,nargs='+')

args = parser.parse_args()
datasetname = args.directory
joblist = []
for j in args.jobs:
	joblist.append(int(j)) 
###############



def parseHeader():
    #get header
    global datasetname
    fsub = open('./'+datasetname+'/src/submit.sh')
    header=[]
    proc_lines = []
    proc_order = {}
    #parse from beginning until '######'
    for line in fsub:
        if 'queue' in line:
        	break
        if "output = ." in line:
            proc_lines.append(line)
            continue
        if "error =" in line:
            proc_lines.append(line)
            continue
        if "transfer_output_files =" in line:
            proc_lines.append(line)
            continue
        if "transfer_output_remaps =" in line:
            proc_lines.append(line)
            continue
        header.append(line)
        if "EXIT" in line:
            header.append("request_memory=4096\n")
    #transfer output_remaps needs to be at the beginning
    transf = proc_lines.pop()
    proc_lines.insert(0,transf)
    return header, proc_lines

def parseJobNum( jobnumber ):
	#get the specified job
	global datasetname
	fsub = open('./'+datasetname+'/src/submit.sh')
	job=[]
	jobstr = '###### job'+str(jobnumber)
	jobnp1 = '###### job'+str(jobnumber+1)
	inJobArgs=False
	for line in fsub:
		if jobstr in line:
			inJobArgs = True
		if jobnp1 in line:
			inJobArgs = False
			break
		if inJobArgs:
			job.append(line)
	
	return job
def copyAndUpdateJob( job, jobnumber , proc_lines):
    updatedJob = []
    for line in proc_lines:
        line = line[line.find(" = ")+3:]
        job.append(line)
    for line in job:
    	if '$(Process)' in line:
    		updatedJob.append( line.replace('$(Process)', str(jobnumber)) )
    	else:
    		updatedJob.append(line)
    return updatedJob



header, proclines = parseHeader()	
#print(proclines)
fresub = open('./'+datasetname+'/src/resubmit_test2.sh','w')
for line in header:
	fresub.write(line)
fresub.write("queue transfer_output_remaps, output, error, transfer_output_files, Arguments from (")
for jobnumber in joblist:
    updatedJob = copyAndUpdateJob( parseJobNum(jobnumber), jobnumber, proclines)
    args = updatedJob.pop(1)
    updatedJob.append(args)
    fresub.write("\n"+updatedJob[0])
    line = ""
    for l in updatedJob[1:]:
        print("job line",l[:-1])
        line += l[:-1]+","
    #print("full line",line)
    fresub.write(line)
fresub.write(")\n")
fresub.close()

print("Resubmission script ready! Launch with:")
print("condor_submit ./"+datasetname+"/src/resubmit_test2.sh")
