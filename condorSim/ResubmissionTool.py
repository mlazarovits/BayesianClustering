
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
    #parse from beginning until '######'
    for line in fsub:
        if '######' in line:
        	break
        header.append(line)
        if "EXIT" in line:
            header.append("request_memory=4096\n")
    return header

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
def copyAndUpdateJob( job, jobnumber ):
#	print('updating',job)
	updatedJob = []
	for line in job:
#		print('line',line)
		if '$(Process)' in line:
			updatedJob.append( line.replace('$(Process)', str(jobnumber)) )
		else:
			updatedJob.append(line)
	return updatedJob



header = parseHeader()	
#print(header)
fresub = open('./'+datasetname+'/src/resubmit.sh','w')
updatedJobHeader = copyAndUpdateJob( header, joblist[0])
for line in updatedJobHeader:
	fresub.write(line)
for jobnumber in joblist:
	updatedJob = copyAndUpdateJob( parseJobNum(jobnumber), jobnumber)
	for line in updatedJob:
		fresub.write(line)
fresub.write(")")
fresub.close()

print("Resubmission script ready! Launch with:")
print("condor_submit ./"+datasetname+"/src/resubmit.sh")
