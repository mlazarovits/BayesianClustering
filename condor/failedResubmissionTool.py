import glob
import os
import shutil

#print(glob.glob("/home/adam/*"))



#test on this BF dir
#BF_dir = "BF_B135_TChipmWW_MCstats"
###job configuration
#cpu="request_cpus = 4"
#disk="request_disk = 15000000 KB"
#mem="request_memory = 4000 MB"


#BF_dir = "BF_B135_bugfix16_TChiWZ_MCstats"
#BF_dir = "BF_B135_bugfix16_HinoN2C1_MCstats"
#cpu="request_cpus = 4"
#disk="request_disk = 2000000 KB"
#mem="request_memory = 4000 MB"

#BF_dir = "BF_B135_bugfix16_T2tt_MCstats"
#BF_dir = "BF_B135_bugfix16_TChipmWW_MCstats"
BF_dir1 = "test/SMS-T2bW_X05_dM-90to170_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X"
BF_dir2 = "Output/SMS-T2bW_X05_dM-90to170_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X"
cpu="request_cpus = 4"
disk="request_disk = 2400000 KB"
mem="request_memory = 4000 MB"




#subs = glob.glob(BF_dir+"/src/submit.sh")
#print(subs)
#for i,sub in enumerate(subs):
#	print(i, sub)
#	jN = sub.split("/")[-1]
#	jN = jN[:-3]
#	jN = jN.split("_")[-1]
#	subs[i] = jN
f = open(BF_dir1+"/src/submit.sh")
subs = []
for l in f:
	if "job" not in l:
		continue
	if "$(Process)" in l:
		continue
	if "log" in l:
		continue
	sub = l[l.find("job")+3:l.find(" #")]
	subs.append(sub)
print(len(subs),"subs")
logs = glob.glob(BF_dir2+"/log/*.out")
for i,log in enumerate(logs):
	jN = log.split("/")[-1]
        jN = jN[jN.find(".")+1:jN.rfind(".")] 
	logs[i] = jN

print(len(logs),"logs")
subSet = set(subs)
logSet = set(logs)

subDiff = subSet.difference(logSet)
print("Analyzing:", BF_dir2)
print("found", len(subDiff), "missing jobs")
print("generating resubmission scripts for this job list:")
#jobs = ""
#for i in subDiff:
#	jobs += i+" "
print(subDiff)

'''
resub_directory = BF_dir+"/src" 

if os.path.exists(resub_directory):
    shutil.rmtree(resub_directory)
os.mkdir(resub_directory) 
print("Directory '% s' created" % resub_directory) 

print("Using these configurations:")
print(cpu)
print(disk)
print(mem)
print("Generating resubmission scripts...")
for key in subDiff:
#	print("procesing #",key)
	subfile = BF_dir+"/src/submit_"+key+".sh"
	resubfile = BF_dir+"/resub_src/submit_"+key+".sh"
#	print(subfile)
	subIn = open(subfile,"r")
	subOut = open(resubfile,"w")
	lines = subIn.readlines()
	for line in lines:
		if( "request_memory" in line):
			subOut.write(cpu+"\n")
			subOut.write(disk+"\n")
			subOut.write(mem+"\n")

		else:
			subOut.write(line)
	subIn.close()
	subOut.close()

print("Generating submission aggregation script")
resub_agg_file = BF_dir+"/condor_resubmit.sh"
aggfile = open(resub_agg_file, "w")
for key in subDiff:
	line = "condor_submit "+BF_dir+"/resub_src/submit_"+key+".sh\n"
	aggfile.write(line)

aggfile.close()
os.chmod(resub_agg_file, 777)
print("done, to launch datacards do this:")
print("./"+BF_dir+"/condor_resubmit.sh")
	
	
#in here
#vim condor_submit.sh
#condor_submit BF_B135_bugfix16_HinoN2C1_MCstats/src/submit_0.sh

#BF_B135_TChiWZ_MCstats/log/job_114.log.out
#how many sub src files are there?
#ls BF_dir/src
'''
