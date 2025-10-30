import glob
import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--directory","-d",help="Directory with src/submit.sh to resubmit",required=True)

args = parser.parse_args()
BF_dir = args.directory





f = open(BF_dir+"/src/submit.sh")
subs = []
for l in f:
    if "job" not in l:
    	continue
    if "$(Process)" in l:
    	continue
    if "log" in l:
    	continue
    sub = l[l.find(" job")+4:]
    sub = sub[:sub.find("#")]
    subs.append(sub)
print(len(subs),"subs")
logs = glob.glob(BF_dir+"/log/*.out")
#logs = glob.glob(BF_dir+"/out/*.csv")
for i,log in enumerate(logs):
    jN = log.split("/")[-1]
    jN = jN[jN.find(".")+1:jN.rfind(".")]
    logs[i] = jN

print(len(logs),"logs")
csvs = glob.glob(BF_dir+"/out/*.csv")
for i,csv in enumerate(csvs):
    jN = csv.split("/")[-1]
    jN = jN[jN.find(".")+1:jN.rfind(".")]
    csvs[i] = jN

print(len(csvs),"csvs")

subSet = set(subs)
logSet = set(logs)
csvSet = set(csvs)

subDiff = subSet.difference(logSet)
print("Analyzing:", BF_dir)
print("found", len(subDiff), "missing logs")
jobs = ""
for i in subDiff:
	jobs += i+" "
if(len(csvs) > 0):
    subDiff = subSet.difference(csvSet)
    print("found", len(subDiff), "missing csvs")
    for i in subDiff:
    	jobs += i+" "
print(jobs)
