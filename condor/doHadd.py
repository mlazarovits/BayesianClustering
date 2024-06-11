import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir","-d",help="top directory that has directories of root files to hadd i.e. Output/GMSB/GMSB_AOD_v13_GMSB_L-150TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix/jets/",required=True)
    parser.add_argument("--force",help='force remake of outfile',action='store_true')
    parser.add_argument("--big",help='run routine for larger rootfiles',action='store_true')
    parser.add_argument("--dryRun",help="dry run, print commands without running",action="store_true")
    args = parser.parse_args()

	
    cmdHadd = "hadd -d /uscmst1b_scratch/lpc1/3DayLifetime/mlazarov/ -j 4"
    for d in os.scandir(args.dir):
        if not os.path.exists(d.path+"/out"):
            continue
        oname = d.name
        proc = d.path[d.path.find("/")+1:]
        proc = proc[proc.find("/")+1:]
        proc = proc[:proc.find("/")]
        proc = proc[:proc.rfind("_AOD")]

        bashfilename = "doHadd_"+proc+"_"+oname+".sh"
        bashfile = open(bashfilename,"w")
        #write cmds to bash script
        if(args.big):
            for i in range(10):
                oname = d.name+"_"+proc+"_"+str(i)+".root"
                oname = "condor_"+oname
                oname = d.path+"/"+oname
                cmd = ""
                #check if file exists
                if os.path.exists(oname):
                    if(args.force):
                        cmd = cmdHadd+" -f"
                    else:
                        print(oname+" exists ")
                        continue
                else:
                    cmd = cmdHadd
                #print(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root")
                bashfile.write(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root\n")
                #if not args.dryRun:
                #    os.system(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root")	
                #print("Wrote to "+oname)
            oname = d.name+"_"+proc+".root"
            oname = "condor_"+oname
            oname = d.path+"/"+oname
            #check if file exists
            if os.path.exists(oname):
            	if(args.force):
            		cmd = cmdHadd+" -f"
            	else:
            		print(oname+" exists ")
            #print(cmd+" "+oname+" "+d.path+"/condor_*.root")	
            bashfile.write(cmd+" "+oname+" "+d.path+"/condor_*.root\n")	
            #if not args.dryRun:
            #    os.system(cmd+" "+oname+" "+d.path+"/condor_*.root")	
            #print("Wrote to "+oname)
            bashfile.close()
        else:
            oname += "_"+proc+".root"
            oname = "condor_"+oname
            oname = d.path+"/"+oname
            cmd = cmdHadd
            #check if file exists
            if os.path.exists(oname):
            	if(args.force):
            		cmd += " -f"
            	else:
            		print(oname+" exists ")
            		continue
            #print(cmd+" "+oname+" "+d.path+"/out/*.root")	
            bashfile.write(cmd+" "+oname+" "+d.path+"/out/*.root\n")	
            #if not args.dryRun:
            #    os.system(cmd+" "+oname+" "+d.path+"/out/*.root")
            #print("Wrote to "+oname)
            bashfile.close()
        print("To hadd run:")
        print("source",bashfilename)

if __name__ == "__main__":
	main()

