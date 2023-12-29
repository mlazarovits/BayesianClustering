import os
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--dir","-d",help="top directory that has directories of root files to hadd i.e. Output/GMSB/GMSB_AOD_v13_GMSB_L-150TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix/jets/",required=True)
	parser.add_argument("--force",help='force remake of outfile',action='store_true')
	args = parser.parse_args()
	print("force? "+str(args.force))

	cmd = "hadd -d /uscmst1b_scratch/lpc1/3DayLifetime/mlazarov/ -j 4"
	for d in os.scandir(args.dir):
		oname = d.name
		proc = d.path[d.path.find("/")+1:]
		proc = proc[proc.find("/")+1:]
		proc = proc[:proc.find("/")]
		
		oname += "_"+proc+".root"
		#check if file exists
		if os.path.exists(oname):
			if(args.force):
				cmd += " -f"
			else:
				continue
		if not os.path.exists(d.path+"/out/"):
			continue
		print(cmd+" "+d.path+"/"+oname+" "+d.path+"/out/*.root")	
		os.system(cmd+" "+d.path+"/"+oname+" "+d.path+"/out/*.root")
		print("Wrote to "+d.path+"/"+oname)


if __name__ == "__main__":
	main()

