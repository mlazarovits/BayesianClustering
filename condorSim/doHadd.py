import os
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--dir","-d",help="top directory that has directories of root files to hadd i.e. Output/[dir]",required=True)
	parser.add_argument("--force",help='force remake of outfile',action='store_true')
	args = parser.parse_args()

	cmd = "hadd -d /uscmst1b_scratch/lpc1/3DayLifetime/mlazarov/ -j 4"
	for d in os.scandir(args.dir):
		if "/out" not in d.path:
			continue
		proc = d.path.split("/")[1]
		opath = d.path[:d.path.find(d.name)]
		oname = "condorSim_"+proc+".root"
		oname = opath+oname
		#check if file exists
		if os.path.exists(oname):
			if(args.force):
				cmd += " -f"
			else:
				print(oname+" exists ")
				continue
		print(cmd+" "+oname+" "+d.path+"/*.root")	
		os.system(cmd+" "+oname+" "+d.path+"/*.root")
		print("Wrote to "+oname)
if __name__ == "__main__":
	main()

