import os
import argparse
import ROOT

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--dir","-d",help="top directory that has directories of root files to hadd i.e. Output/GMSB/GMSB_AOD_v13_GMSB_L-150TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix/jets/",required=True)
	args = parser.parse_args()

	cmd = "hadd -d /uscmst1b_scratch/lpc1/3DayLifetime/mlazarov/ -j 4"
	for d in os.scandir(args.dir):
		oname = d.name
		if "GMSB" in d.path:
			oname += "_GMSB"
		elif "JetHT" in d.path:
			oname += "_JetHT"
		if "jets" in d.path:
			oname += "_jets"
		elif "photons" in d.path:
			oname += "_photons"
		oname += ".root"
		#os.system(cmd+" "+d.path+"/"+oname+" "+d.path+"/out/*.root")
		print(cmd+" "+d.path+"/"+oname+" "+d.path+"/out/*.root")	



if __name__ == "__main__":
	main()

