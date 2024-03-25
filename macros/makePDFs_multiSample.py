import subprocess
import argparse
from shlex import split

parser = argparse.ArgumentParser()
parser.add_argument("--files","-f",help="formatted ('_formatted.root') files to run over",required=True,nargs='+')
parser.add_argument("--observable","-obs",help="observable to plot (string match)", required=True)
parser.add_argument("--outdir","-odir",help="output directory of plots (per file)",required=True)
parser.add_argument("--log",action='store_true')
args = parser.parse_args()

for f in args.files:
	if f.find("formatted") == -1:
		print("Skipping file",f,"not formatted (run macros/HistFormat.C on this sample)")
		continue
	if args.log:
	        strlog = "true"
	else:
	        strlog = "false"
	cmd = 'root -l -b -q \'macros/MakePDFs.C("'+f+'","'+args.outdir+'","'+args.observable+'","'+strlog+'")\''
	print("Executing",cmd)
	subprocess.run(split(cmd))
