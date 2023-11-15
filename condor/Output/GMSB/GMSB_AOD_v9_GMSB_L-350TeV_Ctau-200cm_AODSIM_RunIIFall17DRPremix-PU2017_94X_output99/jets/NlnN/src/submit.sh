universe = vanilla
executable = execute_script.sh
output = ./Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/log/job.$(Process).out
error = ./Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/log/job.$(Process).err
log = ./Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/log/job.log
transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2, config.tgz, /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skim.$(Process).root
transfer_output_remaps = "Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skim.$(Process).root=/uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skim.$(Process).root"



###### job0######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 0 --evtLast 9 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job1######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 10 --evtLast 19 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job2######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 20 --evtLast 29 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job3######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 30 --evtLast 39 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job4######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 40 --evtLast 49 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job5######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 50 --evtLast 59 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job6######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 60 --evtLast 69 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job7######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 70 --evtLast 79 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job8######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 80 --evtLast 89 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job9######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 90 --evtLast 99 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job10######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 100 --evtLast 109 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job11######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 110 --evtLast 119 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job12######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 120 --evtLast 129 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job13######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 130 --evtLast 139 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job14######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 140 --evtLast 149 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job15######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 150 --evtLast 159 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job16######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 160 --evtLast 169 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job17######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 170 --evtLast 179 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job18######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 180 --evtLast 189 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job19######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 190 --evtLast 199 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job20######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 200 --evtLast 209 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job21######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 210 --evtLast 219 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job22######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 220 --evtLast 229 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job23######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 230 --evtLast 239 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job24######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 240 --evtLast 249 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job25######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 250 --evtLast 259 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job26######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 260 --evtLast 269 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job27######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 270 --evtLast 279 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job28######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 280 --evtLast 289 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job29######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 290 --evtLast 299 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job30######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 300 --evtLast 309 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job31######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 310 --evtLast 319 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job32######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 320 --evtLast 329 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job33######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 330 --evtLast 339 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job34######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 340 --evtLast 349 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job35######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 350 --evtLast 359 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job36######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 360 --evtLast 369 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job37######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 370 --evtLast 379 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job38######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 380 --evtLast 389 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job39######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 390 --evtLast 399 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job40######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 400 --evtLast 409 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job41######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 410 --evtLast 419 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job42######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 420 --evtLast 429 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job43######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 430 --evtLast 439 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job44######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 440 --evtLast 449 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job45######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 450 --evtLast 459 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job46######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 460 --evtLast 469 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job47######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 470 --evtLast 479 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job48######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 480 --evtLast 489 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job49######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 490 --evtLast 499 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job50######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 500 --evtLast 509 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job51######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 510 --evtLast 519 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job52######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 520 --evtLast 529 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job53######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 530 --evtLast 539 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job54######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 540 --evtLast 549 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job55######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 550 --evtLast 559 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job56######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 560 --evtLast 569 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job57######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 570 --evtLast 579 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job58######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 580 --evtLast 589 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job59######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 590 --evtLast 599 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job60######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 600 --evtLast 609 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job61######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 610 --evtLast 619 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job62######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 620 --evtLast 629 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job63######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 630 --evtLast 639 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job64######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 640 --evtLast 649 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job65######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 650 --evtLast 659 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job66######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 660 --evtLast 669 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job67######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 670 --evtLast 679 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job68######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 680 --evtLast 689 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job69######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 690 --evtLast 699 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job70######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 700 --evtLast 709 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job71######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 710 --evtLast 718 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job72######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 719 --evtLast 727 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job73######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 728 --evtLast 736 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job74######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 737 --evtLast 745 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job75######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 746 --evtLast 754 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job76######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 755 --evtLast 763 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job77######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 764 --evtLast 772 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job78######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 773 --evtLast 781 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job79######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 782 --evtLast 790 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job80######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 791 --evtLast 799 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job81######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 800 --evtLast 808 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job82######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 809 --evtLast 817 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job83######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 818 --evtLast 826 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job84######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 827 --evtLast 835 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job85######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 836 --evtLast 844 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job86######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 845 --evtLast 853 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job87######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 854 --evtLast 862 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job88######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 863 --evtLast 871 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job89######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 872 --evtLast 880 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job90######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 881 --evtLast 889 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job91######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 890 --evtLast 898 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job92######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 899 --evtLast 907 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job93######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 908 --evtLast 916 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job94######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 917 --evtLast 925 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job95######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 926 --evtLast 934 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job96######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 935 --evtLast 943 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job97######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 944 --evtLast 952 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job98######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 953 --evtLast 961 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue



###### job99######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 962 --evtLast 970 -o Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/skims
Queue
