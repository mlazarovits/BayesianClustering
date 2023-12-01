universe = vanilla
executable = execute_script.sh
output = ./Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/GMMonly/log/job.$(Process).out
error = ./Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/GMMonly/log/job.$(Process).err
log = ./Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/GMMonly/log/job.log
transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2, config.tgz, /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = condorSkim.$(Process).root
transfer_output_remaps = "condorSkim.$(Process).root=/uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/GMMonly/out/condorSkim.$(Process).root"



###### job0######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 0 --evtLast 19 -o condorSkim.$(Process)
Queue



###### job1######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 19 --evtLast 39 -o condorSkim.$(Process)
Queue



###### job2######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 39 --evtLast 59 -o condorSkim.$(Process)
Queue



###### job3######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 59 --evtLast 79 -o condorSkim.$(Process)
Queue



###### job4######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 79 --evtLast 99 -o condorSkim.$(Process)
Queue



###### job5######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 99 --evtLast 119 -o condorSkim.$(Process)
Queue



###### job6######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 119 --evtLast 139 -o condorSkim.$(Process)
Queue



###### job7######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 139 --evtLast 159 -o condorSkim.$(Process)
Queue



###### job8######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 159 --evtLast 179 -o condorSkim.$(Process)
Queue



###### job9######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 179 --evtLast 199 -o condorSkim.$(Process)
Queue



###### job10######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 199 --evtLast 219 -o condorSkim.$(Process)
Queue



###### job11######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 219 --evtLast 239 -o condorSkim.$(Process)
Queue



###### job12######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 239 --evtLast 259 -o condorSkim.$(Process)
Queue



###### job13######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 259 --evtLast 279 -o condorSkim.$(Process)
Queue



###### job14######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 279 --evtLast 299 -o condorSkim.$(Process)
Queue



###### job15######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 299 --evtLast 319 -o condorSkim.$(Process)
Queue



###### job16######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 319 --evtLast 339 -o condorSkim.$(Process)
Queue



###### job17######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 339 --evtLast 359 -o condorSkim.$(Process)
Queue



###### job18######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 359 --evtLast 379 -o condorSkim.$(Process)
Queue



###### job19######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 379 --evtLast 399 -o condorSkim.$(Process)
Queue



###### job20######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 399 --evtLast 419 -o condorSkim.$(Process)
Queue



###### job21######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 419 --evtLast 438 -o condorSkim.$(Process)
Queue



###### job22######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 438 --evtLast 457 -o condorSkim.$(Process)
Queue



###### job23######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 457 --evtLast 476 -o condorSkim.$(Process)
Queue



###### job24######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 476 --evtLast 495 -o condorSkim.$(Process)
Queue



###### job25######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 495 --evtLast 514 -o condorSkim.$(Process)
Queue



###### job26######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 514 --evtLast 533 -o condorSkim.$(Process)
Queue



###### job27######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 533 --evtLast 552 -o condorSkim.$(Process)
Queue



###### job28######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 552 --evtLast 571 -o condorSkim.$(Process)
Queue



###### job29######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 571 --evtLast 590 -o condorSkim.$(Process)
Queue



###### job30######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 590 --evtLast 609 -o condorSkim.$(Process)
Queue



###### job31######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 609 --evtLast 628 -o condorSkim.$(Process)
Queue



###### job32######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 628 --evtLast 647 -o condorSkim.$(Process)
Queue



###### job33######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 647 --evtLast 666 -o condorSkim.$(Process)
Queue



###### job34######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 666 --evtLast 685 -o condorSkim.$(Process)
Queue



###### job35######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 685 --evtLast 704 -o condorSkim.$(Process)
Queue



###### job36######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 704 --evtLast 723 -o condorSkim.$(Process)
Queue



###### job37######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 723 --evtLast 742 -o condorSkim.$(Process)
Queue



###### job38######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 742 --evtLast 761 -o condorSkim.$(Process)
Queue



###### job39######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 761 --evtLast 780 -o condorSkim.$(Process)
Queue



###### job40######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 780 --evtLast 799 -o condorSkim.$(Process)
Queue



###### job41######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 799 --evtLast 818 -o condorSkim.$(Process)
Queue



###### job42######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 818 --evtLast 837 -o condorSkim.$(Process)
Queue



###### job43######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 837 --evtLast 856 -o condorSkim.$(Process)
Queue



###### job44######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 856 --evtLast 875 -o condorSkim.$(Process)
Queue



###### job45######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 875 --evtLast 894 -o condorSkim.$(Process)
Queue



###### job46######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 894 --evtLast 913 -o condorSkim.$(Process)
Queue



###### job47######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 913 --evtLast 932 -o condorSkim.$(Process)
Queue



###### job48######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 932 --evtLast 951 -o condorSkim.$(Process)
Queue



###### job49######
Arguments = -i GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 2 --evtFirst 951 --evtLast 970 -o condorSkim.$(Process)
Queue
