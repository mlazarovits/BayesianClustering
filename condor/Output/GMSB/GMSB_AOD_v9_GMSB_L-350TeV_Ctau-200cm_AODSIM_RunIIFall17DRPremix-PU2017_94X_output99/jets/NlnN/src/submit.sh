universe = vanilla
executable = execute_script.sh
output = ./condor/Output/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/log/job.$(Process).out
error = ./condor/Output/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/log/job.$(Process).err
log = ./condor/Output/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/log/job.log
transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN.$(Process).root
transfer_output_remaps = "Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN.$(Process).root=/uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/Output/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN/out/Output/GMSB/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99/jets/NlnN.$(Process).root"



###### job0######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 0 --evtLast 97
Queue



###### job1######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 98 --evtLast 194
Queue



###### job2######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 195 --evtLast 291
Queue



###### job3######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 292 --evtLast 388
Queue



###### job4######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 389 --evtLast 485
Queue



###### job5######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 486 --evtLast 582
Queue



###### job6######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 583 --evtLast 679
Queue



###### job7######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 680 --evtLast 776
Queue



###### job8######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 777 --evtLast 873
Queue



###### job9######
Arguments = -i /uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condor/../rootfiles/GMSB_AOD_v9_GMSB_L-350TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix-PU2017_94X_output99.root  --skim --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 -s 0 --evtFirst 874 --evtLast 970
Queue
