universe = vanilla
executable = execute_script.sh
output = ./Output/QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100/log/job.$(Process).out
error = ./Output/QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100/log/job.$(Process).err
log = ./Output/QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100/log/job.log
transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2, configSim.tgz, 
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process).root
transfer_output_remaps = "condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process).root=/uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condorSim/Output/QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100/out/condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process).root"



###### job0######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 0 --evtLast 100 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job1######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 100 --evtLast 200 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job2######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 200 --evtLast 300 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job3######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 300 --evtLast 400 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job4######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 400 --evtLast 500 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job5######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 500 --evtLast 600 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job6######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 600 --evtLast 700 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job7######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 700 --evtLast 800 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job8######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 800 --evtLast 900 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue



###### job9######
Arguments = -i root://cmseos.fnal.gov//store/user/mlazarov/KUCMSNtuples/condorSimNtuples_QCD.root --alpha 0.1 --EMalpha 0.5 -v 0 -t 1.0 --gev 0.1 --minpt 10.0 --minNrhs 2 --minemE 0 --minRhE 0.5 -s 0 --evtFirst 900 --evtLast 1000 -o condorSim_QCD_defaultv1_bhcAlpha0p100_emAlpha0p500_thresh1p0_NperGeV0p100.$(Process)
Queue
