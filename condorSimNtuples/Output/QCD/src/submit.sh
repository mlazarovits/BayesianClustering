universe = vanilla
executable = execute_script.sh
output = ./Output/QCD/log/job.$(Process).out
error = ./Output/QCD/log/job.$(Process).err
log = ./Output/QCD/log/job.log
transfer_input_files = /uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2, configSim.tgz, 
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = condorSimNtuples_QCD.$(Process).root
transfer_output_remaps = "condorSimNtuples_QCD.$(Process).root=/uscms/home/mlazarov/nobackup/CMSSW_10_6_5/src/BayesianClustering/condorSimNtuples/Output/QCD/out/condorSimNtuples_QCD.$(Process).root"



###### job0######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 0 --evtLast 10 -o condorSimNtuples_QCD.$(Process)
Queue



###### job1######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 10 --evtLast 20 -o condorSimNtuples_QCD.$(Process)
Queue



###### job2######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 20 --evtLast 30 -o condorSimNtuples_QCD.$(Process)
Queue



###### job3######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 30 --evtLast 40 -o condorSimNtuples_QCD.$(Process)
Queue



###### job4######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 40 --evtLast 50 -o condorSimNtuples_QCD.$(Process)
Queue



###### job5######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 50 --evtLast 60 -o condorSimNtuples_QCD.$(Process)
Queue



###### job6######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 60 --evtLast 70 -o condorSimNtuples_QCD.$(Process)
Queue



###### job7######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 70 --evtLast 80 -o condorSimNtuples_QCD.$(Process)
Queue



###### job8######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 80 --evtLast 90 -o condorSimNtuples_QCD.$(Process)
Queue



###### job9######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 90 --evtLast 100 -o condorSimNtuples_QCD.$(Process)
Queue



###### job10######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 100 --evtLast 110 -o condorSimNtuples_QCD.$(Process)
Queue



###### job11######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 110 --evtLast 120 -o condorSimNtuples_QCD.$(Process)
Queue



###### job12######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 120 --evtLast 130 -o condorSimNtuples_QCD.$(Process)
Queue



###### job13######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 130 --evtLast 140 -o condorSimNtuples_QCD.$(Process)
Queue



###### job14######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 140 --evtLast 150 -o condorSimNtuples_QCD.$(Process)
Queue



###### job15######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 150 --evtLast 160 -o condorSimNtuples_QCD.$(Process)
Queue



###### job16######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 160 --evtLast 170 -o condorSimNtuples_QCD.$(Process)
Queue



###### job17######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 170 --evtLast 180 -o condorSimNtuples_QCD.$(Process)
Queue



###### job18######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 180 --evtLast 190 -o condorSimNtuples_QCD.$(Process)
Queue



###### job19######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 190 --evtLast 200 -o condorSimNtuples_QCD.$(Process)
Queue



###### job20######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 200 --evtLast 210 -o condorSimNtuples_QCD.$(Process)
Queue



###### job21######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 210 --evtLast 220 -o condorSimNtuples_QCD.$(Process)
Queue



###### job22######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 220 --evtLast 230 -o condorSimNtuples_QCD.$(Process)
Queue



###### job23######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 230 --evtLast 240 -o condorSimNtuples_QCD.$(Process)
Queue



###### job24######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 240 --evtLast 250 -o condorSimNtuples_QCD.$(Process)
Queue



###### job25######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 250 --evtLast 260 -o condorSimNtuples_QCD.$(Process)
Queue



###### job26######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 260 --evtLast 270 -o condorSimNtuples_QCD.$(Process)
Queue



###### job27######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 270 --evtLast 280 -o condorSimNtuples_QCD.$(Process)
Queue



###### job28######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 280 --evtLast 290 -o condorSimNtuples_QCD.$(Process)
Queue



###### job29######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 290 --evtLast 300 -o condorSimNtuples_QCD.$(Process)
Queue



###### job30######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 300 --evtLast 310 -o condorSimNtuples_QCD.$(Process)
Queue



###### job31######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 310 --evtLast 320 -o condorSimNtuples_QCD.$(Process)
Queue



###### job32######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 320 --evtLast 330 -o condorSimNtuples_QCD.$(Process)
Queue



###### job33######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 330 --evtLast 340 -o condorSimNtuples_QCD.$(Process)
Queue



###### job34######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 340 --evtLast 350 -o condorSimNtuples_QCD.$(Process)
Queue



###### job35######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 350 --evtLast 360 -o condorSimNtuples_QCD.$(Process)
Queue



###### job36######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 360 --evtLast 370 -o condorSimNtuples_QCD.$(Process)
Queue



###### job37######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 370 --evtLast 380 -o condorSimNtuples_QCD.$(Process)
Queue



###### job38######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 380 --evtLast 390 -o condorSimNtuples_QCD.$(Process)
Queue



###### job39######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 390 --evtLast 400 -o condorSimNtuples_QCD.$(Process)
Queue



###### job40######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 400 --evtLast 410 -o condorSimNtuples_QCD.$(Process)
Queue



###### job41######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 410 --evtLast 420 -o condorSimNtuples_QCD.$(Process)
Queue



###### job42######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 420 --evtLast 430 -o condorSimNtuples_QCD.$(Process)
Queue



###### job43######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 430 --evtLast 440 -o condorSimNtuples_QCD.$(Process)
Queue



###### job44######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 440 --evtLast 450 -o condorSimNtuples_QCD.$(Process)
Queue



###### job45######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 450 --evtLast 460 -o condorSimNtuples_QCD.$(Process)
Queue



###### job46######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 460 --evtLast 470 -o condorSimNtuples_QCD.$(Process)
Queue



###### job47######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 470 --evtLast 480 -o condorSimNtuples_QCD.$(Process)
Queue



###### job48######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 480 --evtLast 490 -o condorSimNtuples_QCD.$(Process)
Queue



###### job49######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 490 --evtLast 500 -o condorSimNtuples_QCD.$(Process)
Queue



###### job50######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 500 --evtLast 510 -o condorSimNtuples_QCD.$(Process)
Queue



###### job51######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 510 --evtLast 520 -o condorSimNtuples_QCD.$(Process)
Queue



###### job52######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 520 --evtLast 530 -o condorSimNtuples_QCD.$(Process)
Queue



###### job53######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 530 --evtLast 540 -o condorSimNtuples_QCD.$(Process)
Queue



###### job54######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 540 --evtLast 550 -o condorSimNtuples_QCD.$(Process)
Queue



###### job55######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 550 --evtLast 560 -o condorSimNtuples_QCD.$(Process)
Queue



###### job56######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 560 --evtLast 570 -o condorSimNtuples_QCD.$(Process)
Queue



###### job57######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 570 --evtLast 580 -o condorSimNtuples_QCD.$(Process)
Queue



###### job58######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 580 --evtLast 590 -o condorSimNtuples_QCD.$(Process)
Queue



###### job59######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 590 --evtLast 600 -o condorSimNtuples_QCD.$(Process)
Queue



###### job60######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 600 --evtLast 610 -o condorSimNtuples_QCD.$(Process)
Queue



###### job61######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 610 --evtLast 620 -o condorSimNtuples_QCD.$(Process)
Queue



###### job62######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 620 --evtLast 630 -o condorSimNtuples_QCD.$(Process)
Queue



###### job63######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 630 --evtLast 640 -o condorSimNtuples_QCD.$(Process)
Queue



###### job64######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 640 --evtLast 650 -o condorSimNtuples_QCD.$(Process)
Queue



###### job65######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 650 --evtLast 660 -o condorSimNtuples_QCD.$(Process)
Queue



###### job66######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 660 --evtLast 670 -o condorSimNtuples_QCD.$(Process)
Queue



###### job67######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 670 --evtLast 680 -o condorSimNtuples_QCD.$(Process)
Queue



###### job68######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 680 --evtLast 690 -o condorSimNtuples_QCD.$(Process)
Queue



###### job69######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 690 --evtLast 700 -o condorSimNtuples_QCD.$(Process)
Queue



###### job70######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 700 --evtLast 710 -o condorSimNtuples_QCD.$(Process)
Queue



###### job71######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 710 --evtLast 720 -o condorSimNtuples_QCD.$(Process)
Queue



###### job72######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 720 --evtLast 730 -o condorSimNtuples_QCD.$(Process)
Queue



###### job73######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 730 --evtLast 740 -o condorSimNtuples_QCD.$(Process)
Queue



###### job74######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 740 --evtLast 750 -o condorSimNtuples_QCD.$(Process)
Queue



###### job75######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 750 --evtLast 760 -o condorSimNtuples_QCD.$(Process)
Queue



###### job76######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 760 --evtLast 770 -o condorSimNtuples_QCD.$(Process)
Queue



###### job77######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 770 --evtLast 780 -o condorSimNtuples_QCD.$(Process)
Queue



###### job78######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 780 --evtLast 790 -o condorSimNtuples_QCD.$(Process)
Queue



###### job79######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 790 --evtLast 800 -o condorSimNtuples_QCD.$(Process)
Queue



###### job80######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 800 --evtLast 810 -o condorSimNtuples_QCD.$(Process)
Queue



###### job81######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 810 --evtLast 820 -o condorSimNtuples_QCD.$(Process)
Queue



###### job82######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 820 --evtLast 830 -o condorSimNtuples_QCD.$(Process)
Queue



###### job83######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 830 --evtLast 840 -o condorSimNtuples_QCD.$(Process)
Queue



###### job84######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 840 --evtLast 850 -o condorSimNtuples_QCD.$(Process)
Queue



###### job85######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 850 --evtLast 860 -o condorSimNtuples_QCD.$(Process)
Queue



###### job86######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 860 --evtLast 870 -o condorSimNtuples_QCD.$(Process)
Queue



###### job87######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 870 --evtLast 880 -o condorSimNtuples_QCD.$(Process)
Queue



###### job88######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 880 --evtLast 890 -o condorSimNtuples_QCD.$(Process)
Queue



###### job89######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 890 --evtLast 900 -o condorSimNtuples_QCD.$(Process)
Queue



###### job90######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 900 --evtLast 910 -o condorSimNtuples_QCD.$(Process)
Queue



###### job91######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 910 --evtLast 920 -o condorSimNtuples_QCD.$(Process)
Queue



###### job92######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 920 --evtLast 930 -o condorSimNtuples_QCD.$(Process)
Queue



###### job93######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 930 --evtLast 940 -o condorSimNtuples_QCD.$(Process)
Queue



###### job94######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 940 --evtLast 950 -o condorSimNtuples_QCD.$(Process)
Queue



###### job95######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 950 --evtLast 960 -o condorSimNtuples_QCD.$(Process)
Queue



###### job96######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 960 --evtLast 970 -o condorSimNtuples_QCD.$(Process)
Queue



###### job97######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 970 --evtLast 980 -o condorSimNtuples_QCD.$(Process)
Queue



###### job98######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 980 --evtLast 990 -o condorSimNtuples_QCD.$(Process)
Queue



###### job99######
Arguments = -v 0 --nevts 1000 --spikeProb 0 --QCD --evtFirst 990 --evtLast 1000 -o condorSimNtuples_QCD.$(Process)
Queue
