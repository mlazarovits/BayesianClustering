if [ -z "$1" ]
then
        echo "Please pass [process] [selection] [year] [HTbin] to make corresponding file list"
        return
fi
PD=$1
SEL=""
NAME=""
YR=""
if [ -z "$3" ]
then
	echo "default year" 
	YR=R17
else
	YR=$3
fi
if [ ! -z "$2" ]
then
	SEL=${YR}_$2_v24
fi
if [ $PD = "GJETS" ]
then
	NAME=GJets_HT-$4_TuneCP5_13TeV-madgraphMLM-pythia8
elif [ $PD = "QCD" ]
then
	if [ $4 = "50to100" ]
	then
		NAME=${PD}_HT$4_TuneCP5_13TeV-madgraphMLM-pythia8
	else
		NAME=${PD}_HT$4_TuneCP5_13TeV-madgraph-pythia8
	fi
elif [ $PD = "MET" ]
then
	NAME=$PD
elif [ $PD = "JetHT" ]
then
	NAME=$PD
elif [ $PD = "DEG" ]
then
	NAME=DoubleEG
elif [ $PD = "SMS-GlGl" ]
then
	NAME=CRAB_UserFiles
	SEL=$2
elif [ $PD = "gogoG_Sig" ]
then
	NAME=CRAB_UserFiles
	SEL=$2
else
	NAME=${PD}_HT-$4_TuneCP5_13TeV-madgraph-pythia8
fi
PREFIX=root://cmseos.fnal.gov/

#if multiple runs in sample, do all runs
FILE1=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_${SEL}/${NAME}/)
components1=($FILE1)
for dir in "${components1[@]}"
do
	DIR1=$dir
	
	FILE2=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_${SEL}/${NAME}/${DIR1})
	IFS="/" read -ra components2 <<< "$FILE2"
	DIR2="${components2[-1]}"
	
	FILE3=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_${SEL}/${NAME}/${DIR1}/${DIR2})
	
	outfile=filelists/${DIR1}_list.txt
	if [ -f ${outfile} ]; then
		echo "Recreating file" $outfile
		rm $outfile #start with a clean slate
	else
		echo "Creating file" $outfile
	fi
	
	for index in ${FILE3[@]}; do
		#eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_R17_${SEL}_v24/${NAME}/${DIR1}/${DIR2}/${index}/ >> ${outfile}
		xrdfs ${PREFIX} ls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_${SEL}/${NAME}/${DIR1}/${DIR2}/${index}/ >> ${outfile}
	done
	sed -i 's/^/root:\/\/cmseos.fnal.gov\//' ${outfile}
	echo Wrote file list to $outfile
done

