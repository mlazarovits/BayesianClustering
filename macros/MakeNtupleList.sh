if [ -z "$1" ]
then
        echo "Please pass [process] [selection] [HTbin] to make corresponding file list"
        return
fi
PD=$1
SEL=""
NAME=""
if [ $2 = "AL1IsoPho" ]
then
	SEL=R17_$2_v24
fi
if [ $PD = "GJETS" ]
then
	NAME=GJets_HT-$3_TuneCP5_13TeV-madgraphMLM-pythia8
elif [ $PD = "MET" ]
then
	NAME=$PD
elif [ $PD = "DEG" ]
then
	NAME=DoubleEG
elif [ $PD = "SMS_GlGl" ]
then
	NAME=CRAB_UserFiles
	SEL=v23
else
	NAME=${PD}_HT-$3_TuneCP5_13TeV-madgraphMLM-pythia8
fi
PREFIX=root://cmseos.fnal.gov/

FILE1=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_${SEL}/${NAME}/)
IFS="/" read -ra components1 <<< "$FILE1"
DIR1="${components1[-1]}"

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
sed -i 's/^/root:\/\/cmseos.fnal.gov\//' ${outfile}
done
echo Wrote file list to $outfile

