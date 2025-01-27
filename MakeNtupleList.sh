PD=GJETS
SEL=AL1IsoPho
NAME=GJets_HT-$1_TuneCP5_13TeV-madgraphMLM-pythia8
outfile=filelists/${NAME}_list.txt
#PD=$1 for first cmd line arg
rm $outfile #start with a clean slate

PREFIX=root://cmseos.fnal.gov/

FILE1=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_R17_${SEL}_v24/${NAME}/)
IFS="/" read -ra components1 <<< "$FILE1"
DIR1="${components1[-1]}"

FILE2=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_R17_${SEL}_v24/${NAME}/${DIR1})
IFS="/" read -ra components2 <<< "$FILE2"
DIR2="${components2[-1]}"

FILE3=$(eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_R17_${SEL}_v24/${NAME}/${DIR1}/${DIR2})
for index in ${FILE3[@]}; do
	#eosls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_R17_${SEL}_v24/${NAME}/${DIR1}/${DIR2}/${index}/ >> ${outfile}
	xrdfs ${PREFIX} ls /store/user/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_${PD}_R17_${SEL}_v24/${NAME}/${DIR1}/${DIR2}/${index}/ >> ${outfile}
sed -i 's/^/root:\/\/cmseos.fnal.gov\//' ${outfile}
done
echo Wrote file list to $outfile

