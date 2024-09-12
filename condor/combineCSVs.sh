if [ -z "$1" ]
then
	echo "Please pass directory with out/*.csv files to concatenate"
	return
fi
for i in "$@"
do
	echo "Combining .csv files in $i"
	part1=$(dirname "$i")
	part2=$(basename "$part1")
	part3=$(basename "$i")
	
	outname=$i"/"$part2"_"$part3
	head -n 1 $i/out/*.0.csv > $outname.csv
	tail -n+2 -q $i/out/*.csv >> $outname.csv
	
	echo "Wrote to" $outname".csv"
	echo ""
done
