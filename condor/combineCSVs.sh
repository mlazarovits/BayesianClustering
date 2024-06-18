echo "Combining .csv files in" $1
part1=$(dirname "$1")
part2=$(basename "$part1")
part3=$(basename "$1")

outname=$1$part2"_"$part3
head -n 1 $1/out/*.0.csv > $outname.csv
tail -n+2 -q $1/out/*.csv >> $outname.csv

echo "Wrote to" $outname".csv"
