#!/bin/bash
#PhiX and UniVec Summary Stats

# Set file names
for i in ./QC/step_5/*.stats; do
        F=`basename $i .stats`;

	NAME="$(head -n 1 ./QC/step_5/"$F".stats | awk -F "[/\t]" '{print$5}')";

	grep -v "#" ./QC/step_5/"$F".stats | cut -f1,3 | sed "1icontaminant\t$NAME" | sed 's/%//g'> ./QC/step_5/"$F"_edited;

	find ./QC/step_5/ -type f -size 0 -exec rm -f {} \;
done

# Merge table
contaminant_stats.R

# Clean up
mkdir -p ./QC/stats/;
cp ./QC/step_5/contaminant.table ./QC/stats;
rm ./QC/step_5/*_edited;
