#!/bin/bash
# Summarize bacterial mappings from step 9 of contaminant_remaval.sh

# Set file names
for i in ./QC/step_9/*_scafstats.txt; do
        F=`basename $i _scafstats.txt`;

# Create table
	cut -f1 ./QC/step_9/"$F"_scafstats.txt | tail -n +2 | awk -F " " '{ print$2 }' > ./QC/step_9/"$F"_taxon.names;
	cut -f2 ./QC/step_9/"$F"_scafstats.txt | tail -n +2 > ./QC/step_9/"$F"_per_unambiguousReads;
	paste -d "\t" ./QC/step_9/"$F"_taxon.names ./QC/step_9/"$F"_per_unambiguousReads > ./QC/step_9/"$F"_bac_contaminants.txt;

# Collapse table as one entry per bacterial taxa, summing the counts
	awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]}}' ./QC/step_9/"$F"_bac_contaminants.txt > ./QC/step_9/"$F"_bac_contaminants_collapsed.txt

# Label columns
	sed -i 's/ /\t/g' ./QC/step_9/"$F"_bac_contaminants_collapsed.txt;
	sed -i "1ibacteria\t$F" ./QC/step_9/"$F"_bac_contaminants_collapsed.txt;
	
# Remove any 0 byte files
	find ./QC/step_9/ -type f -size 0 -exec rm -f {} \;

done

# Merge table
bacterial_contaminants.R

# Clean up
mkdir -p ./QC/stats/;
cp ./QC/step_9/bacterial_contaminant.table ./QC/stats;
rm ./QC/step_9/*_taxon.names;
rm ./QC/step_9/*_per_unambiguousReads;
rm ./QC/step_9/*_contaminants.txt;
rm ./QC/step_9/*_contaminants_collapsed.txt;
