#!/bin/bash
# Tabuluates primer detection statistics
IN=../STATS/step_01
mkdir -p ../STATS/primer_stats
OUT=../STATS/primer_stats

for i in $IN/*.tsv; do

	NAME=`basename $i .tsv`;

	grep -v "#" $IN/"$NAME".tsv | cut -f1,2 | sed "1iprimerB\t$NAME" > $OUT/"$NAME".primer_stats_num_reads.out;

	grep -v "#" $IN/"$NAME".tsv | cut -f1,3 | sed "1iprimerB\t$NAME" | sed 's/%//' > $OUT/"$NAME".primer_stats_percent.out;

done

# Join all output tables
./primer_stats.R

# Clean up
rm $OUT/*.out;
