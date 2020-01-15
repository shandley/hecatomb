#!/bin/bash

# Prep output directories
mkdir -p ./assembly/megahit_contigs/annotation

# Set Variables
ANNOTATEDB=/mnt/data1/databases/CAT/CAT_prepare_20190108
IN=./assembly/megahit_contigs
OUT=./assembly/megahit_contigs/annotation
MIN=1000

# Set file names
for i in $IN/*.mh.contigs.fa; do
        F=`basename $i .mh.contigs.fa`;

# Filter out sequences shorter than length MIN

	seqkit seq $IN/"$F".mh.contigs.fa -m $MIN -j 24 -o $OUT/"$F"_n1000.contigs.fa;

# CAT Annotations

	CAT contigs -c $OUT/"$F"_n1000.contigs.fa -d $ANNOTATEDB/2019-05-06_CAT_database -t $ANNOTATEDB/2019-05-06_taxonomy \
	--out_prefix $OUT/$F --index_chunks 1;

	CAT add_names -i $OUT/"$F".contig2classification.txt -o $OUT/"$F".contig2classification.official_names.txt -t $ANNOTATEDB/2019-05-06_taxonomy --only_official;

	CAT summarise -c $OUT/"$F"_n1000.contigs.fa -i $OUT/"$F".contig2classification.official_names.txt -o $OUT/"$F"_CAT.summary.txt;

done
