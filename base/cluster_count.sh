#!/bin/bash
# Script to dereplicate and count sequences
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Summary:
	# Step 1: Remove exact duplicates
        # Step 2: Deduplicate
        # Step 3: Reformat and prepare for sequence table generation

# Prep output directory
mkdir -p ./QC/step_8/clustered;
IN=./QC/step_8
OUT=./QC/step_8/clustered;

# Set file names
for i in $IN/*_viral_amb.fastq; do
        F=`basename $i _viral_amb.fastq`;

# Step 1: Remove exact duplicates
dedupe.sh in=$IN/"$F"_viral_amb.fastq ow=t \
	out=$OUT/"$F"_R1.s8.deduped.out.fastq \
	ac=f \
	-Xmx128g;

# Step 2: Dereplicate
dedupe.sh in=$OUT/"$F"_R1.s8.deduped.out.fastq ow=t \
	s=4 \
	rnc=t pbr=t \
	csf=$OUT/"$F"_stats.txt out=$OUT/"$F"_best.fasta \
	-Xmx128g;

# Step 3: Extract sequences and counts for seqtable (count table)
# Convert to fasta
reformat.sh in=$OUT/"$F"_best.fasta out=$OUT/"$F"_reformated.fasta \
	deleteinput=t fastawrap=0 \
	ow=t \
	-Xmx128g;

# Parse and combine stats and contig files
# Extract sequences
grep -v '>' $OUT/"$F"_reformated.fasta | sed '1i sequence' > $OUT/"$F"_seqs.txt;

# Extract sequence IDs
# grep '>' $OUT/"$F"_reformated.fasta | sed 's|>Cluster_||' | awk -F "," '{ print$1 }' | sort -n | sed '1i contig_ids' > $OUT/"$F"_contig_ids.txt;

# Extract counts
cut -f 2 $OUT/"$F"_stats.txt | sed "1s/size/$F/" > $OUT/"$F"_counts.txt;

# Create sequence table
paste $OUT/"$F"_seqs.txt $OUT/"$F"_counts.txt > $OUT/"$F"_seqtable.txt;

done

