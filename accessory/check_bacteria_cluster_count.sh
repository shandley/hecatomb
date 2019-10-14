#!/bin/bash
# Script to deduplicate, count and assign taxonomy to sequences assigned as BACTERIAL
# During Step 9 (bbsplit) of the contaminant_removal.sh script
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/
	# mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq
        # SeqKit: https://bioinf.shenwei.me/seqkit/

# Prep output directory
mkdir -p ./results/bacterial_check
IN=./QC/step_9
OUT=./results/bacterial_check

# Set file names
for i in $IN/*_bacterial.fastq; do
        F=`basename $i _bacterial.fastq`;

# Summary:
	# Step 1: Remove exact duplicates (consider to be PCR artifacts)
	# Step 2: Deduplicate
	# Step 3: Reformat and create sequence table

# Remove exact duplicates
dedupe.sh in=$IN/"$F"_bacterial.fastq ow=t \
	out=$OUT/"$F"_bacterial.deduped.out.fastq \
	ac=f \
	-Xmx64g;

# Deduplication 
dedupe.sh in=$OUT/"$F"_bacterial.deduped.out.fastq ow=t \
	s=4 \
	rnc=t pbr=t \
	csf=$OUT/"$F"_stats.txt out=$OUT/"$F"_best.fasta \
	-Xmx64g;

# Convert to fasta
reformat.sh in=$OUT/"$F"_best.fasta out=$OUT/"$F"_reformated.fasta \
	deleteinput=t fastawrap=0 \
	ow=t \
	-Xmx64g;

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
