#!/bin/bash

# cluster_count_htcf_modified.sbatch

#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-49
#SBATCH -J qv_cnt_rm
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=ALL

set -ue

module load bbtools

# Script to dereplicate and count sequences
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Summary:
	# Step 1: Remove exact duplicates (consider to be PCR artifacts)
        # Step 2: Deduplicate
        # Step 3: Reformat and prepare for sequence table generation

# Prep output directory
mkdir -p ./QC/step_9/clustered
IN=./QC/step_9
OUT=./QC/step_9/clustered

read1=$(sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt)
sample=$(echo ${read1} | sed s/_R1_001.fastq.gz//g)

# Set file names
#for i in $IN/*_viral_amb.fastq; do
        #F=`basename $i _viral_amb.fastq`;

# Step 1: Remove exact duplicates
echo
echo "dedupe step 1"
echo

dedupe.sh \
	in=$IN/${sample}_viral_amb.fastq \
	ow=t \
	out=$OUT/${sample}_R1.s9.deduped.out.fastq \
	ac=f \
	-Xmx64g \
	t=${SLURM_CPUS_PER_TASK}

# Step 2: Dereplicate
echo
echo "dedupe step 2"
echo

dedupe.sh \
	in=$OUT/${sample}_R1.s9.deduped.out.fastq \
	ow=t \
	s=4 \
	rnc=t \
	pbr=t \
	csf=$OUT/${sample}_stats.txt \
	out=$OUT/${sample}_best.fasta \
	-Xmx64g \
	t=${SLURM_CPUS_PER_TASK}

# Step 3: Extract sequences and counts for seqtable (count table)
# Convert to fasta
echo
echo "dedupe step 3"
echo

reformat.sh \
	in=$OUT/${sample}_best.fasta \
	out=$OUT/${sample}_reformated.fasta \
	deleteinput=t \
	fastawrap=0 \
	ow=t \
	-Xmx64g

# Parse and combine stats and contig files
# Extract sequences
echo "extract sequences"

grep -v '>' $OUT/${sample}_reformated.fasta | sed '1i sequence' > $OUT/${sample}_seqs.txt

# Extract sequence IDs
# grep '>' $OUT/"$F"_reformated.fasta | sed 's|>Cluster_||' | awk -F "," '{ print$1 }' | sort -n | sed '1i contig_ids' > $OUT/"$F"_contig_ids.txt;

# Extract counts
echo
echo "extract counts"

cut -f 2 $OUT/${sample}_stats.txt | sed "1s/size/${sample}/" > $OUT/${sample}_counts.txt

# Create sequence table
echo
echo "create sequence table"

paste $OUT/${sample}_seqs.txt $OUT/${sample}_counts.txt > $OUT/${sample}_seqtable.txt

