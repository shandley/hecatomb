#!/bin/bash

# fetch_viral_fastas_by_family_htcf_mod.sbatch

# This script will create fasta entries for every family in viruse_tax_table.tsv in the ./results/output directory

#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH -J viral_fastas
#SBATCH --mail-user=mihindu
#SBATCH --mail-type=ALL

# Immediately exit on errors
set -ue

# Modules
module load seqkit
module load pullseq

# Create directory to hold your fastas
mkdir -p ./results/family_reads
OUT=./results/family_reads

# Create list of families present in viruses_tax_table.sh
tail -n +2 ./results/viruses_tax_table.tsv | cut -f6 | sort | uniq > ./results/family.list;

# Create a temporary file list of each family to make for-do-looping simple
cd $OUT;
xargs touch <../family.list;
for f in *; do
	mv "$f" $(echo "$f".id);
done;
