#!/bin/bash

# Script to query target DNA sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

# References:
# Heavy reliance on:
        # mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq
        # SeqKit: https://bioinf.shenwei.me/seqkit/

# REQUIRES that targetDB has already been indexed
# If it has not been index then run the following script in the directory of your choice: create_mmseqs_viral_nt_taxaDB.sh (found in /accessory)
# Note: mmseqs2 taxonomy is currently most useful if you have UniProt formatted fasta databases
# more details about database building can be found at: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy

# Set targetDB
# This nt DB is all refseq viral gene sequences and nearest neighbors (exact duplicates removed)
DB=/mnt/data1/databases/refseq_viral_genomes/refseq_virus_nn_nt/masked/refseq_virus_nt_UniVec_masked/nt.fnaDB

# Create output directory
mkdir -p ./results/mmseqs_nt_out;
OUT=./results/mmseqs_nt_out;

# Create Query databases
mmseqs createdb ./results/mmseqs_aa_out/pviral_aa_unclassified_seqs.fasta $OUT/seqtable_queryDB --dbtype 2;

# mmseqs search
mmseqs search $OUT/seqtable_queryDB $DB $OUT/resultDB $OUT/tmp_nt -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95;

# extract top-hit
mmseqs filterdb $OUT/resultDB $OUT/resultDB.firsthit --extract-lines 1;
mmseqs convertalis $OUT/seqtable_queryDB $DB $OUT/resultDB.firsthit $OUT/resultDB.firsthit.m8;

# Annotate
mmseqs_pviral_nt_annotate.R
