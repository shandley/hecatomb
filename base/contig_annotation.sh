#!/bin/bash

# Set variables
DIR=./assembly/contig_dictionary

# CAT annotation
# GitHub: https://github.com/dutilh/CAT

# Main CAT pipeline
CAT contigs -c $DIR/contig_dictionary.fasta -o $DIR/out.CAT \
        --force \
        --sensitive \
        -d /mnt/data1/databases/CAT/CAT_prepare_20200304/2020-03-04_CAT_database \
        -t /mnt/data1/databases/CAT/CAT_prepare_20200304/2020-03-04_taxonomy;

# Add taxonomic names
CAT add_names -i $DIR/out.CAT.contig2classification.txt -o $DIR/out.CAT.taxonomy \
        -t /mnt/data1/databases/CAT/CAT_prepare_20200304/2020-03-04_taxonomy \
        --only_official --exclude_scores --force;

# Summarize results
CAT summarise -c $DIR/contig_dictionary.fasta -i $DIR/out.CAT.taxonomy \
        -o $DIR/out.CAT.summary;

# Creat contig taxonomy table
# Seperate names from CAT scores and lineage information
        cut -f1,6-12 $DIR/out.CAT.taxonomy > $DIR/contig.taxonomy;
        cut -f1-5 $DIR/out.CAT.taxonomy > $DIR/CAT.scores;


