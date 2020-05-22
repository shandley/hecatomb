#!/bin/bash

# Takes quality controlled reads (through contaminant_remval.sh Step 7) and performs a metagenomic assembly

# Steps:
# 0) Kmer stats file generation
# 1) Digital normalization (bbnorm)
# 2) Individual sample assembly (megahit)
# 3) Contig assembly to create contig dictionary (flye)
# 4) Individual sample mapping to contig dictionary

# Prep output directories
mkdir -p ./assembly
mkdir -p ./assembly/stats
mkdir -p ./assembly/megahit_contigs
mkdir -p ./assembly/contig_dictionary

# Set variables
IN=./QC/step_7
OUT=./assembly
MIN=1000 # Minimum contig length

# Set file names
for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

# Step 0: Tabulate some kmer statistics
bbcountunique.sh in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq \
	interval=2500 \
	out=$OUT/stats/"$F"_uniq_kmer_stats.txt \
	ow=t;

# Step 1: Digital Normalization
bbnorm.sh in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq extra=$IN/"$F"_singletons_R1.out.fastq,$IN/"$F"_singletons_R2.out.fastq \
	out=$OUT/"$F"_R1.norm.out.fastq out2=$OUT/"$F"_R2.norm.out.fastq outt=$OUT/"$F"_tossed.norm.fastq \
	target=20 mindepth=2 \
	hist=$OUT/"$F"_norm.hist \
	ow=t;

# Step 2: Assembly
megahit -1 $OUT/"$F"_R1.norm.out.fastq -2 $OUT/"$F"_R2.norm.out.fastq -r $IN/"$F"_singletons_R1.out.fastq,$IN/"$F"_singletons_R2.out.fastq -o $OUT/"$F"_megahit_out \
	--out-prefix "$F".mh;

done;

find . -type f -name "*.mh.contigs.fa" -exec cp -rfp {} $OUT/megahit_contigs \;

# Step 3: Contig assembly to create contig dictionary (flye)
cat $OUT/megahit_contigs/*.mh.contigs.fa > $OUT/contig_dictionary/all.mh.contigs.fa;

reformat.sh in=$OUT/contig_dictionary/all.mh.contigs.fa out=$OUT/contig_dictionary/all.mh.contigs_trunc.fa \
	ml=$MIN \
	ow=t;

rename.sh in=$OUT/contig_dictionary/all.mh.contigs_trunc.fa out=$OUT/contig_dictionary/all.mh.contigs_for_flye.fa

flye --subassemblies $OUT/contig_dictionary/all.mh.contigs_for_flye.fa -t 24 --meta --plasmids -o $OUT/contig_dictionary -g 1g

reformat.sh in=$OUT/contig_dictionary/assembly.fasta out=$OUT/contig_dictionary/contig_dictionary.fasta \
        ml=$MIN \
        ow=t;

# Step 3: Quantification by mapping
# Set file names
for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

bbmap.sh ref=$OUT/contig_dictionary/contig_dictionary.fasta in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq \
	out=$OUT/contig_dictionary/"$F".aln.sam.gz \
        kfilter=22 subfilter=15 maxindel=80 \
        ow=t;

	# Calculate coverage
	pileup.sh in=$OUT/contig_dictionary/"$F".aln.sam.gz out=$OUT/contig_dictionary/"$F"_cov.txt;

	# Output mappeed/unmapped reads
	reformat.sh in=$OUT/contig_dictionary/"$F".aln.sam.gz out=$OUT/contig_dictionary/"$F"_mapped.fastq mappedonly;
     	reformat.sh in=$OUT/contig_dictionary/"$F".aln.sam.gz out=$OUT/contig_dictionary/"$F"_unmapped.fastq unmappedonly;

done

# Step 4: Contig Annotation with CAT
CAT contigs -c $OUT/contig_dictionary/contig_dictionary.fasta -o $OUT/contig_dictionary/out.CAT \
	--force \
	--sensitive \
	-d /mnt/data1/databases/CAT/CAT_prepare_20200304/2020-03-04_CAT_database \
	-t /mnt/data1/databases/CAT/CAT_prepare_20200304/2020-03-04_taxonomy;

CAT add_names -i $OUT/contig_dictionary/out.CAT.contig2classification.txt -o $OUT/contig_dictionary/out.CAT.taxonomy \
	-t /mnt/data1/databases/CAT/CAT_prepare_20200304/2020-03-04_taxonomy \
	--only_official --exclude_scores --force;

CAT summarise -c $OUT/contig_dictionary/contig_dictionary.fasta -i $OUT/contig_dictionary/out.CAT.taxonomy \
	-o $OUT/contig_dictionary/out.CAT.summary;

# Step 5: Create contig count table
DIR=./assembly/contig_dictionary

# Extract the average coverage from individual sample mappings
for i in $DIR/*_cov.txt; do

        F=`basename $i _cov.txt`;

        # Extract contig IDs and average fold coverage
        tail -n+2 $DIR/"$F"_cov.txt | cut -f2 > $DIR/"$F"_cov.avg;

        sed -i "1i $F" $DIR/"$F"_cov.avg;

done

# Combine individual average coverage files
paste $DIR/*_cov.avg > $DIR/average.coverage;

# Extract contig IDs and join with average coverages
cut -f1 $DIR/out.CAT.taxonomy > $DIR/contig.ids;
paste $DIR/contig.ids $DIR/average.coverage > $DIR/covtable.all;

# Remove temporary files
rm $DIR/*_cov.avg;
rm $DIR/contig.ids;
