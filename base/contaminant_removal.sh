#!/bin/bash bash

# Set bash options for proper failure recognition
set -o errext
set -o pipefail

### DESCRIPTION ###
# Script to remove non-biological sequences (primers, adapters), low-quality bases, host-sequences and obiovus (100% ID)  bacterial sequences from virome sequenced libraries
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Set Variables
CONPATH=/mnt/data1/databases/contaminants # adapter/primers
HOSTPATH=/mnt/data1/databases/human_masked # host
BACPATH=/mnt/data1/databases/masked_dbs/genome-collection-combos/reformated/bac_giant # masked bacterial and giant virus genome

# Prep output directories
mkdir -p ./clumped
mkdir -p ./QC/step_{1..9};

# Set file names
for i in *_R1.fastq.gz; do
        F=`basename $i _R1.fastq.gz`;

# Summary:
	# Step 0: Clumpify reads (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/)
	# Step 1: Remove 5' amplification primer
	# Step 2: Remove 3' read through contaminant (Reverse complement of amplification primer + 6 bases of the adapter)
	# Step 3: Remove primer free adapter (both orientations)
	# Step 4: Remove adapter free primer (both orientations)
	# Step 5: PhiX Removal and vector contamination removal
	# Step 6: Host-removal
	# Step 7: Trim low-quality bases
	# Step 8: Read correction, merging and extension
	# Step 9: Remove bacterial contaminants reserving viral and aambiguous sequences

# Step 0: Clumpify and deduplicate reads
clumpify.sh in="$F"_R1.fastq.gz in2="$F"_R2.fastq.gz \
	out=./clumped/"$F"_R1.fastq.clumped.gz out2=./clumped/"$F"_R2.fastq.clumped.gz \
	reorder=a;

# Step 1: Remove leftmost primerB. Not the reverse complements
bbduk.sh in=./clumped/"$F"_R1.fastq.clumped.gz in2=./clumped/"$F"_R2.fastq.clumped.gz \
	ref=$CONPATH/primerB.fa \
	out=./QC/step_1/"$F"_R1.s1.out.fastq out2=./QC/step_1/"$F"_R2.s1.out.fastq \
	stats=./QC/step_1/"$F".s1.stats \
	k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
	removeifeitherbad=f \
	trimpolya=10 \
	ordered=t \
	rcomp=f \
	ow=t;

# Step 2: Remove 3' read through contaminant
bbduk.sh in=./QC/step_1/"$F"_R1.s1.out.fastq in2=./QC/step_1/"$F"_R2.s1.out.fastq \
	ref=$CONPATH/rc_primerB_ad6.fa \
	out=./QC/step_2/"$F"_R1.s2.out.fastq out2=./QC/step_2/"$F"_R2.s2.out.fastq \
	stats=./QC/step_2/"$F".s2.stats	\
	k=16 hdist=1 mink=11 ktrim=r \
        removeifeitherbad=f \
	ordered=t \
	rcomp=f \
	ow=t;

# Step 3: Remove primer free adapter (both orientations)
bbduk.sh in=./QC/step_2/"$F"_R1.s2.out.fastq in2=./QC/step_2/"$F"_R2.s2.out.fastq \
        ref=$CONPATH/nebnext_adapters.fa \
        out=./QC/step_3/"$F"_R1.s3.out.fastq out2=./QC/step_3/"$F"_R2.s3.out.fastq \
	stats=./QC/step_3/"$F".s3.stats	\
        k=16 hdist=1 mink=11 ktrim=r \
        removeifeitherbad=f \
	ordered=t \
        rcomp=t \
	ow=t;

# Step 4: Remove adapter free primer (both orientations)
bbduk.sh in=./QC/step_3/"$F"_R1.s3.out.fastq in2=./QC/step_3/"$F"_R2.s3.out.fastq \
	ref=$CONPATH/rc_primerB_ad6.fa \
	out=./QC/step_4/"$F"_R1.s4.out.fastq out2=./QC/step_4/"$F"_R2.s4.out.fastq \
	stats=./QC/step_4/"$F".s4.stats	\
	k=16 hdist=0 \
	removeifeitherbad=f \
	ordered=t \
	rcomp=t \
	ow=t;

# Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
bbduk.sh in=./QC/step_4/"$F"_R1.s4.out.fastq in2=./QC/step_4/"$F"_R2.s4.out.fastq \
        ref=$CONPATH/vector_contaminats.fa.gz \
        out=./QC/step_5/"$F"_R1.s5.out.fastq out2=./QC/step_5/"$F"_R2.s5.out.fastq \
	stats=./QC/step_5/"$F".s5.stats \
        k=31 hammingdistance=1 \
	ordered=t \
	ow=t;

# Step 6: Host removal
bbmap.sh in=./QC/step_5/"$F"_R1.s5.out.fastq in2=./QC/step_5/"$F"_R2.s5.out.fastq \
	outu=./QC/step_6/"$F"_unmapped.s6.out.fastq outm=./QC/step_6/"$F"_hostmapped.s6.out.fastq \
	semiperfectmode=t \
	quickmatch fast \
	ordered=t \
	path=$HOSTPATH \
	ow=t;

# Step 7: Trim low-quality bases
bbduk.sh in=./QC/step_6/"$F"_unmapped.s6.out.fastq \
	out=./QC/step_7/"$F"_R1.s7.out.fastq out2=./QC/step_7/"$F"_R2.s7.out.fastq outs=./QC/step_7/"$F"_singletons.s7.out.fastq \
	stats=./QC/step_7/"$F".s7.stats \
	qtrim=r trimq=30 \
	maxns=2 minlength=50 \
	ordered=t \
	ow=t;

# Step 8: Merge forward (R1) and reverse (R2) reads
bbmerge.sh in1=./QC/step_7/"$F"_R1.s7.out.fastq in2=./QC/step_7/"$F"_R2.s7.out.fastq \
	out=./QC/step_8/"$F"_merged.fastq outu1=./QC/step_8/"$F"_R1.unmerged.fastq outu2=./QC/step_8/"$F"_R2.unmerged.fastq \
	rem k=62 extend2=50 ecct vstrict=t \
	ordered=t \
	-Xmx128g \
	ow=t;

	# Split singletons and combine R1 and R2 files
	grep -A 3 '1:N:' ./QC/step_7/"$F"_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/"$F"_singletons_R1.out.fastq;
	grep -A	3 '2:N:' ./QC/step_7/"$F"_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/"$F"_singletons_R2.out.fastq;

	cat ./QC/step_8/"$F"_merged.fastq ./QC/step_8/"$F"_R1.unmerged.fastq ./QC/step_7/"$F"_singletons_R1.out.fastq > ./QC/step_8/"$F"_R1.s8.out.fastq;
	cat ./QC/step_8/"$F"_merged.fastq ./QC/step_8/"$F"_R2.unmerged.fastq ./QC/step_7/"$F"_singletons_R2.out.fastq > ./QC/step_8/"$F"_R2.s8.out.fastq;

# Step 9: Remove bacterial contaminants reserving viral and ambiguous sequences
bbmap.sh in=./QC/step_8/"$F"_R1.s8.out.fastq \
	path=$BACPATH \
	outm=./QC/step_9/"$F"_bacterial.fastq outu=./QC/step_9/"$F"_viral_amb.fastq \
	semiperfectmode=t \
	quickmatch fast \
	refstats=./QC/step_9/"$F"_refstats.txt scafstats=./QC/step_9/"$F"_scafstats.txt \
	ordered=t \
	ow=t;

done
