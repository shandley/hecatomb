#!/bin/bash

#########################################################################
### Important: You may need to set the HOSTPATH. The default is human ###
#########################################################################

### DESCRIPTION ###
# Script to remove non-biological sequences (primers, adapters), low-quality bases, host-sequences and obiovus (100% ID) bacterial sequences from virome sequenced libraries
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Set Variables
CONPATH=/mnt/data1/databases/hecatomb/contaminants # non-biological sequences
HOSTPATH=/mnt/data1/databases/hecatomb/host/human # host sequence database
BACPATH=/mnt/data1/databases/hecatomb/bacteria # masked bacterial and giant virus genome

# Prep output directories
mkdir -p ./clumped
mkdir -p ./QC/step_{1..8};

# Begin time-log
echo -e "S0\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8" > contaminant_removal_runtimes.txt

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
	# Step 8: Remove bacterial contaminants reserving viral and aambiguous sequences

# Step 0: Clumpify and deduplicate reads
T0=$SECONDS
clumpify.sh in="$F"_R1.fastq.gz in2="$F"_R2.fastq.gz \
	out=./clumped/"$F"_R1.fastq.clumped.gz out2=./clumped/"$F"_R2.fastq.clumped.gz \
	reorder=a \
	ow=t;
RUNTIME_0=$(($SECONDS - $T0))

# Step 1: Remove leftmost primerB. Not the reverse complements
T1=$SECONDS
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
RUNTIME_1=$(($SECONDS -	$T1))

# Step 2: Remove 3' read through contaminant
T2=$SECONDS
bbduk.sh in=./QC/step_1/"$F"_R1.s1.out.fastq in2=./QC/step_1/"$F"_R2.s1.out.fastq \
	ref=$CONPATH/rc_primerB_ad6.fa \
	out=./QC/step_2/"$F"_R1.s2.out.fastq out2=./QC/step_2/"$F"_R2.s2.out.fastq \
	stats=./QC/step_2/"$F".s2.stats	\
	k=16 hdist=1 mink=11 ktrim=r \
        removeifeitherbad=f \
	ordered=t \
	rcomp=f \
	ow=t;
RUNTIME_2=$(($SECONDS -	$T2))

# Step 3: Remove primer free adapter (both orientations)
T3=$SECONDS
bbduk.sh in=./QC/step_2/"$F"_R1.s2.out.fastq in2=./QC/step_2/"$F"_R2.s2.out.fastq \
        ref=$CONPATH/nebnext_adapters.fa \
        out=./QC/step_3/"$F"_R1.s3.out.fastq out2=./QC/step_3/"$F"_R2.s3.out.fastq \
	stats=./QC/step_3/"$F".s3.stats	\
        k=16 hdist=1 mink=10 ktrim=r \
        removeifeitherbad=f \
	ordered=t \
        rcomp=t \
	ow=t;
RUNTIME_3=$(($SECONDS -	$T3))

# Step 4: Remove adapter free primer (both orientations)
T4=$SECONDS
bbduk.sh in=./QC/step_3/"$F"_R1.s3.out.fastq in2=./QC/step_3/"$F"_R2.s3.out.fastq \
	ref=$CONPATH/rc_primerB_ad6.fa \
	out=./QC/step_4/"$F"_R1.s4.out.fastq out2=./QC/step_4/"$F"_R2.s4.out.fastq \
	stats=./QC/step_4/"$F".s4.stats	\
	k=16 hdist=0 \
	removeifeitherbad=f \
	ordered=t \
	rcomp=t \
	ow=t;
RUNTIME_4=$(($SECONDS -	$T4))

# Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
T5=$SECONDS
bbduk.sh in=./QC/step_4/"$F"_R1.s4.out.fastq in2=./QC/step_4/"$F"_R2.s4.out.fastq \
        ref=$CONPATH/vector_contaminats.fa.gz \
        out=./QC/step_5/"$F"_R1.s5.out.fastq out2=./QC/step_5/"$F"_R2.s5.out.fastq \
	stats=./QC/step_5/"$F".s5.stats \
        k=31 hammingdistance=1 \
	ordered=t \
	ow=t;
RUNTIME_5=$(($SECONDS -	$T5))

# Step 6: Host removal
T6=$SECONDS
bbmap.sh in=./QC/step_5/"$F"_R1.s5.out.fastq in2=./QC/step_5/"$F"_R2.s5.out.fastq \
	outu=./QC/step_6/"$F"_unmapped.s6.out.fastq outm=./QC/step_6/"$F"_hostmapped.s6.out.fastq \
	semiperfectmode=t \
	quickmatch fast \
	ordered=t \
	path=$HOSTPATH \
	-Xmx48g;
	ow=t;

repair.sh in=./QC/step_6/"$F"_unmapped.s6.out.fastq \
	out=./QC/step_6/"$F"_R1.s6.out.fastq out2=./QC/step_6/"$F"_R2.s6.out.fastq \
	ow=t;

RUNTIME_6=$(($SECONDS - $T6))

# Step 7: Trim low-quality bases
T7=$SECONDS
bbduk.sh in=./QC/step_6/"$F"_R1.s6.out.fastq in2=./QC/step_6/"$F"_R2.s6.out.fastq \
	out=./QC/step_7/"$F"_R1.s7.out.fastq out2=./QC/step_7/"$F"_R2.s7.out.fastq outs=./QC/step_7/"$F"_singletons.s7.out.fastq \
	stats=./QC/step_7/"$F".s7.stats \
	qtrim=r trimq=20 \
	maxns=2 minlength=50 \
	ordered=t \
	ow=t;

	# Split singletons and combine R1 and R2 files
	grep -A 3 '1:N:' ./QC/step_7/"$F"_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/"$F"_singletons_R1.out.fastq;
	grep -A 3 '2:N:' ./QC/step_7/"$F"_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/"$F"_singletons_R2.out.fastq;

	cat ./QC/step_7/"$F"_R1.s7.out.fastq ./QC/step_7/"$F"_singletons_R1.out.fastq > ./QC/step_7/"$F"_R1.s7.combined.out.fastq;
        cat ./QC/step_7/"$F"_R2.s7.out.fastq ./QC/step_7/"$F"_singletons_R2.out.fastq > ./QC/step_7/"$F"_R2.s7.combined.out.fastq;

RUNTIME_7=$(($SECONDS - $T7))

# Step 8: Remove bacterial contaminants reserving viral and ambiguous sequences
T8=$SECONDS
bbmap.sh in=./QC/step_7/"$F"_R1.s7.combined.out.fastq \
	path=$BACPATH \
	outm=./QC/step_8/"$F"_bacterial.fastq outu=./QC/step_8/"$F"_viral_amb.fastq \
	semiperfectmode=t \
	quickmatch fast \
	refstats=./QC/step_8/"$F"_refstats.txt scafstats=./QC/step_8/"$F"_scafstats.txt \
	ordered=t \
	-Xmx48g \
	ow=t;
RUNTIME_8=$(($SECONDS - $T8))

echo -e "$RUNTIME_0\t$RUNTIME_1\t$RUNTIME_2\t$RUNTIME_3\t$RUNTIME_4\t$RUNTIME_5\t$RUNTIME_6\t$RUNTIME_7\t$RUNTIME_8" >> contaminant_removal_runtimes.txt

done

# Add filenames to time table

ls -1 *_R1.fastq.gz > filenames;
sed -i '1ifilenames' filenames;
paste filenames contaminant_removal_runtimes.txt >RUNTIMES_contaminant_removal.log;
