#!/bin/bash

# contaminant_removal_htcf_modified.sbatch

#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-49
#SBATCH -J qv_cnt_rm
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=ALL

set -ue

module load bbtools

#########################################################################
### Important: You may need to set the HOSTPATH. The default is human ###
#########################################################################

### DESCRIPTION ###
# Script to remove non-biological sequences (primers, adapters), low-quality bases, host-sequences and obiovus (100% ID) bacterial sequences from virome sequenced libraries
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Set Variables
CONPATH=/scratch/ref/hecatomb_databases/contaminants # non-biological sequences
HOSTPATH=/scratch/ref/hecatomb_databases/human_masked # host sequence database
BACPATH=/scratch/ref/hecatomb_databases/bac_giant/unique_species # masked bacterial and giant virus genome

# Prep output directories
mkdir -p ./clumped
mkdir -p ./QC/step_{1..9}
mkdir -p ./runtimes

read1=$(sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt)
read2=$(echo ${read1} | sed s/R1/R2/g)

sample=$(echo ${read1} | sed s/_R1_001.fastq.gz//g)

# Begin time-log
echo -e "S0\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\ts9" > ./runtimes/${sample}_contaminant_removal_runtimes.txt

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
echo
echo "step 0"
echo

T0=$SECONDS

clumpify.sh \
	in=${read1} \
	in2=${read2} \
	out=./clumped/${sample}_R1.fastq.clumped.gz \
	out2=./clumped/${sample}_R2.fastq.clumped.gz \
	reorder=a \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_0=$(($SECONDS - $T0))

# Step 1: Remove leftmost primerB. Not the reverse complements
echo
echo "step 1"
echo

T1=$SECONDS

bbduk.sh \
	in=./clumped/${sample}_R1.fastq.clumped.gz \
	in2=./clumped/${sample}_R2.fastq.clumped.gz \
	ref=$CONPATH/primerB.fa \
	out=./QC/step_1/${sample}_R1.s1.out.fastq \
	out2=./QC/step_1/${sample}_R2.s1.out.fastq \
	stats=./QC/step_1/${sample}.s1.stats \
	k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
	removeifeitherbad=f \
	trimpolya=10 \
	ordered=t \
	rcomp=f \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_1=$(($SECONDS -	$T1))

# Step 2: Remove 3' read through contaminant
echo
echo "step 2"
echo

T2=$SECONDS

bbduk.sh \
	in=./QC/step_1/${sample}_R1.s1.out.fastq \
	in2=./QC/step_1/${sample}_R2.s1.out.fastq \
	ref=$CONPATH/rc_primerB_ad6.fa \
	out=./QC/step_2/${sample}_R1.s2.out.fastq \
	out2=./QC/step_2/${sample}_R2.s2.out.fastq \
	stats=./QC/step_2/${sample}.s2.stats	\
	k=16 hdist=1 mink=11 ktrim=r \
	removeifeitherbad=f \
	ordered=t \
	rcomp=f \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_2=$(($SECONDS -	$T2))

# Step 3: Remove primer free adapter (both orientations)
echo
echo "step 3"
echo

T3=$SECONDS

bbduk.sh \
	in=./QC/step_2/${sample}_R1.s2.out.fastq \
	in2=./QC/step_2/${sample}_R2.s2.out.fastq \
    ref=$CONPATH/nebnext_adapters.fa \
    out=./QC/step_3/${sample}_R1.s3.out.fastq \
	out2=./QC/step_3/${sample}_R2.s3.out.fastq \
	stats=./QC/step_3/${sample}.s3.stats \
    k=16 hdist=1 mink=10 ktrim=r \
    removeifeitherbad=f \
	ordered=t \
    rcomp=t \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_3=$(($SECONDS -	$T3))

# Step 4: Remove adapter free primer (both orientations)
echo
echo "step 4"
echo

T4=$SECONDS

bbduk.sh \
	in=./QC/step_3/${sample}_R1.s3.out.fastq \
	in2=./QC/step_3/${sample}_R2.s3.out.fastq \
	ref=$CONPATH/rc_primerB_ad6.fa \
	out=./QC/step_4/${sample}_R1.s4.out.fastq \
	out2=./QC/step_4/${sample}_R2.s4.out.fastq \
	stats=./QC/step_4/${sample}.s4.stats \
	k=16 hdist=0 \
	removeifeitherbad=f \
	ordered=t \
	rcomp=t \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_4=$(($SECONDS -	$T4))

# Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
echo
echo "step 5"
echo

T5=$SECONDS

bbduk.sh \
	in=./QC/step_4/${sample}_R1.s4.out.fastq \
	in2=./QC/step_4/${sample}_R2.s4.out.fastq \
    ref=$CONPATH/vector_contaminats.fa.gz \
    out=./QC/step_5/${sample}_R1.s5.out.fastq \
	out2=./QC/step_5/${sample}_R2.s5.out.fastq \
	stats=./QC/step_5/${sample}.s5.stats \
    k=31 hammingdistance=1 \
	ordered=t \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_5=$(($SECONDS -	$T5))

# Step 6: Host removal
echo
echo "step 6"
echo

T6=$SECONDS

bbmap.sh \
	in=./QC/step_5/${sample}_R1.s5.out.fastq \
	in2=./QC/step_5/${sample}_R2.s5.out.fastq \
	outu=./QC/step_6/${sample}_unmapped.s6.out.fastq \
	outm=./QC/step_6/${sample}_hostmapped.s6.out.fastq \
	semiperfectmode=t \
	quickmatch fast \
	ordered=t \
	path=$HOSTPATH \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

repair.sh \
	in=./QC/step_6/${sample}_unmapped.s6.out.fastq \
	out=./QC/step_6/${sample}_R1.s6.out.fastq \
	out2=./QC/step_6/${sample}_R2.s6.out.fastq \
	ow=t

RUNTIME_6=$(($SECONDS - $T6))

# Step 7: Trim low-quality bases
echo 
echo "step 7"
echo

T7=$SECONDS

bbduk.sh \
	in=./QC/step_6/${sample}_R1.s6.out.fastq \
	in2=./QC/step_6/${sample}_R2.s6.out.fastq \
	out=./QC/step_7/${sample}_R1.s7.out.fastq \
	out2=./QC/step_7/${sample}_R2.s7.out.fastq \
	outs=./QC/step_7/${sample}_singletons.s7.out.fastq \
	stats=./QC/step_7/${sample}.s7.stats \
	qtrim=r trimq=20 \
	maxns=2 minlength=50 \
	ordered=t \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_7=$(($SECONDS - $T7))

# Step 8: Merge forward (R1) and reverse (R2) reads
echo
echo "step 8"
echo

T8=$SECONDS

bbmerge.sh \
	in1=./QC/step_7/${sample}_R1.s7.out.fastq \
	in2=./QC/step_7/${sample}_R2.s7.out.fastq \
	out=./QC/step_8/${sample}_merged.fastq \
	outu1=./QC/step_8/${sample}_R1.unmerged.fastq \
	outu2=./QC/step_8/${sample}_R2.unmerged.fastq \
	rem k=62 extend2=50 ecct vstrict=t \
	ordered=t \
	-Xmx64g \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

# Split singletons and combine R1 and R2 files
grep -A 3 '1:N:' ./QC/step_7/${sample}_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/${sample}_singletons_R1.out.fastq
grep -A	3 '2:N:' ./QC/step_7/${sample}_singletons.s7.out.fastq | sed '/^--$/d' > ./QC/step_7/${sample}_singletons_R2.out.fastq

cat ./QC/step_8/${sample}_merged.fastq ./QC/step_8/${sample}_R1.unmerged.fastq ./QC/step_7/${sample}_singletons_R1.out.fastq > ./QC/step_8/${sample}_R1.s8.out.fastq;
cat ./QC/step_8/${sample}_merged.fastq ./QC/step_8/${sample}_R2.unmerged.fastq ./QC/step_7/${sample}_singletons_R2.out.fastq > ./QC/step_8/${sample}_R2.s8.out.fastq;

RUNTIME_8=$(($SECONDS - $T8))

# Step 9: Remove bacterial contaminants reserving viral and ambiguous sequences
echo
echo "step 9"
echo

T9=$SECONDS

bbmap.sh \
	in=./QC/step_8/${sample}_R1.s8.out.fastq \
	path=$BACPATH \
	outm=./QC/step_9/${sample}_bacterial.fastq \
	outu=./QC/step_9/${sample}_viral_amb.fastq \
	semiperfectmode=t \
	quickmatch fast \
	refstats=./QC/step_9/${sample}_refstats.txt \
	scafstats=./QC/step_9/${sample}_scafstats.txt \
	ordered=t \
	ow=t \
	t=${SLURM_CPUS_PER_TASK}

RUNTIME_9=$(($SECONDS - $T9))

echo -e "$RUNTIME_0\t$RUNTIME_1\t$RUNTIME_2\t$RUNTIME_3\t$RUNTIME_4\t$RUNTIME_5\t$RUNTIME_6\t$RUNTIME_7\t$RUNTIME_8\t$RUNTIME_9" >> ./runtimes/${sample}_contaminant_removal_runtimes.txt

# Add filenames to time table

#ls -1 *_R1.fastq.gz > filenames;
#sed -i '1ifilenames' filenames;
#paste filenames contaminant_removal_runtimes.txt >RUNTIMES_contaminant_removal.log;

