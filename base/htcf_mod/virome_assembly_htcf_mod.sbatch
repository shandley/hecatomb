#!/bin/bash

# virome_assembly_htcf_mod.sbatch

# Takes quality controlled reads (though contaminant_remval.sh Step 7) and performs a metagenomic assembly

# Steps:
# 0) Some stats file generation
# 1) Digital normalization
# 2) Assembly
# 3) Assembly decontamination
# ?) Sketch comparison
# ?) Virfinder
# ?) Bracken

#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-100%20 # change 100 to total number of samples
#SBATCH -J vir_asmbly
#SBATCH --mail-user=mihindu
#SBATCH --mail-type=ALL

# Immediately exit on errors
set -ue

# Modules
module load bbtools
module load megahit

# Prep output directories
mkdir -p ./assembly
mkdir -p ./assembly/stats
mkdir -p ./assembly/megahit_contigs
mkdir -p ./assembly/contig_graphs

# Set variables
IN=./QC/step_7
OUT=./assembly

read1=$(sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt)
sample=$(echo ${read1} | sed s/_R1.fastq.gz//g)

# Set file names
#for i in $IN/*_R1.s7.out.fastq
#do
#	F=`basename $i _R1.s7.out.fastq`

	# Step 0: Tabulate some kmer statistics
	echo
	echo "Viral Assembly: step_0. Tabulate kmer stats"
	echo
	bbcountunique.sh \
		in=$IN/${sample}_R1.s7.out.fastq \
		in2=$IN/${sample}_R2.s7.out.fastq \
		interval=2500 \
		out=$OUT/stats/${sample}_uniq_kmer_stats.txt \
		ow=t \
		-Xmx64g

	# Step 1: Digital Normalization
	echo
	echo "Viral Assembly: step_1. Digital Norm"
	echo
	bbnorm.sh \
		in=$IN/${sample}_R1.s7.out.fastq \
		in2=$IN/${sample}_R2.s7.out.fastq \
		extra=$IN/${sample}_singletons.s7.out.fastq \
		out=$OUT/${sample}_R1.norm.out.fastq \
		out2=$OUT/${sample}_R2.norm.out.fastq \
		outt=$OUT/${sample}_tossed.norm.fastq \
		target=20 \
		mindepth=2 \
		hist=$OUT/${sample}_norm.hist \
		ow=t \
		-Xmx64g
		#t=${SLURM_CPUS_PER_TASK}

	# Step 2: Assembly
	echo
	echo "Viral Assembly: step_2. Megahit"
	echo
	megahit \
		-1 $OUT/${sample}_R1.norm.out.fastq \
		-2 $OUT/${sample}_R2.norm.out.fastq \
		-o $OUT/${sample}_megahit_out \
		--out-prefix ${sample}.mh \
		-t ${SLURM_CPUS_PER_TASK}

	cp $OUT/${sample}_megahit_out/*.mh.contigs.fa $OUT/megahit_contigs
	#cp $OUT/${sample}_megahit_out/k99_${sample}.fastq $OUT/contig_graphs

	# Step 3: Quantification by mapping
	echo
	echo "Viral Assembly: step_3. Quant by mapping"
	echo
	bbmap.sh \
		ref=$OUT/megahit_contigs/${sample}.mh.contigs.fa \
		in=$IN/${sample}_R1.s7.out.fastq \
		in2=$IN/${sample}_R2.s7.out.fastq \
		out=$OUT/megahit_contigs/${sample}.aln.sam.gz \
		kfilter=22 subfilter=15 maxindel=80 \
		ow=t \
		t=${SLURM_CPUS_PER_TASK} \
		-Xmx64g

	# Calculate coverage
	pileup.sh in=$OUT/megahit_contigs/${sample}.aln.sam.gz out=$OUT/megahit_contigs/${sample}_cov.txt

	# Output mappeed/unmapped reads
	reformat.sh in=$OUT/megahit_contigs/${sample}.aln.sam.gz out=$OUT/megahit_contigs/${sample}_mapped.fastq mappedonly
     	reformat.sh in=$OUT/megahit_contigs/${sample}.aln.sam.gz out=$OUT/megahit_contigs/${sample}_unmapped.fastq unmappedonly
#done
