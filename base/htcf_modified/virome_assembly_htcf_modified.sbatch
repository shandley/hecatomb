#!/bin/bash

# virome_assembly_htcf_modified.sbatch

# Takes quality controlled reads (though contaminant_remval.sh Step 7) and performe a metagenomic assembly

# Steps:
# 0) Some stats file generation
# 1) Digital normalization
# 2) Assembly
# 3) Assembly decontamination
# ?) Sketch comparison
# ?) Virfinder
# ?) Bracken

# Modules
module load bbtools
module load megahit

# Immediately exit on errors
set -ue

#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -J vir_asmbly
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=ALL

# Prep output directories
mkdir -p ./assembly
mkdir -p ./assembly/stats
mkdir -p ./assembly/megahit_contigs
mkdir -p ./assembly/contig_graphs

# Set variables
IN=./QC/step_7
OUT=./assembly

# Set file names
for i in $IN/*_R1.s7.out.fastq
do
	F=`basename $i _R1.s7.out.fastq`

	# Step 0: Tabulate some kmer statistics
	echo
	echo "0. Tabulate kmer stats"
	echo
	bbcountunique.sh \
		in=$IN/"$F"_R1.s7.out.fastq \
		in2=$IN/"$F"_R2.s7.out.fastq \
		interval=2500 \
		out=$OUT/stats/"$F"_uniq_kmer_stats.txt \
		ow=t \
		-Xmx64g

	# Step 1: Digital Normalization
	echo
	echo "1. Digital Norm"
	echo
	bbnorm.sh \
		in=$IN/"$F"_R1.s7.out.fastq \
		in2=$IN/"$F"_R2.s7.out.fastq \
		extra=$IN/"$F"_singletons.s7.out.fastq \
		out=$OUT/"$F"_R1.norm.out.fastq \
		out2=$OUT/"$F"_R2.norm.out.fastq \
		outt=$OUT/"$F"_tossed.norm.fastq \
		target=20 \
		mindepth=2 \
		hist=$OUT/"$F"_norm.hist \
		ow=t \
		-Xmx64g
		#t=${SLURM_CPUS_PER_TASK}

	# Step 2: Assembly
	echo
	echo "2. Megahit"
	echo
	megahit \
		-1 $OUT/"$F"_R1.norm.out.fastq \
		-2 $OUT/"$F"_R2.norm.out.fastq \
		-o $OUT/"$F"_megahit_out \
		--out-prefix "$F".mh \
		-t ${SLURM_CPUS_PER_TASK}

	cp $OUT/"$F"_megahit_out/*.mh.contigs.fa $OUT/megahit_contigs
	cp $OUT/"$F"_megahit_out/k99_"$F".fastg $OUT/contig_graphs

	# Step 3: Quantification by mapping
	echo
	echo "3. Quant by mapping"
	echo
	bbmap.sh \
		ref=$OUT/megahit_contigs/"$F".mh.contigs.fa \
		in=$IN/"$F"_R1.s7.out.fastq \
		in2=$IN/"$F"_R2.s7.out.fastq \
		out=$OUT/megahit_contigs/"$F".aln.sam.gz \
		kfilter=22 subfilter=15 maxindel=80 \
		ow=t \
		t=${SLURM_CPUS_PER_TASK} \
		-Xmx64g

	# Calculate coverage
	pileup.sh in=$OUT/megahit_contigs/"$F".aln.sam.gz out=$OUT/megahit_contigs/"$F"_cov.txt

	# Output mappeed/unmapped reads
	reformat.sh in=$OUT/megahit_contigs/"$F".aln.sam.gz out=$OUT/megahit_contigs/"$F"_mapped.fastq mappedonly
     	reformat.sh in=$OUT/megahit_contigs/"$F".aln.sam.gz out=$OUT/megahit_contigs/"$F"_unmapped.fastq unmappedonly
done 
