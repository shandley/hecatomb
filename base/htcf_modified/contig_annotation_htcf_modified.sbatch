#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=1-12
#SBATCH -J ctg_anno
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=all

set -ue

module load seqkit
module load diamond/0.9.21
module load prodigal

# Prep output directories
mkdir -p ./assembly/megahit_contigs/annotation

# Set Variables
ANNOTATEDB=/scratch/mtblab/DYLAN/200210_unmapepd_bifid
IN=./assembly/megahit_contigs
OUT=./assembly/megahit_contigs/annotation
MIN=1000

contigFile=$(sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt)
sample=$(echo ${contigFile} | cut -d'.' -f 1)

# Set file names
#for i in $IN/*.mh.contigs.fa; do
        #F=`basename $i .mh.contigs.fa`;

# Filter out sequences shorter than length MIN
seqkit seq $IN/${contigFile} -m $MIN -j 24 -o $OUT/${sample}_n1000.contigs.fa

# CAT Annotations

/home/rachel.rodgers/software/CAT/CAT_pack/CAT contigs \
	-c $OUT/${sample}_n1000.contigs.fa \
	-d $ANNOTATEDB/2020-02-10_CAT_database \
	-t $ANNOTATEDB/2020-02-10_taxonomy \
	-n 4 \
	--out_prefix $OUT/${sample} \
	--index_chunks 1

/home/rachel.rodgers/software/CAT/CAT_pack/CAT add_names \
	-i $OUT/${sample}.contig2classification.txt \
	-o $OUT/${sample}.contig2classification.official_names.txt \
	-t $ANNOTATEDB/2020-02-10_taxonomy \
	--only_official

/home/rachel.rodgers/software/CAT/CAT_pack/CAT summarise \
	-c $OUT/${sample}_n1000.contigs.fa \
	-i $OUT/${sample}.contig2classification.official_names.txt \
	-o $OUT/${sample}_CAT.summary.txt
