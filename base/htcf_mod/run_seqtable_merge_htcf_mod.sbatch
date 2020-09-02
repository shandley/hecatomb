#!/bin/bash

# run_seqtable_merge_htcf_mod.sbatch

#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -J sqMrg
#SBATCH --mail-user=mihindu
#SBATCH --mail-type=ALL

module load r/3.5.3-python-3.6.5-java-11

Rscript /scratch/sahlab/2020_03_20_IBS_virome/seqtable_merge.R
