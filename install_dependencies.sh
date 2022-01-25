#!/bin/bash

wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
echo "Add me to your PATH by editing your ~/.bashrc or ~/.bash_profile"

echo "installing snakemake and mamba"
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
