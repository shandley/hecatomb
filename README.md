![](https://anaconda.org/beardymcjohnface/hecatomb/badges/platforms.svg)
[![](https://anaconda.org/beardymcjohnface/hecatomb/badges/license.svg)](https://opensource.org/licenses/MIT)
[![](https://anaconda.org/beardymcjohnface/hecatomb/badges/installer/conda.svg)](https://anaconda.org/beardymcjohnface/hecatomb)
![](https://anaconda.org/beardymcjohnface/hecatomb/badges/downloads.svg)
[![Documentation Status](https://readthedocs.org/projects/hecatomb/badge/?version=latest)](https://hecatomb.readthedocs.io/en/latest/?badge=latest)

# Hecatomb

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to *'sacrifice'* false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

**For detailed pipeline overview, installation, usage and customisation instructions,
[please refer to the documentation hosted at Read the Docs](https://hecatomb.readthedocs.io).**

## Quick start guide

### Install

```bash
# create conda env and install
conda create -n hecatomb -c conda-forge -c bioconda -c beardymcjohnface hecatomb

# activate conda env
conda activate hecatomb

# check the installation
hecatomb -h
```

### Running on HPC

Hecatomb is powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/#) and greatly benefits from the use of 
Snakemake profiles for HPC Clusters.
[More information and example for setting up Snakemake profiles for Hecatomb in the documentation](https://hecatomb.readthedocs.io/en/latest/advanced/#profiles-for-hpc-clusters).

### Run the test dataset

```bash
# locally: requires 32 threads and 64 GB RAM
hecatomb run --test --assembly

# HPC: using a profile named 'slurm'
hecatomb run --test --assembly --profile slurm
```

### Current limitations

Hecatomb is currently designed to only work with paired-end reads. 
We have considered making a branch for single-end reads, but that is not currently avaialbe. 
The workflow will crash if you do not supply paired-end reads at this time.

Hecatomb expects paired sequencing reads in the format sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

If your files don't follow this convention then you will need to rename them before running the pipeline.

### Dependencies

The only dependency you need to get up and running with Hecatomb is [conda](https://docs.conda.io/en/latest/).
Hecatomb relies on [conda](https://docs.conda.io/en/latest/) (and [mamba](https://github.com/mamba-org/mamba))
to ensure portability and ease of installation of its dependencies.
All of Hecatomb's dependencies are installed during installation or runtime, so you don't have to worry about a thing!


