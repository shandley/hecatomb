[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/latest_release_date.svg)](https://anaconda.org/bioconda/hecatomb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/platforms.svg)](https://anaconda.org/bioconda/hecatomb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/license.svg)](https://anaconda.org/bioconda/hecatomb)
[![Documentation Status](https://readthedocs.org/projects/hecatomb/badge/?version=latest)](https://hecatomb.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/downloads.svg)](https://anaconda.org/bioconda/hecatomb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

![](docs/img/hecatombLogo.png)

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to *'sacrifice'* false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

## Contents

- [Documentation](#documentation)
- [Quick Start Guide](#quick-start-guide)
- [Inputs](#inputs)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [Links](#links)

## Documentation

[Complete documentation is hosted at Read the Docs](https://hecatomb.readthedocs.io)

### Citation

[Hecatomb is currently on BioRxiv!](https://www.biorxiv.org/content/10.1101/2022.05.15.492003v1)

## Quick start guide

### Running on HPC

Hecatomb is powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/#) and greatly benefits from the use of 
Snakemake profiles for HPC Clusters.
[More information and example for setting up Snakemake profiles for Hecatomb in the documentation](https://hecatomb.readthedocs.io/en/latest/profiles/).

### Install

```bash
# create conda env and install
conda create -n hecatomb -c conda-forge -c bioconda hecatomb

# activate conda env
conda activate hecatomb

# check the installation
hecatomb --help

# download the databases - you only have to do this once
  # locally: using 8 threads (default is 32 threads)
hecatomb install --threads 8

  # HPC: using a snakemake profile named 'slurm'
hecatomb install --profile slurm
```

### Run the test dataset

```bash
# locally: uses 32 threads and 64 GB RAM by default
hecatomb test

# HPC: using a profile named 'slurm'
hecatomb test --profile slurm
```

## Inputs

Hecatomb can process paired- or single-end short-read sequencing, longread sequencing, 
and paired-end sequencing for round A/B library protocol.

```bash
hecatomb run --preprocessing paired
hecatomb run --preprocessing single
hecatomb run --preprocessing longread
hecatomb run --preprocessing roundAB
```

When you specify a directory of reads with `--reads` for paried-end sequencing, 
Hecatomb expects paired-end sequencing reads in the format sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

When you specify a TSV file with `--reads`, Hecatomb expects a 2- or 3-column tab separated file (depending on 
preprocessing method) with the first column specifying a sample name, and the other columns the relative or full paths 
to the forward (and reverse) read files. e.g.

```text
sample1    /path/to/reads/sample1.1.fastq.gz    /path/to/reads/sample1.2.fastq.gz
sample2    /path/to/reads/sample2.1.fastq.gz    /path/to/reads/sample2.2.fastq.gz
```

## Dependencies

The only dependency you need to get up and running with Hecatomb is [conda](https://docs.conda.io/en/latest/).
Hecatomb relies on [conda](https://docs.conda.io/en/latest/) (and [mamba](https://github.com/mamba-org/mamba))
to ensure portability and ease of installation of its dependencies.
All of Hecatomb's dependencies are installed during installation or runtime, so you don't have to worry about a thing!

## Citation

The Hecatomb preprint is available on BioRxiv:
[https://doi.org/10.1101/2022.05.15.492003](https://doi.org/10.1101/2022.05.15.492003)

## Links

[Hecatomb @ bio.tools](https://bio.tools/hecatomb)

[Hecatomb @ WorkflowHub](https://workflowhub.eu/workflows/235)

