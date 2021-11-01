![](https://anaconda.org/beardymcjohnface/hecatomb/badges/platforms.svg)
[![](https://anaconda.org/beardymcjohnface/hecatomb/badges/license.svg)](https://opensource.org/licenses/MIT)
[![](https://anaconda.org/beardymcjohnface/hecatomb/badges/installer/conda.svg)](https://anaconda.org/beardymcjohnface/hecatomb)
![](https://anaconda.org/beardymcjohnface/hecatomb/badges/downloads.svg)
[![Documentation Status](https://readthedocs.org/projects/hecatomb/badge/?version=latest)](https://hecatomb.readthedocs.io/en/latest/?badge=latest)

![](docs/img/hecatombLogo.png)

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to *'sacrifice'* false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

**For detailed pipeline overview, installation, usage and customisation instructions,
[please refer to the documentation hosted at Read the Docs](https://hecatomb.readthedocs.io).**

## Quick start guide

### Running on HPC

Hecatomb is powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/#) and greatly benefits from the use of 
Snakemake profiles for HPC Clusters.
[More information and example for setting up Snakemake profiles for Hecatomb in the documentation](https://hecatomb.readthedocs.io/en/latest/advanced/#profiles-for-hpc-clusters).

### Install

```bash
# create conda env and install
conda create -n hecatomb -c conda-forge -c bioconda -c beardymcjohnface hecatomb

# activate conda env
conda activate hecatomb

# check the installation
hecatomb -h

# download the databases - you only have to do this once
  # locally: uses 32 threads by default
hecatomb install

  # HPC: using a profile named 'slurm'
hecatomb install --profile slurm
```

### Run the test dataset

```bash
# locally: uses 32 threads and 64 GB RAM by default
hecatomb run --test

# HPC: using a profile named 'slurm'
hecatomb run --test --profile slurm
```

### Current limitations

Hecatomb is currently designed to only work with paired-end reads. 
We have considered making a branch for single-end reads, but that is not currently available.

When you specify a directory of reads with `--reads`, Hecatomb expects paired sequencing reads in the format 
sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

When you specify a TSV file with `--reads`, Hecatomb expects a 3-column tab separated file with the first column
specifying a sample name, and the other columns the relative or full paths to the forward and reverse read files. e.g.

```text
sample1    /path/to/reads/sample1.1.fastq.gz    /path/to/reads/sample1.2.fastq.gz
sample2    /path/to/reads/sample2.1.fastq.gz    /path/to/reads/sample2.2.fastq.gz
```

### Dependencies

The only dependency you need to get up and running with Hecatomb is [conda](https://docs.conda.io/en/latest/).
Hecatomb relies on [conda](https://docs.conda.io/en/latest/) (and [mamba](https://github.com/mamba-org/mamba))
to ensure portability and ease of installation of its dependencies.
All of Hecatomb's dependencies are installed during installation or runtime, so you don't have to worry about a thing!

### Links

[Hecatomb @ bio.tools](https://bio.tools/hecatomb)


