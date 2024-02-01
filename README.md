![](hecatombLogo.png)

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/license.svg)
![Anaconda-Server Badge](https://anaconda.org/bioconda/hecatomb/badges/latest_release_date.svg)
[![Documentation Status](https://readthedocs.org/projects/hecatomb/badge/?version=latest&style=flat-square)](https://hecatomb.readthedocs.io/en/latest/?badge=latest)
[![install with bioconda](https://img.shields.io/badge/Install%20with-conda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/hecatomb/README.html)
![](https://img.shields.io/conda/dn/bioconda/hecatomb?label=Conda%20downloads&style=flat-square)
[![install with PyPI](https://img.shields.io/badge/Install%20with-PyPI-brightgreen.svg?style=flat-square)](https://pypi.org/project/hecatomb/)
[![Unit tests](https://github.com/shandley/hecatomb/actions/workflows/unit-tests.yaml/badge.svg)](https://github.com/shandley/hecatomb/actions/workflows/unit-tests.yaml)
[![Env builds](https://github.com/shandley/hecatomb/actions/workflows/build-hecatomb-envs.yaml/badge.svg)](https://github.com/shandley/hecatomb/actions/workflows/build-hecatomb-envs.yaml)

---

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to *'sacrifice'* false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

## Contents

- [Documentation](#documentation)
- [Citation](#citation)
- [Quick Start Guide](#quick-start-guide)
- [Inputs](#inputs)
- [Dependencies](#dependencies)
- [Links](#links)

## Documentation

[Complete documentation is hosted at Read the Docs](https://hecatomb.readthedocs.io)

## Citation

[Hecatomb is currently on BioRxiv!](https://www.biorxiv.org/content/10.1101/2022.05.15.492003v1)

## Quick start guide

### Snakemake profiles (for running on HPCs)

Hecatomb is powered by [Snakemake](https://snakemake.readthedocs.io/en/stable/#) and greatly benefits from the use of 
Snakemake profiles for HPC Clusters.
[More information and example for setting up Snakemake profiles for Hecatomb in the documentation](https://hecatomb.readthedocs.io/en/latest/profiles/).

### Install Hecatomb

__option 1: PIP__

```bash
# Optional: create a virtual with conda or venv
conda create -n hecatomb python=3.10

# activate
conda activte hecatomb

# Install
pip install hecatomb
```

__option 2: Conda__

```bash
# Create the conda env and install hecatomb in one step
conda create -n hecatomb -c conda-forge -c bioconda hecatomb

# activate
conda activate hecatomb
```

__Check installation__

```bash
hecatomb --help
```

### Install databases

```bash
# locally: using 8 threads (default is 32 threads)
hecatomb install --threads 8

# HPC: using a snakemake profile named 'slurm'
hecatomb install --profile slurm
```

### Run test dataset

```bash
# locally: using 32 threads and 64 GB RAM by default
hecatomb test

# HPC: using a profile named 'slurm'
hecatomb test --profile slurm
```

## Inputs

### Parsing samples with `--reads`

You can pass either a directory of reads or a TSV file to `--reads`. 
Note that Hecatomb expects your read file names to include common R1/R2 tags. 
 - __Directory:__ Hecatomb will infer sample names and \_R1/\_R2 pairs from the filenames.
 - __TSV file:__ Hecatomb expects 2 or 3 columns, with column 1 being the sample name and columns 2 and 3 the reads files.

[More information and examples are available here](https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8#file-readme-md)

### Library preprocessing with `--trim`

Hecatomb uses [Trimnami](https://github.com/beardymcjohnface/Trimnami) for read trimming which supports many different
trimming methods. Current options are `fastp` (default), `prinseq`, `roundAB`, `filtlong` (longreads), 
`cutadapt` (FASTA input), and `notrim` (skip trimming). See Trimnami's documentation for more information.

## Dependencies

The only dependency you need to get up and running with Hecatomb is [conda](https://docs.conda.io/en/latest/) or 
the python package manager [pip](https://pypi.org/project/pip/).
Hecatomb relies on [conda](https://docs.conda.io/en/latest/) (and [mamba](https://github.com/mamba-org/mamba))
to ensure portability and ease of installation of its dependencies.
All of Hecatomb's dependencies are installed during installation or runtime, so you don't have to worry about a thing!

## Links

[Hecatomb @ PyPI](https://pypi.org/project/hecatomb/)

[Hecatomb @ bioconda](https://bioconda.github.io/recipes/hecatomb/README.html)

[Hecatomb @ bio.tools](https://bio.tools/hecatomb)

[Hecatomb @ WorkflowHub](https://workflowhub.eu/workflows/235)

