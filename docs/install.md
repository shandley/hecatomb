## System requirements

We have set some sensible defaults for Hecatomb on the following minimum recommended system requirements:

 - 32 CPUs (or 16 CPUs with hyperthreading)
 - 64 GB of RAM
 - approximately 55 GB HDD space (for the databases)
 - Additional HDD space for temporary and output files (highly dependent on input file sizes)

Larger datasets may require more RAM (and CPUs to speed things along).
As such it is highly recommended to run Hecatomb on a HPC cluster.

## Dependencies

Hecatomb was developed as a set of [Snakemake](https://snakemake.readthedocs.io/en/stable/#) workflows, which are 
controlled by a python launcher for your convenience.
While the pipeline utilises a number of other programs and tools to run, these are all managed by the pipeline using 
[conda](https://docs.conda.io/en/latest/) and [mamba](https://github.com/mamba-org/mamba).

Install Hecatomb via conda:

```bash
# This command will create a new conda env called 'hecatomb' and will install Hecatomb and all of it's dependencies.
conda create -n hecatomb -c conda-forge -c bioconda hecatomb
```

That's it!

```bash
# To use Hecatomb, activate your new conda env.
conda activate hecatomb

# Check that it's installed
hecatomb -h
```

## Customisation

If you're running Hecatomb on a cluster it is highly recommended to use Snakemake profiles 
([Set up a profile for Snakemake](profiles.md)).

You may also want to customise the available resources for your system, whether you're using a cluster or running 
locally ([Advanced configuration](configuration.md)).

## Download the databases

Before running Hecatomb for the first time you will need to download the databases.
You will only need to do this step once (unless we update the databases).
You can rerun this step as much as you like; the pipeline will only download any database files that are missing.

```bash
# Either, run locally
hecatomb install

# Or, for a HPC cluster using a Snakemake profile (replace 'slurm' with your profile name)
hecatomb install --profile slurm
```

## Run the test dataset

Hecatomb comes with a test dataset that you can run which will take a few hours to complete.
Use `--test` in place of specifying your read directory with `--reads`.

```bash
# run locally
hecatomb run --test

# run on cluster using a Snakemake profile
hecatomb run --test --profile slurm
```

## Build the docs

These document pages can be built from the repo like so.

```bash
# install mkdocs
pip install mkdocs

# cd to your install directory
cd /path/to/hecatomb

# build
mkdocs build
```
