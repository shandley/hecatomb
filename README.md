[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards-lab.science)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)                                                            
![GitHub language count](https://img.shields.io/github/languages/count/shandley/hecatomb)
[![Downloads](https://img.shields.io/github/downloads/shandley/hecatomb/total?style=flat-square)](https://github.com/shandley/hecatomb/releases)


# hecatomb

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice. In this case we provide a hecatomb of false-positive viral sequences. Hecatomb is desi

# What is `hecatomb`?

Hecatomb is a pipline for culling spurious sequences from viral metagenomes. 

# Is `hecatomb` for me?

If you have:
- viral metagenomes
- sequences that might be viral
- RoundA/B (DNA and cDNA from RNA viral genomes) sequences

Then hecatomb is for you!

# Installing `hecatomb`

Please see the detailed [installation instructions](INSTALLATION.md). Its really quite simple, you need conda and snakemake and we'll do the rest.

## Running hecatomb

You will need a directory with some gzip-compressed `fastq` files in them. At the moment, `hecatomb` requires paired end reads. The only requirement is that your paired end files contain \_R1 for the first mate pair and \_R2 for the second mate pair. We will figure out the file names by just looking for `_R1`. 

The [test_data/fastq](../test_data/fastq) directory has some example datasets that you can run through the pipeline to see if it is working.

You should be able to `cd` into the GitHub directory and run `snakemake` directly:

```bash
conda install -c bioconda -c conda-forge snakemake
git clone https://github.com/shandley/hecatomb.git
cd hecatomb
snakemake --configfile snakemake/config/sample_config.yaml -s snakemake/workflow/download_databases.smk --cores 4 --use-conda
snakemake --configfile snakemake/config/sample_config.yaml -s snakemake/workflow/Snakefile --cores 4 --use-conda
```

Our installation documentation includes help on setting up snakemake profiles, so that you can run that same command as:

```
cd hecatomb
snakemake --configfile configs/sample_config.json --snakefile snakemake/contaminant_removal.snakefile
```


## Config File

we recommend that you make a copy of the config file in your working directory, and edit it there. That way you can keep track of any changes you've made, and rerun code as ncessary.

```
cd my_working_directory
cp ~/hecatomb/configs/sample_config.json ./config.json
nano config.json
```