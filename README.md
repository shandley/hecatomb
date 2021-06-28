[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)                                                            
![GitHub language count](https://img.shields.io/github/languages/count/shandley/hecatomb)
[![Downloads](https://img.shields.io/github/downloads/shandley/hecatomb/total?style=flat-square)](https://github.com/shandley/hecatomb/releases)


# hecatomb

A [hecatomb](https://en.wiktionary.org/wiki/hecatomb) is a great sacrifice or an extensive loss. Heactomb the software empowers an analyst to make data driven decisions to *'sacrifice'* false-positive viral reads from metagenomes to enrich for true-positive viral reads. This process frequently results in a great loss of suspected viral sequences / contigs.

Hecatomb was developed in response to the challenges associated with the detection of viral sequences in metagenomes. Virus detection or virome profiling is typically performed on samples containing extensive host nucleic acid (e.g. a tissue biopsy) or nucleic acid from a variety of other organisms such as bacteria, fungi and archaea from soil samples or bacteria, fungi, archaea and food from mammalian stool samples. All of these non-viral nucleic acid types constitute the *background* and present a variety of issues that must be evaluated in order to confidently detect viral sequences.

An additional challenge when attempting to detect viral sequences in metagenomes is the propensity of viruses to rapidly acquire mutations and evolve. Thus, it is possible that the viral sequences in your samples are likely relatively distinct to those in any known reference database. Thus, when assigning taxonomy of your unknown sequences by comparing to any reference database one must take advantage of algorithms that are able to detect statistically significant sequence similarity over large evolutionary distances. Hecatomb takes advantage of [MMSeqs2](https://github.com/soedinglab/MMseqs2) and succesive search strategies to enable to detection of evolutionary distinct viral sequences.

Finally, hecatomb utilizes information in individual reads **and** in metagenome assembled contigs. There is frequently valuable information available only at the read level as metagenome assembly is challenging even at its best. Thus hecatomb pairs analysis of both reads and contigs to provide a rich and detailed data structure for deep virome investigation.

# Installing and running `hecatomb`

Hecatomb was developed as a [Snakemake](https://snakemake.readthedocs.io/en/stable/#) workflow. In theory, one only needs to download the Snakefile (and associated files such as the environmental yaml's) and the hecatomb databases, create a configuration run file and everything should work. However, there may be some system specific considerations (e.g. system RAM, available CPUs, cluster environment specifications) that will require some tuning. We are working to provide settings that will work on most environments and will post additional instructions for different systems on the [hecatomb wiki](https://github.com/shandley/hecatomb/wiki), but can not gauruntee some level of system customization may still be required. For more specific assistance please post an [issue]((https://github.com/shandley/hecatomb/issues).

**Current Limitations**

Hecatomb is currently designed to only work with paired-end reads. We have considered making a branch for single-end reads, but that is not currently avaialbe. The workflow will crash if you do not supply paired-end reads at this time.

> _Tip:_ Currently, `hecatomb` expects an R1 file and an R2. Because of a few limitations, it is expected that all the files have the same ending after `_R1` and `_R2`. This is very occassionaly a problem, for example if you have something like `sample1_R1_001.fastq.gz` and `sample2_R1_002.fastq.gz`. In that case, rename the `002` file to `001`. We are working on a general solution to this problem.

**Dependencies**

Hecatomb relies on [conda](https://docs.conda.io/en/latest/) to ensure portability and ease of installation. So the only two dependencies you will need to install are [Snakemake](https://snakemake.readthedocs.io/en/stable/#) and [conda](https://docs.conda.io/en/latest/).

We highly recommend you use [mamba](https://github.com/mamba-org/mamba) as it is _a lot_ faster than the base conda.

```bash
conda install -c conda-forge -c bioconda snakemake mamba
```

Once you have Snakemake and Conda installed you will need to complete the following 4-steps prior to running hecatomb. Of note, Steps 1 and 2 will only need to be run once.

**Step 1:** Clone the hecatomb repository.

```bash
git clone https://github.com/shandley/hecatomb.git
```

**Step 2:** Download the hecatomb databases

TBD

**Step 3:** Modify the config file

- Move into the config file directory
- Make a copy of the provided sample config file (sample_config.yaml) and give it your own name
- Modify the config file with your project specifications

```bash
cd /hecatomb/snakemake/config
cp sample_config.yaml my_config.yaml
nano my_config.yaml
```

There are 4 things you will need to update in your new config.yaml

1. The location of the databases you downloaded in Step 2
2. The location of your unprocessed (raw) fastq files
3. The destination directory of your results
4. The name of the host species your samples came from (e.g. mouse, human, bat)

There are additional notes in the provided sample_config.yaml on each of these topics. There are also several other variables you can modify in your configuration file but these are completely optional tuning paramaters.

4) Run hecatomb

Snakemake enables a good deal of customization from the [command line](https://snakemake.readthedocs.io/en/stable/executing/cli.html). An example launch command is below. In this case the command is launched from the directory where the Snakefile is located (/workflow).

```bash
snakemake --snakefile ./Snakefile --configfile ../config/my_config.yaml --resources mem_mb=100000 --cores 64 --use-conda --conda-frontend mamba
```

With this command we are:

- Specifying the location of the Snakefile (--snakefile ./Snakefile)
- Specifying the location of our config file (--configfile ../config/my_config.yaml)
- Specifying the amount of ram (in mb) and number of cores we want to use to run the entire workflow. In this case we ran this command on a server with 128Gb of ram and 64 cores. --resources mem_mb=100000 --cores 64 will specify that 100Gb of ram and 64 cores will be reserved for the running of the entire workflow
- Tell Snakemake to download and install software using conda

## Using hecatomb on a cluster

You can read several recommendations on running hecatomb on a cluster using `slurm` or `sge` on our wiki

## Test Data

You will need a directory with some gzip-compressed `fastq` files in them. At the moment, `hecatomb` requires paired end reads. The only requirement is that your paired end files contain \_R1 for the first mate pair and \_R2 for the second mate pair. We will figure out the file names by just looking for `_R1`. 

The [test_data/fastq](../test_data/fastq) directory has some example datasets that you can run through the pipeline to see if it is working.

You should be able to `cd` into the GitHub directory and run `snakemake` directly:



