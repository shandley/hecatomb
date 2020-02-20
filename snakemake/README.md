# Snakemake Implementations of hecatomb

These [snakemake](https://snakemake.readthedocs.io/) implementations are designed to be minimally painful to install and run, and take advantage of the parallelizations that `snakemake` provides.

# Installation

1. Install [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
2. Install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
3. Clone this git repo:
```
cd
git clone https://github.com/shandley/hecatomb.git
```

3. Make a copy of the config file and edit it. Note: it is generally useful to copy the config file to the working directory, that way you can rerun it whenever needed!
```
cd my_working_directory
cp ~/hecatomb/configs/sample_config.json ./config.json
nano config.json
```

See the [configfile](#configfile) section below for more details on this file.

You are now read to run `hecatomb`!

# Running hecatomb

You will need a directory with some gzip-compressed `fastq` files in them. At the moment, `hecatomb` requires paired end reads (i.e a file whose name ends `_R1.fastq.gz` and a file whose name ends `_R2.fastq.gz`).

The [test_data/fastq](../test_data/fastq) directory has some example datasets that you can run through the pipeline to see if it is working.

Once you have edited the [config file](#config-file), you can run the pipeline. If you are runing on a single machine, this command should be appropriate:

```
snakemake --configfile config.json ~/hecatomb/snakemake/contaminant_removal.snakefile
```

You can take advantage of all the usual `snakemake` options, including running across a cluster:

```
mkdir sge_out sge_err
snakemake --configfile config.json --snakefile ~/hecatomb/snakemake/contaminant_removal.snakefile --cluster 'qsub -cwd -o sge_out -e sge_err -V -q default' -j 1000 --latency-wait 100
```

*Notes:*
1. I always add a --latency-wait when running across the cluster. It really helps.
2. Your `qsub` options will be slightly different
3. This makes `stdout` and `stderr` files in the `sge_out/` and `sge_err/` directories as appropriate.


# Config File

The `config file` can be in either [json](https://www.json.org/) or [YAML](https://yaml.org/), and we have provided the same file, called sample_config in both [json](sample_config.json) and [yaml](sample_config.yaml) formats. 

We currently need the following inputs:

- **Paths**
    - **Databases**: the path to the directory where the databases reside. If you have not [installed the databases](installing-the-databases) (see below), we will install them for you.
    - **Reads**: the path to the directory of `fastq` files. This is where your data should be.
- **DatabaseFiles**:
    - **bacteria**: The name of the bacterial database file. This is expected to be in the *Databases* path defined above.
    - **host**: The name of the host file. The default host is human, but you can remove any other sequences too.
    - **contaminants**: The name of the directory with potential contaminant sequences. By default we search for `nebnext_adapters.fa`, `primerB.fa`, ` rc_primerB_ad6.fa` and `vector_contaminants.fa`

# Installing the databases

We use several databases to remove contaminants (as shown in the config file section). We provide a compressed tarball that you can easily install yourself, or if you just run snakemake without installing the databases, it will install them for you! (And if they are installed, it will skip that step).

However, we know that sometimes downloading large files over the net causes issues, and so here are the steps to manually install the databases:

1. Create the directory that you want to use and change into that directory. This is the directory used in `sample_config.yaml` but you probably are not user `redwards` and so change that part as appropriate:

```
mkdir /home3/redwards/IBD/CERVAID/databases
cd /home3/redwards/IBD/CERVAID/databases
```

2. Download the databases file. If you have `curl`, I would recommend that:

```
curl -LO https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2
```

or if not, try wget

```
wget https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2
```

3. Extact the databases:

```
tar xf hecatomb.databases.tar.bz2
```

You're done! That was easy!


