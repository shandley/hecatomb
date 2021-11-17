## Commands

* `hecatomb install` - Install the databases (you should only need to do this once)
* `hecatomb run` - Run the pipeline
* `hecatomb listHosts` - List the currently-available host genomes
* `hecatomb addHost` - Add your own host genome
* `hecatomb config` - Copy the default config file to the current directory (for use with `--configfile`)

## Input

Hecatomb is currently designed to only work with paired-end reads.
You can either specify a directory of reads, and Hecatomb will infer the sample names and forward/reverse files, or,
you can specify a TSV file to explicitly assign sample names and point to the corresponding read files.
In either case you just use `--reads` and Hecatomb will figure out if it's a file or directory.

When you specify a directory of reads, e.g. `hecatomb run --reads readDir/`, 
Hecatomb expects paired sequencing reads in the format sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

When you specify a TSV file, e.g. `hecatomb run --reads samples.tsv`, 
Hecatomb expects a 3-column tab separated file with the first column specifying the sample name, 
and the other columns the relative or full paths to the forward and reverse read files. e.g.

```text
sample1    /path/to/reads/sample1.1.fastq.gz    /path/to/reads/sample1.2.fastq.gz
sample2    /path/to/reads/sample2.1.fastq.gz    /path/to/reads/sample2.2.fastq.gz
```

## Read annotation + assembly

By default, Hecatomb will annotate your reads and perform an assembly.
If you have more than 32 threads available, you can increase the threads provided to the pipeline with `--threads`:

```bash
hecatomb run --reads fastq/ --threads 64
```

If you're running on a HPC cluster, you should first set up a 
[Snakemake Profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
[More info and example for Hecatomb here](configuration.md#profiles-for-hpc-clusters).
Then you would specify your profile name when running Hecatomb.
Assuming your profile is called `slurm`:

```bash
hecatomb run --reads fastq/ --profile slurm
```

Running Hecatomb on a HPC with a Snakemake profile is THE BEST WAY to run the pipeline.

## Read annotation only

To optionally skip generating an assembly when running Hecatomb, 
the command is exactly the same as above with the addition of the `--skipAssembly` flag:

```bash
hecatomb run --reads fastq/ --profile slurm --skipAssembly
```

## Quicker read annotation

The pipeline bottleneck is the MMSeqs searches.
Use the `--fast` flag to run Hecatomb with less sensitive settings for MMSeqs.
In limited testing, we find it performs almost as well but with considerable runtime improvements.

```bash
hecatomb run --reads fastq/ --profile slurm --fast
```

## Specifying a host genome

Hecatomb includes a thorough host read removal step which utilises a processed host genome.
You can specify a host, or add your own.

By default, Hecatomb will use the human genome.
If your sample is from a different source you will need to specify the host genome for your sample source.

To see what host genomes are available:

```bash
hecatomb listHosts
```

The following should be available by default: 
bat, mouse, camel, celegans, macaque, rat, dog, cat, tick, mosquito, cow, human

So if you are working with mouse samples you would run:

```bash
hecatomb run --reads fastq/ --host mouse
```

## Add your own host genome

If the genome for the host you're working with isn't included in the available hosts, or you have a reference genome
which you think is better, you can add it with `addHost`.
This script will mask viral-like regions from your genome and add it to your Hecatomb host database.

You will need to specify the host genome FASTA file, as well as a name for this host.
Assuming you want to add the llama genome and the FASTA genome file is called `llama.fasta`:

```bash
hecatomb addHost --host llama --hostfa llama.fasta
```

You will then be able to run Hecatomb with your new host genome:

```bash
hecatomb run --reads fastq/ --host llama
```
