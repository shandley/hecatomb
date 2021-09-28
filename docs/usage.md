# Running Hecatomb

## Commands

* `hecatomb install` - Install the databases (you should only need to do this once)
* `hecatomb run` - Run the pipeline
* `hecatomb listHosts` - List the currently-available host genomes
* `hecatomb addHost` - Add your own host genome

## Run Hecatomb

### Read directory

Hecatomb currently expects paired sequencing reads in the format sampleName_R1/R2.fastq(.gz). e.g. 

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

If your files don't follow this convention then you will need to rename them before running the pipeline.

### Read annotation

To annotate your reads, assuming your reads are in a directory called `fastq/`:

```bash
hecatomb run --reads fastq/
```

If you have more than 32 threads available, you can increase the threads provided to the pipeline with `--threads`:

```bash
hecatomb run --reads fastq/ --threads 64
```

If you're running on a HPC cluster, you should first set up a 
[Snakemake Profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
[More info and example for Hecatomb here](advanced.md#profiles-for-hpc-clusters).
Then you would specify your profile name when running Hecatomb.
Assuming your profile is called `slurm`:

```bash
hecatomb run --reads fastq/ --profile slurm
```

Running Hecatomb on a HPC with a Snakemake profile is THE BEST WAY to run the pipeline.

### Read annotation + assembly

To optionally perform an assembly while running Hecatomb, the command is exactly the same as above with the addition
of the `--assembly` flag:

```bash
hecatomb run --reads fastq/ --profile slurm --assembly
```

That's it!

### Specifying a host genome

Hecatomb includes a thorough host read removal step and by default it will use the human genome.
If your sample is from a different source you will need to specify the host genome for your sample source.

To see what host genomes are available:

```bash
hecatomb listHosts
```

The following should be available by default: 
bat, mouse, camel, celegans, macaque, rat, dog, cat, tick, mosquito, cow, human

If you are working with mouse samples:

```bash
hecatomb run --reads fastq/ --host mouse
```

### Adding your own host genome

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
