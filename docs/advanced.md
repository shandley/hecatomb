# Advanced options for Hecatomb
If you're running Hecatomb on a HPC cluster, we absolutely recommend setting up a 
[Snakemake profile](advanced.md#snakemake-profiles-for-hpc-clusters).

There are many parameters that can be changed to tailor the pipeline to your system, 
or to tweak the pipeline's filtering.

We also recommend reviewing the Snakemake `config.yaml` file in your Hecatomb installation directory.
The config file will be at `hecatomb/snakemake/config/config.yaml` and you can find your installation directory with:

```bash
which hecatomb
```

## Profiles for HPC clusters

Snakemake profiles are a must-have for running Snakemake pipelines on HPC clusters.
While they can be a pain to set up, you only need to do this once and then life is easy.
For more information, check the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), 
or our recent [blog post on Snakemake profiles](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated).

Hecatomb ships with an example profile for the Slurm workload manager in `hecatomb/snakemake/profile/example_slurm/`.
The example _profile_ `config.yaml` file (not to be confused with the _Hecatomb_ config file) contains all the Snakemake 
options for jobs that are submitted to the scheduler.
Hecatomb expects the following in the cluster command:

 - `resources.time` for time in minutes
 - `resources.mem_mb` for requested memory in Mb
 - `threads` for requested CPUs
 
We have tried to use what we believe is the most common nomenclature for these variables in Snakemake pipelines 
in the hopes that Hecatomb is compatible with existing Snakemake profiles.

We recommend redirecting STDERR and STDOUT messages to log files using the Snakemake variables `{rule}` and `{jobid}`, 
for instance like this `--output=logs/{rule}/{jobid}.out`.
You should also prepend the scheduler command with a command to make the log directories in case they don't exists
(it can cause errors for some schedulers), in this example like so: `mkdir -p logs/{rule}/ && sbatch ...`. 
This will make troubleshooting easier for jobs that fail due to scheduler issues.

The example profile includes a 'watcher' script.
Snakemake won't always pick up when a scheduler prematurely terminates a job, which is why we need a watcher script.
This line in the config file tells Snakemake how to check on the status of a job:
`cluster-status: ~/.config/snakemake/slurm/slurm-status.py` (be sure to check and update your file path).
The `slurm-status.py` script will query the scheduler with the jobid and report back to Snakemake on the job's status.

## Hecatomb configuration

The Hecatomb configuration file `hecatomb/snakemake/config/config.yaml` contains settings related to resources and 
cutoffs for various stages of the pipeline. The different config settings are outlined further on.
You can permanently change the behaviour of your Hecatomb installation by modifying the values in the config file.
Alternatively, you can change the behaviour of a single run by specifying the new values on the command line. 
See [below for passing your own Snakemake commands](advanced.md#additional-snakemake-commands).

## Database location

The databases are large (~55 GB) and if your Hecatomb installation is on a partition with limited on space,
you might want to specify a new location to house the database files. 
By default, this config setting is blank and the pipeline will use the install location (`hecatomb/databases/`).
You can specify the directory in the Hecatomb config file (`hecatomb/snakemake/config/config.yaml`) under `Databases: `, 
e.g:

```yaml
Databases: /scratch/HecatombDatabases
```

and run or rerun the installation 

```bash
hecatomb install
```

## Default resources

The Hecatomb config file contains some sensible defaults for resources.
While these should work for most datasets, they may fail for larger ones.
You may also have more CPUs etc at your disposal and want to minimise runtime of the pipeline.
Currently, the slowest steps are the MMSeqs searches and increasing the CPUs and RAM could significantly improve runtime.
The other settings (for Megahit and Minimap2, BBTools, and misc) will probably only show modest improvement.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# Memory for MMSeqs in megabytes (e.g 64GB = 64000, recommend >= 64000)
MMSeqsMem: 64000
# Threads for MMSeqs (recommend >= 16)
MMSeqsCPU: 32
# Max runtime for MMSeqs
MMSeqsTimeMin: 5760  # 4 days

# Memory for BBTools in megabytes (recommend >= 16000)
BBToolsMem: 16000
# CPUs for BBTools (recommend >= 8)
BBToolsCPU: 8

# Memory for Megahit/Flye in megabytes (recommend >= 32000)
MhitMem: 32000
# CPUs for Megahit/Flye in megabytes (recommend >= 16)
MhitCPU: 16

# Memory for slightly RAM-hungry jobs in megabytes (recommend >= 16000)
MiscMem: 16000
# CPUs for slightly RAM-hungry jobs (recommend >= 2)
MiscCPU: 2
```

## Filtering settings

There are many filtering etc. cutoff values that are specified in the Hecatomb config file.
For instance `READ_MINLENGTH: ` specifies the minimum allowed readlength after trimming.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# Preprocessing
QSCORE: 15 # Read quality trimming score (rule remove_low_quality in 01_preprocessing.smk)
READ_MINLENGTH: 90 # Minimum read length during QC steps (rule remove_low_quality in 01_preprocessing.smk)
CONTIG_MINLENGTH: 1000 # Read minimum length (rule contig_reformating_and_stats in 01_preprocessing.smk)
ENTROPY: 0.5 # Read minimum entropy (rule remove_low_quality in 01_preprocessing.smk)
ENTROPYWINDOW: 25 # entropy window for low qual read filter
CLUSTERID: 0.97 # Linclust read clustering percent identity (rule cluster_similar_sequences: in 01_preprocessing.smk)

# addHost settings
Entropy: 0.7
```

There are additional settings further down in the config file for users that are familiar with MMSeqs, 
as well as config settings that you should not alter.

## Additional Snakemake commands

As mentioned, Hecatomb is powered by Snakemake but runs via a launcher for your convenience.
The launcher--called with `hecatomb`--lets you specify the directory with your reads, host genome, where to save the results,
whether to do an assembly, and either specify the number of threads to use or a profile to use.
Snakemake itself has many command line options and the launcher can pass additional commands to Snakemake using the `--snake` option.

One such example is if you're not production ready you might wish to do a 'dry-run', where the run is simulated but no 
jobs are submitted, just to see if everything is configured correctly.
To do that, Snakemake needs the dry run flag (`--dry-run`, `--dryrun`, or `-n`).
In Hecatomb, you can pass this flag like so:

```bash
hecatomb run --reads fasq/ --profile slurm --snake=--dry-run
```

Hecatomb prints the Snakemake command to the terminal window and you should see these additional options added to the 
Snakemake command:

```bash
hecatomb run --test --snake=--dry-run
```
```text
Running Hecatomb
Running snakemake command:
snakemake -j 32 --use-conda --conda-frontend mamba --rerun-incomplete --printshellcmds \
  --nolock --conda-prefix /scratch/hecatomb/snakemake/workflow/conda \
  --dry-run -s /scratch/hecatomb/snakemake/workflow/Hecatomb.smk \
  -C Reads=/scratch/hecatomb/test_data Host=human Output=hecatomb_out Assembly=False
Building DAG of jobs...
```

You can use this to specify new config settings to overwrite Hecatomb's default config.
**NOTE: Wrap your Snake commands in quotes if you want to pass whitespace** (like in this example):

```bash
hecatomb run --test --snake="-C QSCORE=20 READ_MINLENGTH=100 ENTROPY=0.7"
```
```text 
Running Hecatomb
Running snakemake command:
snakemake -j 32 --use-conda --conda-frontend mamba --rerun-incomplete --printshellcmds \
  --nolock --conda-prefix /scratch/hecatomb/snakemake/workflow/conda \
  -C QSCORE=20 READ_MINLENGTH=100 ENTROPY=0.7 -s /scratch/hecatomb/snakemake/workflow/Hecatomb.smk \
  -C Reads=/scratch/hecatomb/test_data Host=human Output=hecatomb_out Assembly=False
Building DAG of jobs...
```

Have a look at the full list of available Snakemake options with `snakemake --help`.
Also note: The launcher will pass anything in `--snake=` verbatim, so use with care.
