## Recommended customisation

If you're running Hecatomb on a HPC cluster, we absolutely recommend setting up a 
[Snakemake profile](profiles.md).

We also recommend reviewing the Snakemake `config.yaml` file in your Hecatomb installation directory.
The config file will be at `hecatomb/snakemake/config/config.yaml` and you can find your installation directory with:

```bash
which hecatomb
```

## Changing the Hecatomb configuration

The Hecatomb configuration file `hecatomb/snakemake/config/config.yaml` contains settings related to resources and 
cutoffs for various stages of the pipeline. The different config settings are outlined further on.
You can permanently change the behaviour of your Hecatomb installation by modifying the values in this config file.
Alternatively, you can specify your own config file. To do this, you first copy the system default config file like so:

```bash
hecatomb config
```

You can then edit your new `hecatomb.config.yaml` file to suit your needs and used it in a Hecatomb run like so:

```bash
hecatomb run --configfile hecatomb.config.yaml
```

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
Currently, the slowest steps are the MMSeqs searches; increasing the CPUs and RAM could significantly improve runtime.
The other settings (for Megahit and Minimap2, BBTools, and misc) will probably only show modest improvement.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# Memory for MMSeqs in megabytes (e.g 64GB = 64000, recommend >= 64000)
MMSeqsMem: 64000
# Threads for MMSeqs (recommend >= 16)
MMSeqsCPU: 32
# Max runtime in minutes for MMSeqs (this is only enforced by the Snakemake profile)
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

# Default memory in megabytes (for use with --profile)
defaultMem: 2000
# Default time in minutes (for use with --profile)
defaultTime: 1440
# Default concurrent jobs (for use with --profile)
defaultJobs: 100
```

## Preprocessing settings

There are many filtering etc. cutoff values that are specified in the Hecatomb config file.
For instance `READ_MINLENGTH: ` specifies the minimum allowed read length after trimming.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# Preprocessing
QSCORE: 15 # Read quality trimming score (rule remove_low_quality in 01_preprocessing.smk)
READ_MINLENGTH: 90 # Minimum read length during QC steps (rule remove_low_quality in 01_preprocessing.smk)
CONTIG_MINLENGTH: 1000 # Read minimum length (rule contig_reformating_and_stats in 01_preprocessing.smk)
ENTROPY: 0.5 # Read minimum entropy (rule remove_low_quality in 01_preprocessing.smk)
ENTROPYWINDOW: 25 # entropy window for low qual read filter

# CLUSTER READS TO SEQTABLE (MMSEQS EASY-LINCLUST)
 # -c = req coverage of target seq
 # --min-seq-id = req identity [0-1] of alignment
linclustParams:
 --kmer-per-seq-scale 0.3
 -c 0.7
 --cov-mode 1
 --min-seq-id 0.95
 --alignment-mode 3
```

There are additional settings further down in the config file for users that are familiar with MMSeqs, 
as well as some settings that you should not alter.

## Alignment filtering

Hecatomb has settings for filtering MMSeqs alignments at each stage of the search strategy.
By default, we use a lenient e-value cutoff to maximise the identification of viral sequences in the primary searches,
and a more stringent e-value cutoff for the multi-kingdom search.
You can lower the evalue cutoffs (`-e`) to improve runtime performance at the cost of lower recall.
The `--min-lenghth` should be the same or lower than the preprocessing cutoffs.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# ALIGNMENT FILTERING CUTOFFS
  # --min-length for AA should be equal or less than 1/3 of READ_MINLENGTH
  # --min-length for NT should be equal or less than READ_MINLENGTH
filtAAprimary:
 --min-length 30
 -e 1e-3
filtAAsecondary:
 --min-length 30
 -e 1e-5
filtNTprimary:
 --min-length 90
 -e 1e-3
filtNTsecondary:
 --min-length 90
 -e 1e-5
```

## Alignment settings

Hecatomb can perform MMSeqs alignments using either sensitive (default) or fast (`--fast`) parameters.
You can tweak the setting in the config file but you should consult the MMSeqs documentation before making any changes.

The relevant section in `hecatomb/snakemake/config/config.yaml` is shown below:

```yaml
# PERFORMANCE SETTINGS - SEE MMSEQS DOCUMENTATION FOR DETAILS
# sensitive AA search
perfAA:
 --start-sens 1
 --sens-steps 3
 -s 7
 --lca-mode 2
 --shuffle 0
# fast AA search
perfAAfast:
 -s 4.0
 --lca-mode 1
 --shuffle 0
# sensitive NT search
perfNT:
 --start-sens 2
 -s 7
 --sens-steps 3
# fast NT search
perfNTfast:
 -s 4.0
```

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

Hecatomb prints the Snakemake command to the terminal window before running and you should see these additional options 
added to the Snakemake command:

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

Have a look at the full list of available Snakemake options with `snakemake --help`. 
The launcher will pass anything in `--snake=` verbatim, so use with care.

**NOTE: Don't use hecatomb's --snake arg to pass additional config settings with Snakemake's -C/--config arg**,
instead, follow the directions above [for changing the config of a Hecatomb run](configuration.md#changing-the-hecatomb-configuration).
