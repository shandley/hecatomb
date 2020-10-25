
# Installing hecatomb


Hecatomb has been written using [snakemake](https://snakemake.readthedocs.io/en/stable/) and [conda](https://docs.conda.io/en/latest/) to ensure portability and ease of installation.

Quick version (TL;DR):

```bash
conda install -c bioconda -c conda-forge snakemake
git clone https://github.com/shandley/hecatomb.git
cd hecatomb
snakemake --configfile snakemake/config/sample_config.yaml -s snakemake/workflow/download_databases.smk --cores 4 --use-conda
snakemake --configfile snakemake/config/sample_config.yaml -s snakemake/workflow/Snakefile --cores 4 --use-conda
```


Our goal is to ensure that you only need to install conda and snakemake, and then everything else should be installed for you. For a discussion about this, please see our [Technical Notes](snakemake/TECHNICAL_NOTES.md). If you find that you need to install something to make hecatomb work please [post an issue](https://github.com/shandley/hecatomb/issues)

# Detailed Installation Instructions.

You only need to install `conda` and `snakemake` and `hecatomb` will do the rest for you. No more software installs!

## Step 1. Install `conda`. 

This is undoubtedly the hardest part, and it is quite easy. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) following their instructions.

## Step 2. Install `snakemake` and `mamba`

We highly recommend you install `mamba` as it is a _lot_ faster than conda. But you don't have to (if you don't, skip the mamba step below)

```
conda install -c conda-forge -c bioconda snakemake mamba
```

### Step 3. Clone this repo

Simply clone this repo anywhere you want to run hecatomb!

```
git clone https://github.com/shandley/hecatomb.git
```

### Step 4. Download the databases

As with any metagenome analysis tool, hecatomb requires a few databases, and we have created a database download snakefile that does this for you. 

You only need to download the databases once, and then you are done. Depending on your internet connection it may take a few hours to do this step.

You may choose to edit the database location in the config file. By default it is `databases` and will be relative to whereever you run the `snakemake` command from:

```bash
snakemake --configfile snakemake/config/sample_config.yaml -s snakemake/workflow/download_databases.smk --cores 4 --use-conda
```

Of course, change the number of cores to use as appropriate. 

> _Hint:_ This command will run on a single node or computer, but if you are using a cluster, see below to set up a general snakemake configuration to allow you to run this on a cluster.

The database download runs in parallel, but the final database is about 72 GB, and so it will take a while to download and process all the data. I recommend runing that overnight!

> _Tip:_ If the download files, we recommend running the `snakemake` command again. Sometimes the download fails mid-stream, but snakemake removes potentially corrupted files and starts the process again.


### Step 5. Run the code

You can either use our test data or edit the config file and point it to your `fastq` files. Change this line:

```yaml
Reads: test_data/fastq
```

and replace `test_data/fastq` with a directory to your `fastq` files.

Then you can run the code with:

```shell script
snakemake --configfile snakemake/config/sample_config.yaml -s snakemake/workflow/Snakefile --cores 4 --use-conda 
```

(Again, change the number of cores as appropriate.)

> _Tip:_ Currently, `hecatomb` expects an R1 file and an R2. Because of a few limitations, it is expected that all the files have the same ending after `_R1` and `_R2`. This is very occassionaly a problem, for example if you have something like `sample1_R1_001.fastq.gz` and `sample2_R1_002.fastq.gz`. In that case, rename the `002` file to `001`. We are working on a general solution to this problem.

# Using hecatomb on a cluster

You can run `hecatomb` on a single machine but it excels when you run it across a cluster. We recommend using `slurm` or `sge`

> _Tip_: `SGE` and `slurm` are two similar cluster queuing systems, and your cluster likely has one or the other. If you routinely use `qsub` to submit jobs, your cluster is probably using `SGE`. If you routinely use `sbatch` to submit jobs, your cluster is probably using `slurm`. Not sure? Ask your systems administrator for advice!

## Running on SGE

If you are using `SGE` or a variant,  you can run `hecatomb` using the snakemake cluster command:

```bash
mkdir sge_out sge_err
snakemake -s snakemake/workflow/Snakefile --configfile config/sample_config.yaml --cluster 'qsub -cwd -o sge_out -e sge_err' --local-cores 6 --cores 600 --latency-wait 60  --default-resources cpus=1 mem_mb=2000 --use-conda --conda-frontend mamba
```

Note that here we have added a couple of additional options that you should vary depending on your cluster configuration:

- `-cwd` means using the current working directory
- `-o sge_out` and `-e sge_err` write the output and error log files to `sge_out` and `sge_err` respectively. Because we make these as directories individual output files are written.
- `--local-cores` is how many cores there are on the master node of the cluster
- `--cores` is how many compute cores there are on the cluster
- `--latency-wait` allows time for the files to be transfered once the job has finished computing.
- `--default-resources` specifies the default memory request and cpu request if it is not overridden by a job
- `--use-conda` is to ensure we use conda for all software installation
- `--conda-frontend mamba` uses the mamba replacement for conda. This is optional, but recommended. 

## Running on Slurm

You can run on slurm with a similar command, just changing the `--cluster` part:

```bash
mkdir logs_slurm
snakemake -s snakemake/workflow/Snakefile --configfile snakemake/config/sample_config.yaml --cluster 'sbatch  --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{jobid} _{jobid}.out -e logs_slurm/{rule}_{jobid}.err' --local-cores 32 --cores 600 --latency-wait 60 --default-resources cpus=1 mem_mb=2000 --use-conda --conda-frontend mamba
```

This puts the output and error files in the directory `logs_slurm`. On my cluster, nothing runs if I forget to make the logs_slurm output directories, though!


### Advanced. Set up your snakemake environment

We recommend using profiles (see these [great](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/) and [great](http://bluegenes.github.io/Using-Snakemake_Profiles/) blogs for more information). 

You can make a profile for a computer with no cluster, sge, or slurm.

In any case, we are going to make a directory and file in your home directory. That way, whenever you run `snakemake` it can access it.


#### SGE profile

For `SGE` we will make a file called `config.yaml` in the sge directory:

```bash
mkdir -p ~/.config/snakemake/sge/
vi ~/.config/snakemake/sge/config.yaml
```

In that file put these contents:

```yaml
# non-slurm settings
jobs: 600
use-conda: True
conda-frontend: mamba
default-resources: [cpus=1, mem_mb=2000]
keep-going: True

# sge settings
cluster: "qsub -cwd -o sge_out -e sge_err -pe smp {resources.cpus} -V "
latency-wait: 60
local-cores: 6
```

> _Tip_: Your cluster may have slightly different parallel environments. To figure out what is avaialble, try the command `qconf -spl` or ask your systems administrator. You may need to change the `smp` part in this file.

> _Tip_: Adjust the `local-cores` and `cores` to reflect your cluster.

Once you have set this up, you can run the same snakemake command as before a lot simpler:

```bash
mkdir sge_out sge_err
snakemake -s snakemake/workflow/Snakefile --configfile config/sample_config.yaml --profile sge
```

#### Slurm profile

For `slurm` we will make a file called `config.yaml` in the slurm directory:

```bash
mkdir -p ~/.config/snakemake/slurm/
vi ~/.config/snakemake/slurm/config.yaml
```

Put this in that file:

```yaml
# non-slurm settings

jobs: 10
use-conda: True
conda-frontend: mamba
default-resources: [cpus=1, mem_mb=2000]
keep-going: True


# slurm settings

cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{jobid}.out -e logs_slurm/{rule}_{jobid}.err "
latency-wait: 60
local-cores: 32
```

Once you have set this up, you can run the same snakemake command as before a lot simpler:

```bash
mkdir logs_slurm
snakemake -s snakemake/workflow/Snakefile --configfile config/sample_config.yaml --profile slurm
```

## The config file

You should make a copy of the config file. We typically make a copy of that file into each directory where we are working. Then if you make any changes to that file they reside with the data. 
There are [example config files](configs/) in both [JSON](configs/sample_config.json) and [YAML](configs/sample_config.yaml), and of course snakemake can use either. (If you are not sure, YAML is probably easier to start with than JSON).

The key things in the config file are:

1. The database file location. You can set that in the config file and then create the database as described below
2. The directory name where your raw reads (`fastq files`) reside. 

You can adjust almost everything else as needed, and the directories will be created when the data is generated.




