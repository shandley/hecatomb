# Installing and testing Hecatomb

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
While the pipeline utilises a number of other programs and tools to run, you only need to install a couple of 
dependencies and Snakemake will take care of the rest.

Hecatomb relies on [conda](https://docs.conda.io/en/latest/) to ensure portability and ease of installation. 
The only other dependencies you will need to install are [Snakemake](https://snakemake.readthedocs.io/en/stable/#), 
[conda](https://docs.conda.io/en/latest/), and the python packages [pysam](https://pysam.readthedocs.io/en/latest/api.html)
and [plotly](https://plotly.com/python/).
We also highly recommend you use [mamba](https://github.com/mamba-org/mamba) as it is _a lot_ faster than the base conda.

After installing conda, install all the other dependencies like so:

```bash
conda install -c conda-forge -c bioconda snakemake mamba pysam plotly
```

To install Hecatomb, first download the GitHub repo:

```bash
git clone https://github.com/shandley/hecatomb.git
```

Then add or link the launcher to your PATH (e.g. if $HOME/bin is already on your PATH):

```bash
cd hecatomb/bin
ln -s $(pwd)/hecatomb $HOME/bin
```

That's it!

## Advanced customisation

If you're running Hecatomb on a cluster it is highly recommended to use Snakemake profiles 
([Running Hecatomb with a Snakemake profile](#)).
You may also want to customise the available resources for your system, whether you're using a cluster or running 
locally ([Advanced configuration](#)).

## Download the Hecatomb databases

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
The test dataset reads are in `test_data/` in the Hecatomb directory.

```bash
# run locally
hecatomb run --reads /path/to/hecatomb/test_data/

# run on cluster using a Snakemake profile
hecatomb run --reads /path/to/hecatomb/test_data/ --profile slurm
```

