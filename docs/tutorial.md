# Hecatomb tutorial

This tutorial will walk through the process of running Hecatomb and performing some preliminary plots and analyses in both R and Python.

## Running Hecatomb

### System information

For this tutorial I'll be running Hecatomb on a 16-core/32-thread workstation with 64 GB of RAM running Ubuntu 18.04 LTS.
While this is fine for smaller datasets it is highly recommended using a HPC cluster or server with more CPUs and RAM for larger datasets.

### New install

```bash
# create new conda env and install hecatomb
conda create -n hecatomb -c conda-forge -c bioconda hecatomb

# activate env
conda activate hecatomb

# install the database files
hecatomb install
```

### Run Hecatomb

We will run hecatomb on the included test dataset, using the fast MMSeqs settings with 32 threads 
(which is the default anyway). This will give us an assembly and some read annotations.

```bash
Hecatomb run --test --threads 32 --fast
```

We should now have all the files we need!

TODO: finish tutorial

### Hecatomb run report

## Analysis and plotting in R and Python

### Dependencies

### Plotting and filtering

### Statistical tests

## Snakemake profiles

We'll walk through two ways to set up a profile for the Slurm workload manager.
If your HPC uses a different workload manager, the process of installing a profile will be similar but different.

### Copy an example profile

We have provided an example profile for the Slurm workload manager that should work for most HPCs using Slurm.
Snakemake will look for profiles in your home directory at:

```text
~/.config/snakemake/
```

First create a directory for your new profile, we'll call it 'slurm':

```bash
mkdir -p ~/.config/snakemake/slurm
```

Now copy the files for the example slurm profile 
(you can view them [here on GitHub](https://github.com/shandley/hecatomb/blob/main/snakemake/profile/example_slurm/)):

```bash
# go to your new profile directory
cd ~/.config/snakemake/slurm/
# copy the files (either from GitHub or from where you installed Hecatomb)
wget https://github.com/shandley/hecatomb/blob/ec1c62cfaddf29ade68cf4f33f4991fa07f9e6e0/snakemake/profile/example_slurm/config.yaml
wget https://github.com/shandley/hecatomb/blob/ec1c62cfaddf29ade68cf4f33f4991fa07f9e6e0/snakemake/profile/example_slurm/slurm-status.py
```

This example includes the necessary `config.yaml` file for profiles and a watcher script called `slurm-status.py`.
Make the watcher script executable:

```bash
chmod +x ~/.config/snakemake/slurm/slurm-status.py
```

Done!
You can now use this profile with hecatomb:

```bash
hecatomb run --test --fast --profile slurm
```

### Create a profile with cookiecutter

Cookiecutter is a nifty tool for creating projects using a template.
For Snakemake profiles, Cookiecutter takes away a lot of the manual configuration steps involved with setting up a profile for a specific scheduler.
There are currently [Cookiecutter templates for Slurm, SGE, PBS, and several other workload managers](https://github.com/Snakemake-Profiles),
and Hecatomb is intended to be compatible with these profiles.

We will walk through installing the [Slurm profile using Cookiecutter](https://github.com/Snakemake-Profiles/slurm).
To begin, create a new directory for your profile:

```bash
mkdir -p ~/.config/snakemake/slurm
```

Move to this directory and run the Cookiecutter command for this profile:

```bash
cd ~/.config/snakemake/slurm/
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```

Follow the prompts and you're done!
For our system, we did not need to specify anything, but you may need to specify account information for billing etc.
There is more detail on the [GitHub page for this profile](https://github.com/Snakemake-Profiles/slurm).

Use the profile with Hecatomb:

```bash
hecatomb run --test --fast --profile slurm
```
