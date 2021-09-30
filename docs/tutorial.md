This tutorial will walk through the process of running Hecatomb and performing some preliminary plots and analyses in both R and Python.

## Run Hecatomb

### System information

For this tutorial I'll be running Hecatomb on a 16-core/32-thread workstation with 64 GB of RAM running Ubuntu 18.04 LTS.
While this is fine for smaller datasets it is highly recommended using a HPC cluster or server with more CPUs and RAM for larger datasets.

### New install

```shell
# create new conda env and install hecatomb
conda create -n hecatomb -c conda-forge -c bioconda -c beardymcjohnface hecatomb

# activate env
conda activate hecatomb

# install the database files
hecatomb install
```

### Run Hecatomb

We will run hecatomb on the included test dataset, using the fast MMSeqs settings.
We'll use all 32 threads (which is the default anyway) and include the `--assembly` flag to make sure Hecatomb performs an assembly as well.

```shell
Hecatomb run --test --threads 32 --assembly --fast
```

We should now have all the files we need!

TODO: finish tutorial

### Hecatomb run report

## Analysis and plotting in R and Python

### Dependencies

### Plotting and filtering

### Statistical tests


