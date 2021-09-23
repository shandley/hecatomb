This tutorial will walk through the process of running Hecatomb and performing some preliminary plots and analyses in both R and Python.

## Run Hecatomb

### System information

For this tutorial I'll be running Hecatomb on a 16-core/32-thread workstation with 64 GB of RAM running Ubuntu 18.04 LTS.
While this is fine for smaller datasets it is highly recommended using a HPC cluster or server with more CPUs and RAM for larger datasets.

### New install

```shell
# create new conda env and install hecatomb
conda create -n hecatomb -c beardymcjohnface hecatomb

# activate env
conda activate hecatomb

# install the database files
hecatomb install
```

### Run Hecatomb

We will run hecatomb on the included test dataset.
Start by linking the test directory (which should be in the installation location).

```shell
ln -s `ls -l $(which hecatomb) | sed 's/.* //' | sed 's/hecatomb/..\/test_data/'` .
```

Now just run Hecatomb!
We'll using all 32 threads and include the `--assembly` flag to make sure Hecatomb performs an assembly as well.

```shell
Hecatomb run --reads test_data/ --threads 32 --assembly
```

### Hecatomb run report

## Analysis and plotting in R and Python

### Dependencies

### Plotting and filtering

### Statistical tests


