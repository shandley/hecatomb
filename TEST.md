# Testing hecatomb

You should be able to change into the `hecatomb` directory and run the test directly. 

You will need:

1. [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
2. [snakemake](https://snakemake.readthedocs.io/)

````
cd hecatomb
snakemake --configfile configs/sample_config.json --snakefile snakemake/contaminant_removal.snakefile
```



