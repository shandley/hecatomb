# Technical Notes and Design Decisions

In implementing `hecatomb` we had to make several design decisions. We keep track of some of those here, for clarity, history, and to explain the choices that we made.


# Databases

We decided to abstract out downloading the databases. This is very often a time-consuming task and we appreciate the unstable network connections make this tricky. We [have described](../databases) how to build the databases without using `snakemake`, but of course provide that as a pipeline too.

We separately remove the LINEs and SINEs from the sequences. The human genome we used is masked to remove those elements (and other viral elements), and so we have a separate screen for LINEs and SINEs. This puts those hits in a separate fastq file and so you can explore them if your virus has a SINE element (for example, the [BeAn 58058 poxvirus](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5334577/) has a SINE element.)

# Snakemake

We chose to use [snakemake](https://snakemake.readthedocs.io/) over other workflow managers largely because of our extensive experience with it. However, there is no reason you could not implement `hecatomb` in nextflow or any one of the other workflow managers.

We have tried to abstract out the configuration from the implementation. In particular, we have kept the input/output directories in the configuration file so that you can easily manipulate those. We strongly recommend that you copy the configuration file to each working directory, and then edit it in that directory (or not, as appropriate, if you like the directory structure we use!).



# MMSEQS

We use [mmseqs](https://github.com/soedinglab/MMseqs2) for a lot of the clustering and searching as in our hands it is at least as robust as [blast](https://blast.ncbi.nlm.nih.gov/), is faster, and includes some convenient taxonomy analysis scripts. 

There are a couple of issues that we have encountered with `mmseqs` that we have worked around.
- `mmseqs` creates one output file per core, and so you don't know, *a priori* how many or which output files are produced (which is required for snakemake rules!). However, per [this issue](https://github.com/soedinglab/MMseqs2/issues/292) you can use the `.dbtype` and `.index` files to track when `mmseqs` has completed.
- `mmseqs` requests a temporary directory for many of its computations, but then sets a symlink to it called "latest". The problem is if you fire off lots of `mmseqs` jobs concurrently they overwrite each other. We use the standard `bash` approach of creating temporary directories using `mktemp`. We start with a single location to hold them (saved as the variable `TMPDIR`), and then use this call in the `shell` to create a new temporary directory: `$(mktemp -d -p {TMPDIR})`. Each directory is unique. 


