This is the current conda build recipe for Hecatomb, currently hosted at 
[https://anaconda.org/beardymcjohnface](https://anaconda.org/beardymcjohnface)

```bash
# from hecatomb/build/
conda build hecatomb/
```

## Dependencies

We make sure we're using Python 3 and we obviously need Snakemake.
We otherwise only need to specify mamba (for Snakemake to build conda envs for each rule), 
and Pysam and Plotly for a few Snakemake 'Run' rules (as you can't specify envs for these).

## Build

The dependencies are installed, and the repo is downloaded and unpacked to the conda env.
We copy the contents of `bin/` and `snakemake/` which are all that is required to run the pipeline.
The `test_data/` files are copied as they are required for running the Hecatomb tests.
We also copy `docs/` and `mkdocs.yaml` if the user wants to build a local copy of the documentation 
(they would need to separately install mkdocs).

## Requirements

Use this to build the conda env that would be generated from a bioconda install. 
Used for checking dependency versions etc for bioconda installation.
