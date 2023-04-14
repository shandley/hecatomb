# Preprocessing modules

__paired_end.smk__

Generic module for preprocessing paired end short reads.

__single_end.smk__

Generic module for preprocessing single end short reads.

__longreads.smk__

Module for processing long reads such as PacBio, Nanopore.

__roundAB.smk__

Module for processing paired end short reads from Round A/B library
prep for viral metagenomics.

__cluster_seqs.smk__

Rules to cluster processed R1 reads--for use with all preprocessing modules.

# Create your own module

Fork this repo and add your own preprocessing module to support a different library type.
You'll need to update the CLI, the config file, and add the rules.

**hecatomb/__main__.py**

Add the new module option to `--preprocess` in `hecatomb/__main__.py`.

**hecatomb/snakemake/config/immutable.yaml**

Update `config/immutable.yaml` and add the appropriate rules for your new module.
Reuse existing files where possible to make life easier.
Lastly, add the outputs (from below) used for assembly as targets by updating 
`config/immutable.yaml`. e.g.:

```python
modules:
    your_new_module:
        targets:
            ["_R1.unmapped.fastq.gz", "_R1.singletons.fastq.gz", "_R1.all.fastq.gz", 
             "_R2.singletons.fastq.gz", "_R2.unmapped.fastq.gz", "_R2.all.fastq.gz"]
        ...
```

Now, create the new file for your preprocessing module.

**hecatomb/snakemake/workflow/rules/preprocessing/your_new_module.smk**

### Inputs

R1 (and R2 for paired) reads accessed from `samples.reads["sample"]["R1/R2"]` e.g.:
```python
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]['R1'],
        r2=lambda wildcards: samples.reads[wildcards.sample]['R2'],
```

### Outputs

__For read annotation__

Reads for clustering via `cluster_seqs.smk`.
For paired reads, these should be the concatenated R1 paired and singleton reads e.g.:

```python
    input:
        os.path.join(dir.out.temp,"p04","{sample}_R1.all.fastq")
```

__For assembly__

QC'd and host removed reads. e.g.:

```python
# paired end
    input:
        r1 = os.path.join(dir.out.assembly, "{sample}_R1.unmapped.fastq.gz"),
        r2 = os.path.join(dir.out.assembly, "{sample}_R2.unmapped.fastq.gz"),
        r1s = os.path.join(dir.out.assembly, "{sample}_R1.singletons.fastq.gz"),
        r2s = os.path.join(dir.out.assembly, "{sample}_R2.singletons.fastq.gz")
```

or

```python
# single end / longreads
    input:
        os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz")
```

### Pull request

Do you think other people will be interested in your new module?
Create a pull request!