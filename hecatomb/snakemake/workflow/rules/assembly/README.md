# Assembly modules

## paired_end.smk

Generic module for individual sample assemblies with Megahit for paired end short reads.

## single_end.smk

Generic module for individual sample assemblies with Megahit for single end short reads.

## longreads.smk

Generic module for individual sample assemblies with Canu for longreads.

## combine_sample_assemblies.smk

Rules to combine all individual sample assemblies into a cross assembly.
Uses FlyE with --subassemblies.

## Create your own sample assembly

### CLI

Assembly method is currently inferred from preprocessing.
Add a new preprocessing option in `hecatomb/__main__.py`.

### Config

Update `config/immutable.yaml` and add the appropriate rules for your new module.
Reuse existing files where possible to make life easier.

### Inputs

QC'd and host removed paired or single reads. e.g.:

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
# single end / longread
    input:
        os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz")
```

### Outputs

Concatenated contigs from all individual sample assemblies, e.g.:

```python
    output:
        os.path.join(dir.out.assembly, "all_sample_contigs.fasta.gz")
```

Make sure all contig names are unique.
The `combine_sample_assembly.smk` rules will generate the cross assembly from this file.
You will need to provide a rule for doing the coverage calculations from `os.path.join(dir.out.assembly, "all_sample_contigs.fasta")`.
Just copy-paste the `rule coverage_calculations` from one of the other files as a template.

### Pull request

Do you think other people will be interested in your new module?
Create a pull request!