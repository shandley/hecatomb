# Assembly modules

__paired_end.smk__

Generic module for individual sample assemblies with Megahit for paired end short reads.

__single_end.smk__

Generic module for individual sample assemblies with Megahit for single end short reads.

__longreads.smk__

Generic module for individual sample assemblies with Canu for longreads.

__combine_sample_assemblies.smk__

Rules to combine all individual sample assemblies into a cross assembly.
Uses FlyE with `--subassemblies`.

# Create your own assembly module

Fork this repo and add your own assembly module to support a different library type.
You'll need to update the CLI, the config file, and add the rules.
Ideally you should add rules for both a cross-assembly and a merged-assembly.

**hecatomb/__main__.py**

Assembly module is currently inferred from `--library`.
Add a new library option in `hecatomb/__main__.py`.

**hecatomb/snakemake/config/immutable.yaml**

Update `config/immutable.yaml` and add the appropriate rules for your new module.
Reuse existing files where possible to make life easier. 
Lasty, create a new file for your assembly module.

**hecatomb/snakemake/workflow/rules/assembly/your_new_module.smk**

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

### Outputs: cross-assembly

Currently, each library type has the option to do a __cross-assembly__ or a __merged-assembly__.
For the merged assembly, the output is the concatenated contigs from all individual sample assemblies, e.g.:

```python
    output:
        os.path.join(dir.out.assembly, "all_sample_contigs.fasta.gz")
```

Make sure all contig names are unique.
The `combine_sample_assembly.smk` rules will generate the merged assembly and assembly graph from this file.

### Outputs: cross-assembly

This will be a one-step assembly for all the samples.
The outputs will the final assembly fasta and the assembly graph.

```python
    output:
        assembly = os.path.join(dir.out.results, "cross_assembly.fasta"),
        graph = os.path.join(dir.out.results, "cross_assembly_graph.gfa"),
```

### Ouputs: coverage calculations

You will need to provide a rule for doing the coverage calculations from the final assembly
`os.path.join(dir.out.results, f"{config.args.assembly}_assembly.fasta")`.
Just copy-paste the `rule coverage_calculations` from one of the other files as a template.
You might not even have to change anything.

### Pull request

Do you think other people will be interested in your new module?
Create a pull request!
