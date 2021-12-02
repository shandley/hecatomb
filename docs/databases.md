A complete list of all database files is available in the Hecatomb `config.yaml` file.

## Contaminants

`databases/contaminants/`

This database consists of a collection of NEBNext and TruSeq adapters, primerB sequences, 
and a collection of cloning vectors that are common sources of false positives. 
These sequences are used in preprocessing to trim reads of non-biological contaminant sequences.

## Host genomes

`databases/host/`

Contamination of viral metagenomes by host DNA can be a significant burden and source of false positive in viral annotation.
We use a host reference genome to filter out any host reads prior to annotation.
However, host genomes typically contain large amounts of virus-like and virus-derived sequence.
This could lead to erroneous removal of true viral reads.
We therefore must process host reference genomes to mask viral-like sequence prior to using it for preprocessing.
Hecatomb comes with several processed host genomes and a tool for users to add their own genomes.
See [Adding your own host genome](usage.md#adding-your-own-host-genome) for more info.

## Amino acid databases

`databases/aa/`

**Primary AA**

The primary AA database is used to quickly identify any reads that are viral-like.
This database is all UniProt viral protein entires clusterd at 99% identity. 
This clustering greatly reduces the size of the database with minimal sacrifice to information content. 
For example, there were 4,635,541 protein entires for viruses in UniProt on October 29, 2020. 
When clustered at 99% identity only 1,840,337 representative sequences remained, 
this is about a 60% reduction in target query space.

**Secondary AA**

The secondary AA database is used to check if any viral-like reads (from the primary AA search) have a better non-viral match.
TODO ...

## Nucleotide databases

`databases/nt/`

**Primary NT**

The primary NT database is used to quickly identify any viral-like reads that were not identified in the primary AA search.
TODO ...

**Secondary NT**

The secondary NT database is used to check if the viral-like reads (from the primary NT search) have a better non-viral match.
TODO ...

## Taxonomy

`databases/tax/`

Hecatomb utilises the NCBI taxonomy database (taxdump) with [TaxonKit](https://github.com/shenwei356/taxonkit) 
for converting taxon IDs into complete lineages for the output bigtable.

## Tables

`databases/tables/`

These are a collection of supplementary tables.
Primarily, the Baltimore classifications are provided here to be merged with the bigtable annotations.
