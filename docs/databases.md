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

## Databases

`databases/aa/` and `databases/nt/`
These are the amino acid (AA) and nucleotide (NT) databases used for sequence annotation of both reads and contigs. 
For each of the AA and NT databases, there is a primary viral database used for identifying reads that match a known virus, 
and a secondary multi-kingdom database which is used for assigning taxonomy to either reads and contigs. 

The primary AA database includes all UniProt viral protein entries clustered at 99% identity. 
The secondary AA database consists of the Uniclust50 database [Mirdita et al. 2017](https://doi.org/10.1093/nar/gkw1081) 
supplemented with the primary AA database. 
The primary NT database consists of all viral sequences in GenBank clustered at 100% identity to remove redundancy. 
The secondary NT database consists of a customised polymicrobial nucleotide database containing representative RefSeq 
genomes from Bacteria (n = 14,933), Archaea (n = 511), Fungi (n = 423), Protozoa (n = 90) and plant (n = 145) genomes. 

## Taxonomy

`databases/tax/`

Hecatomb utilises the NCBI taxonomy database (taxdump) with [TaxonKit](https://github.com/shenwei356/taxonkit) 
for converting taxon IDs into complete lineages for the output bigtable.

## Tables

`databases/tables/`

These are a collection of supplementary tables.
Primarily, the Baltimore classifications are provided here to be merged with the bigtable annotations.
