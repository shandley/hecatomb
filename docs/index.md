# Hecatomb

![](img/hecatombBanner.png)

A hecatomb is a great sacrifice or an extensive loss. 
Heactomb the software empowers an analyst to make data driven decisions to 'sacrifice' false-positive viral reads from 
metagenomes to enrich for true-positive viral reads. 
This process frequently results in a great loss of suspected viral sequences / contigs.

## Motivation

Hecatomb was developed in response to the challenges associated with the detection of viral sequences in metagenomes. 
Virus detection or virome profiling is typically performed on samples containing extensive host nucleic acid (e.g. a 
tissue biopsy) or nucleic acid from a variety of other organisms such as bacteria, fungi and archaea from soil samples 
or bacteria, fungi, archaea and food from mammalian stool samples. 
All of these non-viral nucleic acid types constitute the background and present a variety of issues that must be 
evaluated in order to confidently detect viral sequences.

Detection and quantitative analysis of viral sequences from mixed metagenomic data has other unique challenges. 
Viruses have a propensity to rapidly acquire mutations and evolve to be highly diverged from reference database sequences.
Coding and non-coding regions of viral genomes have disparate challenges (e.g. long-terminal 5' and 3' repeats).
Viruses are further clssified as either single or double stranded, and either DNA or RNA.
False-positives are rampant in many studies due to significant sequence homology between viruses and other domains of life.

Hecatomb utilises a unique search strategy to address these issues and enable researchers to accurately identify and quantify 
viruses in their metagenome samples.

## Pipeline overview

The pipeline consists of three main sections: Read preprocessing, Virus identification, and Assembly.

![](img/hecatombPipeline.png)

### Preprocessing

Heactomb is designed to perform rigorous quality control prior to assembly and taxonomic assignment. 
This rigor is justified primarily by the philosophy of [Garbage In, Garbage out](https://en.wikipedia.org/wiki/Garbage_in,_garbage_out). 
More specifically, the following issues are dealt with to ensure that only non-contaminat biological sequence is used for downstream analysis.

1. Non-biological sequence removal (primers, adapters)
2. Host sequence removal
3. Removal of redundant sequences (clustering)
	- Creation of sequence count table
	- Calculation of sequence properties (e.g. GC content, tetramer frequencies)

Hecatomb comes with a number of host genomes for use with this preprocessing, and bespoke host genomes can be added by the user.

### Virus identification

Hecatomb uses mmseqs2 to assign taxonomy to individual sequences. 
It does so through an iterative search approach, which is a large part of the heart of hecatomb.

The first search takes all query seqeunces (the seqtable from the preprocessing steps) and queries them against a virus protein target database. 
This initial search is a translated (nt-to-aa) search. 
The target database is all UniProt viral protein entires clusterd at 99% identity. 
This clustering greatly reduces the size of the database with minimal sacrifice to information content. 
For example, there were 4,635,541 protein entires for viruses in UniProt on October 29, 2020. 
When clustered at 99% identity only 1,840,337 representative sequences remained, this is about a 60% reduction in target query space.

All of the reads assigned a viral taxonomy in this first search need to be confirmed to be truly viral in origin. 
The issue here is that querying a target database consisting solely of viral proteins does not permit each sequence to see if it really originates from a different domain of life. 
So while a sequence may be classified as viral at this stage, once queried against a more comprehensive sequence database 
(consisting of bacteria, plants, vertebrates, fungi, etc.) those sequences may have higher statistically similar sequences in those domains of life.

It is formally possible to just do this secondary search, skipping the first step querying each sequence against a virus-only protien databse. 
However, this can be computationally very expensive. 
For example, if you have 1e6 reads as your query input (sequences in your seqtable) and you query these against the 
virus-only protien database with 1.8e6 entries that is over 25 orders of magnitude smaller than UniRef50 which has 49,410,134 entries when this was written. 
The primary search against a virus-only database serves as a 'net' to (relatively) quickly capture sequences of interest. 
In our experience, however, the secondary search against the larger database is still required as many of the sequences 
assigned a viral taxonomy when only querying the smaller virus-only sequence database will reveal themselves to be entirely something else upon this secondary search. 
Thus, the iterative approach will expedite sample processing by capturing only a small portion of the input sequences as 'potentially' viral and confirming them in the secondary search as 'likely' viral.

This process wherein a sequence is classified as viral in the first search against a virus-only database and then no 
longer classfied as viral after the second search against a more comprehensive database is where hecatomb gets it's namesake. 
Many investigators, particularly those lacking access to large computing infrastructure will rely on the results of the 
primary search which, due to the relatively small target database size, can be done on commodity hardware. 
Examining the results at this stage will regularly identify a large number of viruses. 
Basically anything with any level of sequence identity to a viral protein will be called viral even if the sequence 
originated from a bacteria, plant, fungi or other organism. 
These sequences are 'sacrificed' through the secondary search when many of the viral calls from the first search are lost.

### Assembly

The preprocessing rule also goes ahead and does assembly as contigs are an important prerequisite for many downstream analysis.
The assembly strategy consists of several steps intended to be resource efficient and to maximise representation of all present species.
First, individual [MEGAHIT](https://github.com/voutcn/megahit) assemblies are produced for each sample.
The reads are then mapped to the combined assemblies and any unmapped reads undergo another round of assembly.
All contigs are then combined and merged into a non-redundant set of contigs.

The assembly contigs are directly annotated with MMSeqs.
The assembly is also subject to a pseudo consensus annotation approach whereby the SeqTable sequences are mapped and their
Taxonomic assignments in the BigTable are combined with the read mapping information.
We find this useful with investigating contigs of interest.

## Snakemake

The Hecatomb pipeline is powered by [Snakemake](https://snakemake.readthedocs.io/).
It comes preconfigured and includes a launcher script to make running the pipeline simple.
We chose to use Snakemake over other workflow managers largely because of our extensive experience with it.
Snakemake does a lot of heavy lifting to make our lives easier, and the pipeline better.
It manages all the pipeline jobs, as well as generates benchmarks and reports for everything, allows the pipeline to 
be naturally reentrant, parallel, portable, robust, containerised, and many other buzzwords!

![](img/hecatombSnakemake.png)

## File outputs

The pipeline outputs a number of files for further analysis and exploration, as well as to provide an overview of the 
read preprocessing and distribution.

### Report

`report.html`

This file is generated by Snakemake and outlines a lot of information relating to the Hecatomb run.
Under the Results tabs are summary files for things like reads for each sample following different preprocessing steps
as well as some summary plots.

### SeqTable

`hecatomb_out/RESULTS/seqtable.fasta`

The SeqTable is the primary output of the read preprocessing and serves as the input for Taxonomic assignment.
It is composed of all the representative sequences from the clustered reads for all samples.
Samples are clustered individually, and the seq IDs for this fasta file follows the format `>sampleID:count:seqNumber`.
Here, `count` is the number of reads in that cluster which is important for statistical exploration.
Sequences are numbered sequentially (`seqNumber`) to ensure unique IDs.

### BigTable

`hecatomb_out/RESULTS/bigtable.tsv`

The BigTable is the main output of Taxonomic assignment and can be directly imported into R or Python.
The BigTable combines the seqtable IDs with their sampleID, counts, normalised counts, alignment information, taxonomic assignments and Baltimore classification.
This file is big, hence the name, but is designed to make merging with sample metadata, plotting, and statistical interrogation as easy as possible.

The header looks like this:

```
seqID  sampleID  count  normCount  alnType  targetID  evalue  pident  fident  nident  mismatches  qcov  tcov  qstart  qend  qlen  tstart  tend  tlen  alnlen  bits  targetName  taxMethod  kingdom  phylum  class  order  family  genus  species  baltimoreType  baltimoreGroup
```

### TaxonLevelCounts

`hecatomb_report/taxonLevelCounts.tsv`

This file is derived from the BigTable and summarises the total sequence counts, for each sample, at all taxonomic levels.
The TaxonLevelCounts combines the sampleID with the taxonomic level for which the counts refer, the full taxonomic path, 
the taxon name, and the total and normalised read counts.
The purpose of this file is to expedite statistical interrogation of your data.
For instance, if you wanted to compare the numbers of say Flaviviridae reads between two groups of samples, 
those counts have already been collected, and you can simply run your analysis and plotting on the relevant slice of the table.  

The file looks something like this:

```
sampleID    taxonLevel  taxonPath                                   taxonName       count   normCount
sample1     Kingdom     k_Bacteria                                  Bacteria        3162    3178.818
sample1     phylum      K_Viruses,p_Phixviricota                    Phixviricota    1216    1222.467
sample1     class       K_Viruses,p_Uroviricota,c_Caudoviricetes    Caudoviricetes  1234    1240.564
etc.
```

### Assembly

`hecatomb_out/RESULTS/assembly.fasta`

These are the contigs generated when running Hecatomb with the `--assembly` flag.
The assembly is used for producing the ContigSeqTable and ContigKrona plots, as well as the direct contig annotations.

### CONTIG ANNOTATIONS

TODO

### ContigSeqTable

`hecatomb_out/RESULTS/contigSeqTable.tsv`

The ContigSeqTable combines the read mapping information for the assembly with the read-based taxonomic assignments.
This file is intended to assist the user in identifying and binning assembly contigs by applying a consensus approach to contig taxonomic assignment.
The file includes the positional mapping information and can also enable investigation of more complex features such as 
chimeric contigs, recombination or horizontal transfer events.

The header looks like this:

```
contigID  seqID  start  stop  len  qual  count  normCount  alnType  taxMethod  kingdom  phylum  class  order  family  genus  species  baltimoreType  baltimoreGroup
```

### Sankey

`hecatomb_report/Sankey.svg`

The sankey diagram shows the fate of all reads throughout the preprocessing and read-based taxonomic assignment steps.
It serves to visualise the read filtering and distribution of taxonomic assignment methods and give you an overall impression of how well things went.
The sankey diagram produced for the test dataset is shown below. 
This dataset is relatively rich in viral sequences and yet the majority of reads are either filtered or non-viral (that we know of).

[![](img/Sankey.svg)](img/Sankey.svg)

### krona.html and contigKrona.html

`hecatomb_report/krona.html`

`hecatomb_report/contigKrona.html`

The Krona plots are to assist in visual exploration of the read annotations.
krona.html is derived from the bigtable and shows the raw distribution of taxon assignments.
contigKrona.html is derived from the contigSeqTable and includeds the taxon assignment method (either tophit or LCA).
The contigKrona plot helps to visualise the distributions of topHit versus LCA assigned reads as well as the 
distributions over contigs of the identified species, and the distribution of taxonomic assignments for each contig.
