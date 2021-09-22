# About Hecatomb

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
more info on the assembly strategy...

## File outputs

### seqtable.fasta

### bigtable.tsv

### assembly.fasta

### contigSeqTable.tsv

### krona.html and contigKrona.html


