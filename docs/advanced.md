## Prefilter the bigtable

If your bigtable is too big, you can remove the stuff you don't want using the linux command line before loading it into R or Python.

**Remove low quality hits**

Use `awk` to apply evalue/bitscore/alignemnt lenght/etcetera cutoffs:

```bash
# copy the header for your new file
head -1 bigtable.tsv > newBigtable.tsv

# filter low quality hits, e.g. using an evalue cutoff of say 1e-20
awk -F '\t' '$7<1e-20' bigtable.tsv >> newBigtable.tsv
```

**Remove irrelevant hits**

If you only plan on analysing viruses or say protein hits, you can remove everything else with `awk`:

```bash
# copy the header
head -1 bigtable.tsv > newBigtable.tsv

# return only rows where "Viruses" is in the "kingdom" column
akw -F '\t' '$24=="Viruses"' bigtable.tsv >> newBigtable.tsv

# alternatively, only return the "aa" protein-based hits
awk -F '\t' '$5=="aa"' bigtable.tsv >> newBigtable.tsv
```

**Remove irrelevant columns**

If you've finished prefiltering, or you only plan on using say the evalues for filtering, 
you can remove all the other alignment metrics, this time using `cut`:

```bash
# remove all the alignment metrics except for evalue
cut --complement -f8-22 bigtable.tsv > newBigtable.tsv
```

## Retrieving sequences of interest

Upon analysing your bigtable, you may have a collection of interesting hits that you want to investigate further.
To retrieve the relevant sequences from your `seqtable.fasta` file you will only need a list of the sequence IDs--the first column in the bigtable.

For instance, get all seqIDs for the Viral family 'Flaviviridae' and save it to a file:

```R
# in R using tidyr/dplyr get the seq IDs
flaviviridaeSeqs = data %>% 
    filter(family=='Flaviviridae') %>% 
    pull(seqID)

# print to a file
lapply(flaviviridaeSeqs, write, "flavSeqIDs.list", append=TRUE, ncolumns=1)
```

Use Samtools to grab all these sequences from the seqtable.fasta file:

```shell
# in bash
samtools faidx seqtable.fasta -r flavSeqIDs.list > flaviviridaeSeqs.fasta
```
