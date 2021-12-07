This section assumes you have finished Tutorials [Part 1](tutorialPt1.md) and [Part 2](tutorialPt2.md)

## Unfiltered taxon counts

The `taxonLevelCounts.tsv` file is intended to make it quick and easy to compare virus hit counts across samples.
Let's take a look at it:

```R
View(taxonCounts)
```

[![](img/taxCountTable.png)](img/taxCountTable.png)

The file contains the counts and normalised counts for each taxonomic level, for each sample, 
based on the raw (unfiltered) hit counts.
If we wanted to compare the number of Adenoviridae hits between each sample we could pull out the Adenoviridae counts and plot them:

```R
# get adenoviridae counts
adenoCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Adenoviridae')

# plot
ggplot(adenoCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
```

[![](img/tuteAdenoBar.png)](img/tuteAdenoBar.png)

We can take this one step further and plot all the viral family normalised counts for each sample, in a stacked bar chart:

```R
# get all viral family counts
viralCounts = taxonCounts %>% 
    filter(taxonLevel=='family',grepl('k_Viruses',taxonPath))

# plot
ggplot(viralCounts) +
    geom_bar(aes(x=sampleID,y=count,fill=taxonName),position='stack',stat='identity') +
    coord_flip()
```

[![](img/tuteViralCounts.png)](img/tuteViralCounts.png)

## Generating taxon counts

The `taxonLevelCounts.tsv` file is convenient for comparing the raw counts,
but you will likely want to generate new counts from your _filtered_ hits.
Recreate the above plot from the filtered hits by first summing the counts
or normalisedCounts, e.g. at the family level:

```R
# Answer for "Challenge: Filter your raw viral hits to only keep protein hits with an evalue < 1e-10"
virusesFiltered = viruses %>% 
    filter(alnType=='aa',evalue<1e-10)

# collect the filtered counts
viralFiltCounts = virusesFiltered %>% 
    group_by(sampleID,family) %>% 
    summarise(n = sum(CPM))
```

Then plot again. 
This time we use `position='fill'` to make it look like 16s data, so we can confuse people.
We'll also add `theme_bw()` because grey is ugly:

```R
ggplot(viralFiltCounts) +
    geom_bar(aes(x=sampleID,y=n,fill=family),position='fill',stat='identity') +
    coord_flip() +
    theme_bw()
```

[![](img/tuteViralFiltCounts.png)](img/tuteViralFiltCounts.png)

These count tables we will use for plotting and some statistical comparisons.

# Challenge

**Make a stacked bar chart of the viral families for the Male and Female monkeys**

![](img/tuteGenderCounts.png)

# Visualising groups

We have a few viral families that are very prominent in our samples.
Let's see if there is a difference in viral loads according to our sample groups.
Collect sample counts for _Microviridae_.
Include the metadata group in `group_by()` so you can use it in the plot.

```R
microCounts = virusesFiltered %>% 
    group_by(family,sampleID,MacGuffinGroup) %>% 
    filter(family=='Microviridae') %>% 
    summarise(n = sum(CPM))
```

And plot. I like jitter plots but boxplots or violin plots might work better if you have hundreds of samples.

```R
ggplot(microCounts) +
    geom_jitter(aes(x=MacGuffinGroup,y=n),width = 0.1) +
    theme_bw()
```

![](img/tuteMicrovirJitter.png)

Let's do the same for _Podoviridae_.

```R
# collect counts
podoCounts = virusesFiltered %>% 
    group_by(family,sampleID,MacGuffinGroup) %>% 
    filter(family=='Podoviridae') %>% 
    summarise(n = sum(CPM))

# plot
ggplot(podoCounts) +
    geom_jitter(aes(x=MacGuffinGroup,y=n),width = 0.1) +
    theme_bw()
```

![](img/tutePodoJitter.png)

# Challenge

**Could gender be a good predictor of viral load for these families?**

While the MacGuffinGroup looks promising for _Podoviridae_, 
we'll need to [move on to Part 4: statistical tests](tutorialPt4.md) to find out for sure. 
