#!/usr/bin/env python3

import os
import atexit
import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

# the index position for each tax level in the bigtable.tsv - off by one because we capture the range as l[21:i]
tsvIndex = {
    "kingdom": 24,
    "phylum": 25,
    "class": 26,
    "order": 27,
    "family": 28,
    "genus": 29,
    "species": 30,
}
idxStart = 23
short = {"23": "k", "24": "p", "25": "c", "26": "o", "27": "f", "28": "g", "29": "s"}

logging.debug(f"Opening {snakemake.output[0]} for writing\n")
out = open(snakemake.output[0], "w")
out.write("sampleID\ttaxonLevel\ttaxonPath\ttaxonName\tcount\tpercent\n")
# re-read the file for each sample to keep the memory happy - this is probably not necessary
for sample in snakemake.params.samples:
    logging.debug(f"parsing {snakemake.input[0]} for sample {sample}\n")
    counts = {}  # counts[taxlevel][taxname] = int
    cpm = {}  # normalised counts, same structure as counts
    infh = open(snakemake.input[0], "r")
    infh.readline()  # skip header
    for line in infh:
        l = line.split("\t")
        if l[1] == sample:
            for t, i in tsvIndex.items():
                try:
                    if len(l[i].strip()) == 0:
                        continue
                except IndexError:
                    continue
                try:
                    counts[t]
                    cpm[t]
                except KeyError:
                    counts[t] = {}
                    cpm[t] = {}
                taxPath = []
                for o in range(idxStart, i):
                    taxPath.append(
                        f"{short[str(o)]}_{l[o]}"
                    )  # taxon path = k_kingName,p_phylName etc.
                outPath = ",".join(taxPath)
                try:
                    counts[t][outPath] += int(l[2])
                    cpm[t][outPath] += float(l[3])
                except KeyError:
                    counts[t][outPath] = int(l[2])
                    cpm[t][outPath] = float(l[3])
    infh.close()
    for taxLevel in counts.keys():
        for taxPath in counts[taxLevel].keys():
            taxName = taxPath.split("_")[-1]
            out.write(
                f"{sample}\t{taxLevel}\t{taxPath}\t{taxName}\t{counts[taxLevel][taxPath]}\t{cpm[taxLevel][taxPath]}\n"
            )
out.close()
