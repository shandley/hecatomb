#!/usr/bin/env python3

import os
import atexit
import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Slurping contig seq table")
counts = {}
with open(snakemake.input[0], "r") as fh:
    for line in fh:
        l = line.strip().split("\t")
        if l[0] == "contigID":
            continue
        t = "\t".join((l[9:] + [l[0]]))
        c = l[1].split(":")
        try:
            counts[t] += int(c[1])
        except KeyError:
            counts[t] = int(c[1])
logging.debug("Sorting and writing contig taxon info")
outFH = open(snakemake.output[0], "w")
for k in sorted(counts.keys()):
    outFH.write(f"{counts[k]}\t{k}\n")
outFH.close()
