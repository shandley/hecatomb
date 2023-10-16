#!/usr/bin/env python3

import os
import pysam
import atexit
import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Slurping taxon info for seq IDs")
tax = {}
with open(snakemake.input.taxon, "r") as inTSV:
    inTSV.readline()  # skip header
    for line in inTSV:
        l = line.strip().split("\t")
        tax[l[0]] = "\t".join((l[2:5] + l[22:]))

smplCounts = {}  # Read in read counts for each sample for normalised counts
with open(snakemake.input.counts, "r") as countfh:
    for line in countfh:
        l = line.strip().split("\t")
        smplCounts[l[0]] = int(l[1])

logging.debug("Parsing mapped reads and pairing read taxon with mapped coords")
with open(snakemake.output[0], "w") as outFH:
    outFH.write("\t".join((snakemake.params.contigTaxonHeader)) + "\n")
    bam = pysam.AlignmentFile(snakemake.input.bam, "rb")
    for read in bam.fetch():
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        infOut = "\t".join(
            [
                str(read.reference_name),
                str(read.query_name),
                str(read.reference_start),
                str(read.reference_end),
                str(read.reference_length),
                str(read.mapping_quality),
            ]
        )
        try:
            taxOut = tax[read.query_name]
        except KeyError:
            c = read.query_name.split(":")
            taxOut = "\t".join(
                ([c[1], str((int(c[1]) / smplCounts[c[0]]) * 1000000)] + ["NA"] * 11)
            )
        outFH.write(f"{infOut}\t{taxOut}\n")
    bam.close()

logging.debug("Done")
