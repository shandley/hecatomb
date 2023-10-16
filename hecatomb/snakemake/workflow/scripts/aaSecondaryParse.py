#!/usr/bin/env python3

import os
import logging
import atexit


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Reading seq ids from AA bigtable")
virSeqs = set()
with open(snakemake.input.bigtable, "r") as btbl:
    for line in btbl:
        l = line.strip().split("\t")
        virSeqs.add(l[0])

logging.debug("Parsing seqtable and writing non AA-hit seqs")
with open(snakemake.input.seqs, "r") as inFa:
    with open(snakemake.output.unclass_seqs, "w") as outFa:
        for line in inFa:
            if line.startswith(">"):
                id = line.strip().replace(">", "")
                seq = inFa.readline().strip()
                if not id in virSeqs:
                    outFa.write(f">{id}\n{seq}\n")
            else:
                logging.error(
                    f"malformed {snakemake.input.seqs} file, or something, complain to Mike"
                )
                exit(1)

logging.debug("Done")
