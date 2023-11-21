#!/usr/bin/env python3

import os
import logging
import atexit


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Reading in seq IDs that were classified in primary AA search")
topHit = set()
with open(snakemake.input.alnsort, "r") as alnfh:
    for line in alnfh:
        l = line.strip().split()
        topHit.add(l[0])

logging.debug("Parsing seqtable and writing seqs that were classified")
with open(snakemake.output.class_seqs, "w") as outClass:
    with open(snakemake.input.seqs, "r") as inFa:
        for line in inFa:
            if line.startswith(">"):
                id = line.strip().replace(">", "")
                seq = inFa.readline().strip()
                if id in topHit:
                    outClass.write(f">{id}\n{seq}\n")
            else:
                sys.stderr.write(
                    f"malformed {snakemake.input.seqs} file, or something, complain to Mike."
                )
                exit(1)

logging.debug("Done")
