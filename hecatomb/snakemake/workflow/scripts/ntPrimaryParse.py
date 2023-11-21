#!/usr/bin/env python3

import os
import logging
import atexit


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Reading hit seq ids from primary NT")
hit = set()
with open(snakemake.input.align, "r") as inAln:
    for line in inAln:
        l = line.strip().split("\t")
        hit.add(l[0])

logging.debug("Parsing AA-unclass seqs and writing NT-hit seqs")
with open(snakemake.output.class_seqs, "w") as outClass:
    with open(snakemake.output.unclass_seqs, "w") as outUnclass:
        with open(snakemake.input.seqs, "r") as inFa:
            for line in inFa:
                if line.startswith(">"):
                    id = line.strip().replace(">", "")
                    seq = inFa.readline().strip()
                    if id in hit:
                        outClass.write(f">{id}\n{seq}\n")
                    else:
                        outUnclass.write(f">{id}\n{seq}\n")
                else:
                    logging.error(
                        f"malformed {snakemake.input.seqs} file, or something, complain to Mike"
                    )
                    exit(1)

logging.debug("Done")
