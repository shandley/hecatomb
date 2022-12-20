#!/usr/bin/env python3

import os
import logging
import atexit


def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None


atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Reading in seq IDs that were classified in primary AA search")
topHit = {}
with open(snakemake.input.alnsort, "r") as alnfh:
    for line in alnfh:
        l = line.strip().split()
        topHit[l[0]] = 1

logging.debug("Parsing seqtable and writing seqs that were classified")
with open(snakemake.output.class_seqs, "w") as outClass:
    with open(snakemake.input.seqs, "r") as inFa:
        for line in inFa:
            if line.startswith(">"):
                id = line.strip().replace(">", "")
                seq = inFa.readline().strip()
                try:
                    topHit[id]
                    outClass.write(f">{id}\n{seq}\n")
                except KeyError:
                    pass
            else:
                sys.stderr.write(
                    f"malformed {input.seqs} file, or something, complain to Mike."
                )
                exit(1)

logging.debug("Done")
