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

logging.debug("Reading hit seq ids from primary NT")
hit = {}
with open(snakemake.input.align, "r") as inAln:
    for line in inAln:
        l = line.strip().split("\t")
        hit[l[0]] = 1

logging.debug("Parsing AA-unclass seqs and writing NT-hit seqs")
with open(snakemake.output.class_seqs, "w") as outClass:
    with open(snakemake.output.unclass_seqs, "w") as outUnclass:
        with open(snakemake.input.seqs, "r") as inFa:
            for line in inFa:
                if line.startswith(">"):
                    id = line.strip().replace(">", "")
                    seq = inFa.readline().strip()
                    try:
                        hit[id]
                        outClass.write(f">{id}\n{seq}\n")
                    except KeyError:
                        outUnclass.write(f">{id}\n{seq}\n")
                else:
                    logging.ERROR(
                        f"malformed {input.seqs} file, or something, complain to Mike"
                    )
                    exit(1)

logging.debug("Done")
