#!/usr/bin/env python3

import os
import atexit
import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Concatenating ORFM ORFs from seqtable.fasta")


def printSeq(seq, fh):
    if seqComp["curID"]:
        fh.write(f'{seq["curID"]}{seq["curSeq"]}\n')


seqComp = {"curID": str(), "CurSeq": str()}
with open(snakemake.input[0], "r") as ifh:
    with open(snakemake.output[0], "w") as ofh:
        for line in ifh:
            if line.startswith(">"):
                if seqComp["curID"] != line:
                    printSeq(seqComp, ofh)
                    seqComp["curID"] = line
                    seqComp["curSeq"] = ifh.readline().strip()
                else:
                    seqComp["curSeq"] = 'XXXXX'.join([seqComp["curSeq"], ifh.readline().strip()])
            else:
                logging.debug(f"Not a seq ID, this shouldn't occur: {line}")
                exit(1)
