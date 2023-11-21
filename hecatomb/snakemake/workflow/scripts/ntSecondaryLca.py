#!/usr/bin/env python3

import os
import logging
import atexit


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

lin = set()
curRead = ""


def dumpRead(lst, fh):
    tOut = ";".join(list(lst))
    fh.write(f"{curRead}\t{tOut}\n")


logging.debug("Reading alignments and extracting taxon IDs")
with open(snakemake.output.top, "w") as outTop:
    topDone = set()
    with open(snakemake.output.lin, "w") as outLin:
        with open(snakemake.input.align, "r") as alnfh:
            for line in alnfh:
                l = line.strip().split("\t")
                if curRead != l[0]:
                    if curRead != "":
                        dumpRead(lin, outLin)
                    curRead = l[0]
                if not l[0] in topDone:
                    outTop.write(f"{l[0]}\t{l[1]}\n")
                    topDone.add(l[0])
                lin.add(l[1])
        dumpRead(lin, outLin)

logging.debug("Done")
