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
                    lin = set()
                t = l[1].split("|")  # e.g. tid|2293023|NZ_QUCO01000049.1
                if not l[0] in topDone:
                    outTop.write(f"{l[0]}\t{t[1]}\n")
                    topDone.add(l[0])
                if not t[1] in lin:
                    lin.add(t[1])
        dumpRead(lin, outLin)

logging.debug("Done")
