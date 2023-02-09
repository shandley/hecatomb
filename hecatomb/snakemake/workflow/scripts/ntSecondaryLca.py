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

lin = []
curRead = ""


def dumpRead(lst, fh):
    tOut = ";".join([str(i) for i in sorted(lst)])
    fh.write(f"{curRead}\t{tOut}\n")


logging.debug("Reading alignments and extracting taxon IDs")
with open(snakemake.output.lin, "w") as out:
    with open(snakemake.input.align, "r") as alnfh:
        for line in alnfh:
            l = line.strip().split("\t")
            if curRead != l[0]:
                if curRead != "":
                    dumpRead(lin, out)
                curRead = l[0]
                lin = []
            t = l[1].split("|")  # e.g. tid|2293023|NZ_QUCO01000049.1
            if not t[1] in lin:
                lin.append(int(t[1]))
    dumpRead(lin, out)

logging.debug("Done")
