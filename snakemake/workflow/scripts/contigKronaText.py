#!/usr/bin/env python3

import os
import atexit
import logging

def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None

atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0], filemode='w', level=logging.DEBUG)

logging.debug('Slurping contig seq table')
counts = {}
for l in stream_tsv(snakemake.input[0]):
    if l[0] == "contigID":
        continue
    t = '\t'.join((l[9:] + [l[0]]))
    c = l[1].split(':')
    try:
        counts[t] += int(c[1])
    except KeyError:
        counts[t] = int(c[1])
logging.debug('Sorting and writing contig taxon info')
outFH = open(snakemake.output[0], 'w')
for k in sorted(counts.keys()):
    outFH.write(f'{counts[k]}\t{k}\n')
outFH.close()