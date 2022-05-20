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

logging.debug('Slurping Tax assignments from bigtable')
counts = {}
for l in stream_tsv(snakemake.input[0]):
    if l[0] == "seqID":
        continue
    t = '\t'.join(l[23:30])
    try:
        counts[t] += int(l[2])
    except KeyError:
        counts[t] = int(l[2])
logging.debug('Sorting, counting, and writing tax assignments')
outFH = open(snakemake.output[0], 'w')
for k in sorted(counts.keys()):
    outFH.write(f'{counts[k]}\t{k}\n')
outFH.close()