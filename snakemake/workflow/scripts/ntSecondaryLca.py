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
logging.basicConfig(filename=snakemake.log[0], filemode='w', level=logging.DEBUG)

logging.debug('Reading alignments and extracting taxon IDs')
lin = {}
with open(snakemake.input.align, 'r') as alnfh:
    for line in alnfh:
        l = line.strip().split('\t')
        t = l[1].split('|')  # e.g. tid|2293023|NZ_QUCO01000049.1
        try:
            if not int(t[1]) in lin[l[0]]:
                lin[l[0]].append(int(t[1]))
        except KeyError:
            lin[l[0]] = [int(t[1])]

logging.debug('Joining tax ids for seqs and writing output')
with open(snakemake.output.lin, 'w') as out:
    for s in lin.keys():
        tOut = ';'.join([str(i) for i in sorted(lin[s])])
        out.write(f'{s}\t{tOut}\n')

logging.debug('Done')
