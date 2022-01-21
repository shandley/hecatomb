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

covStat = {}
logging.debug(f'Reading {snakemake.input.covstats}')
with open(snakemake.input.covstats, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            l = line.split('\t')
            covStat[l[0]] = [l[1] ,l[3] ,l[4] ,l[5] ,l[9]]

total = 0
logging.debug(f'Reading {snakemake.input.rpkm}')
with open(snakemake.input.rpkm, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            l = line.split('\t')
            rpk = int(l[2]) / (int(l[1]) / 1000)      # RPK
            total += rpk / 1000000         # size factor per million

logging.debug(f'Writing to {snakemake.output.count_tbl}')
with open(snakemake.output.count_tbl ,'w') as o:
    o.write('\t'.join([
        '#Sample',
        'Contig',
        'Length',
        'Reads',
        'RPKM',
        'FPKM',
        'CPM',
        'AverageFold',
        'ReferenceGC',
        'CoveragePercentage',
        'coverageBases',
        'MedianFold']
        ) + '\n')
    with open(snakemake.input.rpkm ,'r') as f:
        for line in f:
            if not line.startswith('#'):
                l = line.split('\t')
                rpk = int(l[2]) / (int(l[1]) / 1000)      # RPK
                try:
                    tpm = rpk / total               # TPM
                except ZeroDivisionError:
                    tpm = 0
                o.write( '\t'.join([
                    wildcards.sample,       # sample
                    l[0],                   # contig
                    l[1],                   # length
                    l[4],                   # reads
                    l[5],                   # RPKM
                    l[7],                   # FPKM
                    str(tpm) ]+              # CPM
                    covStat[l[0]]           # aveFole -> medianFold
                ) + '\n')