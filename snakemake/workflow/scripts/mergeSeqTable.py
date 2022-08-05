#!/usr/bin/env python3

import os
import logging
import atexit

def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None

def streamCountTable(file):
    """Parse the seq count table"""
    with open(file, 'r') as counts:
        line = counts.readline() # skip header
        for line in counts:
            l = line.split()
            if len(l) == 2:
                yield l
            else:
                logging.warning(f'Possible malformed sample seqtable {file}')

atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0], filemode='w', level=logging.DEBUG)
logging.debug('Reading in individual sample seqtables')
outFa = open(snakemake.output.fa, 'w')
outTsv = open(snakemake.output.tsv, 'w')
for sample in snakemake.params.samples:
    logging.debug(f'sample {sample}')
    seqId = 0
    seqCounts = 0
    st = os.path.join(snakemake.params.tmpdir, f"{sample}_R1.seqtable")
    for l in streamCountTable(st):
        seqCounts += int(l[1])
    outTsv.write(f'{sample}\t{seqCounts}\n')
    for l in streamCountTable(st):
        percCount = (int(l[1]) / seqCounts) * 100
        id = ':'.join((sample, l[1], str(percCount), str(seqId)))  # fasta header = >sample:count:percentCount:seqId
        seqId = seqId + 1
        outFa.write(f'>{id}\n{l[0]}\n')
outFa.close()
outTsv.close()
logging.debug('Done')
