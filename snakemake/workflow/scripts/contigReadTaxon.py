#!/usr/bin/env python3

import os
import pysam
import atexit
import logging

def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None

atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0], filemode='w', level=logging.DEBUG)

logging.debug('Slurping taxon info for seq IDs')
tax = {}
with open(snakemake.input.taxon, 'r') as inTSV:
    inTSV.readline() # skip header
    for line in inTSV:
        l = line.strip().split('\t')
        tax[l[0]] = '\t'.join((l[2:5] + l[22:]))

logging.debug('Parsing mapped reads and pairing read taxon with mapped coords')
with open(snakemake.output[0], 'w') as outFH:
    outFH.write('\t'.join((
        'contigID',
        'seqID',
        'start',
        'stop',
        'len',
        'qual',
        'count',
        'CPM',
        'alnType',
        'taxMethod',
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'baltimoreType',
        'baltimoreGroup\n'
    )))
    bam = pysam.AlignmentFile(snakemake.input.bam, 'rb')
    for read in bam.fetch():
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        infOut = '\t'.join([str(read.reference_name),
                            str(read.query_name),
                            str(read.reference_start),
                            str(read.reference_end),
                            str(read.reference_length),
                            str(read.mapping_quality)])
        try:
            taxOut = tax[read.query_name]
        except KeyError:
            taxOut = '\t'.join((['NA'] * 8))
        outFH.write(f'{infOut}\t{taxOut}\n')
    bam.close()

logging.debug('Done')
