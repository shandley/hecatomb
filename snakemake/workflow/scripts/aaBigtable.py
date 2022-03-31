#!/usr/bin/env python3

import os
import atexit
import logging
import re
from statistics import median

def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None

atexit.register(exitLogCleanup, snakemake.log[0])
logging.basicConfig(filename=snakemake.log[0],filemode='w',level=logging.DEBUG)

logging.debug('Reading in Baltimore classifications')
balt = {}
with open(snakemake.input.balt, 'r') as baltfh:
    baltfh.readline() # skip header
    for line in baltfh:
        l = line.strip().split('\t')
        if len(l) == 3:
            balt[l[0]] = f'{l[1]}\t{l[2]}'

logging.debug('Reading in reads counts for each sample')
counts = []
smplCounts = {}     # Read in read counts for each sample for normalised counts
with open(snakemake.input.counts, 'r') as countfh:
    for line in countfh:
        l = line.strip().split('\t')
        counts.append(int(l[1]))
        smplCounts[l[0]] = int(l[1])

logging.debug('Reading in Taxon info for LCA seqs')
lcaLin = {}
with open(snakemake.input.lca, 'r') as lcafh:
    for line in lcafh:
        l = line.strip().split('\t')
        # dont use lca lineage for root taxIDs - see config.yaml for IDs
        if not l[1] in snakemake.params.taxIdIgnore:
            lcaLin[l[0]] = '\t'.join((l[2:]))

logging.debug('Reading in Taxon inf for tophit seqs')
topLin = {}
with open(snakemake.input.top, 'r') as topfh:
    for line in topfh:
        l = line.strip().split('\t')
        try:
            lcaLin[l[0]]
        except KeyError:
            topLin[l[0]] = '\t'.join((l[2:]))

logging.debug('Parsing alignments and printing output on the fly')
with open(snakemake.output[0], 'w') as out:
    out.write('\t'.join(('seqID',
                         'sampleID',
                         'count',
                         'CPM',
                         'alnType',     # aa or nt
                         'targetID',
                         'evalue',
                         'pident',
                         'fident',
                         'nident',
                         'mismatches',
                         'qcov',
                         'tcov',
                         'qstart',
                         'qend',
                         'qlen',
                         'tstart',
                         'tend',
                         'tlen',
                         'alnlen',
                         'bits',
                         'targetName',
                         'taxMethod',
                         'kingdom',
                         'phylum',
                         'class',
                         'order',
                         'family',
                         'genus',
                         'species',
                         'baltimoreType',
                         'baltimoreGroup\n')))

    with open(snakemake.input.aln, 'r') as alnfh:
        for line in alnfh:
            l = line.strip().split('\t')
            try:
                taxOut = 'LCA\t' + lcaLin[l[0]]
            except KeyError:
                try:
                    taxOut = 'TopHit\t' + topLin[l[0]]
                except KeyError:
                    taxOut = '\t'.join((['NA'] * 8))
            # seq ID = sample:count:seqNum
            seqInf = l[0].split(':')
            tName =  re.sub('.*\||[a-zA-Z]+=.*','',l[18])
            cpm = str(( int(seqInf[1]) / smplCounts[seqInf[0]] ) * 1000000)
            seqOut = '\t'.join((l[0], seqInf[0], seqInf[1], cpm))
            # convert aa alignment len to equivalent nt alignment len
            l[15] = str(int(l[15]) * 3)
            alnOut = 'aa\t' + '\t'.join((l[1:17]))
            try:
                baltOut = balt[taxOut.split()[5]]
            except KeyError:
                baltOut = 'NA\tNA'
            out.write('\t'.join((seqOut, alnOut, tName, taxOut, baltOut)))
            out.write('\n')

logging.debug('Done')
