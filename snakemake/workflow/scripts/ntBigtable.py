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
logging.basicConfig(filename=snakemake.log[0] ,filemode='w' ,level=logging.DEBUG)

logging.debug('Slurping baltimore classifications')
counts = []
smplCounts = {}
balt = {}
with open(snakemake.input.balt, 'r') as balfh:
    balfh.readline()    # skip header
    for line in balfh:
        l = line.strip().split('\t')
        if len(l) == 3:
            balt[l[0]] = f'{l[1]}\t{l[2]}'

logging.debug('Reading sample total read counts')
with open(snakemake.input.counts, 'r') as cntfh:
    for line in cntfh:
        l = line.strip().split('\t')
        counts.append(int(l[1]))
        smplCounts[l[0]] = int(l[1])

logging.debug('Reading in lca-hit seq IDs')
lcaLin = {}
with open(snakemake.input.lca, 'r') as lcafh:
    for line in lcafh:
        l = line.strip().split('\t')
        # dont use lca lineage for root taxIDs - see config.yaml for IDs
        if not l[1] in snakemake.params.taxIdIgnore:
            lcaLin[l[0]] = '\t'.join((l[2:]))

logging.debug('Reading in tophit seq IDs')
topLin = {}
with open(snakemake.input.top, 'r') as topfh:
    for line in topfh:
        l = line.strip().split('\t')
        # skip if using lca lineage
        try:
            lcaLin[l[0]]
        except KeyError:
            topLin[l[0]] = '\t'.join((l[2:]))

logging.debug('Parsing alignments and writing NT bigtable output')
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
                         'baltimoreGroup')))
    out.write('\n')
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
            cpm = str((int(seqInf[1]) / smplCounts[seqInf[0]]) * 1000000)
            tName =  re.sub('.*\||[a-zA-Z]+=.*','',l[18])
            seqOut = '\t'.join((l[0], seqInf[0], seqInf[1], cpm))
            alnOut = 'nt\t' + '\t'.join((l[1:17]))
            try:
                baltOut = balt[taxOut.split()[5]]
            except KeyError:
                baltOut = 'NA\tNA'
            out.write('\t'.join((seqOut, alnOut, tName, taxOut, baltOut)))
            out.write('\n')

logging.debug('Done')