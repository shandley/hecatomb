#!/usr/bin/env python3

import os
import atexit
import logging
import re


logging.basicConfig(filename=snakemake.log[0] ,filemode='w' ,level=logging.DEBUG)

logging.debug('Slurping baltimore classifications')
balt = {}
with open(snakemake.input.balt, 'r') as balfh:
    balfh.readline()    # skip header
    for line in balfh:
        l = line.strip().split('\t')
        if len(l) == 3:
            balt[l[0]] = f'{l[1]}\t{l[2]}'

logging.debug('Reading in lca-hit seq IDs')
lcaLin = {}
with open(snakemake.input.lca, 'r') as lcafh:
    for line in lcafh:
        l = line.strip().split('\t')
        if not l[1] in snakemake.params.taxIdIgnore:
            lcaLin[l[0]] = l[2:]

logging.debug('Reading in tophit seq IDs')
topLin = {}
with open(snakemake.input.top, 'r') as topfh:
    for line in topfh:
        l = line.strip().split('\t')
        try:
            lcaLin[l[0]]
        except KeyError:
            topLin[l[0]] = l[2:]

logging.debug('Parsing alignments and writing NT bigtable output')
outVir = open(snakemake.output.vir, 'w')
outNonVir = open(snakemake.output.nonvir, 'w')

outVir.write('\t'.join((snakemake.params.bigtableHeader)) + '\n')
outNonVir.write('\t'.join((snakemake.params.bigtableHeader)) + '\n')

prevAln = str()
with open(snakemake.input.aln, 'r') as alnfh:
    for line in alnfh:
        l = line.strip().split('\t')
        if l[0] != prevAln:
            prevAln = l[0]
            try:
                taxOut = ['LCA'] + lcaLin[l[0]]
            except KeyError:
                try:
                    taxOut = ['TopHit'] + topLin[l[0]]
                except KeyError:
                    taxOut = ['NA'] * 8
            taxOutPrint = '\t'.join(taxOut)
            if taxOut[1] == 'Viruses':
                out = outVir
            else:
                out = outNonVir
            seqInf = l[0].split(':')
            seqOut = '\t'.join((l[0], seqInf[0], seqInf[1], seqInf[2]))
            alnOut = 'nt\t' + '\t'.join((l[1:17]))
            try:
                baltOut = balt[taxOut[5]]
            except KeyError:
                baltOut = 'NA\tNA'
            out.write('\t'.join((seqOut, alnOut, l[1], taxOutPrint, baltOut)))
            out.write('\n')

logging.debug('Done')