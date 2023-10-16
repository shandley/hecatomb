#!/usr/bin/env python3

import os
import atexit
import logging
import re


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Reading in Baltimore classifications")
balt = {}
with open(snakemake.input.balt, "r") as baltfh:
    baltfh.readline()  # skip header
    for line in baltfh:
        l = line.strip().split("\t")
        if len(l) == 3:
            balt[l[0]] = f"{l[1]}\t{l[2]}"

logging.debug("Reading in Taxon info for LCA seqs")
lcaLin = {}
with open(snakemake.input.lca, "r") as lcafh:
    for line in lcafh:
        l = line.strip().split("\t")
        if (
            not l[1] in snakemake.params.taxIdIgnore
        ):  # dont use lca lineage for root taxIDs - see config.yaml for IDs
            lcaLin[l[0]] = l[2:]

logging.debug("Reading in Taxon inf for tophit seqs")
topLin = {}
with open(snakemake.input.top, "r") as topfh:
    for line in topfh:
        l = line.strip().split("\t")
        try:
            lcaLin[l[0]]
        except KeyError:
            topLin[l[0]] = l[2:]

logging.debug("Parsing alignments and printing output on the fly")
outVir = open(snakemake.output.vir, "w")
outNonVir = open(snakemake.output.nonvir, "w")

outVir.write("\t".join((snakemake.params.bigtableHeader)) + "\n")
outNonVir.write("\t".join((snakemake.params.bigtableHeader)) + "\n")

prevAln = str()
with open(snakemake.input.aln, "r") as alnfh:
    for line in alnfh:
        l = line.strip().split("\t")
        if l[0] != prevAln:
            prevAln = l[0]
            try:
                taxOut = ["LCA"] + lcaLin[l[0]]
            except KeyError:
                try:
                    taxOut = ["TopHit"] + topLin[l[0]]
                except KeyError:
                    taxOut = ["NA"] * 8
            taxOutPrint = "\t".join(taxOut)
            if taxOut[1] == "Viruses":
                out = outVir
            else:
                out = outNonVir
            seqInf = l[0].split(":")  # seq ID = sample:count:perc:seqNum
            tName = re.sub(".*\||[a-zA-Z]+=.*", "", l[18])
            seqOut = "\t".join((l[0], seqInf[0], seqInf[1], seqInf[2]))
            l[15] = str(
                int(l[15]) * 3
            )  # convert aa alignment len to equivalent nt alignment len
            alnOut = "aa\t" + "\t".join((l[1:17]))
            try:
                baltOut = balt[taxOut[5]]
            except KeyError:
                baltOut = "NA\tNA"
            out.write("\t".join((seqOut, alnOut, tName, taxOutPrint, baltOut)))
            out.write("\n")


logging.debug("Done")
