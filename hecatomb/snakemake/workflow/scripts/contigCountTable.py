#!/usr/bin/env python3

import os
import atexit
import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

covStat = {}
logging.debug(f"Reading {snakemake.input.covstats}")
with open(snakemake.input.covstats, "r") as f:
    for line in f:
        if not line.startswith("#"):
            l = line.strip().split("\t")
            for i in [1, 2, 4, 9]:
                l[i] = str(round(float(l[i]), 2))
            covStat[l[0]] = [l[1], l[3], l[4], l[5], l[9]]

total = 0
logging.debug(f"Reading {snakemake.input.rpkm}")
with open(snakemake.input.rpkm, "r") as f:
    for line in f:
        if not line.startswith("#"):
            l = line.split("\t")
            rpk = int(l[2]) / (int(l[1]) / 1000)  # RPK
            total += rpk / 1000000  # size factor per million

logging.debug(f"Writing to {snakemake.output.count_tbl}")
with open(snakemake.output.count_tbl, "w") as o:
    o.write(
        "\t".join(
            [
                "Sample",
                "Contig",
                "Length",
                "Reads",
                "RPKM",
                "FPKM",
                "SPM",
                "AverageFold",
                "ReferenceGC",
                "CoveragePercentage",
                "coverageBases",
                "MedianFold",
            ]
        )
        + "\n"
    )
    with open(snakemake.input.rpkm, "r") as f:
        for line in f:
            if not line.startswith("#"):
                l = line.strip().split("\t")
                rpk = round(int(l[2]) / (int(l[1]) / 1000), 2)  # RPK
                try:
                    tpm = round(rpk / total, 2)  # TPM
                except ZeroDivisionError:
                    tpm = 0
                o.write(
                    "\t".join(
                        [
                            snakemake.wildcards.sample,  # sample
                            l[0],  # contig
                            l[1],  # length
                            l[4],  # reads
                            str(round(float(l[5]), 2)),  # RPKM
                            str(round(float(l[7]), 2)),  # FPKM
                            str(tpm),
                        ]
                        + covStat[l[0]]  # SPM  # aveFole -> medianFold
                    )
                    + "\n"
                )
