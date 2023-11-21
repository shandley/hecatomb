#!/usr/bin/env python3

import os
import atexit
import logging


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

gc = {}
logging.debug("Reading GC counts")
with open(snakemake.input.gc, "r") as fh:
    for line in fh:
        l = line.strip().split("\t")
        if len(l) == 2:
            gc[l[0]] = l[1]

with open(snakemake.output[0], "w") as out:
    logging.debug("Parsing tet feqs and printing output")
    with open(snakemake.input.tet, "r") as fh:
        for line in fh:
            l = line.strip().split("\t")
            if l[0] == "scaffold":
                out.write("id\tGC\t")
            else:
                out.write(f"{l[0]}\t{gc[l[0]]}\t")
            out.write("\t".join(l[2:]))
            out.write("\n")

logging.debug("Done")
