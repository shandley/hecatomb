#!/usr/bin/env python3

import os
import logging
import plotly.graph_objects as go
import atexit


def sum_counts(fname, R1=False):
    """Collect the sum of all reads for all samples from a summary count file (e.g. from collect_counts)"""
    count = 0
    with open(fname, "r") as infh:
        for line in infh:
            l = line.split()
            if R1:
                if l[2] == "R1":
                    count += float(l[3])
            else:
                count += float(l[3])
    return count


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

logging.debug("Collecting summary counts")
s0 = sum_counts(snakemake.input[0])
s1 = sum_counts(snakemake.input[1])
s2 = sum_counts(snakemake.input[2])
s3 = sum_counts(snakemake.input[3])
s4 = sum_counts(snakemake.input[4])
s5 = sum_counts(snakemake.input[5])
s6 = sum_counts(snakemake.input[6])
s7 = sum_counts(snakemake.input[7])
s7R1 = sum_counts(snakemake.input[7], R1=True)
s8 = sum_counts(snakemake.input[8])
s9 = sum_counts(snakemake.input[9])
s1disc = s0 - s1
s2disc = s1 - s2
s3disc = s2 - s3
s4disc = s3 - s4
s5disc = s4 - s5
s6disc = s5 - s6
s7disc = s6 - s7
s8R2 = s7 - s7R1
s8disc = s7R1 - s8
s9rdnt = s8 - s9
searches = {}
for f in snakemake.input[10:]:
    with open(f, "r") as fh:
        for line in fh:
            l = line.strip().split("\t")
            searches[l[0]] = int(l[1])

paaV = searches["Primary AA viral:"]
paaU = searches["Primary AA non-viral:"]
saaV = searches["Secondary AA viral:"]
saaU = searches["Secondary AA non-viral:"]
pntV = searches["Primary NT viral:"]
pntU = searches["Primary NT non-viral:"]
sntV = searches["Secondary NT viral:"]
sntU = searches["Secondary NT non-viral:"]

labels = [
    "Raw reads",  # 0
    "5' primer",  # 1
    "3' read-through",  # 2
    "Primer-free adapter",  # 3
    "Adapter-free primer",  # 4
    "Vector removal",  # 5
    "Low-qual trim",  # 6
    "Host removal",  # 7
    "Cluster seqs",  # 8
    "Representative seqs",  # 9
    "Redundant seqs",  # 10
    "Redundant (R2 reads)",  # 11
    "Discarded",  # 12
    "Seqtable (+ counts)",  # 13
    "Primary AA viral",  # 14
    "Primary AA non-viral",  # 15
    "Secondary AA viral",  # 16
    "Secondary AA non-viral",  # 17
    "Primary NT viral",  # 18
    "Primary NT non-viral",  # 19
    "Secondary NT viral",  # 20
    "Secondary NT non-viral",  # 21
    "Viruses",  # 22
    "Non-viral (virus-like)",  # 23
    "Non-viral",  # 24
]

source = [
    0,
    0,  # s1
    1,
    1,  # s2
    2,
    2,  # s3
    3,
    3,  # s4
    4,
    4,  # s5
    5,
    5,  # s6
    6,
    6,  # s7
    7,
    7,
    7,  # s8.1
    11,  # s8.2
    8,
    8,  # s9.1
    9,
    10,  # s9.2
    13,
    13,  # pAA
    14,
    14,  # sAA
    16,
    17,  # sAA
    15,
    15,  # pNT
    19,  # pNT
    18,
    18,  # sNT
    20,
    21,  # sNT
]

target = [
    1,
    12,  # s1
    2,
    12,  # s2
    3,
    12,  # s3
    4,
    12,  # s4
    5,
    12,  # s5
    6,
    12,  # s6
    7,
    12,  # s7
    8,
    11,
    12,  # s8.1
    12,  # s8.2
    9,
    10,  # s9.1
    13,
    13,  # s9.2
    14,
    15,  # pAA
    16,
    17,  # sAA
    22,
    23,  # sAA
    18,
    19,  # pNT
    24,  # pNT
    20,
    21,  # sNT
    22,
    23,  # sNT
]

values = [
    s1,
    s1disc,  # s1
    s2,
    s2disc,  # s2
    s3,
    s3disc,  # s3
    s4,
    s4disc,  # s4
    s5,
    s5disc,  # s5
    s6,
    s6disc,  # s6
    s7,
    s7disc,  # s7
    s8,
    s8R2,
    s8disc,  # s8.1
    s8R2,  # s8.2
    s9,
    s9rdnt,  # s9.1
    s9,
    s9rdnt,  # s9.2
    paaV,
    paaU,  # pAA
    saaV,
    saaU,  # sAA
    saaV,
    saaU,  # sAA
    pntV,
    pntU,  # pNT
    pntU,  # pNT
    sntV,
    sntU,  # sNT
    sntV,
    sntU,  # sNT
]

logging.debug("Generating sanky diagram")
link = dict(source=source, target=target, value=values)
node = dict(label=labels, pad=20, thickness=5)
data = go.Sankey(link=link, node=node)
fig = go.Figure(data)
fig.write_image(snakemake.output[0], width=2000, height=1000)
logging.debug("Done")
