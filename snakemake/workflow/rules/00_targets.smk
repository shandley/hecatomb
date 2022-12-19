"""
All target output files for Hecatomb are declared here
"""

import attrmap as ap

targets = ap.AttrMap()


# Preprocessing files (more preprocessing-specific targets specified in 01_preprocessing*.smk files)
targets.preprocessing = [
    os.path.join(dir.out.results, "seqtable.fasta"),
    os.path.join(dir.out.results, "sampleSeqCounts.tsv"),
    os.path.join(dir.out.results, "seqtable.properties.tsv"),
]


# Assembly files
targets.assembly = [
    os.path.join(dir.out.results,"assembly.fasta"),
    os.path.join(dir.out.results,"contig_count_table.tsv"),
    os.path.join(dir.out.results,"assembly.properties.tsv"),
    ]


# Contig annotations
targets.contigAnnotations = [
    os.path.join(dir.out.results,"contigAnnotations.tsv"),
]


# Mapping files
targets.mapping = [
    os.path.join(dir.out.mapping,"assembly.seqtable.bam"),
    os.path.join(dir.out.mapping,"assembly.seqtable.bam.bai"),
    os.path.join(dir.out.results,"contigSeqTable.tsv"),
    os.path.join(dir.out.results,"contigKrona.html")
]


# Secondary AA search files
targets.readAnnotations = [
    os.path.join(dir.out.secondaryAA, "AA_bigtable.tsv"),
    os.path.join(dir.out.secondaryNT, "NT_bigtable.tsv"),
    os.path.join(dir.out.results, "bigtable.tsv"),
]


# Summary files
targets.summary = [
    os.path.join(dir.out.results, 'hecatomb.samples.tsv'),
    os.path.join(dir.out.results, "taxonLevelCounts.tsv"),
    os.path.join(dir.out.results, "krona.html")
]


# optionalSummary = [
#     os.path.join(SUMDIR,"Step00_counts.tsv"),
#     os.path.join(SUMDIR,"Step01_counts.tsv"),
#     os.path.join(SUMDIR,"Step02_counts.tsv"),
#     # os.path.join(SUMDIR,"Step03_counts.tsv"),
#     # os.path.join(SUMDIR,"Step04_counts.tsv"),
#     # os.path.join(SUMDIR,"Step05_counts.tsv"),
#     # os.path.join(SUMDIR,"Step06_counts.tsv"),
#     # os.path.join(SUMDIR,"Step07_counts.tsv"),
#     # os.path.join(SUMDIR,"Step08_counts.tsv"),
#     # os.path.join(SUMDIR,"Step09_counts.tsv"),
#     os.path.join(SUMDIR,"Step10_counts.tsv"),
#     os.path.join(SUMDIR,"Step11_counts.tsv"),
#     os.path.join(SUMDIR,"Step12_counts.tsv"),
#     os.path.join(SUMDIR,"Step13_counts.tsv"),
#     # os.path.join(SUMDIR,"Sankey.svg"),
# ]

