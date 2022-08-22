"""
All target output files for Hecatomb are declared here
"""

# Preprocessing files (more preprocessing-specific targets specified in 01_preprocessing*.smk files)
PreprocessingFiles = [
    os.path.join(RESULTS, "seqtable.fasta"),
    os.path.join(RESULTS, "sampleSeqCounts.tsv"),
    os.path.join(RESULTS, "seqtable.properties.tsv"),
]


# Assembly files
AssemblyFiles = [
    os.path.join(RESULTS,"assembly.fasta"),
    os.path.join(RESULTS,"contig_count_table.tsv"),
    os.path.join(RESULTS,"assembly.properties.tsv"),
    ]


# Contig annotations
ContigAnnotFiles = [
    os.path.join(RESULTS,"contigAnnotations.tsv"),
]


# Mapping files
MappingFiles = [
    os.path.join(MAPPING,"assembly.seqtable.bam"),
    os.path.join(MAPPING,"assembly.seqtable.bam.bai"),
    os.path.join(RESULTS,"contigSeqTable.tsv"),
    os.path.join(SUMDIR,"contigKrona.html")
]


# Secondary AA search files
ReadAnnotationFiles = [
    os.path.join(SECONDARY_AA_OUT, "AA_bigtable.tsv"),
    os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv"),
    os.path.join(RESULTS, "bigtable.tsv"),
]


# Summary files
SummaryFiles = [
    optionalSummary,
    os.path.join(RESULTS, 'hecatomb.samples.tsv'),
    os.path.join(RESULTS, "taxonLevelCounts.tsv"),
    os.path.join(RESULTS, "krona.html")
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

