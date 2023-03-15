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
    os.path.join(dir.out.results, 'hecatomb.samples.tsv'),
    expand(
        os.path.join(dir.out.assembly,"{sample}{file}"),
        sample=samples.names,
        file=config.modules[config.args.library]["targets"]
    )
]


# Assembly files
targets.assembly = [
    os.path.join(dir.out.results,f"{config.args.assembly}_assembly.fasta"),
    os.path.join(dir.out.results,f"{config.args.assembly}_assembly_graph.gfa"),
    os.path.join(dir.out.results,"contig_count_table.tsv"),
    os.path.join(dir.out.results,f"{config.args.assembly}_assembly.properties.tsv"),
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
    os.path.join(dir.out.results, "taxonLevelCounts.tsv"),
    os.path.join(dir.out.results, "krona.html"),
    os.path.join(dir.out.results, "seqtable.unclassified.fasta")
]


# Summary files
targets.summary = [
    os.path.join(dir.out.results, "summary.tsv"),
]


