import attrmap as ap
import attrmap.utils as au

targets = ap.AttrMap()

### READ TRIMMING TARGETS
    # Also need some python logic for MEGAHIT input params
targets.trimnami = []
targets.cross.r1 = []
targets.cross.r2 = []
targets.cross.s = []
targets.unmapped.r1 = []
targets.unmapped.r2 = []
targets.unmapped.s = []

if config.args.host.lower() == "none":
    config.args.hostStr = ""
else:
    config.args.hostStr = ".host_rm"

for sample_name in samples.names:
    if samples.reads[sample_name]["R2"] is not None:
        samples.trimmed[sample_name]["R1"] = os.path.join(dir.out.trim, sample_name + config.args.hostStr + ".paired.R1.fastq.gz")
        samples.trimmed[sample_name]["R2"] = os.path.join(dir.out.trim, sample_name + config.args.hostStr + ".paired.R2.fastq.gz")
        samples.trimmed[sample_name]["S"] = os.path.join(dir.out.trim, sample_name + config.args.hostStr + ".paired.S.fastq.gz")

        targets.cross.r1.append(samples.trimmed[sample_name]["R1"])
        targets.cross.r2.append(samples.trimmed[sample_name]["R2"])
        targets.cross.s.append(samples.trimmed[sample_name]["S"])

        targets.unmapped.r1.append(os.path.join(dir.out.assembly, sample_name, sample_name + ".assemblyUnmapped_R1.fastq"))
        targets.unmapped.r2.append(os.path.join(dir.out.assembly, sample_name, sample_name + ".assemblyUnmapped_R2.fastq"))
        targets.unmapped.s.append(os.path.join(dir.out.assembly, sample_name, sample_name + ".assemblyUnmapped.s.fastq"))

        targets.trimnami += [
            samples.trimmed[sample_name]["R1"],
            samples.trimmed[sample_name]["R2"],
            samples.trimmed[sample_name]["S"],
        ]
    else:
        samples.trimmed[sample_name]["R1"] = os.path.join(dir.out.trim, sample_name + config.args.hostStr + ".single.fastq.gz")
        samples.trimmed[sample_name]["R2"] = None
        samples.trimmed[sample_name]["S"] = None
        targets.cross.s.append(samples.trimmed[sample_name]["R1"])
        targets.unmapped.s.append(os.path.join(dir.out.assembly,sample_name,sample_name + ".assemblyUnmapped.s.fastq"))
        targets.trimnami.append(samples.trimmed[sample_name]["R1"])

targets.cross.r1 = "-1 " + ",".join(targets.cross.r1)
targets.cross.r2 = "-2 " + ",".join(targets.cross.r2)
targets.cross.s  = "-r " + ",".join(targets.cross.s)

samples = au.convert_state(samples, read_only=True)


### PREPROCESSING
targets.preprocessing = targets.trimnami


### ASSEMBLY
targets.assembly = [
    os.path.join(dir.out.results, f"{config.args.assembly}_assembly.fasta"),
    os.path.join(dir.out.results, f"{config.args.assembly}_assembly.gfa"),
    os.path.join(dir.out.results, "sample_coverage.tsv"),
    os.path.join(dir.out.results, "all_coverage.tsv"),
    ]


### CONTIG ANNOTATIONS
targets.contigAnnotations = [
    os.path.join(dir.out.results, "contigAnnotations.tsv"),
]

targets.mapping = [
    os.path.join(dir.out.mapping, "assembly.seqtable.bam"),
    os.path.join(dir.out.mapping, "assembly.seqtable.bam.bai"),
    os.path.join(dir.out.results, "contigSeqTable.tsv"),
    os.path.join(dir.out.results, "contigKrona.html")
]


### READ ANNOTATIONS
targets.readAnnotations = [
    os.path.join(dir.out.results,"seqtable.fasta"),
    os.path.join(dir.out.results,"sampleSeqCounts.tsv"),
    os.path.join(dir.out.results,"seqtable.properties.tsv"),
    os.path.join(dir.out.secondaryAA, "AA_bigtable.tsv"),
    os.path.join(dir.out.secondaryNT, "NT_bigtable.tsv"),
    os.path.join(dir.out.results, "bigtable.tsv"),
    os.path.join(dir.out.results, "taxonLevelCounts.tsv"),
    os.path.join(dir.out.results, "krona.html"),
    os.path.join(dir.out.results, "seqtable.unclassified.fasta")
]


### SUMMARY REPORTS
targets.summary = [
    os.path.join(dir.out.results, "summary.tsv"),
]
