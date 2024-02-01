
targets = dict()


### READ TRIMMING TARGETS
    # Also need some python logic for MEGAHIT input params
targets["trimnami"] = []
targets["cross"] = dict()
targets["cross"]["r1"] = []
targets["cross"]["r2"] = []
targets["cross"]["s"] = []
targets["unmapped"] = dict()
targets["unmapped"]["r1"] = []
targets["unmapped"]["r2"] = []
targets["unmapped"]["s"] = []
targets["trimmed"] = dict()
samples["trimmed"] = dict()


if config["args"]["host"].lower() == "none":
    config["args"]["hostStr"] = ""
else:
    config["args"]["hostStr"] = ".host_rm"


for sample_name in samples["names"]:
    targets["trimmed"][sample_name] = dict()
    samples["trimmed"][sample_name] = dict()
    if samples["reads"][sample_name]["R2"] is not None:
        samples["trimmed"][sample_name]["R1"] = os.path.join(dir["out"]["trim"], sample_name + "_R1" + config["args"]["hostStr"] + ".fastq.gz")
        samples["trimmed"][sample_name]["R2"] = os.path.join(dir["out"]["trim"], sample_name + "_R2" + config["args"]["hostStr"] + ".fastq.gz")
        samples["trimmed"][sample_name]["S"] = os.path.join(dir["out"]["trim"], sample_name + "_RS" + config["args"]["hostStr"] + ".fastq.gz")

        targets["cross"]["r1"].append(samples["trimmed"][sample_name]["R1"])
        targets["cross"]["r2"].append(samples["trimmed"][sample_name]["R2"])
        targets["cross"]["s"].append(samples["trimmed"][sample_name]["S"])

        targets["unmapped"]["r1"].append(os.path.join(dir["out"]["assembly"], sample_name, sample_name + ".assemblyUnmapped_R1.fastq"))
        targets["unmapped"]["r2"].append(os.path.join(dir["out"]["assembly"], sample_name, sample_name + ".assemblyUnmapped_R2.fastq"))
        targets["unmapped"]["s"].append(os.path.join(dir["out"]["assembly"], sample_name, sample_name + ".assemblyUnmapped_RS.fastq"))

        targets["trimnami"] += [
            samples["trimmed"][sample_name]["R1"],
            samples["trimmed"][sample_name]["R2"],
            samples["trimmed"][sample_name]["S"],
        ]
    else:
        samples["trimmed"][sample_name]["R1"] = os.path.join(dir["out"]["trim"], sample_name + "_S" + config["args"]["hostStr"] + ".fastq.gz")
        samples["trimmed"][sample_name]["R2"] = None
        samples["trimmed"][sample_name]["S"] = None
        targets["cross"]["s"].append(samples["trimmed"][sample_name]["R1"])
        targets["unmapped"]["s"].append(os.path.join(dir["out"]["assembly"],sample_name,sample_name + ".assemblyUnmapped_S.fastq"))
        targets["trimnami"].append(samples["trimmed"][sample_name]["R1"])

targets["cross"]["r1"] = "-1 " + ",".join(targets["cross"]["r1"])
targets["cross"]["r2"] = "-2 " + ",".join(targets["cross"]["r2"])
# targets["cross"]["s"]  = "-r " + ",".join(targets["cross"]["s"])


### PREPROCESSING
targets["preprocessing"] = targets["trimnami"]


### ASSEMBLY
targets["assembly"] = [
    os.path.join(dir["out"]["results"], config["args"]["assembly"] + '_assembly.fasta'),
    os.path.join(dir["out"]["results"], config["args"]["assembly"] + '_assembly.gfa'),
    os.path.join(dir["out"]["results"], "sample_coverage.tsv"),
    os.path.join(dir["out"]["results"], "all_coverage.tsv"),
]


### CONTIG ANNOTATIONS
targets["contigAnnotations"] = [
    os.path.join(dir["out"]["results"], "contigAnnotations.tsv"),
]

targets["mapping"] = [
    os.path.join(dir["out"]["mapping"], "assembly.seqtable.bam"),
    os.path.join(dir["out"]["mapping"], "assembly.seqtable.bam.bai"),
    os.path.join(dir["out"]["results"], "contigSeqTable.tsv"),
    os.path.join(dir["out"]["results"], "contigKrona.html")
]


### READ ANNOTATIONS
targets["readAnnotations"] = [
    os.path.join(dir["out"]["results"],"seqtable.fasta"),
    os.path.join(dir["out"]["results"],"sampleSeqCounts.tsv"),
    os.path.join(dir["out"]["results"],"seqtable.properties.tsv"),
    os.path.join(dir["out"]["secondaryAA"], "AA_bigtable.tsv"),
    os.path.join(dir["out"]["secondaryNT"], "NT_bigtable.tsv"),
    os.path.join(dir["out"]["results"], "bigtable.tsv"),
    os.path.join(dir["out"]["results"], "taxonLevelCounts.tsv"),
    os.path.join(dir["out"]["results"], "krona.html"),
    os.path.join(dir["out"]["results"], "seqtable.unclassified.fasta")
]


### SUMMARY REPORTS
targets["summary"] = [
    os.path.join(dir["out"]["results"], "summary.tsv"),
]


# BUILD ENV TARGETS
targets["envs"] = []

for filename in os.listdir(dir["env"]):
    if filename.endswith(".yaml") or filename.endswith(".yml"):
        targets["envs"].append(os.path.join(dir["out"]["temp"], filename + ".done"))