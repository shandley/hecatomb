
targets = dict()


# NEW HOST TARGET
config["hecatomb"]["outFasta"] = os.path.join(
    config["hecatomb"]["args"]["database_paths"]["hostBase"], config["hecatomb"]["args"]["host"], "masked_ref.fa.gz"
)
config["hecatomb"]["virRefSeq"] = os.path.join(config["hecatomb"]["args"]["database_paths"]["hostBase"], config["hecatomb"]["dbvirRefSeq"]["file"])


targets["cross"] = dict()
targets["cross"]["r1"] = []
targets["cross"]["r2"] = []
targets["cross"]["s"] = []

for sample in config["trimnami"]["trimmed"]:
    if "R1" in config["trimnami"]["trimmed"][sample]:
        targets["cross"]["r1"].append(config["trimnami"]["trimmed"][sample]["R1"])
    if "R2" in config["trimnami"]["trimmed"][sample]:
        targets["cross"]["r2"].append(config["trimnami"]["trimmed"][sample]["R2"])
    if "S" in config["trimnami"]["trimmed"][sample]:
        targets["cross"]["s"].append(config["trimnami"]["trimmed"][sample]["S"])

targets["cross"]["r1"] = "-1 " + ",".join(targets["cross"]["r1"]) if targets["cross"]["r1"] else ""
targets["cross"]["r2"] = "-2 " + ",".join(targets["cross"]["r2"]) if targets["cross"]["r2"] else ""


### PREPROCESSING
targets["preprocessing"] = config["trimnami"]["targets"]["reads"]
targets["preprocessing"].append(os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"seqtable.fasta"))


### ASSEMBLY
targets["assembly"] = [
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], config["hecatomb"]["args"]["assembly"] + '_assembly.fasta'),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], config["hecatomb"]["args"]["assembly"] + '_assembly.gfa'),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "sample_coverage.tsv"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "all_coverage.tsv"),
]


### CONTIG ANNOTATIONS
targets["contigAnnotations"] = [
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "contigAnnotations.tsv"),
]

targets["mapping"] = [
    os.path.join(config["hecatomb"]["args"]["temp_paths"]["mapping"], "assembly.seqtable.bam"),
    os.path.join(config["hecatomb"]["args"]["temp_paths"]["mapping"], "assembly.seqtable.bam.bai"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "contigSeqTable.tsv"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "contigKrona.html")
]


### READ ANNOTATIONS
targets["readAnnotations"] = [
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"seqtable.fasta"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"sampleSeqCounts.tsv"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"seqtable.properties.tsv"),
    os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "AA_bigtable.tsv"),
    os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "NT_bigtable.tsv"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "bigtable.tsv"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "taxonLevelCounts.tsv"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "krona.html"),
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "seqtable.unclassified.fasta")
]


### SUMMARY REPORTS
targets["summary"] = [
    os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "summary.tsv"),
]


# BUILD ENV TARGETS
targets["envs"] = []
# Envs in the envs/ directory
for filename in os.listdir(os.path.join(workflow.basedir, "envs")):
    if filename.endswith(".yaml") or filename.endswith(".yml"):
        targets["envs"].append(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], filename + ".done"))
# Envs build by sub-snaketools
#targets["envs"].append(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "subenvs.trimnami"))
#targets["envs"].append(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "subenvs.koverage"))
