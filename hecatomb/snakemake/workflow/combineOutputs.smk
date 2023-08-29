
### CONFIG
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
resources = config["resources"]
config = config["hecatomb"]

# folders to combine
if len(config["args"]["combineRuns"]) < 2:
    sys.stderr.write(
        "\nError: Please specify two or more Hecatomb directories to combine\n\n"
    )
    sys.exit(1)


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


# check for sample duplications and assembly files
allDirSmplLen = {}  # dir["out"]["results"]:sample:len
allSamples = []
assemblyFiles = True
for f in config["args"]["combineRuns"]:
    for assemblyFile in [
        os.path.join(f, "results", "assembly.fasta"),
        os.path.join(f, "results", "contigSeqTable.tsv"),
    ]:
        if not is_non_zero_file(assemblyFile):
            sys.stderr.write(
                f"No/missing assembly files for {f}, skipping assembly files.\n"
            )
            assemblyFiles = False
    with open(os.path.join(f, "results", "sampleSeqCounts.tsv"), "r") as t:
        for line in t:
            l = line.strip().split()
            if not l[0] in allSamples:
                allSamples.append(l[0])
                try:
                    allDirSmplLen[f][l[0]] = l[1]
                except KeyError:
                    allDirSmplLen[f] = {}
                    allDirSmplLen[f][l[0]] = l[1]
            else:
                sys.stderr.write(f"Ignoring duplicated sample {l[0]} in {f}\n")


# hijack contig_mapping.smk for remaking the contigSeqTable
include: "rules/preflight/directories.smk"
include: "rules/preflight/functions.smk"
include: "rules/preflight/contig_mapping.smk"


# TARGETS
targets = [
    os.path.join(dir["out"]["results"], "bigtable.tsv"),
    os.path.join(dir["out"]["results"], "seqtable.fasta"),
    os.path.join(dir["out"]["results"], "seqtable.properties.tsv"),
    os.path.join(dir["out"]["results"], "assembly.fasta"),
    os.path.join(dir["out"]["results"], "sampleSeqCounts.tsv"),
    os.path.join(dir["out"]["results"], "contigSeqTable.tsv"),
]


def combineResultDirOutput(
    outFile, inFileName, sampleCol=1, header=True, dirs=list(allDirSmplLen.keys())
):
    with open(outFile, "w") as o:
        if header:
            with open(os.path.join(dirs[0], "results", inFileName), "r") as f:
                o.write(f.readline())  # print header
        for d in dirs:
            with open(os.path.join(d, "results", inFileName), "r") as f:
                if header:
                    f.readline()  # skip header
                for line in f:
                    l = line.strip().split("\t")
                    try:
                        allDirSmplLen[d][l[sampleCol]]  # check if need to skip sample
                        o.write(line)
                    except KeyError:
                        pass


# make targets
rule all:
    input:
        targets


# rules
rule combineSampleSeqCounts:
    """Combine all sampleSeqCounts.tsv files"""
    input:
        expand(os.path.join('{dir}','results','sampleSeqCounts.tsv'), dir=config["args"]["combineRuns"])
    output:
        os.path.join(dir["out"]["results"], 'sampleSeqCounts.tsv')
    run:
        with open(output[0],'w') as o:
            for oDir in allDirSmplLen.keys():
                for smpl in allDirSmplLen[oDir].keys():
                    o.write(f'{smpl}\t{allDirSmplLen[oDir][smpl]}\n')


rule combineBigtables:
    """Combine all bigtable.tsv files"""
    input:
        expand(os.path.join('{dir}','results','bigtable.tsv'),dir=config["args"]["combineRuns"])
    output:
        os.path.join(dir["out"]["results"], 'bigtable.tsv')
    run:
        combineResultDirOutput(output[0],'bigtable.tsv')


rule combineSeqtableProperties:
    """Combine all seqtable.properties.tsv files"""
    input:
        expand(os.path.join('{dir}','results','seqtable.properties.tsv'),dir=config["args"]["combineRuns"])
    output:
        os.path.join(dir["out"]["results"], 'seqtable.properties.tsv')
    run:
        combineResultDirOutput(output[0],'seqtable.properties.tsv',sampleCol=0)


rule combineSeqTables:
    """Combine all seqtable.fasta files"""
    input:
        expand(os.path.join('{dir}','results','seqtable.fasta'),dir=config["args"]["combineRuns"])
    output:
        os.path.join(dir["out"]["results"], 'seqtable.fasta')
    run:
        with open(output[0],'w') as o:
            for oDir in allDirSmplLen.keys():
                with open(os.path.join(oDir,'results','seqtable.fasta'),'r') as f:
                    p=True
                    for line in f:
                        if line.startswith('>'):
                            s = line.replace('>','').split(':')
                            try:
                                allDirSmplLen[oDir][s[0]]
                                p=True
                                o.write(line)
                            except KeyError:
                                p=False
                        else:
                            if p:
                                o.write(line)


rule combineAssemblies:
    "Combine assemblies by running flye with subassemblies"
    input:
        expand(os.path.join('{dir}','results','assembly.fasta'),dir=config["args"]["combineRuns"])
    output:
        assembly = os.path.join(dir["out"]["results"], 'assembly.fasta'),
        contigs = temp(os.path.join(dir["out"]["results"], 'allContigs.fa')),
        flye = temp(directory(os.path.join(dir["out"]["results"], 'Flye')))
    params:
        os.path.join(dir["out"]["results"], 'Flye', 'assembly.fasta')
    threads:
        config["resources"]["med"]["cpu"]
    resources:
        mem_mb = config["resources"]["med"]["mem"]
    conda:
        os.path.join('envs', 'metaflye.yaml')
    shell:
        """
        n=0
        for i in {input}; do
          cat $i | sed "s/>/>$n/"
          n=0$n
        done > {output.contigs}
        flye --subassemblies {output.contigs} -t {threads} --plasmids -o {output.flye} -g 1g
        mv {params} {output.assembly}
        """
