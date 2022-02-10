
# folders to combine
allOutputDir = []
for f in config['outDirs']:
    allDbFiles.append(allOutputDir)

def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

# check for sample duplications and assembly files
allDirSmplLen = {}     # outDir:sample:len
allSamples = []
assemblyFiles = True
for f in allOutputDir:
    for assemblyFile in [os.path.join(f, 'RESULTS', 'assembly.fasta'), os.path.join(f, 'RESULTS', 'contigSeqTable.tsv')]:
        if not is_non_zero_file(assemblyFile):
            sys.stderr.write(f'No/missing assembly files for {f}, skipping assembly files.\n')
            assemblyFiles=False
    with open(os.path.join(f,'RESULTS','sampleSeqCounts.tsv'),'r') as t:
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
                sys.stderr.write(f'Ignoring duplicated sample {l[0]} in {f}\n')

# combined output directory
outDir = 'combinedHecatombOutput'

# TARGETS
targets = [
    os.path.join(outDir, 'bigtable.tsv'),
    os.path.join(outDir, 'seqtable.fasta'),
    os.path.join(outDir, 'seqtable.properties.tsv'),
    os.path.join(outDir, 'assembly.fasta'),
    os.path.join(outDir, 'sampleSeqCounts.tsv'),
    os.path.join(outDir, 'contigSeqTable.tsv'),
]


def combineResultDirOutput(outFile, inFileName, sampleCol=1, header=True, dirs=allDirSmplLen.keys()):
    with open(outFile,'w') as o:
        if header:
            with open(os.path.join(dirs[0],'RESULTS',inFileName),'r') as f:
                o.write(f.readline())   # print header
        for d in dirs:
            with open(os.path.join(d,'RESULTS',inFileName),'r') as f:
                if header:
                    f.readline()  # skip header
                for line in f:
                    l = line.strip().split('\t')
                    try:
                        allDirSmplLen[dir][l[sampleCol]]    # check if need to skip sample
                        o.write(line)
                    except KeyError:
                        pass

# make targets
rule all:
    input:
        targets

# hijack 03_mapping.smk for remaking the contigSeqTable
RESULTS = outDir
MAPPING = os.path.join(outDir, 'tmp')
include: "rules/00_functions.smk"
include: "rules/03_mapping.smk"

# rules
rule combineSampleSeqCounts:
    """Combine all sampleSeqCounts.tsv files"""
    input:
        expand(os.path.join('{dir}','RESULTS','sampleSeqCounts.tsv'), dir=allOutputDir)
    output:
        os.path.join(outDir, 'sampleSeqCounts.tsv')
    run:
        with open(output[0],'w') as o:
            for oDir in allDirSmplLen.keys():
                for smpl in allDirSmplLen[oDir].keys():
                    o.write(f'{smpl}\t{allDirSmplLen[oDir][smpl]}\n')

rule combineBigtables:
    """Combine all bigtable.tsv files"""
    input:
        expand(os.path.join('{dir}','RESULTS','bigtable.tsv'),dir=allOutputDir)
    output:
        os.path.join(outDir, 'bigtable.tsv')
    run:
        combineResultDirOutput(output[0],'bigtable.tsv')

rule combineSeqtableProperties:
    """Combine all seqtable.properties.tsv files"""
    input:
        expand(os.path.join('{dir}','RESULTS','seqtable.properties.tsv'),dir=allOutputDir)
    output:
        os.path.join(outDir, 'seqtable.properties.tsv')
    run:
        combineResultDirOutput(output[0],'seqtable.properties.tsv',sampleCol=0)

rule combineSeqTables:
    """Combine all seqtable.fasta files"""
    input:
        expand(os.path.join('{dir}','RESULTS','seqtable.fasta'),dir=allOutputDir)
    output:
        os.path.join(outDir, 'seqtable.fasta')
    run:
        with open(output[0],'w') as o:
            for oDir in allDirSmplLen.keys():
                with open(os.path.join(oDir,'RESULTS','seqtable.fasta'),'r') as f:
                    p=True
                    for line in f:
                        if line.startswith('>'):
                            s = line.replace('>','').split(':')
                            try:
                                allDirSmplLen[oDir][s]
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
        expand(os.path.join('{dir}','RESULTS','assembly.fasta'),dir=allOutputDir)
    output:
        assembly = os.path.join(outDir, 'assembly.fasta'),
        contigs = temp(os.path.join(outDir, 'allContigs.fa')),
        flye = temp(directory(os.path.join(outDir, 'Flye')))
    params:
        os.path.join(outDir, 'Flye', 'assembly.fasta')
    threads:
        config['MediumJobCpu']
    resources:
        mem_mb = config['MediumJobMem']
    shell:
        """
        cat {input} > {output.contigs}
        flye --subassemblies {output.contigs} -t {threads} --plasmids -o {output.flye} -g 1g
        mv {params} {output.assembly}
        """
