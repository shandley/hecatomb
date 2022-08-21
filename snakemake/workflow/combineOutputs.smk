
# folders to combine
allOutputDir = config['outDirs']
if len(allOutputDir) < 2:
    sys.stderr.write('Error: Please specify two or more Hecatomb directories to combine\n')
    sys.exit(1)

def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

# check for sample duplications and assembly files
allDirSmplLen = {}     # RESULTS:sample:len
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


MMSeqsMem = config['BigJobMem']
MMSeqsCPU = config['BigJobCpu']
MMSeqsMemSplit = str(int(0.75 * int(MMSeqsMem))) + 'M'
MMSeqsTimeMin = config['BigJobTimeMin']
MhitMem = config['MediumJobMem']
MhitCPU = config['MediumJobCpu']
BBToolsMem = config['SmallJobMem']
BBToolsCPU = config['SmallJobCpu']
MiscMem = config['MoreRamMem']
MiscCPU = config['MoreRamCpu']
FastpMem = config['MediumJobMem']
FastpCPU = config['MediumJobCpu']


# hijack 03_mapping.smk for remaking the contigSeqTable
include: "rules/00_directories.smk"
include: "rules/00_functions.smk"
include: "rules/03_mapping.smk"


# TARGETS
targets = [
    os.path.join(RESULTS, 'bigtable.tsv'),
    os.path.join(RESULTS, 'seqtable.fasta'),
    os.path.join(RESULTS, 'seqtable.properties.tsv'),
    os.path.join(RESULTS, 'assembly.fasta'),
    os.path.join(RESULTS, 'sampleSeqCounts.tsv'),
    os.path.join(RESULTS, 'contigSeqTable.tsv'),
]


def combineResultDirOutput(outFile, inFileName, sampleCol=1, header=True, dirs=list(allDirSmplLen.keys())):
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
                        allDirSmplLen[d][l[sampleCol]]    # check if need to skip sample
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
        expand(os.path.join('{dir}','RESULTS','sampleSeqCounts.tsv'), dir=allOutputDir)
    output:
        os.path.join(RESULTS, 'sampleSeqCounts.tsv')
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
        os.path.join(RESULTS, 'bigtable.tsv')
    run:
        combineResultDirOutput(output[0],'bigtable.tsv')

rule combineSeqtableProperties:
    """Combine all seqtable.properties.tsv files"""
    input:
        expand(os.path.join('{dir}','RESULTS','seqtable.properties.tsv'),dir=allOutputDir)
    output:
        os.path.join(RESULTS, 'seqtable.properties.tsv')
    run:
        combineResultDirOutput(output[0],'seqtable.properties.tsv',sampleCol=0)

rule combineSeqTables:
    """Combine all seqtable.fasta files"""
    input:
        expand(os.path.join('{dir}','RESULTS','seqtable.fasta'),dir=allOutputDir)
    output:
        os.path.join(RESULTS, 'seqtable.fasta')
    run:
        with open(output[0],'w') as o:
            for oDir in allDirSmplLen.keys():
                with open(os.path.join(oDir,'RESULTS','seqtable.fasta'),'r') as f:
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
        expand(os.path.join('{dir}','RESULTS','assembly.fasta'),dir=allOutputDir)
    output:
        assembly = os.path.join(RESULTS, 'assembly.fasta'),
        contigs = temp(os.path.join(RESULTS, 'allContigs.fa')),
        flye = temp(directory(os.path.join(RESULTS, 'Flye')))
    params:
        os.path.join(RESULTS, 'Flye', 'assembly.fasta')
    threads:
        config['MediumJobCpu']
    resources:
        mem_mb = config['MediumJobMem']
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
