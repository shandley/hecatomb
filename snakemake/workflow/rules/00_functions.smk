"""
Generic rules and functions for Hecatomb
"""

### PYTHON FUNCTIONS
def exitLogCleanup(*args):
    """Cleanup the logging file(s) prior to exiting"""
    for logFile in args:
        os.unlink(logFile)
    return None

def stream_tsv(tsvFile):
    """Read a file line-by-line and split by whitespace"""
    with open(tsvFile, 'r') as filehandle:
        for line in filehandle:
            line = line.strip()
            l = line.split('\t')
            yield l

def file_len(fname):
    """Return the number of lines in a file"""
    if fname.endswith('.gz'):
        import gzip
        f = gzip.open(fname, 'rb')
    else:
        f = open(fname, 'r')
    for i, l in enumerate(f):
        pass
    f.close()
    return i + 1

def fasta_clust_counts(fname):
    """Return the sum of seq count values in a fasta file of clustered seqs"""
    count = 0
    with open(fname, 'r') as filehandle:
        for line in filehandle:
            if line.startswith('>'):
                count += int(line.split(':')[1])     # fasta id = >sample:count:seqID
    return count

def collect_start_counts(sampleReads, outFile):
    """Collect the read counts of the raw input files"""
    with open(outFile, 'w') as outfh:
        for sample in SAMPLES:
            R1c = file_len(sampleReads[sample]['R1']) / 4
            R2c = file_len(sampleReads[sample]['R2']) / 4
            outfh.write(f'{sample}\tStep_00\tR1\t{R1c}\n')
            outfh.write(f'{sample}\tStep_00\tR2\t{R2c}\n')
    return None

def collect_counts(inPrefix, inSuffix, stepName, outFile):
    """Collect fastq counts for all samples and print to file"""
    with open(outFile,'w') as outfh:
        for sample in SAMPLES:
            R1c = file_len(os.path.join(inPrefix, f'{sample}_R1{inSuffix}')) / 4
            R2c = file_len(os.path.join(inPrefix, f'{sample}_R2{inSuffix}')) / 4
            outfh.write(f'{sample}\t{stepName}\tR1\t{R1c}\n')
            outfh.write(f'{sample}\t{stepName}\tR2\t{R2c}\n')
    return None

def sum_counts(fname, R1=False):
    """Collect the sum of all reads for all samples from a summary count file (e.g. from collect_counts)"""
    count = 0
    with open(fname, 'r') as infh:
        for line in infh:
            l = line.split()
            if R1:
                if l[2] == "R1":
                    count += float(l[3])
            else:
                count += float(l[3])
    return count


### RECIPES
rule fasta_index:
    """Index a .fasta file for rapid access with samtools faidx."""
    input:
        "{file}.fasta"
    output:
        "{file}.fasta.fai"
    log:
        "{file}.samtools.stderr"
    conda:
        os.path.join('..','envs','samtools.yaml')
    resources:
        mem_mb=MiscMem
    shell:
        "samtools faidx {input[0]} > {output[0]} 2> {log}"

rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    log:
        "{file}.samtools.stderr"
    conda:
        os.path.join('..','envs','samtools.yaml')
    threads:
        MiscCPU
    resources:
        mem_mb=MiscMem
    shell:
        "samtools index -@ {threads} {input[0]} {output[0]} 2> {log}"

rule calculate_gc:
    """Calculate GC content for sequences"""
    input:
        os.path.join(RESULTS, "{file}.fasta")
    output:
        temp(os.path.join(RESULTS, "{file}.properties.gc"))
    benchmark:
        os.path.join(BENCH, "calculate_gc.{file}.txt")
    log:
        os.path.join(STDERR, "calculate_gc.{file}.log")
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    threads:
        MiscCPU
    resources:
        mem_mb = MiscMem
    shell:
        """
        countgc.sh in={input} format=2 ow=t > {output} 2> {log}
        rm {log}
        """

rule calculate_tet_freq:
    """Calculate tetramer frequency

    The tail commands trims the first line which is junk that should be printed to stdout.
    """
    input:
        os.path.join(RESULTS, "{file}.fasta")
    output:
        temp(os.path.join(RESULTS, "{file}.properties.tetramer"))
    benchmark:
        os.path.join(BENCH, "calculate_tet_freq.{file}.txt")
    log:
        os.path.join(STDERR, "calculate_tet_freq.{file}.log")
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    threads:
        MiscCPU
    resources:
        mem_mb = MiscMem
    shell:
        """
        tetramerfreq.sh in={input} w=0 ow=t \
            | tail -n+2 \
            > {output} 2> {log}
        rm {log}
        """

rule seq_properties_table:
    """Combine GC and tet freq tables"""
    input:
        gc=os.path.join(RESULTS, "{file}.properties.gc"),
        tet=os.path.join(RESULTS, "{file}.properties.tetramer")
    output:
        os.path.join(RESULTS, "{file}.properties.tsv")
    benchmark:
        os.path.join(BENCH, "seq_properties_table.{file}.txt")
    threads:
        MiscCPU
    resources:
        mem_mb = MiscMem
    log:
        os.path.join(STDERR, '{file}.seq_properties_table.log')
    run:
        import logging
        import atexit
        atexit.register(exitLogCleanup,log[0])
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        gc = {}
        logging.debug('Reading GC counts')
        for l in stream_tsv(input.gc):
            if len(l) == 2:
                gc[l[0]] = l[1]
        out = open(output[0],'w')
        logging.debug('Parsing tet feqs and printing output')
        for l in stream_tsv(input.tet):
            if l[0] == 'scaffold':
                out.write('id\tGC\t')
            else:
                out.write(f'{l[0]}\t{gc[l[0]]}\t')
            out.write('\t'.join(l[2:]))
            out.write('\n')
        out.close()
        logging.debug('Done')

