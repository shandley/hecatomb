"""
Generic rules and functions for Hecatomb
"""

import re


### PYTHON FUNCTIONS
def stream_tsv(tsvFile, skip_header=False):
    """Read a file line-by-line and split by whitespace"""
    with open(tsvFile, "r") as filehandle:
        if skip_header:
            filehandle.readline()
        for line in filehandle:
            line = line.strip()
            l = line.split("\t")
            yield l


def file_len(fname):
    """Return the number of lines in a file"""
    if fname.endswith(".gz"):
        import gzip
        f = gzip.open(fname, "rb")
    else:
        f = open(fname, "r")
    if re.search("fastq(\.gz)?$", fname):
        for i, l in enumerate(f):
            pass
        f.close()
        return int((i + 1) / 4)
    else:
        n=0
        for line in f:
            if line.startswith('>'):
                n+=1
        f.close()
        return n


## GENERIC RECIPES
rule fasta_index:
    """Index a .fasta file for rapid access with samtools faidx."""
    input:
        "{file}.fasta"
    output:
        "{file}.fasta.fai"
    log:
        "{file}.samtools.stderr"
    conda:
        os.path.join(dir.env,   "samtools.yaml")
    resources:
        mem_mb = config.resources.ram.mem
    shell:
        "samtools faidx {input} > {output} 2> {log} && rm {log}"


rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    log:
        "{file}.samtools.stderr"
    conda:
        os.path.join(dir.env,   "samtools.yaml")
    threads:
        config.resources.ram.cpu
    resources:
        mem_mb = config.resources.ram.mem
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log} && rm {log}"


rule calculate_gc:
    """Calculate GC content for sequences"""
    input:
        os.path.join(dir.out.results, "{file}.fasta")
    output:
        temp(os.path.join(dir.out.results, "{file}.properties.gc"))
    benchmark:
        os.path.join(dir.out.bench, "calculate_gc.{file}.txt")
    log:
        os.path.join(dir.out.stderr, "calculate_gc.{file}.log")
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    threads:
        config.resources.ram.cpu
    resources:
        mem_mb = config.resources.ram.mem
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
        os.path.join(dir.out.results, "{file}.fasta")
    output:
        temp(os.path.join(dir.out.results, "{file}.properties.tetramer"))
    benchmark:
        os.path.join(dir.out.bench, "calculate_tet_freq.{file}.txt")
    log:
        os.path.join(dir.out.stderr, "calculate_tet_freq.{file}.log")
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    threads:
        config.resources.ram.cpu
    resources:
        mem_mb = config.resources.ram.mem
    shell:
        """
        {{
        tetramerfreq.sh in={input} w=0 ow=t -Xmx{resources.mem_mb}m \
            | tail -n+2;
        }} > {output} 2>> {log}
        rm {log}
        """


rule seq_properties_table:
    """Combine GC and tet freq tables"""
    input:
        gc=os.path.join(dir.out.results, "{file}.properties.gc"),
        tet=os.path.join(dir.out.results, "{file}.properties.tetramer")
    output:
        os.path.join(dir.out.results, "{file}.properties.tsv")
    benchmark:
        os.path.join(dir.out.bench, "seq_properties_table.{file}.txt")
    threads:
        config.resources.ram.cpu
    resources:
        mem_mb = config.resources.ram.mem
    log:
        os.path.join(dir.out.stderr, "{file}.seq_properties_table.log")
    script:
        os.path.join(dir.scripts,  "seqPropertyTable.py")
