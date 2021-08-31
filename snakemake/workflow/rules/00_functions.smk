"""
Generic rules and functions for Hecatomb
"""

### RECIPES
rule fasta_index:
    """Index a .fasta file for rapid access with samtools faidx."""
    input:
        "{file}.fasta"
    output:
        "{file}.fasta.fai"
    log:
        "{file}.samtools.stderr"
    threads:
        MiscCPU
    resources:
        mem_mb=MiscMem
    wrapper:
        "0.77.0/bio/samtools/faidx"

rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    log:
        "{file}.samtools.stderr"
    threads:
        MiscCPU
    resources:
        mem_mb=MiscMem
    wrapper:
        "0.77.0/bio/samtools/index"

rule calculate_gc:
    """Step 13a: Calculate GC content for sequences"""
    input:
        os.path.join(RESULTS,"{file}.fasta")
    output:
        temp(os.path.join(RESULTS,"{file}.properties.gc"))
    benchmark:
        os.path.join(BENCH,"calculate_gc.{file}.txt")
    log:
        os.path.join(STDERR,"calculate_gc.{file}.log")
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    threads:
        MiscCPU
    resources:
        mem_mb=MiscMem
    shell:
        """
        countgc.sh in={input} format=2 ow=t > {output} 2> {log}
        """

rule calculate_tet_freq:
    """Step 13b: Calculate tetramer frequency

    The tail commands trims the first line which is junk that should be printed to stdout.
    """
    input:
        os.path.join(RESULTS,"{file}.fasta")
    output:
        temp(os.path.join(RESULTS,"{file}.properties.tetramer"))
    benchmark:
        os.path.join(BENCH,"calculate_tet.{file}.txt")
    log:
        os.path.join(STDERR,"calculate_tet.{file}.log")
    conda:
        os.path.join("..", "envs", "bbmap.yaml")
    threads:
        MiscCPU
    resources:
        mem_mb=MiscMem
    shell:
        """
        tetramerfreq.sh in={input} w=0 ow=t \
            | tail -n+2 \
            > {output} 2> {log}
        """

rule seq_properties_table:
    """Step 13c: Combine GC and tet freq tables"""
    input:
        gc=os.path.join(RESULTS,"{file}.properties.gc"),
        tet=os.path.join(RESULTS,"{file}.properties.tetramer")
    output:
        os.path.join(RESULTS,"{file}.properties.tsv")
    benchmark:
        os.path.join(BENCH,"seq_properties_table.{file}.txt")
    threads:
        MiscCPU
    resources:
        mem_mb=MiscMem
    run:
        gc = {}
        for l in stream_tsv(input.gc):
            if len(l) == 2:
                gc[l[0]] = l[1]
        out = open(output[0],'w')
        for l in stream_tsv(input.tet):
            if l[0] == 'scaffold':
                out.write('id\tGC\t')
            else:
                out.write(f'{l[0]}\t{gc[l[0]]}\t')
            out.write('\t'.join(l[2:]))
            out.write('\n')
        out.close()

# rule zip:
#     input:
#         '{file}'
#     output:
#         '{file}.gz'
#     wildcard_constraints:
#         file=".*\.(?!gz$)\w*$"
#     shell:
#         "gzip {input}"
#
# rule zstd_decomp:
#     input:
#         '{file}.zst'
#     output:
#         '{file}'
#     wildcard_constraints:
#         file=".*\.(?!zst$)\w*$"
#     conda:
#         "../envs/curl.yaml"
#     shell:
#         'zstd -qd {input}'

### FUNCTIONS
def stream_tsv(tsvFile):
    """Read a file line-by-line and split by whitespace"""
    filehandle = open(tsvFile, 'r')
    for line in filehandle:
        line = line.strip()
        l = line.split('\t')
        yield l
    filehandle.close()
