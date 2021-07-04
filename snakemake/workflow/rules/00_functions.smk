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
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input} -o {output}
        """

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
        "../envs/bbmap.yaml"
    shell:
        "countgc.sh in={input} format=2 ow=t > {output} 2> {log}"

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
        "../envs/bbmap.yaml"
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

rule zip:
    input:
        '{file}'
    output:
        '{file}.gz'
    shell:
        "gzip {input}"

### FUNCTIONS
def stream_tsv(tsvFile):
    """Read a file line-by-line and split by whitespace"""
    filehandle = open(tsvFile, 'r')
    for line in filehandle:
        l = line.split()
        yield l
    filehandle.close()
