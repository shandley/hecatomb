rule map_seq_table:
    """Mapping step 01: Map the seq table to the population assembly.
    
    The seqtable reads are mapped to the assembly and the taxon information for the reads is used to help identify
    the organism to which each contig belongs.
    """
    input:
        assembly = os.path.join(RESULTS, "assembly.fasta"),
        seqtable = os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(MAPPING, "assembly.seqtable.bam")
    log:
        mm2 = os.path.join(STDERR, "map_seq_table.mm2.log"),
        stool = os.path.join(STDERR, "map_seq_table.stools.log")
    conda:
        os.path.join('../', 'envs', 'minimap2.yaml')
    benchmark:
        os.path.join(BENCH, 'map_seq_table.txt')
    threads:
        MhitCPU
    resources:
        mem_mb = MhitMem
    shell:
        """
        minimap2 -ax sr --secondary=no -t {threads} {input.assembly} {input.seqtable} 2> {log.mm2} \
            | samtools sort -m 1G -o {output} 2> {log.stool}
        rm {log.mm2} {log.stool}
        """

rule contig_read_taxonomy:
    """Mapping step 02: Create a table of seq mapping info and their taxonomic information
    
    1. read the tax information from bigtable.tsv to a dictionary
    2. parse the aligned sequences and print a row for each primary alignment
    output format:
    contigID\tseqID\tstart\tstop\tlen\tqual\ttaxMethod\tK\tP\tC\tO\tF\tG\tS\n
    """
    input:
        bam = os.path.join(MAPPING, "assembly.seqtable.bam"),
        bai = os.path.join(MAPPING, "assembly.seqtable.bam.bai"),
        taxon = os.path.join(RESULTS, "bigtable.tsv")
    output:
        os.path.join(RESULTS, "contigSeqTable.tsv")
    threads:
        MiscCPU
    resources:
        mem_mb = MiscMem
    benchmark:
        os.path.join(BENCH, 'contig_read_taxonomy.txt')
    log:
        os.path.join(STDERR, 'contig_read_taxonomy.log')
    conda:
        os.path.join('..', 'envs', 'pysam.yaml')
    script:
        os.path.join('..', 'scripts', 'contigReadTaxon.py')
