rule map_seq_table:
    """Mapping step 01: Map the seq table to the population assembly.
    
    The seqtable reads are mapped to the assembly and the taxon information for the reads is used to help identify
    the organism to which each contig belongs.
    """
    input:
        assembly = os.path.join(dir["out"]["results"], config["args"]["assembly"] + '_assembly.fasta'),
        seqtable = os.path.join(dir["out"]["results"], "seqtable.fasta")
    output:
        os.path.join(dir["out"]["mapping"], "assembly.seqtable.bam")
    log:
        mm2 = os.path.join(dir["out"]["stderr"], "map_seq_table.mm2.log"),
        stool = os.path.join(dir["out"]["stderr"], "map_seq_table.stools.log")
    conda:
        os.path.join(dir["env"], "minimap2.yaml")
    benchmark:
        os.path.join(dir["out"]["bench"], "map_seq_table.txt")
    threads:
        resources["lrg"]["cpu"]
    resources:
        mem_mb = resources["lrg"]["mem"],
        mem = str(resources["lrg"]["mem"]) + "MB",
        time = resources["lrg"]["time"]
    group:
        "contigmap"
    shell:
        "minimap2 -ax sr --secondary=no -t {threads} {input.assembly} {input.seqtable} 2> {log.mm2} "
            "| samtools sort -@ {threads} -o {output} 2> {log.stool}; "


rule contig_read_taxonomy:
    """Mapping step 02: Create a table of seq mapping info and their taxonomic information
    
    1. read the tax information from bigtable.tsv to a dictionary
    2. parse the aligned sequences and print a row for each primary alignment
    output format:
    contigID\tseqID\tstart\tstop\tlen\tqual\tcount\tpercent\ttaxMethod\tK\tP\tC\tO\tF\tG\tS\n
    """
    input:
        bam = os.path.join(dir["out"]["mapping"], "assembly.seqtable.bam"),
        bai = os.path.join(dir["out"]["mapping"], "assembly.seqtable.bam.bai"),
        taxon = os.path.join(dir["out"]["results"], "bigtable.tsv"),
        counts = os.path.join(dir["out"]["results"], "sampleSeqCounts.tsv"),
    output:
        os.path.join(dir["out"]["results"], "contigSeqTable.tsv")
    params:
        contigTaxonHeader = config["immutable"]["contigTaxonHeader"]
    threads:
        resources["lrg"]["cpu"]
    resources:
        mem_mb = resources["lrg"]["mem"],
        mem = str(resources["lrg"]["mem"]) + "MB",
        time = resources["lrg"]["time"]
    benchmark:
        os.path.join(dir["out"]["bench"], "contig_read_taxonomy.txt")
    log:
        os.path.join(dir["out"]["stderr"], "contig_read_taxonomy.log")
    conda:
        os.path.join(dir["env"], "pysam.yaml")
    group:
        "contigmap"
    script:
        os.path.join(dir["scripts"],  "contigReadTaxon.py")
