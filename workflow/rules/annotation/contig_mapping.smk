rule contig_mapping_map_seq_table:
    """Mapping step 01: Map the seq table to the population assembly.
    
    The seqtable reads are mapped to the assembly and the taxon information for the reads is used to help identify
    the organism to which each contig belongs.
    """
    input:
        assembly = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], config["hecatomb"]["args"]["assembly"] + '_assembly.fasta'),
        seqtable = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "seqtable.fasta")
    output:
        os.path.join(config["hecatomb"]["args"]["temp_paths"]["mapping"], "assembly.seqtable.bam")
    log:
        mm2 = os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "map_seq_table.mm2.log"),
        stool = os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "map_seq_table.stools.log")
    conda:
        os.path.join("..", "..", "envs", "minimap2_samtools.yaml")
    container:
        config["hecatomb"]["container"]["minimap2_samtools"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "map_seq_table.txt")
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    shell:
        "minimap2 -ax sr --secondary=no -t {threads} {input.assembly} {input.seqtable} 2> {log.mm2} "
            "| samtools sort -@ {threads} -o {output} 2> {log.stool}; "


rule contig_mapping_read_taxonomy:
    """Mapping step 02: Create a table of seq mapping info and their taxonomic information
    
    1. read the tax information from bigtable.tsv to a dictionary
    2. parse the aligned sequences and print a row for each primary alignment
    output format:
    contigID\tseqID\tstart\tstop\tlen\tqual\tcount\tpercent\ttaxMethod\tK\tP\tC\tO\tF\tG\tS\n
    """
    input:
        bam = os.path.join(config["hecatomb"]["args"]["temp_paths"]["mapping"], "assembly.seqtable.bam"),
        bai = os.path.join(config["hecatomb"]["args"]["temp_paths"]["mapping"], "assembly.seqtable.bam.bai"),
        taxon = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "bigtable.tsv"),
        counts = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "sampleSeqCounts.tsv"),
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "contigSeqTable.tsv")
    params:
        contigTaxonHeader = config["hecatomb"]["immutable"]["contigTaxonHeader"]
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "contig_read_taxonomy.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "contig_read_taxonomy.log")
    conda:
        os.path.join("..", "..", "envs", "krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    script:
        os.path.join("..", "..", "scripts",  "contigReadTaxon.py")
