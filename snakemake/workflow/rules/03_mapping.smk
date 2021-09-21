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
        mm2 = os.path.join(STDERR, "m01.map_seq_table.mm2.log"),
        stool = os.path.join(STDERR, "m01.map_seq_table.stools.log")
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
        os.path.join(BENCH, 'm02.contig_read_taxtable.txt')
    log:
        os.path.join(STDERR, 'm02.contig_read_taxtable.log')
    run:
        import logging
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        import pysam
        logging.debug('Slurping taxon info for seq IDs')
        tax = {}
        inTSV = open(input.taxon, 'r')
        inTSV.readline() # skip header
        for line in inTSV:
            l = line.strip().split('\t')
            tax[l[0]] = '\t'.join((l[3:5] + l[23:]))
        inTSV.close()
        logging.debug('Parsing mapped reads and pairing read taxon with mapped coords')
        outFH = open(output[0], 'w')
        outFH.write('contigID\tseqID\tstart\tstop\tlen\tqual\tcount\talnType\ttaxMethod\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n')
        bam = pysam.AlignmentFile(input.bam, 'rb')
        for read in bam.fetch():
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            infOut = '\t'.join([str(read.reference_name),
                                str(read.query_name),
                                str(read.reference_start),
                                str(read.reference_end),
                                str(read.reference_length),
                                str(read.mapping_quality)])
            try:
                taxOut = tax[read.query_name]
            except KeyError:
                taxOut = '\t'.join((['NA'] * 8))
            outFH.write(f'{infOut}\t{taxOut}\n')
        bam.close()
        outFH.close()
        logging.debug('Done')

rule contig_krona_text_format:
    input:
        os.path.join(RESULTS, "contigSeqTable.tsv")
    output:
        os.path.join(RESULTS, "contigKrona.txt")
    log:
        os.path.join(STDERR, 'm03.kronaText.log')
    run:
        import logging
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        logging.debug('Slurping contig seq table')
        counts = {}
        for l in stream_tsv(input[0]):
            if l[0] == "contigID":
                continue
            t = '\t'.join((l[9:] + [l[0]]))
            c = l[1].split(':')
            try:
                counts[t] += int(c[1])
            except KeyError:
                counts[t] = int(c[1])
        logging.debug('Sorting and writing contig taxon info')
        outFH = open(output[0],'w')
        for k in sorted(counts.keys()):
            outFH.write(f'{counts[k]}\t{k}\n')
        outFH.close()

rule contig_krona_plot:
    input:
        os.path.join(RESULTS, "contigKrona.txt")
    output:
        os.path.join(RESULTS, "contigKrona.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    log:
        os.path.join(STDERR, 'm04.kronaPlotContig.log')
    shell:
        """
        ktImportText {input} -o {output} &> {log}
        """