if makeReport:
    include: '04_summaries_optional.smk'
else:
    rule touchSummCounts:
        """dummy rule to 'skip' summary counting"""
        output:
            touch(optionalSummary)

rule tax_level_counts:
    """Generate a long table of hit counts at different taxon levels
    
    (excluding at the species level as that is essentially the bigtable.tsv file)
    """
    input:
        os.path.join(RESULTS, "bigtable.tsv")
    output:
        report(os.path.join(SUMDIR, "taxonLevelCounts.tsv"),
            caption = "../report/tax_level_counts.rst",
            category = "Output")
    params:
        samples = list(SAMPLES)
    log:
        os.path.join(STDERR, 'tax_level_counts.log')
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    script:
        os.path.join('..', 'scripts', 'taxLevelCounts.py')


rule dumpSamplesTsv:
    output:
        os.path.join(SUMDIR, 'hecatomb.samples.tsv')
    run:
        writeSamplesTsv(sampleReads, output[0])


rule krona_text_format:
    """Taxon step 18: Text format summary of bigtable for a krona plot"""
    input:
        os.path.join(RESULTS, "bigtable.tsv")
    output:
        os.path.join(SUMDIR, "krona.txt")
    benchmark:
        os.path.join(BENCH, "krona_text_format.txt")
    log:
        os.path.join(STDERR, 'krona_text_format.log')
    run:
        import logging
        import atexit
        atexit.register(exitLogCleanup,log[0])
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        logging.debug('Slurping Tax assignments from bigtable')
        counts = {}
        for l in stream_tsv(input[0]):
            if l[0]=="seqID":
                continue
            t = '\t'.join(l[23:])
            try:
                counts[t] += int(l[2])
            except KeyError:
                counts[t] = int(l[2])
        logging.debug('Sorting, counting, and writing tax assignments')
        outFH = open(output[0], 'w')
        for k in sorted(counts.keys()):
            outFH.write(f'{counts[k]}\t{k}\n')
        outFH.close()


rule krona_plot:
    """Taxon step 19: Krona plot of bigtable"""
    input:
        os.path.join(SUMDIR, "krona.txt")
    output:
        os.path.join(SUMDIR, "krona.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    benchmark:
        os.path.join(BENCH, "krona_plot.txt")
    log:
        os.path.join(STDERR, 'krona_plot.log')
    shell:
        """
        ktImportText {input} -o {output} &> {log}
        rm {log}
        """


rule contig_krona_text_format:
    input:
        os.path.join(RESULTS, "contigSeqTable.tsv")
    output:
        os.path.join(SUMDIR, "contigKrona.txt")
    log:
        os.path.join(STDERR, 'contig_krona_text_format.log')
    run:
        import logging
        import atexit
        atexit.register(exitLogCleanup,log[0])
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
        os.path.join(SUMDIR, "contigKrona.txt")
    output:
        os.path.join(SUMDIR, "contigKrona.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    log:
        os.path.join(STDERR, 'contig_krona_plot.log')
    shell:
        """
        ktImportText {input} -o {output} &> {log}
        rm {log}
        """
