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
    script:
        os.path.join('..','scripts','kronaText.py')


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
    script:
        os.path.join('..','scripts','contigKronaText.py')


rule contig_krona_plot:
    input:
        os.path.join(SUMDIR, "contigKrona.txt")
    output:
        os.path.join(SUMDIR, "contigKrona.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    log:
        os.path.join(STDERR, 'contig_krona_plot.log')
    resources:
        mem_mb = MiscMem
    shell:
        """
        ktImportText {input} -o {output} &> {log}
        rm {log}
        """
