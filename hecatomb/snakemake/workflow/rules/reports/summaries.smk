
rule tax_level_counts:
    """Generate a long table of hit counts at different taxon levels
    
    (excluding at the species level as that is essentially the bigtable.tsv file)
    """
    input:
        os.path.join(dir["out"]["results"], "bigtable.tsv")
    output:
        report(os.path.join(dir["out"]["results"], "taxonLevelCounts.tsv"),
            caption = "../report/tax_level_counts.rst",
            category = "Output")
    params:
        samples = samples["names"]
    log:
        os.path.join(dir["out"]["stderr"], "tax_level_counts.log")
    resources:
        mem_mb=resources["ram"]["mem"],
        mem=str(resources["ram"]["mem"]) + "MB",
    threads:
        resources["ram"]["cpu"]
    script:
        os.path.join(dir["scripts"],  "taxLevelCounts.py")


rule dumpSamplesTsv:
    output:
        os.path.join(dir["out"]["results"],"hecatomb.samples.tsv")
    params:
        samples["reads"]
    localrule:
        True
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0], output[0])


rule krona_text_format:
    """Taxon step 18: Text format summary of bigtable for a krona plot"""
    input:
        os.path.join(dir["out"]["results"], "bigtable.tsv")
    output:
        os.path.join(dir["out"]["temp"], "krona.txt")
    benchmark:
        os.path.join(dir["out"]["bench"], "krona_text_format.txt")
    group:
        "krona"
    log:
        os.path.join(dir["out"]["stderr"], "krona_text_format.log")
    script:
        os.path.join(dir["scripts"],  "kronaText.py")


rule krona_plot:
    """Taxon step 19: Krona plot of bigtable"""
    input:
        os.path.join(dir["out"]["temp"], "krona.txt")
    output:
        os.path.join(dir["out"]["results"], "krona.html")
    conda:
        os.path.join(dir["env"], "krona.yaml")
    benchmark:
        os.path.join(dir["out"]["bench"], "krona_plot.txt")
    group:
        "krona"
    log:
        os.path.join(dir["out"]["stderr"], "krona_plot.log")
    resources:
        mem_mb=resources["ram"]["mem"],
        mem=str(resources["ram"]["mem"]) + "MB",
    shell:
        "ktImportText {input} -o {output} &> {log} "


rule contig_krona_text_format:
    input:
        os.path.join(dir["out"]["results"], "contigSeqTable.tsv")
    output:
        os.path.join(dir["out"]["temp"], "contigKrona.txt")
    group:
        "contig_krona"
    log:
        os.path.join(dir["out"]["stderr"], "contig_krona_text_format.log")
    script:
        os.path.join(dir["scripts"],  "contigKronaText.py")


rule contig_krona_plot:
    input:
        os.path.join(dir["out"]["temp"], "contigKrona.txt")
    output:
        os.path.join(dir["out"]["results"], "contigKrona.html")
    conda:
        os.path.join(dir["env"], "krona.yaml")
    group:
        "contig_krona"
    log:
        os.path.join(dir["out"]["stderr"], "contig_krona_plot.log")
    resources:
        mem_mb = resources["ram"]["mem"],
        mem = str(resources["ram"]["mem"]) + "MB",
    shell:
        "ktImportText {input} -o {output} &> {log} "
