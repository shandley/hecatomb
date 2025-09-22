
rule summary_tax_level_counts:
    """Generate a long table of hit counts at different taxon levels
    
    (excluding at the species level as that is essentially the bigtable.tsv file)
    """
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "bigtable.tsv")
    output:
        report(os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "taxonLevelCounts.tsv"),
            caption = "../report/tax_level_counts.rst",
            category = "Output")
    params:
        samples = samples["names"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "tax_level_counts.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    script:
        os.path.join("..", "..", "scripts",  "taxLevelCounts.py")


rule summary_dump_samples_tsv:
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"hecatomb.samples.tsv")
    params:
        samples["reads"]
    localrule:
        True
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0], output[0])


rule summary_krona_text_format:
    """Taxon step 18: Text format summary of bigtable for a krona plot"""
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "bigtable.tsv")
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "krona.txt")
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "krona_text_format.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "krona_text_format.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    script:
        os.path.join("..", "..", "scripts",  "kronaText.py")


rule summary_krona_plot:
    """Taxon step 19: Krona plot of bigtable"""
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "krona.txt")
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "krona.html")
    conda:
        os.path.join("..", "..", "envs", "krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "krona_plot.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "krona_plot.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    shell:
        "ktImportText {input} -o {output} &> {log} "


rule summary_contig_krona_text_format:
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "contigSeqTable.tsv")
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "contigKrona.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "contig_krona_text_format.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    script:
        os.path.join("..", "..", "scripts",  "contigKronaText.py")


rule summary_contig_krona_plot:
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "contigKrona.txt")
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "contigKrona.html")
    conda:
        os.path.join("..", "..", "envs", "krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "contig_krona_plot.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    shell:
        "ktImportText {input} -o {output} &> {log} "
