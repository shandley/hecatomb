rule population_assembly:
    """Assembly step 05: Create 'contig dictionary' of all unique contigs present in the study (aka: population assembly)"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_samples_contigs_size_selected.fasta")
    output:
        assembly = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")),
        stats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary.stats")
    params:
        flye_out = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE")
    benchmark:
        os.path.join(BENCH, "population_assembly.txt")
    log:
        log1 = os.path.join(STDERR, "population_assembly.flye.log"),
        log2 = os.path.join(STDERR, "population_assembly.stats.log")
    resources:
        mem_mb = MhitMem
    threads:
        MhitCPU
    conda:
        "../envs/metaflye.yaml"
    shell:
        """
        flye --subassemblies {input} -t {threads} --plasmids -o {params.flye_out} -g 1g &>> {log.log1}
        rm {log.log1}
        statswrapper.sh in={output.assembly} out={output.stats} \
            format=2 \
            ow=t 2> {log.log2}
        rm {log.log2}
        """

rule link_assembly:
    """Assembly step 06: Link the final assembly to the results directory; not really a step."""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")
    output:
        os.path.join(RESULTS, "assembly.fasta")
    run:
        os.rename(os.path.abspath(input[0]), os.path.abspath(output[0]))

# rule coverage_calculations

rule create_contig_count_table:
    """Assembly step 08: Transcript Per Million (TPM) calculator

    Useful resource: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/"""
    input:
        rpkm = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.rpkm"),
        covstats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.cov_stats")
    output:
        count_tbl = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_contig_counts.tsv"))
    benchmark:
        os.path.join(BENCH, "create_contig_count_table.{sample}.txt")
    log:
        os.path.join(STDERR, "create_contig_count_table.{sample}.log")
    script:
        os.path.join('..', 'scripts', 'contigCountTable.py')

rule concatentate_contig_count_tables:
    """Assembly step 09: Concatenate contig count tables

    Note: this is done as a separate rule due to how snakemake handles i/o files. It does not work well in Step 20b as 
    the i/o PATTERNS are different.
    """
    input:
        expand(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_contig_counts.tsv"), sample=SAMPLES)
    output:
        os.path.join(RESULTS, "contig_count_table.tsv")
    benchmark:
        os.path.join(BENCH, "concatentate_contig_count_tables.txt")
    log:
        os.path.join(STDERR, "concatentate_contig_count_tables.log")
    shell:
        """
        {{ cat {input} > {output};
            sed -i '1i sample_id\tcontig_id\tlength\treads\tRPKM\tFPKM\tTPM\tavg_fold_cov\tcontig_GC\tcov_perc\tcov_bases\tmedian_fold_cov' {output};
        }} 2> {log}
        rm {log}
        """
