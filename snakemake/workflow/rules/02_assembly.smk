"""
What is accomplished with these rules?
    - Assembly
        - Sample assembly
        - Population assembly
        - Contig abundance esitmation
"""

rule assembly_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(TMPDIR, "p07", "{sample}_R1.unmapped.fastq"),
        r2 = os.path.join(TMPDIR, "p07", "{sample}_R2.unmapped.fastq"),
        r1s = os.path.join(TMPDIR, "p07", "{sample}_R1.singletons.fastq"),
        r2s = os.path.join(TMPDIR, "p07", "{sample}_R2.singletons.fastq")
    output:
        r1_norm = os.path.join(ASSEMBLY, "{sample}_R1.norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, "{sample}_R2.norm.fastq")
    benchmark:
        os.path.join(BENCH, "a01.kmer_normalization_{sample}.txt")
    log:
        os.path.join(STDERR, "a01.kmer_norm_{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} \
            extra={input.r1s},{input.r2s} \
            out={output.r1_norm} out2={output.r2_norm} \
            target=100 \
            ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        """

rule individual_sample_assembly:
    """Assembly step 02: Individual sample assemblies

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        r1_norm = os.path.join(ASSEMBLY, "{sample}_R1.norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, "{sample}_R2.norm.fastq"),
        r1s = os.path.join(TMPDIR, "p07", "{sample}_R1.singletons.fastq"),
        r2s = os.path.join(TMPDIR, "p07", "{sample}_R2.singletons.fastq")
    output:
        contigs = os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa")
    params:
        mh_dir = directory(os.path.join(ASSEMBLY, '{sample}')),
        contig_dic = directory(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY"))
    benchmark:
        os.path.join(BENCH, "a02.megahit_{sample}.txt")
    log:
        os.path.join(STDERR, "a02.megahit_{sample}.log")
    resources:
        mem_mb = MhitMem
    threads:
        MhitCPU
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        rmdir {params.mh_dir};

        megahit -1 {input.r1_norm} -2 {input.r2_norm} -r {input.r1s},{input.r2s} \
            -o {params.mh_dir} --out-prefix {wildcards.sample} \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log};
        """

rule concatenate_contigs:
    """Assembly step 03: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        expand(os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa"), sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    benchmark:
        os.path.join(BENCH, "a03.concatenate_assemblies.txt")
    log:
        os.path.join(STDERR, "a03.cat.log")
    shell:
        "cat {input} > {output} 2> {log}"

rule contig_reformating_and_stats:
    """Assembly step 04: Remove short contigs (Default: 1000). Defined in config[CONTIG_MINLENGTH]"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    output:
        rename = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.renamed.fasta")),
        size = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.fasta"),
        stats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.stats")
    benchmark:
        os.path.join(BENCH, "a04.contig_reformating.txt")
    log:
        log1 = os.path.join(STDERR, "a04.rename.log"),
        log2 = os.path.join(STDERR, "a04.reformat.log"),
        log3 = os.path.join(STDERR, "a04.stats.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        rename.sh in={input} out={output.rename} \
            prefix=contig_ \
            ow=t \
            -Xmx{resources.mem_mb}m 2> {log.log1};

        reformat.sh in={output.rename} out={output.size} \
            ml={config[CONTIG_MINLENGTH]} \
            ow=t \
            -Xmx{resources.mem_mb}m 2> {log.log2};

        statswrapper.sh in={input} out={output.stats} \
            format=2 \
            ow=t 2> {log.log3};
        """

rule population_assembly:
    """Assembly step 05: Create 'contig dictionary' of all unique contigs present in the study (aka: population assembly)"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.fasta")
    output:
        assembly = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta"),
        stats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary.stats")
    params:
        flye_out = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE")
    benchmark:
        os.path.join(BENCH, "a05.population_assembly.txt")
    log:
        log1 = os.path.join(STDERR, "a05.flye.log"),
        log2 = os.path.join(STDERR, "a05.stats.log")
    resources:
        mem_mb = MhitMem
    threads:
        MhitCPU
    conda:
        "../envs/metaflye.yaml"
    shell:
        """
        flye --subassemblies {input} -t {threads} --plasmids -o {params.flye_out} -g 1g &>> {log.log1};
        statswrapper.sh in={output.assembly} out={output.stats} \
            format=2 \
            ow=t 2> {log.log2};
        """

rule link_assembly:
    """Assembly step 06: Link the final assembly to the results directory; not really a step."""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")
    output:
        os.path.join(RESULTS, "assembly.fasta")
    run:
        os.symlink(os.path.abspath(input[0]), os.path.abspath(output[0]))

rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(TMPDIR, "p07", "{sample}_R1.all.fastq"),
        r2 = os.path.join(TMPDIR, "p07", "{sample}_R2.all.fastq"),
        ref = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")
    output:
        sam = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.aln.sam.gz"),
        unmap = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.unmapped.fastq"),
        covstats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.cov_stats"),
        rpkm = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.rpkm"),
        statsfile = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.statsfile"),
        scafstats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.scafstats")
    benchmark:
        os.path.join(BENCH, "a07.coverage_calculations_{sample}.txt")
    log:
        os.path.join(STDERR, "a07.bbmap_{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbmap.sh ref={input.ref} in={input.r1} in2={input.r2} \
            nodisk \
            out={output.sam} \
            outu={output.unmap} \
            ambiguous=random \
            slow=t \
            physcov=t \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            statsfile={output.statsfile} \
            scafstats={output.scafstats} \
            maxindel=100 minid=90 \
            ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log};
        """

rule create_contig_count_table:
    """Assembly step 08: Transcript Per Million (TPM) calculator

    Useful resource: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/"""
    input:
        rpkm = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.rpkm"),
        covstats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.cov_stats")
    output:
        counts_tmp = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_counts.tmp")),
        TPM_tmp = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_TPM.tmp")),
        TPM = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_TPM")),
        TPM_final = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_TPM.final")),
        cov_temp = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_cov.tmp")),
        count_tbl = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_contig_counts.tsv")
    benchmark:
        os.path.join(BENCH, "a08.tpm_caluclator_{sample}.txt")
    log:
        os.path.join(STDERR, "a08.tmp_calc_{sample}.log")
    shell:
        """
        ## TPM Calculator
        # Prepare table & calculate RPK
        {{ tail -n+6 {input.rpkm} | \
            cut -f1,2,5,6,8 | \
            awk 'BEGIN{{ FS=OFS="\t" }} {{ print $0, $3/($2/1000) }}' > {output.counts_tmp};

        # Calculate size factor (per million)
        sizef=$(awk 'BEGIN{{ total=0 }} {{ total=total+$6/1000000 }} END{{ printf total }}' {output.counts_tmp});

        # Calculate TPM
        awk -v awkvar="$sizef" 'BEGIN{{ FS=OFS="\t" }} {{ print $0, $6/awkvar }}' < {output.counts_tmp} > {output.TPM_tmp};

        # Add sample name
        awk -v awkvar="{wildcards.sample}" 'BEGIN{{FS=OFS="\t"}} {{ print awkvar, $0 }}' < {output.TPM_tmp} > {output.TPM};

        # Remove RPK
        cut --complement -f 7 {output.TPM} > {output.TPM_final};

        ## Coverage stats modifications
        tail -n+2 {input.covstats} | cut -f2,4,5,6,9 > {output.cov_temp};

        ## Combine tables
        paste {output.TPM_final} {output.cov_temp} > {output.count_tbl};
        }} 2> {log}
        """

rule concatentate_contig_count_tables:
    """Assembly step 09: Concatenate contig count tables

    Note: this is done as a separate rule due to how snakemake handles i/o files. It does not work well in Step 20b as 
    the i/o PATTERNS are different.
    """
    input:
        expand(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}_contig_counts.tsv"), sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "contig_count_table.tsv")
    benchmark:
        os.path.join(BENCH, "a09.concat_contig_count_tables.txt")
    log:
        os.path.join(STDERR, "a09.concat_tables.log")
    shell:
        """
        {{ cat {input} > {output};
        sed -i '1i sample_id\tcontig_id\tlength\treads\tRPKM\tFPKM\tTPM\tavg_fold_cov\tcontig_GC\tcov_perc\tcov_bases\tmedian_fold_cov' {output}; \
        }} 2> {log}
        """
