"""
What is accomplished with these rules?
    - Assembly
        - Sample assembly
        - Population assembly
        - Contig abundance esitmation
"""

rule assembly_kmer_normalization:
    """Step 14: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1=os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq"),
        r2=os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq"),
        r1s=os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s=os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        r1_norm=os.path.join(ASSEMBLY, PATTERN_R1 + ".norm.fastq"),
        r2_norm=os.path.join(ASSEMBLY, PATTERN_R2 + ".norm.fastq")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s14.kmer_normalization_{sample}.txt")
    log:
        log=os.path.join(STDERR, "step_14", "s14_{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 15: Individual sample assemblies

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        r1_norm=os.path.join(ASSEMBLY, PATTERN_R1 + ".norm.fastq"),
        r2_norm=os.path.join(ASSEMBLY, PATTERN_R2 + ".norm.fastq"),
        r1s=os.path.join(QC,"HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s=os.path.join(QC,"HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        contigs=os.path.join(ASSEMBLY, PATTERN, PATTERN + ".contigs.fa")
    params:
        mh_dir=directory(os.path.join(ASSEMBLY,'{sample}')),
        contig_dic=directory(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY"))
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s15.megahit_{sample}.txt")
    log:
        log=os.path.join(STDERR, "step_15", "s15_{sample}.log")
    resources:
        mem_mb=MhitMem
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
    """Step 16: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        lambda wildcards: expand(os.path.join(ASSEMBLY, PATTERN, PATTERN + ".contigs.fa"), sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s16.concatenate_assemblies.txt")
    log:
        log=os.path.join(STDERR, "step_16", "s16.log")
    shell:
        "cat {input} > {output}"

rule contig_reformating_and_stats:
    """Step 17: Remove short contigs (Default: 1000). Defined in config[CONTIG_MINLENGTH]"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    output:
        rename=temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.renamed.fasta")),
        size=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.fasta"),
        stats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.stats")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s17.contig_reformating.txt")
    log:
        log1=os.path.join(STDERR, "step_17", "s17.rename.log"),
        log2=os.path.join(STDERR, "step_17", "s17.reformat.log"),
        log3=os.path.join(STDERR, "step_17", "s17.stats.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 18: Create 'contig dictionary' of all unique contigs present in the study (aka: population assembly)"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.fasta")
    output:
        assembly=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta"),
        stats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary.stats")
    params:
        flye_out=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s18.population_assembly.txt")
    log:
        log1=os.path.join(STDERR, "step_18", "s18.flye.log"),
        log2=os.path.join(STDERR, "step_18", "s18.stats.log")
    resources:
        mem_mb=MhitMem
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
    """Link the final assembly to the RESULTS"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY","FLYE","assembly.fasta")
    output:
        os.path.join(RESULTS, "assembly.fasta")
    run:
        os.symlink(os.path.abspath(input[0]), os.path.abspath(output[0]))

rule coverage_calculations:
    """Step 19: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1=os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq"),
        r2=os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".all.fastq"),
        ref=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")
    output:
        sam=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".aln.sam.gz"),
        unmap=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".unmapped.fastq"),
        covstats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".cov_stats"),
        rpkm=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".rpkm"),
        statsfile=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".statsfile"),
        scafstats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".scafstats")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s19.coverage_calculations_{sample}.txt")
    log:
        log=os.path.join(STDERR,"step_16","s16_{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 20: Transcript Per Million (TPM) calculator

    Useful resource: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/"""
    input:
        rpkm=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + ".rpkm"),
        covstats=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + ".cov_stats")
    output:
        counts_tmp=temporary(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_counts.tmp")),
        TPM_tmp=temporary(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_TPM.tmp")),
        TPM=temporary(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_TPM")),
        TPM_final=temporary(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_TPM.final")),
        cov_temp=temporary(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_cov.tmp")),
        count_tbl=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_contig_counts.tsv")
    benchmark:
        os.path.join(BENCH,"PREPROCESSING","s20.tpm_caluclator_{sample}.txt")
    log:
        log=os.path.join(STDERR,"step_20","s20_{sample}.log")
    shell:
        """
        ## TPM Calculator
        # Prepare table & calculate RPK
        tail -n+6 {input.rpkm} | \
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

        """

rule concatentate_contig_count_tables:
    """Step 21: Concatenate contig count tables

    Note: this is done as a separate rule due to how snakemake handles i/o files. It does not work well in Step 20b as 
    the i/o PATTERNS are different.
    """
    input:
        lambda wildcards: expand(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",PATTERN + "_contig_counts.tsv"),
            sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING","contig_count_table.tsv")
    benchmark:
        os.path.join(BENCH,"PREPROCESSING","s21.concat_contig_count_tables.txt")
    log:
        log=os.path.join(STDERR,"step_21","s21.log")
    shell:
        """
        cat {input} > {output};
        sed -i '1i sample_id\tcontig_id\tlength\treads\tRPKM\tFPKM\tTPM\tavg_fold_cov\tcontig_GC\tcov_perc\tcov_bases\tmedian_fold_cov' {output};
        """

# rule calculate_contig_dictionary_properties:
#     """Step 22: Calculate contig sequence properties properties (ie. GC-content, tetramer frequencies) per sequence"""
#     input:
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","assembly.fasta")
#     output:
#         gc=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","contig_dictionary_properties.gc"),
#         tetramer=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","contig_dictionary_properties.tetramer"),
#         seq_properties=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","contig_dictionary_properties.tsv")
#     benchmark:
#         os.path.join(BENCH,"PREPROCESSING","s22.calculate_contig_dictionary_properties.txt")
#     log:
#         log1=os.path.join(LOGS,"step_20","s20.gc.log"),
#         log2=os.path.join(LOGS,"step_20","s20.tetramer.log")
#     resources:
#         mem_mb=100000,
#         cpus=64
#     conda:
#         "../envs/bbmap.yaml"
#     shell:
#         """
#         # Calcualate per sequence GC content
#         countgc.sh in={input} format=2 ow=t | awk 'NF' > {output.gc};
#         sed -i '1i id\tGC' {output.gc} 2> {log.log1};
#
#         # Calculate per sequence tetramer frequency
#         tetramerfreq.sh in={input} w=0 ow=t | \
#         tail -n+2 | \
#         cut --complement -f2 > {output.tetramer} 2> {log.log2};
#
#         sed -i 's/scaffold/id/' {output.tetramer};
#
#         # Combine
#         csvtk join -f 1 {output.gc} {output.tetramer} -t -T > {output.seq_properties};
#
#         """