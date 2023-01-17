"""
Per-sample assemblies for short paired reads
    Take all trimmed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in combine_sample_assemblies.smk
"""

rule assembly_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(dir.out.assembly, "{sample}_R1.unmapped.fastq.gz"),
        r1s = os.path.join(dir.out.assembly, "{sample}_R1.singletons.fastq.gz")
    output:
        r1_norm = temp(os.path.join(dir.out.assembly, "{sample}_R1.norm.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "kmer_normalization_{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "kmer_norm_{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    shell:
        """
        bbnorm.sh in={input.r1} \
            extra={input.r1s} \
            out={output.r1_norm} \
            target=100 \
            ow=t \
            threads={threads} -Xmx{resources.javaAlloc}m 2> {log}
        rm {log}
        """


rule individual_sample_assembly:
    """Assembly step 02: Individual sample assemblies

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        r1_norm = os.path.join(dir.out.assembly, "{sample}_R1.unmapped.fastq.gz"),
        r1s = os.path.join(dir.out.assembly, "{sample}_R1.singletons.fastq.gz"),
    output:
        contigs = os.path.join(dir.out.assembly, "{sample}", "{sample}.contigs.fa"),
        renamed = os.path.join(dir.out.assembly, "{sample}", "{sample}.rename.contigs.fa"),
        tar = os.path.join(dir.out.assembly,'{sample}.tar.zst')
    params:
        mh_dir=lambda w, output: os.path.split(output.contigs)[0],
        minlen=config.qc.contigMinLen
    benchmark:
        os.path.join(dir.out.bench, "megahit_{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "megahit_{sample}.log")
    resources:
        mem_mb = config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "megahit.yaml")
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -r {input.r1_norm},{input.r1s} \
            -o {params.mh_dir} --out-prefix {wildcards.sample} --min-contig-len {params.minlen} \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log}
        sed 's/>/>{wildcards.sample}/' {output.contigs} > {output.renamed}
        tar cf - {params.mh_dir} | zstd -T8 -9 > {output.tar} 2> {log}
        rm {log}
        """


rule mapSampleAssemblyPairedReads:
    """Map the sample paired reads to the sample assembly"""
    input:
        r1 = os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq.gz"),
        contigs = os.path.join(dir.out.assembly, "{sample}", "{sample}.contigs.fa"),
    output:
        temp(os.path.join(dir.out.assembly, '{sample}', '{sample}.se.bam'))
    conda:
        os.path.join(dir.env, 'minimap2.yaml')
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem
    log:
        os.path.join(dir.out.stderr, 'sampleAssemblyMapSe.{sample}.log')
    benchmark:
        os.path.join(dir.out.bench, 'sampleAssemblyMapSe.{sample}.txt')
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1} | \
            samtools sort -n -o {output[0]};
        }} 2> {log}
        rm {log}
        """


rule mapSampleAssemblyUnpairedReads:
    """Map the sample unpaired reads to the sample assembly"""
    input:
        r1s = os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq.gz"),
        contigs = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.fa")
    output:
        temp(os.path.join(dir.out.assembly,'{sample}','{sample}.assemblyUnmapped.s.fastq'))
    conda:
        os.path.join(dir.env, 'minimap2.yaml')
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem
    log:
        os.path.join(dir.out.stderr, 'sampleAssemblyMapS.{sample}.log')
    benchmark:
        os.path.join(dir.out.bench, 'sampleAssemblyMapS.{sample}.txt')
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1s} | \
            samtools sort -n | \
            samtools fastq -f 4 > {output[0]};
        }} 2> {log}
        rm {log}
        """


rule pullPairedUnmappedReads:
    """Grab the paired unmapped reads (neither pair mapped)"""
    input:
        os.path.join(dir.out.assembly,'{sample}','{sample}.se.bam')
    output:
        r1 = temp(os.path.join(dir.out.assembly, '{sample}', '{sample}.assemblyUnmapped_R1.fastq')),
    conda:
        os.path.join(dir.env,   'samtools.yaml')
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem
    log:
        os.path.join(dir.out.stderr, 'pullPairedUnmappedReads.{sample}.log')
    benchmark:
        os.path.join(dir.out.bench, 'pullPairedUnmappedReads.{sample}.txt')
    shell:
        """
        samtools fastq -f 77 {input} > {output.r1} 2> {log}
        rm {log}
        """


rule pullPairedUnmappedReadsMateMapped:
    """Grab the paired unmapped reads (mate is mapped)"""
    input:
        os.path.join(dir.out.assembly,'{sample}','{sample}.se.bam')
    output:
        temp(os.path.join(dir.out.assembly, '{sample}', '{sample}.assemblyUnmapped.se.s.fastq'))
    conda:
        os.path.join(dir.env,   'samtools.yaml')
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem
    log:
        os.path.join(dir.out.stderr, 'pullPairedUnmappedReads.{sample}.log')
    benchmark:
        os.path.join(dir.out.bench, 'pullPairedUnmappedReads.{sample}.txt')
    shell:
        """
        samtools fastq -f5 -F8 {input[0]} > {output[0]} 2> {log}
        rm {log}
        """


rule poolR1Unmapped:
    """Concatenate the unmapped, paired R1 reads for all samples"""
    input:
        expand(os.path.join(dir.out.assembly, '{sample}', '{sample}.assemblyUnmapped_R1.fastq'), sample=samples.names)
    output:
        temp(os.path.join(dir.out.assembly, 'unmapRescue_R1.fastq'))
    shell:
        """cat {input} > {output}"""


rule poolUnpairedUnmapped:
    """Concatenate the unmapped, unpaired reads for all samples"""
    input:
        fq = expand(os.path.join(dir.out.assembly,'{sample}','{sample}.assemblyUnmapped.{sSe}.fastq'), sample=samples.names, sSe = ['s','se.s']),
    output:
        temp(os.path.join(dir.out.assembly, 'unmapRescue.s.fastq'))
    shell:
        """
        cat {input.fq} > {output}
        """


rule rescue_read_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(dir.out.assembly, 'unmapRescue_R1.fastq'),
        s = os.path.join(dir.out.assembly, 'unmapRescue.s.fastq')
    output:
        r1_norm = temp(os.path.join(dir.out.assembly, 'unmapRescueNorm_R1.fastq'))
    benchmark:
        os.path.join(dir.out.bench, "rescue_read_kmer_normalization.txt")
    log:
        os.path.join(dir.out.stderr, "rescue_read_kmer_normalization.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    shell:
        """
        bbnorm.sh in={input.r1} \
            extra={input.s} \
            out={output.r1_norm} \
            target=100 \
            ow=t \
            threads={threads} -Xmx{resources.javaAlloc}m 2> {log}
        rm {log}
        """


rule unmapped_read_rescue_assembly:
    """Assemble the unmapped reads from all samples

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        r1_norm = os.path.join(dir.out.assembly, 'unmapRescue_R1.fastq'),
        s = os.path.join(dir.out.assembly, 'unmapRescue.s.fastq')
    output:
        contigs = os.path.join(dir.out.assembly, 'rescue', "rescue.contigs.fa"),
        renamed = os.path.join(dir.out.assembly, 'rescue', 'rescue.rename.contigs.fa')
    params:
        mh_dir=lambda w, output: os.path.split(output.contigs)[0],
        minlen=config.qc.contigMinLen
    benchmark:
        os.path.join(dir.out.bench, "unmapped_read_rescue_assembly.txt")
    log:
        os.path.join(dir.out.stderr, "unmapped_read_rescue_assembly.log")
    resources:
        mem_mb = config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "megahit.yaml")
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -r {input.r1_norm},{input.s} \
            -o {params.mh_dir} --out-prefix rescue --min-contig-len {params.minlen} \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log}
        sed 's/>/>rescue/' {output.contigs} > {output.renamed}
        rm {log}
        """


rule concatenate_contigs:
    """Assembly step 03: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        expand(os.path.join(dir.out.assembly, "{sample}", "{sample}.rename.contigs.fa"), sample=samples.names),
        os.path.join(dir.out.assembly,'rescue',"rescue.rename.contigs.fa")
    output:
        temp(os.path.join(dir.out.assembly, "all_sample_contigs.fasta"))
    params:
        expand(os.path.join(dir.out.assembly,'{sample}'), sample=samples.names),
    benchmark:
        os.path.join(dir.out.bench, "concatenate_assemblies.txt")
    log:
        os.path.join(dir.out.stderr, "concatenate_assemblies.log")
    shell:
        """
        cat {input} > {output} 2> {log} && rm {log}
        rm -rf {params}
        """


rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(dir.out.assembly, "{sample}_R1.all.fastq.gz"),
        ref = os.path.join(dir.out.results, "assembly.fasta")
    output:
        sam = temp(os.path.join(dir.out.mapping, "{sample}.aln.sam.gz")),
        unmap = temp(os.path.join(dir.out.mapping, "{sample}.unmapped.fastq")),
        covstats = temp(os.path.join(dir.out.mapping, "{sample}.cov_stats")),
        rpkm = temp(os.path.join(dir.out.mapping, "{sample}.rpkm")),
        statsfile = temp(os.path.join(dir.out.mapping, "{sample}.statsfile")),
        scafstats = temp(os.path.join(dir.out.mapping, "{sample}.scafstats"))
    benchmark:
        os.path.join(dir.out.bench, "coverage_calculations.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "coverage_calculations.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    shell:
        """
        bbmap.sh ref={input.ref} in={input.r1} \
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
            threads={threads} -Xmx{resources.javaAlloc}m 2> {log}
        rm {log}
        """