"""
Per-sample assemblies for short paired reads
    Take all trimmed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in combine_sample_assemblies.smk
"""


rule co_assembly:
    """Alternative to cross assembly; assemble everything together in one hit.

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        expand(
            os.path.join(dir.out.assembly,"{sample}_R1.{unmapsingle}.fastq.gz"),
            sample=samples.names,
            unmapsingle=["unmapped", "singletons"]
        )
    output:
        assembly=os.path.join(dir.out.results,"co_assembly.fasta"),
        graph=os.path.join(dir.out.results,"co_assembly_graph.gfa"),
        tmp=temp(os.path.join(dir.out.assembly,"coAssembly","co_assembly_graph.fastg")),
        tar=os.path.join(dir.out.assembly,"coAssembly.tar.zst")
    params:
        r=','.join(expand(
                os.path.join(dir.out.assembly,"{sample}_R1.{unmapsingle}.fastq.gz"),
                sample=samples.names,
                unmapsingle=["unmapped","singletons"]
            )),
        mh_dir=os.path.join(dir.out.assembly,"coAssembly"),
        mh_int=os.path.join(dir.out.assembly,"coAssembly","intermediate_contigs"),
        params=config.assembly.megahit,
        assembly=os.path.join(dir.out.assembly,"coAssembly","assembly.fasta"),
        graph=os.path.join(dir.out.assembly,"coAssembly","assembly_graph.gfa"),
    benchmark:
        os.path.join(dir.out.bench,"megahit_coassembly.txt")
    log:
        os.path.join(dir.out.stderr,"megahit_coassembly.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env,"megahit.yaml")
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -r {params.r} -o {params.mh_dir} -t {threads} {params.params} &> {log}
        kctg=$(ls -t {params.mh_int}/*.contigs.fa | grep -v final | head -1)
        kmax=$([[ $kctg =~ ([0-9]+) ]] && echo "${{BASH_REMATCH[1]}}")
        megahit_toolkit contig2fastg $kmax {params.mh_int}/$kctg > {output.tmp}
        Bandage reduce {output.tmp} {output.graph}
        cp {params.assembly} {output.assembly}
        tar cf - {params.mh_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
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
        tar = os.path.join(dir.out.assembly,"{sample}.tar.zst")
    params:
        mh_dir=lambda w, output: os.path.split(output.contigs)[0],
        params=config.assembly.megahit
    benchmark:
        os.path.join(dir.out.bench, "megahit_{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "megahit_{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "megahit.yaml")
    group:
        "assembly"
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -r {input.r1_norm},{input.r1s} \
            -o {params.mh_dir} --out-prefix {wildcards.sample} -t {threads} \
            {params.params} &>> {log}
        sed 's/>/>{wildcards.sample}/' {output.contigs} > {output.renamed}
        tar cf - {params.mh_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        rm {log}
        """


rule mapSampleAssemblyPairedReads:
    """Map the sample paired reads to the sample assembly"""
    input:
        r1 = os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq.gz"),
        contigs = os.path.join(dir.out.assembly, "{sample}", "{sample}.contigs.fa"),
    output:
        temp(os.path.join(dir.out.assembly, "{sample}", "{sample}.se.bam"))
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    log:
        os.path.join(dir.out.stderr, "sampleAssemblyMapSe.{sample}.log")
    benchmark:
        os.path.join(dir.out.bench, "sampleAssemblyMapSe.{sample}.txt")
    group:
        "assembly"
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1} | \
            samtools sort -n -o {output};
        }} 2> {log}
        rm {log}
        """


rule mapSampleAssemblyUnpairedReads:
    """Map the sample unpaired reads to the sample assembly"""
    input:
        r1s = os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq.gz"),
        contigs = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.fa")
    output:
        temp(os.path.join(dir.out.assembly,"{sample}","{sample}.assemblyUnmapped.s.fastq"))
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    benchmark:
        os.path.join(dir.out.bench, "sampleAssemblyMapS.{sample}.txt")
    group:
        "assembly"
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1s} | \
            samtools sort -n | \
            samtools fastq -f 4 > {output};
        }}
        """


rule pullPairedUnmappedReads:
    """Grab the paired unmapped reads (neither pair mapped)"""
    input:
        os.path.join(dir.out.assembly,"{sample}","{sample}.se.bam")
    output:
        temp(os.path.join(dir.out.assembly, "{sample}", "{sample}.assemblyUnmapped_R1.fastq")),
    conda:
        os.path.join(dir.env,   "samtools.yaml")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "pullPairedUnmappedReads.{sample}.txt")
    group:
        "assembly"
    shell:
        """
        samtools fastq -f 77 {input} > {output}
        """


rule pullPairedUnmappedReadsMateMapped:
    """Grab the paired unmapped reads (mate is mapped)"""
    input:
        os.path.join(dir.out.assembly,"{sample}","{sample}.se.bam")
    output:
        temp(os.path.join(dir.out.assembly, "{sample}", "{sample}.assemblyUnmapped.se.s.fastq"))
    conda:
        os.path.join(dir.env,   "samtools.yaml")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "pullPairedUnmappedReads.{sample}.txt")
    group:
        "assembly"
    shell:
        """
        samtools fastq -f5 -F8 {input} > {output}
        """


rule poolR1Unmapped:
    """Concatenate the unmapped, paired R1 reads for all samples"""
    input:
        expand(os.path.join(dir.out.assembly, "{sample}", "{sample}.assemblyUnmapped_R1.fastq"), sample=samples.names)
    output:
        os.path.join(dir.out.assembly, "rescue_R1.unmapped.fastq.gz")
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    conda:
        os.path.join(dir.env, "pigz.yaml")
    group:
        "assemblyRescue"
    shell:
        """
        cat {input} | pigz -p {threads} -c > {output}
        """


rule poolUnpairedUnmapped:
    """Concatenate the unmapped, unpaired reads for all samples"""
    input:
        fq = expand(os.path.join(dir.out.assembly,"{sample}","{sample}.assemblyUnmapped.{sSe}.fastq"), sample=samples.names, sSe = ["s","se.s"]),
    output:
        temp(os.path.join(dir.out.assembly, "rescue_R1.singletons.fastq.gz"))
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    conda:
        os.path.join(dir.env, "pigz.yaml")
    group:
        "assemblyRescue"
    shell:
        """
        cat {input} | pigz -p {threads} -c > {output}
        """


rule concatenate_contigs:
    """Assembly step 03: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        expand(os.path.join(dir.out.assembly, "{sample}", "{sample}.rename.contigs.fa"), sample=samples.names),
        os.path.join(dir.out.assembly,"rescue","rescue.rename.contigs.fa")
    output:
        os.path.join(dir.out.assembly,"all_sample_contigs.fasta.gz")
    params:
        dirs = expand(os.path.join(dir.out.assembly,"{sample}"), sample=samples.names + ["rescue"]),
        compression = '-' + str(config.qc.compression)
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    conda:
        os.path.join(dir.env, "pigz.yaml")
    group:
        "assembly"
    shell:
        """
        cat {input} | pigz -p {threads} {params.compression} -c > {output}
        rm -rf {params.dirs}
        """


rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(dir.out.assembly, "{sample}_R1.all.fastq.gz"),
        ref = os.path.join(dir.out.results, f"{config.args.assembly}_assembly.fasta")
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
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = config.resources.med.time
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
