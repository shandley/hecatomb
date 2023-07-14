"""
Per-sample assemblies for longreads
    Take all host-removed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in coverage.smk
"""
import os

rule lr_cross_assembly:
    """Alternative to merged assembly; assemble everything together in one hit."""
    input:
        targets.trimnami
    output:
        assembly = os.path.join(dir.out.results, "cross_assembly.fasta"),
        graph = os.path.join(dir.out.results, "cross_assembly.gfa"),
        tar = os.path.join(dir.out.assembly,"crossAssembly.tar.zst")
    params:
        params = config.assembly.metaflye,
        dir = os.path.join(dir.out.assembly, "crossAssembly"),
        assembly = os.path.join(dir.out.assembly, "crossAssembly", "assembly.fasta"),
        graph = os.path.join(dir.out.assembly, "crossAssembly", "assembly_graph.gfa"),
    resources:
        mem_mb = resources.big.mem,
        mem = resources.big.mem + "MB",
        time = resources.big.time
    threads:
        resources.big.cpu
    log:
        os.path.join(dir.out.stderr, "canu_cross_assembly.log")
    benchmark:
        os.path.join(dir.out.bench, "canu_cross_assembly.txt")
    conda:
        os.path.join(dir.env, "metaflye.yaml")
    shell:
        """
        flye -o {params.dir} -t {threads} {params.params} {input} 2> {log}
        mv {params.assembly} {output.assembly}
        mv {params.graph} {output.graph}
        tar cf - {params.dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        rm {log}
        """

# rule cross_assembly_minimap:
#     """Run all-v-all minimap2 on longreads for minisam"""
#     input:
#         targets.trimnami
#     output:
#         reads = temp(os.path.join(dir.out.assembly, "reads.fastq.gz")),
#         paf = temp(os.path.join(dir.out.assembly, "reads.paf.gz"))
#     resources:
#         mem_mb = resources.med.mem,
#         mem = resources.med.mem + "MB",
#         time = resources.med.time
#     threads:
#         resources.med.cpu
#     log:
#         os.path.join(dir.out.stderr, "cross_assembly_minimap.log")
#     benchmark:
#         os.path.join(dir.out.bench, "cross_assembly_minimap.txt")
#     conda:
#         os.path.join(dir.env, "minimap2.yaml")
#     shell:
#         """
#         cat {input} > {output.reads}
#
#         minimap2 \
#             -x ava-ont \
#             -t{threads} \
#             {output.reads} \
#             {output.reads} \
#             2> {log} \
#             | gzip -1 > {output.paf}
#         """
#
#
# rule cross_assembly_miniasm:
#     """Run miniasm on the all-v-all aligned reads"""
#     input:
#         reads = os.path.join(dir.out.assembly, "reads.fastq.gz"),
#         paf = os.path.join(dir.out.assembly, "reads.paf.gz")
#     output:
#         graph = os.path.join(dir.out.results,"cross_assembly.gfa"),
#         unitigs = os.path.join(dir.out.results, "cross_assembly.fasta"),
#     resources:
#         mem_mb = resources.med.mem,
#         mem = resources.med.mem + "MB",
#         time = resources.med.time
#     threads:
#         resources.med.cpu
#     log:
#         os.path.join(dir.out.stderr, "cross_assembly_miniasm.log")
#     benchmark:
#         os.path.join(dir.out.bench, "cross_assembly_miniasm.txt")
#     conda:
#         os.path.join(dir.env, "miniasm.yaml")
#     shell:
#         """
#         miniasm \
#             -f {input.reads} \
#             {input.paf} \
#             2> {log} \
#             > {output.graph}
#
#         awk \
#             '/^S/{{print ">"$2"\\n"$3}}' \
#             {output.graph} \
#             > {output.unitigs} \
#             2>> {log}
#         """


rule canu_sample_assembly:
    """Per-sample assembly with canu; also works for unmapped rescue reads"""
    input:
        os.path.join(dir.out.trim, "{sample}" + config.args.hostStr + ".single.fastq.gz")
    output:
        ctg = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.fasta"),
        ctgq = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.uniq.fasta"),
        un = os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.fasta"),
        unq = os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.uniq.fasta"),
        tar = os.path.join(dir.out.assembly,"{sample}.tar.zst")
    params:
        settings = config.assembly.canu,
        canu_dir = lambda w, output: os.path.split(output.ctg)[0]
    resources:
        mem_mb = resources.med.mem,
        mem = resources.med.mem + "MB",
        time = resources.med.time
    threads:
        resources.med.cpu
    log:
        os.path.join(dir.out.stderr, "canu_sample_assembly.{sample}.log")
    conda:
        os.path.join(dir.env, "canu.yaml")
    shell:
        """
        canu {params.settings} {input} \
            batThreads={threads} \
            batMemory={resources.mem_mb}M \
            -p {wildcards.sample} \
            -d {params.canu_dir} \
            &>> {log}
        sed 's/>tig/>{wildcards.sample}./' {output.ctg} > {output.ctgq}
        sed 's/>tig/>{wildcards.sample}./' {output.un} > {output.unq}
        tar cf - {params.canu_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        rm {log}
        """


rule combine_canu_unassembled:
    """Combine the unassembled reads from all canu assemblies"""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.uniq.fasta"), sample=samples.names)
    output:
        temp(os.path.join(dir.out.assembly,"unmappedRescue_R1.all.fasta.gz"))
    resources:
        time = resources.sml.time
    group:
        "assembly"
    shell:
        """cat {input} > {output}"""


rule combine_canu_contigs:
    """Combine contigs from all samples plus unmapped rescue assembly"""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.uniq.fasta"), sample=samples.names)
    output:
        os.path.join(dir.out.assembly,"all_sample_contigs.fasta.gz")
    threads:
        resources.med.cpu
    resources:
        time = resources.sml.time
    params:
        compression = "-" + str(config.qc.compression)
    conda:
        os.path.join(dir.env, "pigz.yaml")
    group:
        "assembly"
    shell:
        """cat {input} | pigz -p {threads} {params.compression} -c > {output}"""
