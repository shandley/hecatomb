
# load default config
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
res = config["resources"]
config = config["hecatomb"]


# directories
include: os.path.join("rules", "preflight", "directories.smk")


# host files
dir["dbs"]["host"]["fasta"] = os.path.join(
    dir["dbs"]["host"]["base"], config["args"]["hostName"], "masked_ref.fa.gz"
)
dir["dbs"]["host"]["virRefSeq"] = os.path.join(dir["dbs"]["host"]["base"], "viral.1.1.genomic.fna.gz")


rule all:
    input:
        dir["dbs"]["host"]["fasta"]


rule dl_refseq_viral:
    """Download the refseq viral genomic file needed"""
    output:
        dir["dbs"]["host"]["virRefSeq"]
    params:
        url = config["dbvirRefSeq"]["url"],
        file = config["dbvirRefSeq"]["file"]
    conda:
        os.path.join(dir["env"], "curl.yaml")
    shell:
        "curl {params.url} -o {params.file}"


rule minimap_viral_refseq:
    """Map the viral refseq db to the host"""
    input:
        ref = config["args"]["hostFa"],
        vir = dir["dbs"]["host"]["virRefSeq"]
    output:
        temp(os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.bed'))
    params:
        config["addHost"]["minViralAlnLen"]
    conda:
        os.path.join(dir["env"], "minimap.yaml")
    resources:
        mem_mb = res["med"]["mem"],
        mem = str(res["med"]["mem"]) + "MB",
        time = res["med"]["time"]
    threads:
        res["med"]["cpu"]
    log:
        os.path.join(dir["out"]["stderr"], 'minimap_viral_refseq.log')
    benchmark:
        os.path.join(dir["out"]["bench"], "minimap_viral_refseq.log")
    shell:
        """
        {{ minimap2 -t {threads} {input.ref} {input.vir}  \
            | awk '$11>{params}' \
            | cut -f6,8,9 \
            | sort -k1,1 -k2,2n -k3,3n; \
        }} > {output} 2> {log}
        """


rule mask_fasta:
    """Mask the host genome using bedtools"""
    input:
        fa = config["args"]["hostFa"],
        aln = os.path.join(dir["out"]["temp"],f'{config["args"]["hostName"]}.bed')
    output:
        mask = (dir["out"]["temp"],f'{config["args"]["hostName"]}.mask.bed'),
        fa = dir["dbs"]["host"]["fasta"]
    conda:
        os.path.join(dir["env"], "bedtools.yaml")
    log:
        os.path.join(dir["out"]["stderr"], "mask_fasta.log")
    benchmark:
        os.path.join(dir["out"]["bench"], "mask_fasta.log")
    shell:
        """
        bedtools merge -i {input.aln} > {output.mask}
        bedtools maskfasta -fi {input.fa} -bed {output.mask} \
            | gzip -1 - > {output.mask}
        """


# rule map_shred_seq:
#     """Add host step 01: Map the virus shred seqs to the input host"""
#     input:
#         ref = config["args"]["hostFa"],
#         shred = dir["dbs"]["host"]["virShred"]
#     output:
#         temp(os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.sam.gz'))
#     conda:
#         os.path.join(dir["env"], "bbmap.yaml")
#     resources:
#         mem_mb = config["resources"]["med"]["mem"],
#     threads:
#         config["resources"]["med"]["cpu"]
#     log:
#         os.path.join(dir["out"]["stderr"], 'map_shred_seq.log')
#     shell:
#         """
#         bbmap.sh ref={input.ref} in={input.shred} \
#             outm={output} path=tmp/ \
#             minid=0.90 maxindel=2 ow=t \
#             threads={threads} -Xmx{resources.mem_mb}m &> {log}
#         """
#
#
# rule mask_host:
#     """Add host step 02: Mask the host"""
#     input:
#         ref = config["args"]["hostFa"],
#         sam = os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.sam.gz')
#     output:
#         fa = temp(os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.processed.fasta')),
#         gz = dir["dbs"]["host"]["fasta"]
#     conda:
#         os.path.join(dir.env, "bbmap.yaml")
#     resources:
#         mem_mb = config["resources"]["med"]["mem"],
#     threads:
#         config["resources"]["med"]["cpu"]
#     params:
#         entropy = config["qc"]["entropy"]
#     log:
#         os.path.join(dir["out"]["stderr"], 'mask_host.log')
#     shell:
#         """
#         bbmask.sh in={input.ref} out={output.fa} \
#             entropy={params.entropy} sam={input.sam} ow=t \
#             threads={threads} -Xmx{resources.mem_mb}m &> {log}
#         gzip -c {output.fa} > {output.gz}
#         """
