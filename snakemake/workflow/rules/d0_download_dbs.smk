"""
Rules for downloading and installing the necessary databases. Imported by DownloadDB.smk
"""

rule download_hecatomb_db:
    """
    Download and unpack the hecatomb databases
    """
    output:
        os.path.join(BACPATH,config['bacteria']),
        os.path.join(CONPATH,config['contaminants']),
        os.path.join(HOSTPATH,config['human']),
        os.path.join(CONPATH,"nebnext_adapters.fa"),
        os.path.join(CONPATH,"primerB.fa"),
        os.path.join(CONPATH,"rc_primerB_ad6.fa")
    conda:
        "../envs/curl.yaml"
    params:
        wd=DBDIR,
        tar=config["url"]["hecatomb"]["tar"],
        md5=config["url"]["hecatomb"]["md5"],
        filename=config["url"]["hecatomb"]["filename"]
    shell:
        """
        cd {params.wd}; \
        curl -Lgo {params.filename} '{params.tar}'; \
        curl -Lgo {params.filename}.md5 '{params.md5}'; \
        md5sum -c {params.filename}.md5; \
        tar -Izstd -xf {params.filename}; \
        rm {params.filename}.md5
        """

rule download_nucleotide_databases:
    """
    Download and unpack the hecatomb nucleotide databases
    """
    output:
        os.path.join(NUCLPATH,"refseq_virus_nt_UniVec_masked/nt.fnaDB.dbtype"),
        os.path.join(NUCLPATH,"bac_virus_masked/nt.fnaDB.dbtype"),
        os.path.join(NUCLPATH,"refseq_virus_nt_UniVec_masked/nt.fnaDB.index"),
        os.path.join(NUCLPATH,"bac_virus_masked/nt.fnaDB.index")
    conda:
        "../envs/curl.yaml"
    params:
        wd=DBDIR,
        filename=config["url"]["hecatomb_nucl"]["filename"],
        tar=config["url"]["hecatomb_nucl"]["tar"],
        md5=config["url"]["hecatomb_nucl"]["md5"]
    shell:
        """
        cd {params.wd}; \
        curl -Lgo {params.filename} "{params.tar}"; \
        curl -Lgo {params.filename}.md5 "{params.md5}"; \
        md5sum -c {params.filename}.md5; \
        rm {params.filename}.md5; \
        tar -Izstd -xf {params.filename}
        """

# rule extract_nucleotide_databases:
#     """
#     Extract the nucleotide databases
#     """
#     input:
#         os.path.join(DBDIR, "hecatomb.nucleotide.databases.tar.bz2")
#     output:
#         os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked/nt.fnaDB.dbtype"),
#         os.path.join(NUCLPATH, "bac_virus_masked/nt.fnaDB.dbtype"),
#         os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked/nt.fnaDB.index"),
#         os.path.join(NUCLPATH, "bac_virus_masked/nt.fnaDB.index"),
#     shell:
#         "tar -C {DBDIR} -xf {input}"

# rule make_bac_databases:
#     """
#     Generate bacteria bbmap index
#     """
#     input:
#         os.path.join(BACPATH, config['DatabaseFiles']['bacteria'])
#     output:
#         directory(os.path.join(BACPATH, "ref"))
#     params:
#         wd = BACPATH,
#         fa = config['DatabaseFiles']['bacteria']
#     resources:
#         mem_mb=50000,
#         cpus=16
#     conda:
#         "../envs/bbmap.yaml"
#     shell:
#         "cd {params.wd} && bbmap.sh -Xmx{resources.mem_mb}m threads={resources.cpus} ref={params.fa}"

# rule make_bac_bt_idx:
#     """
#     Generate bacteria bowtie2 index
#     """
#     input:
#         os.path.join(BACPATH,config['bacteria'])
#     output:
#         expand(os.path.join(BACPATH,"bac_uniquespecies_giant.masked_Ns_removed.{n}.bt2l"),n=[1, 2, 3, 4])
#     benchmark:
#         "benchmarks/make_bac_bt_idx.txt"
#     resources:
#         mem_mb=50000,
#         cpus=16
#     params:
#         wd=BACPATH,
#         fa=config['bacteria'],
#         bt="bac_uniquespecies_giant.masked_Ns_removed"
#     conda:
#         "../envs/bowtie2.yaml"
#     shell:
#         "cd {params.wd} && bowtie2-build --threads {resources.cpus} --large-index {params.fa} {params.bt}"

# rule make_host_databases:
#     """
#     Generate host bbmap index
#     """
#     input:
#         os.path.join(HOSTPATH, config['DatabaseFiles']['host'])
#     output:
#         directory(os.path.join(HOSTPATH, "ref"))
#     params:
#         wd = HOSTPATH,
#         fa = config['DatabaseFiles']['host']
#     resources:
#         mem_mb=50000,
#         cpus=16
#     conda:
#         "../envs/bbmap.yaml"
#     shell:
#         "cd {params.wd} && bbmap.sh -Xmx{resources.mem_mb}m threads={resources.cpus} ref={params.fa}"

# rule make_host_bt_idx:
#     """
#     Generate host bowtie2 index
#     """
#     input:
#         os.path.join(HOSTPATH,config['human'])
#     output:
#         expand(os.path.join(HOSTPATH,"human_virus_masked.{n}.bt2l"),n=[1, 2, 3, 4])
#     benchmark:
#         "benchmarks/make_host_bt_idx.txt"
#     resources:
#         mem_mb=50000,
#         cpus=16
#     params:
#         wd=HOSTPATH,
#         btidx="human_virus_masked",
#         fa=config['human']
#     conda:
#         "../envs/bowtie2.yaml"
#     shell:
#         "cd {params.wd} && bowtie2-build --threads {resources.cpus} --large-index {params.fa} {params.btidx}"

# REMOVED: Download the processed mapping file instead
# rule download_id_taxonomy_mapping:
#     output:
#         os.path.join(TAXPATH, "idmapping.dat.gz")
#     conda: "envs/curl.yaml"
#     shell:
#         """
#         cd {TAXPATH};
#         curl -LO "{id_mapping_url}"
#         """
#
# rule uniprot_to_ncbi_mapping:
#     input:
#         os.path.join(TAXPATH, "idmapping.dat.gz")
#     output:
#         os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat")
#     shell:
#         """
#         zcat {input} | awk '$2 == "NCBI_TaxID" {{print $1"\t"$3 }}' > {output}
#         """

rule download_uniprot_ncbi_mapping:
    """
    Download and unpack the processed uniprot ncbi mappings
    """
    output:
        os.path.join(TAXPATH,"uniprot_ncbi_mapping.dat")
    params:
        wd=TAXPATH,
        filename=config["url"]["id_map"]["filename"],
        zst=config["url"]["id_map"]["zst"],
        md5=config["url"]["id_map"]["md5"]
    shell:
        """
        cd {params.wd}; \
        curl -Lo {params.filename} "{params.zst}"; \
        curl -Lo {params.filename}.md5 "{params.md5}"; \
        md5sum -c {params.filename}.md5; \
        zstd -d {params.filename} --rm; \
        rm {params.filename}.md5
        """

rule download_ncbi_taxonomy:
    """
    Download taxdump
    """
    output:
        temp(os.path.join(TAXPATH,"taxdump.tar.gz"))
    conda: "../envs/curl.yaml"
    params:
        wd=TAXPATH,
        filename=config["url"]["taxdump"]["filename"],
        tar=config["url"]["taxdump"]["tar"],
        md5=config["url"]["taxdump"]["md5"]
    shell:
        """
        cd {params.wd}; \
        curl -Lo {params.filename} "{params.tar}"; \
        curl -Lo {params.filename}.md5 "{params.md5}"; \
        md5sum -c {params.filename}.md5; \
        rm {params.filename}.md5
        """

rule extract_ncbi_taxonomy:
    """
    Unpack taxdump
    """
    input:
        os.path.join(TAXPATH,"taxdump.tar.gz")
    output:
        os.path.join(TAXPATH,"citations.dmp"),
        os.path.join(TAXPATH,"delnodes.dmp"),
        os.path.join(TAXPATH,"division.dmp"),
        os.path.join(TAXPATH,"gc.prt"),
        os.path.join(TAXPATH,"gencode.dmp"),
        os.path.join(TAXPATH,"merged.dmp"),
        os.path.join(TAXPATH,"readme.txt"),
        os.path.join(TAXPATH,"names.dmp"),
        os.path.join(TAXPATH,"nodes.dmp"),
    params:
        wd=TAXPATH
    shell:
        "cd {params.wd} && tar xf taxdump.tar.gz"

rule download_accession_to_tax:
    """
    Download accession to tax id
    """
    output:
        os.path.join(TAXPATH,"nucl_gb.accession2taxid.gz")
    conda: "../envs/curl.yaml"
    params:
        wd=TAXPATH,
        filename=config["url"]["ntacc2tax"]["filename"],
        gz=config["url"]["ntacc2tax"]["gz"],
        md5=config["url"]["ntacc2tax"]["md5"]
    shell:
        """
        cd {params.wd}; \
        curl -Lgo {params.filename} "{params.gz}"; \
        curl -Lgo {params.filename}.md5 "{params.md5}"; \
        md5sum -c {params.filename}.md5; \
        rm {params.filename}.md5
        """

rule extract_accession_to_tax:
    """
    Unzip accession to tax id
    """
    input:
        os.path.join(TAXPATH,"nucl_gb.accession2taxid.gz")
    output:
        os.path.join(TAXPATH,"nucl_gb.accession2taxid")
    conda:
        "../envs/pigz.yaml"
    shell:
        "unpigz {input}"

# rule create_nt_tax_table:
#     input:
#         nt = os.path.join(NUCLPATH, "nt.fna"),
#         tx = os.path.join(TAXPATH, "nucl_gb.accession2taxid")
#     output:
#         tx = os.path.join(NUCLPATH, "nt.tax")
#     shell:
#         """
#         grep '^>' {input.nt} | cut -f 1 -d ' ' | sed -e 's/^>//' | \
#         xargs -n 10 | sed -e 's/ /|/g' | \
#         xargs -i egrep {} {input.tx} | cut -f 2,3 > {output.tx}
#         """

# REMOVED: instead download the pre-processed uniprot virus db
# rule download_uniprot_viruses:
#     output:
#         os.path.join(PROTPATH, "uniprot_virus.faa")
#     conda: "envs/curl.yaml"
#     shell:
#         """
#         mkdir -p {PROTPATH} && curl -Lgo {output} "{uniprot_virus_url}"
#         """
#
# rule cluster_uniprot:
#     input:
#         os.path.join(PROTPATH, "uniprot_virus.faa")
#     output:
#         db = os.path.join(PROTPATH, "uniprot_virus_c99.faa"),
#         cl = os.path.join(PROTPATH, "uniprot_virus_c99.faa.clstr")
#     benchmark:
#         "benchmarks/cluster_uniprot.txt"
#     resources:
#         mem_mb=20000,
#         cpus=8
#     conda:
#         "envs/cdhit.yaml"
#     shell:
#         "cd-hit -i {input} -o {output.db} -d 0 -c 0.99 -M 0 -T 0"
#
# rule mmseqs_uniprot_clusters:
#     """
#     REQUIRED
#     """
#     input:
#         os.path.join(PROTPATH, "uniprot_virus_c99.faa")
#     output:
#         os.path.join(PROTPATH, "uniprot_virus_c99.db")
#     benchmark:
#         "benchmarks/mmseqs_uniprot_clusters.txt"
#     resources:
#         mem_mb=20000,
#         cpus=8
#     conda:
#         "envs/mmseqs2.yaml"
#     shell:
#         "mmseqs createdb {input} {output}"

rule download_uniprot_clusters:
    """
    Download and unpack the uniprot_virus clusters
    """
    output:
        os.path.join(PROTPATH,"uniprot_virus_c99.db"),
        os.path.join(PROTPATH,"uniprot_virus_c99.faa")
    conda:
        "../envs/curl.yaml"
    params:
        wd=PROTPATH,
        filename=config["url"]["uniprot_virus"]["filename"],
        tar=config["url"]["uniprot_virus"]["tar"],
        md5=config["url"]["uniprot_virus"]["md5"]
    shell:
        """
        cd {params.wd}; \
        curl -Lgo {params.filename} "{params.tar}"; \
        curl -Lgo {params.filename}.md5 "{params.md5}"; \
        md5sum -c {params.filename}.md5; \
        rm {params.filename}.md5; \
        tar -Izstd -xf {params.filename}
        """

rule mmseqs_uniprot_taxdb:
    input:
        vdb=os.path.join(PROTPATH,"uniprot_virus_c99.db"),
        idm=os.path.join(TAXPATH,"uniprot_ncbi_mapping.dat"),
        nms=os.path.join(TAXPATH,"names.dmp"),
        nds=os.path.join(TAXPATH,"nodes.dmp"),
        mgd=os.path.join(TAXPATH,"merged.dmp"),
        dln=os.path.join(TAXPATH,"delnodes.dmp"),
    params:
        tax=TAXPATH
    benchmark:
        "benchmarks/mmseqs_uniprot_taxdb.txt"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/mmseqs2.yaml"
    output:
        multiext(os.path.join(PROTPATH,"uniprot_virus_c99"),
            ".db_mapping",".db_names.dmp",".db_nodes.dmp",
            ".db_merged.dmp",".db_delnodes.dmp"),
        tmp=temp(directory(os.path.join(TMPDIR,'mmseqs_uniprot_taxdb')))
    shell:
        """
        mmseqs createtaxdb --ncbi-tax-dump {params.tax} --tax-mapping-file {input.idm} {input.vdb} $(mkdir -p {output.tmp})
        """

rule download_uniref50:
    """
    Download the uniprot uniref50 database
    """
    output:
        os.path.join(PROTPATH,"uniref50.fasta.gz")
    conda: "../envs/curl.yaml"
    params: config["url"]["uniref50"]
    resources:
        time_min=480
    shell:
        """
        curl -Lgo {output} "{params}"
        """

rule uniref_plus_viruses:
    input:
        ur=os.path.join(PROTPATH,"uniref50.fasta.gz"),
        vr=os.path.join(PROTPATH,"uniprot_virus_c99.faa"),
    output:
        temp(os.path.join(UNIREF50VIR,"uniref50_virus.fasta"))
    conda:
        "../envs/pigz.yaml"
    shell:
        """
        unpigz -c {input.ur} | cat - {input.vr} > {output}
        """

rule mmseqs_urv:
    input:
        os.path.join(UNIREF50VIR,"uniref50_virus.fasta")
    output:
        os.path.join(UNIREF50VIR,"uniref50_virus.db")
    benchmark:
        "benchmarks/mmseqs_urv.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        "mmseqs createdb {input} {output}"

rule mmseqs_urv_taxonomy:
    input:
        vdb=os.path.join(UNIREF50VIR,"uniref50_virus.db"),
        idm=os.path.join(TAXPATH,"uniprot_ncbi_mapping.dat"),
        nms=os.path.join(TAXPATH,"names.dmp"),
        nds=os.path.join(TAXPATH,"nodes.dmp"),
        mgd=os.path.join(TAXPATH,"merged.dmp"),
        dln=os.path.join(TAXPATH,"delnodes.dmp"),
    params:
        tax=TAXPATH
    output:
        multiext(os.path.join(UNIREF50VIR,"uniref50_virus"),".db_mapping",".db_names.dmp",".db_nodes.dmp",".db_merged.dmp",
            ".db_delnodes.dmp"),
        tmp=temp(directory(os.path.join(TMPDIR,'mmseqs_urv_taxonomy')))
    benchmark:
        "benchmarks/mmseqs_urv_taxonomy.txt"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createtaxdb --ncbi-tax-dump {params.tax} --tax-mapping-file {input.idm} {input.vdb} \
        $(mkdir -p {output.tmp}) --threads {resources.cpus}
        """

# rule mmseqs_nt_db:
#     input:
#         nt = os.path.join(NUCLPATH, "nt.fna")
#     output:
#         idx = os.path.join(NUCLPATH, "ntDB.index"),
#         dbt = os.path.join(NUCLPATH, "ntDB.dbtype")
#     benchmark:
#         "benchmarks/mmseqs_nt_db.txt"
#     resources:
#         mem_mb=20000,
#         cpus=8
#     conda:
#         "../envs/mmseqs2.yaml"
#     params:
#         db = os.path.join(NUCLPATH, "ntDB")
#     shell:
#         """
#         mmseqs createdb {input} {params.db} --dbtype 2 --shuffle 0
#         """

# rule line_sine_download:
#     """
#     A database of LINES and SINES that we screen against to
#     remove contaminants.
#
#     LINES: "http://sines.eimb.ru/banks/SINEs.bnk"
#     SINES: "http://sines.eimb.ru/banks/LINEs.bnk"
#
#     Note that the current version has two ids with the same name
#     but different sequences, so we used seqtk to rename them
#     """
#     output:
#         os.path.join(CONPATH,"line_sine.fasta")
#     conda:
#         "../envs/curl.yaml"
#     params:
#         wd=CONPATH,
#         filename=config["url"]["lineSine"]["filename"],
#         zst=config["url"]["lineSine"]["zst"],
#         md5=config["url"]["lineSine"]["md5"]
#     shell:
#         """
#         cd {params.wd};
#         curl -Lgo {params.filename} "{params.zst}";
#         curl -Lgo {params.filename}.md5 "{params.md5}";
#         md5sum -c {params.filename}.md5;
#         rm {params.filename}.md5;
#         zstd -d {params.filename} --rm
#         """

# rule line_sine_format:
#     """
#     Build the bowtie2 indices
#     """
#     input:
#         os.path.join(CONPATH,"line_sine.fasta")
#     output:
#         expand(os.path.join(CONPATH,"line_sine.{n}.bt2"),n=[1, 2, 3, 4]),
#         expand(os.path.join(CONPATH,"line_sine.rev.{m}.bt2"),m=[1, 2])
#     params:
#         idx=os.path.join(CONPATH,"line_sine")
#     resources:
#         mem_mb=100000,
#         cpus=8
#     conda:
#         "../envs/bowtie2.yaml"
#     shell:
#         """
#         bowtie2-build {input} {params.idx}
#         """

rule download_taxonomizr:
    """
    Download and extract taxonomiser
    """
    output:
        os.path.join(TAXPATH,"taxonomizr_accessionTaxa.sql")
    conda:
        "../envs/curl.yaml"
    params:
        wd=TAXPATH,
        filename=config["url"]["taxonomizer"]["filename"],
        zst=config["url"]["taxonomizer"]["zst"],
        md5=config["url"]["taxonomizer"]["md5"]
    shell:
        """
        cd {params.wd}; \
        curl -Lgo {params.filename} "{params.zst}"; \
        curl -Lgo {params.filename}.md5 "{params.md5}"; \
        md5sum -c {params.filename}.md5; \
        rm {params.filename}.md5; \
        zstd -d {params.filename} --rm
        """
