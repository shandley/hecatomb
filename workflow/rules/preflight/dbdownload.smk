# DOWNLOAD RULES
rule download_db_file:
    """Download a Hecatomb-maintained DB file."""
    output:
        os.path.join(config["hecatomb"]["args"]["databases"], "{path}","{file}")
    wildcard_constraints:
        path=r"virus_.*_aa|virus_.*_nt|host|contaminants|tables"
    run:
        import urllib.request
        import urllib.parse
        import shutil
        dlUrl1 = urllib.parse.urljoin(config["hecatomb"]["dbs"]["mirror1"], os.path.join(wildcards.path, wildcards.file))
        dlUrl2 = urllib.parse.urljoin(config["hecatomb"]["dbs"]["mirror2"], os.path.join(wildcards.path, wildcards.file))
        try:
            with urllib.request.urlopen(dlUrl1) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
        except:
            with urllib.request.urlopen(dlUrl2) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)


rule unzip_test_db:
    """Unzip Hecatomb test DB files"""
    input:
        os.path.join(os.path.join(workflow.basedir, "..", "hecatomb-testdb"), "{path}","{file}.zst")
    output:
        os.path.join(config["hecatomb"]["args"]["databases"], "{path}","{file}")
    wildcard_constraints:
        path=r"virus_.*_testing",
        file=r"(?!.*zst$).*"
    conda:
        os.path.join("..", "..", "envs", "krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    shell:
        """
        zstd -d {input} -o {output}
        """


rule download_taxdump:
    """Download the current NCBI taxdump."""
    output:
        expand(os.path.join(config["hecatomb"]["args"]["databases"], "{file}"), file=config["hecatomb"]["dbtax"]["files"]),
        temp(config["hecatomb"]["dbtax"]["tar"])
    params:
        url = config["hecatomb"]["dbtax"]["url"],
        tar = config["hecatomb"]["dbtax"]["tar"],
        md5 = config["hecatomb"]["dbtax"]["md5"],
        dir = os.path.join(config["hecatomb"]["args"]["databases"], "taxonomy")
    conda:
        os.path.join("..", "..", "envs", "krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    shell:
        """
        curl {params.url} -o {params.tar}
        curl {params.md5} | md5sum -c
        mkdir -p {params.dir}
        tar xvf {params.tar} -C {params.dir}
        """
