"""

Just download the databases. You should only need to do this once (unless we update the databases!)

Rob Edwards, Feb 2020
"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()

DBDIR = config['Paths']['Databases']
if not os.path.exists(DBDIR):
    os.mkdir(DBDIR)

TMPDIR = config['Paths']['Temp']
if not os.path.exists(TMPDIR):
    os.mkdir(TMPDIR)

# paths for our databases
BACPATH  = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "human_masked")
CONPATH  = os.path.join(DBDIR, "contaminants")
PROTPATH = os.path.join(DBDIR, "proteins")
TAXPATH  = os.path.join(DBDIR, "taxonomy")

# URLs where we download the data from
id_mapping_url    = "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
hecatomb_db_url   = "https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2"
taxdump_url       = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
uniprot_virus_url = "https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&format=fasta&&sort=score&fil=reviewed:no"


rule all:
    input:
        # the database directories
        os.path.join(BACPATH, "ref"),
        os.path.join(HOSTPATH, "ref"),
        os.path.join(PROTPATH, "uniprot_virus.faa"),
        multiext(os.path.join(PROTPATH, "uniprot_virus_c99"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp")


"""
These rules set up the databases. If we don't have them already installed,
we download them and index them using bbmap.sh.

If the databases already exist, these rules will be ignored
"""

rule download_databases:
    """
    If we don't have the database directories, download them
    """
    output:
        os.path.join(BACPATH, config['DatabaseFiles']['bacteria']),
        os.path.join(CONPATH, config['DatabaseFiles']['contaminants']),
        os.path.join(HOSTPATH, config['DatabaseFiles']['host']),
        os.path.join(CONPATH, "nebnext_adapters.fa"),
        os.path.join(CONPATH, "primerB.fa"),
        os.path.join(CONPATH, "rc_primerB_ad6.fa")
    shell:
        "cd {DBDIR} && curl -LO {hecatomb_db_url} && tar xf hecatomb.databases.tar.bz2"

rule make_bac_databases:
    input:
        os.path.join(BACPATH, config['DatabaseFiles']['bacteria'])
    output:
        directory(os.path.join(BACPATH, "ref"))
    params:
        wd = BACPATH,
        fa = config['DatabaseFiles']['bacteria']
    shell:
        "cd {params.wd} && bbmap.sh ref={params.fa}"

rule make_host_databases:
    input:
        os.path.join(HOSTPATH, config['DatabaseFiles']['host'])
    output:
        directory(os.path.join(HOSTPATH, "ref"))
    params:
        wd = HOSTPATH,
        fa = config['DatabaseFiles']['host']
    shell:
        "cd {params.wd} && bbmap.sh ref={params.fa}"


"""

Download the NCBI taxonomy databases and uniprot databases.

We run the downloads through snakemake rather than through bbmap as
snakemake will easily parallize everything.

This section is adapted from `accessory/uniprot_viral_DB_build.sh`

"""


rule download_uniprot_viruses:
    output:
        os.path.join(PROTPATH, "uniprot_virus.faa")
    shell:
        """
        mkdir -p {PROTPATH} && curl -Lo {output} {uniprot_virus_url}"
        """

rule download_ncbi_taxonomy:
    output:
        os.path.join(TAXPATH, "taxdump.tar.gz")
    shell:
        "curl -Lo {output} {taxdump_url}"


rule download_id_taxonomy_mapping:
    output:
        os.path.join(TAXPATH, "idmapping.dat.gz")
    shell:
        "cd {TAXPATH} && curl -LO {id_mapping_url}"

rule extract_ncbi_taxonomy:
    input:
        os.path.join(TAXPATH, "taxdump.tar.gz")
    params:
        path = TAXPATH
    output:
        os.path.join(TAXPATH, "citations.dmp"),
        os.path.join(TAXPATH, "delnodes.dmp"),
        os.path.join(TAXPATH, "division.dmp"),
        os.path.join(TAXPATH, "gc.prt"),
        os.path.join(TAXPATH, "gencode.dmp"),
        os.path.join(TAXPATH, "merged.dmp"),
        os.path.join(TAXPATH, "names.dmp"),
        os.path.join(TAXPATH, "nodes.dmp"),
        os.path.join(TAXPATH, "readme.txt")
    shell:
        "cd {params.path} && tar xf taxdump.tar.gz"

rule cluster_uniprot:
    input:
        os.path.join(PROTPATH, "uniprot_virus.faa")
    output:
        db = os.path.join(PROTPATH, "uniprot_virus_c99.faa"),
        cl = os.path.join(PROTPATH, "uniprot_virus_c99.faa.clstr")
    shell:
        "cd-hit -i {input} -o {output.db} -d 0 -c 0.99 -M 0 -T 0"

rule mmseqs_uniprot_clusters:
    input:
        os.path.join(PROTPATH, "uniprot_virus_c99.faa")
    output:
        os.path.join(PROTPATH, "uniprot_virus_c99.db")
    shell:
        "mmseqs createdb {input} {output}"

rule mmseqs_uniprot_taxdb:
    input:
        vdb = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        tax = os.path.join(TAXPATH, "nodes.dmp"),
        idm = os.path.join(TAXPATH, "idmapping.dat.gz")
    params:
        tax = TAXPATH
    output:
        multiext(os.path.join(PROTPATH, "uniprot_virus_c99"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp")
    shell:
        """
        mmseqs createtaxdb --ncbi-tax-dump {params.tax} --tax-mapping-file {input.idm} {input.vdb} $(mktemp -d -p {TMPDIR})
        """

