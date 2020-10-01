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
NUCLPATH = os.path.join(DBDIR, "nucleotide")

# these are just derivied from above
URVPATH = os.path.join(PROTPATH, "uniref_plus_virus") # uniref50 + viruses


# URLs where we download the data from
id_mapping_url    = "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
hecatomb_db_url   = "https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2"
taxdump_url       = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
uniprot_virus_url = "https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&format=fasta&&sort=score&fil=reviewed:no"
ntacc2tax         = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"


rule all:
    input:
        # the database directories
        os.path.join(BACPATH, "ref"),
        os.path.join(HOSTPATH, "ref"),
        os.path.join(CONPATH, "line_sine.fasta"),
        os.path.join(PROTPATH, "uniprot_virus.faa"),
        os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat"),
        multiext(os.path.join(PROTPATH, "uniprot_virus_c99"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp"),
        multiext(os.path.join(URVPATH, "uniref50_virus"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp")



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
        mkdir -p {PROTPATH} && curl -Lo {output} '{uniprot_virus_url}'
        """

rule download_uniref50:
    output:
        os.path.join(PROTPATH, "uniref50.fasta.gz")
    shell:
        """
        mkdir -p {PROTPATH} && curl -Lo {output} "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
        """

rule download_ncbi_taxonomy:
    output:
        os.path.join(TAXPATH, "taxdump.tar.gz")
    shell:
        """
        curl -Lo {output} "{taxdump_url}"
        """


rule download_id_taxonomy_mapping:
    output:
        os.path.join(TAXPATH, "idmapping.dat.gz")
    shell:
        """
        cd {TAXPATH};
        curl -LO "{id_mapping_url}"
        """


rule uniprot_to_ncbi_mapping:
    input:
        os.path.join(TAXPATH, "idmapping.dat.gz")
    output:
        os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat")
    shell:
        """
        zcat {input} | awk '$2 == "NCBI_TaxID" {{print $1"\t"$3 }}' > {output}
        """

rule extract_ncbi_taxonomy:
    input:
        os.path.join(TAXPATH, "taxdump.tar.gz")
    params:
        path = TAXPATH
    output:
        temp(os.path.join(TAXPATH, "citations.dmp")),
        temp(os.path.join(TAXPATH, "delnodes.dmp")),
        temp(os.path.join(TAXPATH, "division.dmp")),
        temp(os.path.join(TAXPATH, "gc.prt")),
        temp(os.path.join(TAXPATH, "gencode.dmp")),
        temp(os.path.join(TAXPATH, "merged.dmp")),
        temp(os.path.join(TAXPATH, "readme.txt")),
        temp(os.path.join(TAXPATH, "names.dmp")),
        temp(os.path.join(TAXPATH, "nodes.dmp")),
    shell:
        "cd {params.path} && tar xf taxdump.tar.gz"

rule download_accession_to_tax:
    output:
        os.path.join(TAXPATH, "nucl_gb.accession2taxid.gz")
    shell:
        """
        cd {TAXPATH};
        curl -LO "{ntacc2tax}"
        """

rule extract_accession_to_tax:
    input:
        os.path.join(TAXPATH, "nucl_gb.accession2taxid.gz")
    output:
        os.path.join(TAXPATH, "nucl_gb.accession2taxid")
    shell:
        "unpigz {input}"

rule create_nt_tax_table:
    input:
        nt = os.path.join(NUCLPATH, "nt.fna"),
        tx = os.path.join(TAXPATH, "nucl_gb.accession2taxid")
    output:
        tx = os.path.join(NUCLPATH, "nt.tax")
    shell:
        """
        grep '^>' {input.nt} | cut -f 1 -d ' ' | sed -e 's/^>//' | \
        xargs -n 10 | sed -e 's/ /|/g' | \
        xargs -i egrep {} {input.tx} | cut -f 2,3 > {output.tx}
        """

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
        idm = os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat"),
        nms = os.path.join(TAXPATH, "names.dmp"),
        nds = os.path.join(TAXPATH, "nodes.dmp"),
        mgd = os.path.join(TAXPATH, "merged.dmp"),
        dln = os.path.join(TAXPATH, "delnodes.dmp"),
    params:
        tax = TAXPATH
    output:
        multiext(os.path.join(PROTPATH, "uniprot_virus_c99"), 
                 ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", 
                 ".db_merged.dmp", ".db_delnodes.dmp")
    shell:
        """
        mmseqs createtaxdb --ncbi-tax-dump {params.tax} --tax-mapping-file {input.idm} {input.vdb} $(mktemp -d -p {TMPDIR})
        """

rule uniref_plus_viruses:
    input:
        ur = os.path.join(PROTPATH, "uniref50.fasta.gz"),
        vr = os.path.join(PROTPATH, "uniprot_virus_c99.faa"),
    output:
        temp(os.path.join(URVPATH, "uniref50_virus.fasta"))
    shell:
        """
        unpigz -c {input.ur} | cat - {input.vr} > {output}
        """

rule mmseqs_urv:
    input:
        os.path.join(URVPATH, "uniref50_virus.fasta")
    output:
        os.path.join(URVPATH, "uniref50_virus.db")
    shell:
        "mmseqs createdb {input} {output}"

rule mmseqs_urv_taxonomy:
    input:
        vdb = os.path.join(URVPATH, "uniref50_virus.db"),
        idm = os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat"),
        nms = os.path.join(TAXPATH, "names.dmp"),
        nds = os.path.join(TAXPATH, "nodes.dmp"),
        mgd = os.path.join(TAXPATH, "merged.dmp"),
        dln = os.path.join(TAXPATH, "delnodes.dmp"),
    params:
        tax = TAXPATH
    output:
        multiext(os.path.join(URVPATH, "uniref50_virus"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp")
    shell:
        """
        mmseqs createtaxdb --ncbi-tax-dump {params.tax} --tax-mapping-file {input.idm} {input.vdb} $(mktemp -d -p {TMPDIR})
        """

rule mmseqs_nt_db:
    input:
        nt = os.path.join(NUCLPATH, "nt.fna")
    output:
        idx = os.path.join(NUCLPATH, "ntDB.index"),
        dbt = os.path.join(NUCLPATH, "ntDB.dbtype")
    params:
        db = os.path.join(NUCLPATH, "ntDB")
    shell:
        """
        mmseqs createdb {input} {params.db} --dbtype 2 --shuffle 0
        """

rule line_sine_download:
    """
    A database of LINES and SINES that we screen against to 
    remove contaminants.
    """
    output:
        os.path.join(CONPATH, "line_sine.fasta")
    shell:
        """
        curl -L http://sines.eimb.ru/banks/SINEs.bnk > {output} && \
        curl -L http://sines.eimb.ru/banks/LINEs.bnk >> {output}
        """

