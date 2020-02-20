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

# paths for our databases
BACPATH = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "human_masked")
CONPATH = os.path.join(DBDIR, "contaminants")


rule all:
    input:
        # the database directories
        os.path.join(BACPATH, "ref"),
        os.path.join(HOSTPATH, "ref"),

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
    params:
        db = DBDIR
    shell:
        "cd {params.db} && curl -LO https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2 && tar xf hecatomb.databases.tar.bz2"

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
