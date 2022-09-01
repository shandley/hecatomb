"""
Database and output locations for Hecatomb

Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""


### DATABASE BASE DIRECTORY
try:
    if config['Databases'] is None:
        DBDIR = os.path.join(workflow.basedir, '..', '..', 'databases')
    else:
        DBDIR = config['Databases']
except KeyError:
    DBDIR = os.path.join(workflow.basedir,'..','..','databases')


### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'hecatomb_out'
else:
    OUTPUT = config['Output']


### DATABASE SUBDIRs
CONPATH = os.path.join(DBDIR, "contaminants")
TAX = os.path.join(DBDIR, "tax", "taxonomy")
TABLES = os.path.join(DBDIR, "tables")
HOSTPATH = os.path.join(DBDIR, "host")


### MMSEQS DBs
UNIVIRDB = os.path.join(DBDIR, "aa", "virus_primary_aa")
UNIREF50VIR = os.path.join(DBDIR, "aa", "virus_secondary_aa")
NCBIVIRDB = os.path.join(DBDIR, "nt", "virus_primary_nt")
POLYMICRODB = os.path.join(DBDIR, "nt", "virus_secondary_nt")


### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'results')
WORKDIR = os.path.join(OUTPUT, 'processing')
TMPDIR = os.path.join(WORKDIR, 'temp')
STDERR = os.path.join(OUTPUT, 'stderr')
BENCH = os.path.join(OUTPUT, 'benchmarks')
# SUMDIR = os.path.join('hecatomb_report')
ASSEMBLY = os.path.join(WORKDIR, 'assembly')
MAPPING = os.path.join(WORKDIR, 'mapping')
STATS = os.path.join(WORKDIR, 'stats')


# MMSEQS OUTPUT DIRs
PRIMARY_AA_OUT = os.path.join(WORKDIR, "mmseqs_aa_primary")
SECONDARY_AA_OUT = os.path.join(WORKDIR, "mmseqs_aa_secondary")
PRIMARY_NT_OUT = os.path.join(WORKDIR, "mmseqs_nt_primary")
SECONDARY_NT_OUT = os.path.join(WORKDIR, "mmseqs_nt_secondary")

