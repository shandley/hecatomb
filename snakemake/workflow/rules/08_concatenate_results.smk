"""
concatenate_results.sh

take all the results and make some cool data sets from them!
"""


### Ignore this header, it comes directly from hecatomb.snakefile



import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


###################################################################
#                                                                 #
# The output directories where results and analyses are written   #
#                                                                 #
###################################################################

# paths for our data. This is where we will read and put things
READDIR = config['Paths']['Reads']
CLUMPED = config['Output']["Clumped"]
QC = config['Output']['QC']
RESULTS = config['Output']['Results']

###################################################################
#                                                                 #
# We create subdirectories in this temp directory as needed       #
#                                                                 #
###################################################################

TMPDIR = config['Paths']['Temp']
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR, exist_ok=True)


###################################################################
#                                                                 #
# The database and path structure. Note that below we define the  #
# databases that we are actaully looking for, and test            #
# to see if each of them exist. You should be able                #
# to automatically download them.                                 #
#                                                                 #
###################################################################

DBDIR = config['Paths']['Databases']
BACBT2 = os.path.join(DBDIR, "bac_giant_unique_species", "bac_uniquespecies_giant.masked_Ns_removed")
HOSTBT2 = os.path.join(DBDIR, "human_masked", "human_virus_masked")
CONPATH = os.path.join(DBDIR, "contaminants")
PROTPATH = os.path.join(DBDIR, "proteins")
NUCLPATH = os.path.join(DBDIR, "nucleotides")
TAXPATH  = os.path.join(DBDIR, "taxonomy")

###################################################################
#                                                                 #
# The bacterial and host bowtie2 indexes                          #
#                                                                 #
###################################################################

for bti in [BACBT2, HOSTBT2]:
    if not os.path.exists(f"{bti}.1.bt2l") and not os.path.exists(f"{bti}.1.bt2"):
        sys.stderr.write(f"FATAL: You do not appear to have the bowtie2 indexes for {bti}\n")
        sys.stderr.write("This version of hecatomb uses bowtie to remove host, bacteria, and contaminants\n")
        sys.stderr.write("Please re-run the alt version of download databases\n")
        sys.exit()

if not os.path.exists(os.path.join(CONPATH, "line_sine.1.bt2")):
    sys.stderr.write("FATAL: You do not appear to have the bowtie2 indexes for the line/sine database\n")
    sys.stderr.write("This version of hecatomb uses bowtie to remove contaminants\n")
    sys.stderr.write("Please re-run the alt version of download databases\n")
    sys.exit()

if not os.path.exists(PROTPATH):
    sys.stderr.write("FATAL: You appear not to have the protein databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()

###################################################################
#                                                                 #
# Amino acid definitions and paths for searches                   #
#                                                                 #
###################################################################

AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    os.makedirs(AA_OUT, exist_ok=True)

AA_OUT_CHECKED  = os.path.join(RESULTS, "mmseqs_aa_checked_out")
if not os.path.exists(AA_OUT_CHECKED):
    os.makedirs(AA_OUT_CHECKED, exist_ok=True)

###################################################################
#                                                                 #
# nucleotide definitions and paths for searches                   #
#                                                                 #
###################################################################

NTDB = os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked", "nt.fnaDB")
if not os.path.exists(NTDB):
    sys.stderr.write(f"FATAL: You appear not to have the nucleotide ")
    sys.stderr.write(f"database {NTDB} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()


NT_OUT = os.path.join(RESULTS, "mmseqs_nt_out")
if not os.path.exists(NT_OUT):
    os.makedirs(NT_OUT)

NT_CHECKED_OUT = os.path.join(RESULTS, "mmseqs_nt_checked_out")
if not os.path.exists(NT_CHECKED_OUT):
    os.makedirs(NT_CHECKED_OUT)

# note that we have two databases called "nt.fnaDB". Sorry.
BVMDB = os.path.join(NUCLPATH, "bac_virus_masked", "nt.fnaDB")
if not os.path.exists(BVMDB):
    sys.stderr.write(f"FATAL: You appear not to have the nucleotide ")
    sys.stderr.write(f"database {BVMDB} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()



###################################################################
#                                                                 #
# Taxonomy databases and related information                      #
#                                                                 #
###################################################################

TAXTAX = os.path.join(TAXPATH, "taxonomizr_accessionTaxa.sql")
if not os.path.exists(TAXTAX):
    sys.stderr.write(f"FATAL: You appear not to have the taxonomizr ")
    sys.stderr.write(f"database {TAXTAX} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()


PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
if not os.path.exists(PHAGE_LINEAGES):
    sys.stderr.write("FATAL: phages/phage_taxonomic_lineages.txt not ")
    sys.stderr.write("found in the databases directory. Please check ")
    sys.stderr.write("you have the latest version of the databases\n")
    sys.exit()

###################################################################
#                                                                 #
# Uniprot databases and related information                       #
#                                                                 #
###################################################################

URVPATH = os.path.join(PROTPATH, "uniref_plus_virus")
URVDB = os.path.join(URVPATH, "uniref50_virus.db") # uniref50 + viruses database
if not os.path.exists(URVDB):
    sys.stderr.write("FATAL: {URVDB} not found.\n")
    sys.stderr.write("Please make sure that you have run ")
    sys.stderr.write("download_databases.snakefile before commencing\n")
    sys.exit()

VIRDB = os.path.join(PROTPATH, "uniprot_virus_c99.db")
if not os.path.exists(VIRDB):
    sys.stderr.write(f"FATAL: {VIRDB} does not exist. Please ensure you")
    sys.stderr.write(" have installed the databases\n")
    sys.exit()


###################################################################
#                                                                 #
# Read the sequence files and parse the file names.               #
#                                                                 #
###################################################################

SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extentions}'))

if not EXTENSIONS:
    sys.stderr.write("""
        FATAL: We could not parse the sequence file names.
        We are expecting {sample}_R1{extension}, and so your files
        should contain the characters '_R1' in the fwd reads
        and '_R2' in the rev reads
        """)
    sys.exit()
# we just get the generic extension. This is changed in Step 1

file_extension = EXTENSIONS[0]
# a convenience so we don't need to use '{sample}_R1' all the time
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'

if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write("You should complain to Rob\n")
    sys.exit()

rule concatenate_results_first:
    input:
        os.path.join(RESULTS, "viruses_tax_table.tsv"),
        os.path.join(RESULTS, "phage_tax_table.tsv"),
        os.path.join(RESULTS, "aa.aln.m8"),
        os.path.join(RESULTS, "nt.aln.m8"),
        directory("family_reads")

rule jive_aa_annotation:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv")
    shell:
        """
        sed 's/uc_//g' {input} > {output}
        """

rule jive_nt_annotation:
    input:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage_edited.tsv")
    shell:
        """
        tail -n +2 {input} | sed 's/NA/unknown/g; s/uc_//g' > {output}
        """

# this is the name of the rule in concatenate_results.sh and it is too cute to change
rule happily_marry_outputs:
    input:
        aa = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv"),
        nt = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage_edited.tsv")
    output:
        temporary(os.path.join(RESULTS, "viruses_tax_table_tmp.tsv"))
    shell:
        """
        cat {input.aa} {input.nt} | sort -n -k 1 > {output}
        """

rule add_crown_to_marriage:
    # not really, its a title
    input:
        os.path.join(RESULTS, "viruses_tax_table_tmp.tsv")
    output:
        os.path.join(RESULTS, "viruses_tax_table.tsv")
    shell:
        """
        sed -e '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}
        """

rule fix_phage_names:
    input:
        tsv = os.path.join(AA_OUT, "phage_table.tsv")
    output:
        temporary(os.path.join(RESULTS, "phage_tax_table_temp.tsv"))
    shell:
        """
        sed 's/uc_//g' {input} > {output}
        """

rule fix_phage_cols:
    input:
        os.path.join(RESULTS, "phage_tax_table_temp.tsv")
    output:
        temporary(os.path.join(RESULTS, "phage_tax_table_temp2.tsv"))
    shell:
        """
        cut -f1-8 {input} > {output}
        """
            
rule add_phage_title:
    input:
        os.path.join(RESULTS, "phage_tax_table_temp2.tsv")
    output:
        os.path.join(RESULTS, "phage_tax_table.tsv")
    shell:
        """
        sed -e '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}
        """

rule add_aa_tax_header:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    output:
        os.path.join(RESULTS, "aa.aln.m8")
    shell:
        """
        sed -e '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' \
                {input} > {output}
        """

rule add_nt_tax_header:
    input:
        os.path.join(NT_OUT, "resultDB.firsthit.m8")
    output:
        os.path.join(RESULTS, "nt.aln.m8")
    shell:
        """
        sed -e '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' \
                {input} > {output}
        """

"""
I am not very happy with this solution and looking for some alternative ideas
Basically this is very linear, and I can't figure out how to make it non-linear!
"""

rule make_fasta:
    input:
        os.path.join(RESULTS, "viruses_tax_table.tsv")
    output:
        directory("family_reads")
    shell:
        """
        mkdir -p {output} && 
        for FAM in $(tail -n +2 {input} | cut -f6 | awk '!s[$0]++');
        do 
            for TID in $(grep $FAM {input} | cut -f 1);
            do
                grep -A1 -Fw $TID results/seqtable.fasta >> {output}/$FAM.fasta
            done
        done
        """




