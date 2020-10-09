"""

The hecatomb snakefile has all the parts of the hecatomb pipeline, however
you will need to download the databases before we can begin!

Rob Edwards, Jan 2020
"""

import os
import sys
import socket

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


DBDIR = config['Paths']['Databases']

# paths for our databases
BACPATH = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "human_masked")
CONPATH = os.path.join(DBDIR, "contaminants")


if not os.path.exists(os.path.join(HOSTPATH, "ref")):
    sys.stderr.write("FATAL: You appear not to have the host databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()

# paths for our data. This is where we will read and put things
READDIR = config['Paths']['Reads']
CLUMPED = config['Output']["Clumped"]
QC = config['Output']['QC']
RESULTS = config['Output']['Results']

# how much memory we have
XMX = config['System']['Memory']

SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}_R1.fastq.gz'))
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'

# Summary:
	# Step 0: Clumpify reads (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/)
	# Step 1: Remove 5' amplification primer
	# Step 2: Remove 3' read through contaminant (Reverse complement of amplification primer + 6 bases of the adapter)
	# Step 3: Remove primer free adapter (both orientations)
	# Step 4: Remove adapter free primer (both orientations)
	# Step 5: PhiX Removal and vector contamination removal
	# Step 6: Host-removal
	# Step 7: Trim low-quality bases
	# Step 8: Read correction, merging and extension
	# Step 9: Remove bacterial contaminants reserving viral and aambiguous sequences
    # Step 10: Remove exact duplicates
    # Step 11: Dereplicate
    # Step 12: Extract sequences and counts for seqtable (count table)
    # Step 13. Parse and combine stats into a seqeunce count table

rule all:
    input:
        # step 9 output
        # expand(os.path.join(QC, "step_9", "{sample}.viral_amb.fastq"), sample=SAMPLES)
        # process all the data into a sequence table
        expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES),
        #convert the sequence table into a single outpupt
        os.path.join(RESULTS, "seqtable_all.tsv"),
        os.path.join(RESULTS, "seqtable.tab2fx")




"""
Process the data. 

These are the steps to do the work!
"""


rule clumpify:
    """
    Step 0: Clumpify and deduplicate reads
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + ".fastq.gz"),
        r2 = os.path.join(READDIR, PATTERN_R2 + ".fastq.gz")
    output:
        r1 = os.path.join(CLUMPED, PATTERN_R1 + ".clumped.fastq.gz"),
        r2 = os.path.join(CLUMPED, PATTERN_R2 + ".clumped.fastq.gz")
    shell:
        """
        clumpify.sh {XMX} in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            reorder=a \
            ow=t;
        """

rule remove_leftmost_primerB:
    """
    Step 1: Remove leftmost primerB. Not the reverse complements
    """
    input:
        r1 = os.path.join(CLUMPED, PATTERN_R1 + ".clumped.fastq.gz"),
        r2 = os.path.join(CLUMPED, PATTERN_R2 + ".clumped.fastq.gz"),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        stats = os.path.join(QC, "step_1", "{sample}.s1.stats.txt")
    shell:
        """
        bbduk.sh {XMX} in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t rcomp=f ow=t
        """

rule remove_3prime_contaminant:
    """
    Step 2: Remove 3' read through contaminant
    """
    input:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = os.path.join(QC, "step_2", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(QC, "step_2", PATTERN_R2 + ".s2.out.fastq"),
        stats = os.path.join(QC, "step_2", "{sample}.s2.stats.txt")
    shell:
        """
        bbduk.sh {XMX} in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t;
        """

rule remove_primer_free_adapter:
    """
    Step 3: Remove primer free adapter (both orientations)
    """
    input:
        r1 = os.path.join(QC, "step_2", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(QC, "step_2", PATTERN_R2 + ".s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = os.path.join(QC, "step_3", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(QC, "step_3", PATTERN_R2 + ".s3.out.fastq"),
        stats = os.path.join(QC, "step_3", "{sample}.s3.stats.txt")
    shell:
        """
        bbduk.sh {XMX} in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t;
        """

rule remove_adapter_free_primer:
    """
    Step 4: Remove adapter free primer (both orientations)
    """
    input:
        r1 = os.path.join(QC, "step_3", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(QC, "step_3", PATTERN_R2 + ".s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = os.path.join(QC, "step_4", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(QC, "step_4", PATTERN_R2 + ".s4.out.fastq"),
        stats = os.path.join(QC, "step_4", "{sample}.s4.stats.txt")
    shell:
        """
        bbduk.sh {XMX} in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t;
        """

rule remove_vector_contamination:
    """
    Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
    """
    input:
        r1 = os.path.join(QC, "step_4", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(QC, "step_4", PATTERN_R2 + ".s4.out.fastq"),
        primers = os.path.join(CONPATH, config['DatabaseFiles']['contaminants'])
    output:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
        stats = os.path.join(QC, "step_5", "{sample}.s5.stats.txt")
    shell:
        """
        bbduk.sh {XMX} in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t;
        """

rule host_removal:
    """
    Step 6: Host removal
    """
    input:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
        refpath = os.path.join(HOSTPATH, "ref")
    output:
        unmapped = os.path.join(QC, "step_6", "{sample}.unmapped.s6.out.fastq"),
        mapped = os.path.join(QC, "step_6", "{sample}.hostmapped.s6.out.fastq")
    params:
        hostpath = HOSTPATH 
    shell:
        """
        bbmap.sh {XMX} in={input.r1} in2={input.r2} \
            outu={output.unmapped} outm={output.mapped} \
            path={params.hostpath} \
            semiperfectmode=t quickmatch fast ordered=t ow=t;
        """

rule repair:
    """
    Step 6b. Repair the paired ends
    """
    input:
        unmapped = os.path.join(QC, "step_6", "{sample}.unmapped.s6.out.fastq"),
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq")
    shell:
        """   
        repair.sh {XMX} in={input.unmapped} \
            out={output.r1} out2={output.r2} \
            ow=t;
        """

rule trim_low_quality:
    """
    Step 7: Trim low-quality bases
    """
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq")
    output:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + ".s7.out.fastq"),
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq"),
        stats = os.path.join(QC, "step_7", "{sample}.s7.stats.txt")
    shell:
        """
        bbduk.sh {XMX} in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} outs={output.singletons} \
            stats={output.stats} \
            qtrim=r trimq=20 maxns=2 minlength=50 ordered=t ow=t;
        """

rule merge_reads:
    """
    Step 8: Merge forward (R1) and reverse (R2) reads
    """
    input:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + ".s7.out.fastq"),
    output:
        merged = os.path.join(QC, "step_8", "{sample}.merged.fastq"),
        r1 = os.path.join(QC, "step_8", PATTERN_R1 + ".unmerged.fastq"),
        r2 = os.path.join(QC, "step_8", PATTERN_R2 + ".unmerged.fastq"),
    shell:
        """
        bbmerge.sh in1={input.r1} in2={input.r2} \
            out={output.merged} outu1={output.r1} outu2={output.r2} \
            rem k=62 extend2=50 ecct vstrict=t ordered=t {XMX} ow=t;
        """

"""
In the original contaminant_removal.sh pipeline these were sequential. 
However, snakemake can easily run these in parallel as they
are independent. So I split them into a few separate rules
"""

rule get_r1_singletons:
    """
    Step 8b.i Split R1 singletons
    """
    input:
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq")
    output:
        r1singletons = os.path.join(QC, "step_8", "{sample}.singletons.R1.out.fastq"),
    shell:
        """
        grep --no-group-separator -A 3 '1:N:' {input.singletons} > {output.r1singletons};
        """
rule get_r2_singletons:
    """
    Step 8b.ii Split R2 singletons
    """
    input:
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq")
    output:
        r2singletons = os.path.join(QC, "step_8", "{sample}.singletons.R2.out.fastq")
    shell:
        """
        grep --no-group-separator -A 3 '2:N:' {input.singletons} > {output.r2singletons};
        """

rule concat_r1:
    """
    Step 8c.i Concatenate reads
    """
    input:
        merged = os.path.join(QC, "step_8", "{sample}.merged.fastq"),
        r1 = os.path.join(QC, "step_8", PATTERN_R1 + ".unmerged.fastq"),
        r1singletons = os.path.join(QC, "step_8", "{sample}.singletons.R1.out.fastq")
    output:
        r1combo = os.path.join(QC, "step_8", PATTERN_R1 + ".s8.out.fastq"),
    shell:
        """
        cat {input.merged} {input.r1} {input.r1singletons} > {output.r1combo};
        """

rule concat_r2:
    """
    Step 8c.ii Concatenate reads
    """
    input:
        merged = os.path.join(QC, "step_8", "{sample}.merged.fastq"),
        r2 = os.path.join(QC, "step_8", PATTERN_R2 + ".unmerged.fastq"),
        r2singletons = os.path.join(QC, "step_8", "{sample}.singletons.R2.out.fastq")
    output:
        r2combo = os.path.join(QC, "step_8", PATTERN_R2 + ".s8.out.fastq")
    shell:
        """
        cat {input.merged} {input.r2} {input.r2singletons} > {output.r2combo};
        """

rule remove_bacteria:
    """
    Step 9: Remove bacterial contaminants reserving viral and ambiguous sequences
    """
    input:
        r1 = os.path.join(QC, "step_8", PATTERN_R1 + ".s8.out.fastq"),
        r2 = os.path.join(QC, "step_8", PATTERN_R2 + ".s8.out.fastq"),
        bacpath = os.path.join(BACPATH, "ref")
    output:
        mapped = os.path.join(QC, "step_9", "{sample}.bacterial.fastq"),
        unmapped = os.path.join(QC, "step_9", "{sample}.viral_amb.fastq"),
        scafstats = os.path.join(QC, "step_9", "{sample}.scafstats.txt")
    params:
        bacpath = BACPATH
    shell:
        """
        bbmap.sh {XMX} in={input.r1} \
            path={params.bacpath} \
            outm={output.mapped} outu={output.unmapped} \
            scafstats={output.scafstats} \
            semiperfectmode=t quickmatch fast ordered=t ow=t;
       """


rule remove_exact_dups:
    """
    Step 10: Remove exact duplicates
    """
    input:
        os.path.join(QC, "step_9", "{sample}.viral_amb.fastq")
    output:
        os.path.join(QC, "step_10", "{sample}_R1.s9.deduped.out.fastq")
    shell:
        """
        dedupe.sh in={input} \
                out={output} \
                ac=f  ow=t {XMX};
        """

rule deduplicate:
    """
    Step 11: Dereplicate
    """
    input:
        os.path.join(QC, "step_10", "{sample}_R1.s9.deduped.out.fastq")
    output:
        fa = os.path.join(QC, "step_11", "{sample}_R1.best.fasta"),
        stats = os.path.join(QC, "step_11", "{sample}_stats.txt")
    shell:
        """
        dedupe.sh in={input} \
            csf={output.stats} out={output.fa} \
            ow=t s=4 rnc=t pbr=t {XMX};
        """

rule extract_seq_counts:
    """
    Step 12: Extract sequences and counts for seqtable (count table)
    """
    input:
        os.path.join(QC, "step_11", "{sample}_R1.best.fasta")
    output:
        os.path.join(QC, "step_12", "{sample}_reformated.fasta")
    shell:
        """
        reformat.sh in={input} out={output} \
            deleteinput=t fastawrap=0 \
            ow=t \
            {XMX};
        """

rule extract_counts:
    """
    Step 13. Parse and combine stats and contig files
    """
    input:
        os.path.join(QC, "step_12", "{sample}_reformated.fasta") 
    output:
        os.path.join(QC, "counts", "{sample}_seqs.txt")
    shell:
        """
        grep -v '>' {input} | sed '1i sequence' > {output}
        """

rule extract_counts_ids:
    """
    Step 13b. Extract sequence IDs
    """
    input:
        os.path.join(QC, "step_12", "{sample}_reformated.fasta")
    output:
        os.path.join(QC, "counts", "{sample}_contig_ids.txt")
    shell:
        """
        grep '>' {input} | sed 's|>Cluster_||' | awk -F "," '{{ print$1 }}' | sort -n | sed '1i contig_ids' > {output}
        """

rule exract_count_stats:
    """
    Step 13c. Extract counts
    """
    input:
        os.path.join(QC, "step_11", "{sample}_stats.txt")
    output:
        os.path.join(QC, "counts", "{sample}_counts.txt")
    params:
        s = "{sample}"
    shell:
        """
        cut -f 2 {input} | sed "1s/size/{params.s}/" > {output}
        """

rule create_seq_table:
    """
    Step 13d. Create sequence table
    """
    input:
        seq = os.path.join(QC, "counts", "{sample}_seqs.txt"),
        cnt = os.path.join(QC, "counts", "{sample}_counts.txt")
    output:
        os.path.join(QC, "counts", "{sample}_seqtable.txt")
    shell:
        """
        paste {input.seq} {input.cnt} > {output}
        """

rule merge_seq_table:
    """
    Merge seq counts
    """
    input:
        files = expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES)
    output:
        seqtable = os.path.join(RESULTS, "seqtable_all.tsv"),
        tab2fx = os.path.join(RESULTS, "seqtable.tab2fx")
    params:
        resultsdir = directory(RESULTS),
    script:
        "scripts/seqtable_merge.R"



