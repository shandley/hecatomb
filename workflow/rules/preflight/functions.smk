
import re


### PYTHON FUNCTIONS
def stream_tsv(tsvFile, skip_header=False):
    """Read a file line-by-line and split by whitespace"""
    with open(tsvFile, "r") as filehandle:
        if skip_header:
            filehandle.readline()
        for line in filehandle:
            line = line.strip()
            l = line.split("\t")
            yield l


def file_len(fname):
    """Return the number of lines in a file"""
    if fname.endswith(".gz"):
        import gzip
        f = gzip.open(fname, "rb")
    else:
        f = open(fname, "r")
    if re.search(r"fastq(\.gz)?$", fname):
        for i, l in enumerate(f):
            pass
        f.close()
        return int((i + 1) / 4)
    else:
        n=0
        for line in f:
            if line.startswith(">"):
                n+=1
        f.close()
        return n


## GENERIC RECIPES
if config["hecatomb"]["args"]["custom_aa"]:
    rule functions_create_custom_primary_aa:
        """Create a custom primary aa database from a FASTA"""
        input:
            config["hecatomb"]["args"]["custom_aa"]
        output:
            os.path.join(config["hecatomb"]["args"]["databases"], "virus_primary_aa", "sequenceDB")
        log:
            os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "create_custom_primary_aa.stderr")
        conda:
            os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
        container:
            config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
        shell:
            "mmseqs createdb {input} {output} --dbtype 1 &> {log}"


if config["hecatomb"]["args"]["custom_nt"]:
    rule functions_create_custom_primary_nt:
        """Create a custom primary nt database from a FASTA"""
        input:
            config["hecatomb"]["args"]["custom_nt"]
        output:
            os.path.join(config["hecatomb"]["args"]["databases"], "virus_primary_nt", "sequenceDB")
        log:
            os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "create_custom_primary_nt.stderr")
        conda:
            os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
        container:
            config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
        shell:
            "mmseqs createdb {input} {output} --dbtype 2 &> {log}"


rule functions_fasta_index:
    """Index a .fasta file for rapid access with samtools faidx."""
    input:
        "{file}.fasta"
    output:
        "{file}.fasta.fai"
    log:
        "{file}.samtools.stderr"
    conda:
        os.path.join("..", "..", "envs", "minimap2_samtools.yaml")
    container:
        config["hecatomb"]["container"]["minimap2_samtools"]
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    shell:
        "samtools faidx {input} > {output} 2> {log}"


rule functions_bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    log:
        "{file}.samtools.stderr"
    conda:
        os.path.join("..", "..", "envs", "minimap2_samtools.yaml")
    container:
        config["hecatomb"]["container"]["minimap2_samtools"]
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"


rule functions_calculate_gc:
    """Calculate GC content for sequences"""
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.fasta")
    output:
        temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.properties.gc"))
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "calculate_gc.{file}.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "calculate_gc.{file}.log")
    conda:
        os.path.join("..", "..", "envs", "bbmap_bedtools_pigz_flye.yaml")
    container:
        config["hecatomb"]["container"]["bbmap_bedtools_pigz_flye"]
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    shell:
        "countgc.sh in={input} format=2 ow=t > {output} 2> {log}"


rule functions_calculate_tet_freq:
    """Calculate tetramer frequency

    The tail commands trims the first line which is junk that should be printed to stdout.
    """
    input:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.fasta")
    output:
        temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.properties.tetramer"))
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "calculate_tet_freq.{file}.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "calculate_tet_freq.{file}.log")
    conda:
        os.path.join("..", "..", "envs", "bbmap_bedtools_pigz_flye.yaml")
    container:
        config["hecatomb"]["container"]["bbmap_bedtools_pigz_flye"]
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    shell:
        "{{ "
        "tetramerfreq.sh in={input} w=0 ow=t -Xmx{resources.mem_mb}m "
            "| tail -n+2; "
        "}} > {output} 2>> {log}; "


rule functions_seq_properties_table:
    """Combine GC and tet freq tables"""
    input:
        gc=os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.properties.gc"),
        tet=os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.properties.tetramer")
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "{file}.properties.tsv")
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "seq_properties_table.{file}.txt")
    resources:
        **config["resources"]["med"]
    threads:
        config["resources"]["med"]["cpu"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "{file}.seq_properties_table.log")
    script:
        os.path.join("..", "..", "scripts",  "seqPropertyTable.py")
