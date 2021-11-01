rule tax_level_counts:
    """Generate a long table of hit counts at different taxon levels
    
    (excluding at the species level as that is essentially the bigtable.tsv file)
    """
    input:
        os.path.join(RESULTS, "bigtable.tsv")
    output:
        report(os.path.join(SUMDIR, "taxonLevelCounts.tsv"),
            caption = "../report/tax_level_counts.rst",
            category = "Output")
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    run:
        # the index position for each tax level in the bigtable.tsv - off by one because we capture the range as l[21:i]
        tsvIndex = {'kingdom': 24, 'phylum': 25, 'class': 26, 'order': 27, 'family': 28, 'genus': 29, 'species': 30}
        idxStart = 23
        short = {'23':'k', '24':'p', '25':'c', '26':'o', '27':'f', '28':'g', '29':'s'}
        out = open(output[0], 'w')
        out.write('sampleID\ttaxonLevel\ttaxonPath\ttaxonName\tcount\tnormCount\n')
        # re-read the file for each sample to keep the memory happy - this is probably not necessary
        for sample in SAMPLES:
            counts = {}         # counts[taxlevel][taxname] = int
            normCounts = {}     # normalised counts, same structure as counts
            infh = open(input[0], 'r')
            infh.readline()     # skip header
            for line in infh:
                l = line.split('\t')
                if l[1] == sample:
                    for t,i in tsvIndex.items():
                        try:
                            if len(l[i].strip()) == 0:
                                continue
                        except IndexError:
                            continue
                        try:
                            counts[t]
                            normCounts[t]
                        except KeyError:
                            counts[t] = {}
                            normCounts[t] = {}
                        taxPath = []
                        for o in range(idxStart,i):
                            taxPath.append(f'{short[str(o)]}_{l[o]}')   # taxon path = k_kingName,p_phylName etc.
                        outPath = ','.join(taxPath)
                        try:
                            counts[t][outPath] += int(l[2])
                            normCounts[t][outPath] += float(l[3])
                        except KeyError:
                            counts[t][outPath] = int(l[2])
                            normCounts[t][outPath] = float(l[3])
            infh.close()
            for taxLevel in counts.keys():
                for taxPath in counts[taxLevel].keys():
                    taxName = taxPath.split('_')[-1]
                    out.write(f'{sample}\t{taxLevel}\t{taxPath}\t{taxName}\t{counts[taxLevel][taxPath]}\t{normCounts[taxLevel][taxPath]}\n')
        out.close()


### Collect read counts following each preprocessing step
rule dumpSamplesTsv:
    output:
        os.path.join(SUMDIR, 'hecatomb.samples.tsv')
    run:
        writeSamplesTsv(sampleReads, output[0])

rule step00_counts:
    output:
        report(os.path.join(SUMDIR,"Step00_counts.tsv"),
            caption = "../report/step00.rst",
            category = "Preprocessing")
    run:
        collect_start_counts(sampleReads, output[0])

rule step01_counts:
    input:
        expand(os.path.join(TMPDIR, "p01", "{sample}_{rn}.s1.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step01_counts.tsv"),
            caption = "../report/step01.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p01"), ".s1.out.fastq", "Step_01", output[0])

rule step02_counts:
    input:
        expand(os.path.join(TMPDIR, "p02", "{sample}_{rn}.s2.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step02_counts.tsv"),
            caption = "../report/step02.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p02"), ".s2.out.fastq", "Step_02", output[0])

rule step03_counts:
    input:
        expand(os.path.join(TMPDIR, "p03", "{sample}_{rn}.s3.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step03_counts.tsv"),
            caption = "../report/step03.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p03"), ".s3.out.fastq", "Step_03", output[0])

rule step04_counts:
    input:
        expand(os.path.join(TMPDIR, "p04", "{sample}_{rn}.s4.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step04_counts.tsv"),
            caption = "../report/step04.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p04"), ".s4.out.fastq", "Step_04", output[0])

rule step05_counts:
    input:
        expand(os.path.join(TMPDIR, "p05", "{sample}_{rn}.s5.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step05_counts.tsv"),
            caption = "../report/step05.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p05"), ".s5.out.fastq", "Step_05", output[0])

rule step06_counts:
    input:
        expand(os.path.join(TMPDIR, "p06", "{sample}_{rn}.s6.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step06_counts.tsv"),
            caption = "../report/step06.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p06"), ".s6.out.fastq", "Step_06", output[0])

rule step07_counts:
    input:
        expand(os.path.join(TMPDIR, "p07", "{sample}_{rn}.all.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step07_counts.tsv"),
            caption = "../report/step07.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "p07"), ".all.fastq", "Step_07", output[0])

rule step08_counts:
    input:
        expand(os.path.join(TMPDIR, "p08", "{sample}_R1.deduped.out.fastq"), sample=SAMPLES)
    output:
        report(os.path.join(SUMDIR, "Step08_counts.tsv"),
            caption = "../report/step08.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        for sample in SAMPLES:
            R1c = file_len(os.path.join(TMPDIR, "p08", f"{sample}_R1.deduped.out.fastq")) / 4
            outfh.write(f'{sample}\tStep_08\tR1\t{R1c}\n')
        outfh.close()

rule step09_counts:
    input:
        expand(os.path.join(TMPDIR, "p09", "{sample}_R1_rep_seq.fasta"), sample=SAMPLES)
    output:
        report(os.path.join(SUMDIR, "Step09_counts.tsv"),
            caption = "../report/step09.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        for sample in SAMPLES:
            R1c = file_len(os.path.join(TMPDIR, "p09", f"{sample}_R1_rep_seq.fasta")) / 2
            outfh.write(f'{sample}\tStep_09\tR1\t{R1c}\n')
        outfh.close()


### Collect rep seq counts following each search stage
rule step10_counts:
    input:
        classSeq = os.path.join(PRIMARY_AA_OUT,"MMSEQS_AA_PRIMARY_classified.fasta"),
        unclassSeq = os.path.join(PRIMARY_AA_OUT,"MMSEQS_AA_PRIMARY_unclassified.fasta")
    output:
        report(os.path.join(SUMDIR, "Step10_counts.tsv"),
            caption = "../report/step10.rst",
            category = "Read annotation")
    run:
        classCount = fasta_clust_counts(input.classSeq)
        unclassCount = fasta_clust_counts(input.unclassSeq)
        outfh = open(output[0],'w')
        outfh.write(f"Primary AA viral:\t{classCount}\n")
        outfh.write(f"Primary AA non-viral:\t{unclassCount}\n")
        outfh.close()

rule step11_counts:
    input:
        tsv = os.path.join(SECONDARY_AA_OUT,"AA_bigtable.tsv")
    output:
        report(os.path.join(SUMDIR, "Step11_counts.tsv"),
            caption = "../report/step11.rst",
            category = "Read annotation")
    run:
        virCnts = 0
        nonVirCnts = 0
        infh = open(input.tsv,'r')
        infh.readline() # skip header
        for line in infh:
            l = line.split('\t')
            if l[23] == "Viruses":
                virCnts += int(l[2])
            else:
                nonVirCnts += int(l[2])
        infh.close()
        outfh = open(output[0],'w')
        outfh.write(f"Secondary AA viral:\t{virCnts}\n")
        outfh.write(f"Secondary AA non-viral:\t{nonVirCnts}\n")
        outfh.close()

rule step12_counts:
    input:
        classSeq = os.path.join(PRIMARY_NT_OUT,"classified_seqs.fasta"),
        unclassSeq = os.path.join(PRIMARY_NT_OUT,"unclassified_seqs.fasta")
    output:
        report(os.path.join(SUMDIR,"Step12_counts.tsv"),
            caption = "../report/step12.rst",
            category = "Read annotation")
    run:
        classCount = fasta_clust_counts(input.classSeq)
        unclassCount = fasta_clust_counts(input.unclassSeq)
        outfh = open(output[0],'w')
        outfh.write(f"Primary NT viral:\t{classCount}\n")
        outfh.write(f"Primary NT non-viral:\t{unclassCount}\n")
        outfh.close()

rule step13_counts:
    input:
        os.path.join(SECONDARY_NT_OUT,"NT_bigtable.tsv")
    output:
        report(os.path.join(SUMDIR,"Step13_counts.tsv"),
            caption = "../report/step13.rst",
            category = "Read annotation")
    run:
        virCnts = 0
        nonVirCnts = 0
        unKnwn = 0
        infh = open(input[0],'r')
        infh.readline()     # skip header
        for line in infh:
            l = line.split('\t')
            if l[23] == "Viruses":
                virCnts += int(l[2])
            else:
                nonVirCnts += int(l[2])
        infh.close()
        outfh = open(output[0], 'w')
        outfh.write(f"Secondary NT viral:\t{virCnts}\n")
        outfh.write(f"Secondary NT non-viral:\t{nonVirCnts}\n")
        outfh.close()


rule krona_text_format:
    """Taxon step 18: Text format summary of bigtable for a krona plot"""
    input:
        os.path.join(RESULTS, "bigtable.tsv")
    output:
        os.path.join(SUMDIR, "krona.txt")
    benchmark:
        os.path.join(BENCH, "krona_text_format.txt")
    log:
        os.path.join(STDERR, 'krona_text_format.log')
    run:
        import logging
        import atexit
        atexit.register(exitLogCleanup,log[0])
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        logging.debug('Slurping Tax assignments from bigtable')
        counts = {}
        for l in stream_tsv(input[0]):
            if l[0]=="seqID":
                continue
            t = '\t'.join(l[23:])
            try:
                counts[t] += int(l[2])
            except KeyError:
                counts[t] = int(l[2])
        logging.debug('Sorting, counting, and writing tax assignments')
        outFH = open(output[0], 'w')
        for k in sorted(counts.keys()):
            outFH.write(f'{counts[k]}\t{k}\n')
        outFH.close()

rule krona_plot:
    """Taxon step 19: Krona plot of bigtable"""
    input:
        os.path.join(SUMDIR, "krona.txt")
    output:
        os.path.join(SUMDIR, "krona.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    benchmark:
        os.path.join(BENCH, "krona_plot.txt")
    log:
        os.path.join(STDERR, 'krona_plot.log')
    shell:
        """
        ktImportText {input} -o {output} &> {log}
        rm {log}
        """


rule contig_krona_text_format:
    input:
        os.path.join(RESULTS, "contigSeqTable.tsv")
    output:
        os.path.join(SUMDIR, "contigKrona.txt")
    log:
        os.path.join(STDERR, 'contig_krona_text_format.log')
    run:
        import logging
        import atexit
        atexit.register(exitLogCleanup,snakemake.log[0])
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        logging.debug('Slurping contig seq table')
        counts = {}
        for l in stream_tsv(input[0]):
            if l[0] == "contigID":
                continue
            t = '\t'.join((l[9:] + [l[0]]))
            c = l[1].split(':')
            try:
                counts[t] += int(c[1])
            except KeyError:
                counts[t] = int(c[1])
        logging.debug('Sorting and writing contig taxon info')
        outFH = open(output[0],'w')
        for k in sorted(counts.keys()):
            outFH.write(f'{counts[k]}\t{k}\n')
        outFH.close()

rule contig_krona_plot:
    input:
        os.path.join(SUMDIR, "contigKrona.txt")
    output:
        os.path.join(SUMDIR, "contigKrona.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    log:
        os.path.join(STDERR, 'contig_krona_plot.log')
    shell:
        """
        ktImportText {input} -o {output} &> {log}
        rm {log}
        """


### Generate a sankey diagram of the whole read-annotation pipeline
rule sankey_diagram:
    input:
        os.path.join(SUMDIR,"Step00_counts.tsv"),
        os.path.join(SUMDIR,"Step01_counts.tsv"),
        os.path.join(SUMDIR,"Step02_counts.tsv"),
        os.path.join(SUMDIR,"Step03_counts.tsv"),
        os.path.join(SUMDIR,"Step04_counts.tsv"),
        os.path.join(SUMDIR,"Step05_counts.tsv"),
        os.path.join(SUMDIR,"Step06_counts.tsv"),
        os.path.join(SUMDIR,"Step07_counts.tsv"),
        os.path.join(SUMDIR,"Step08_counts.tsv"),
        os.path.join(SUMDIR,"Step09_counts.tsv"),
        os.path.join(SUMDIR,"Step10_counts.tsv"),
        os.path.join(SUMDIR,"Step11_counts.tsv"),
        os.path.join(SUMDIR,"Step12_counts.tsv"),
        os.path.join(SUMDIR,"Step13_counts.tsv")
    output:
        report(os.path.join(SUMDIR, "Sankey.svg"),
            caption = "../report/sankey.rst",
            category = "Summary")
    conda:
        os.path.join('..', 'envs', 'plotly.yaml')
    log:
        os.path.join(STDERR, 'sankey_diagram.log')
    script:
        os.path.join('..', 'scripts', 'sankey.py')
