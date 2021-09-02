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
        tsvIndex = {'kingdom': 23, 'phylum': 24, 'class': 25, 'order': 26, 'family': 27, 'genus': 28, 'species': 29}
        idxStart = 22
        short = {'22':'k', '23':'p', '24':'c', '25':'o', '26':'f', '27':'g', '28':'s'}
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
rule step00_counts:
    input:
        expand(os.path.join(READDIR, "{sample}_{rn}" + file_extension), sample=SAMPLES, rn=['R1','R2'])
    output:
        report(os.path.join(SUMDIR,"Step00_counts.tsv"),
            caption = "../report/step00.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(READDIR), file_extension, "Step_00", output[0])

rule step01_counts:
    input:
        expand(os.path.join(TMPDIR, "step_01", "{sample}_{rn}.s1.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step01_counts.tsv"),
            caption = "../report/step01.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "step_01"), ".s1.out.fastq", "Step_01", output[0])

rule step02_counts:
    input:
        expand(os.path.join(TMPDIR, "step_02", "{sample}_{rn}.s2.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step02_counts.tsv"),
            caption = "../report/step02.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "step_02"), ".s2.out.fastq", "Step_02", output[0])

rule step03_counts:
    input:
        expand(os.path.join(TMPDIR, "step_03", "{sample}_{rn}.s3.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step03_counts.tsv"),
            caption = "../report/step03.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "step_03"), ".s3.out.fastq", "Step_03", output[0])

rule step04_counts:
    input:
        expand(os.path.join(TMPDIR, "step_04", "{sample}_{rn}.s4.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step04_counts.tsv"),
            caption = "../report/step04.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "step_04"), ".s4.out.fastq", "Step_04", output[0])

rule step05_counts:
    input:
        expand(os.path.join(TMPDIR, "step_05", "{sample}_{rn}.s5.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step05_counts.tsv"),
            caption = "../report/step05.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "step_05"), ".s5.out.fastq", "Step_05", output[0])

rule step06_counts:
    input:
        expand(os.path.join(TMPDIR, "step_06", "{sample}_{rn}.s6.out.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step06_counts.tsv"),
            caption = "../report/step06.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(TMPDIR, "step_06"), ".s6.out.fastq", "Step_06", output[0])

rule step07_counts:
    input:
        expand(os.path.join(QC, "HOST_REMOVED", "{sample}_{rn}.all.fastq"), sample=SAMPLES, rn=['R1', 'R2'])
    output:
        report(os.path.join(SUMDIR, "Step07_counts.tsv"),
            caption = "../report/step07.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(QC, "HOST_REMOVED"), ".all.fastq", "Step_07", output[0])

rule step08_counts:
    input:
        expand(os.path.join(QC, "CLUSTERED", "{sample}_R1.deduped.out.fastq"), sample=SAMPLES)
    output:
        report(os.path.join(SUMDIR, "Step08_counts.tsv"),
            caption = "../report/step08.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        for sample in SAMPLES:
            R1c = file_len(os.path.join(QC, "CLUSTERED", f"{sample}_R1.deduped.out.fastq")) / 4
            outfh.write(f'{sample}\tStep_08\tR1\t{R1c}\n')
        outfh.close()

rule step09_counts:
    input:
        expand(os.path.join(QC, "CLUSTERED", "{sample}_R1_rep_seq.fasta"), sample=SAMPLES)
    output:
        report(os.path.join(SUMDIR, "Step09_counts.tsv"),
            caption = "../report/step09.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        for sample in SAMPLES:
            R1c = file_len(os.path.join(QC, "CLUSTERED", f"{sample}_R1_rep_seq.fasta")) / 2
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
            if l[22] == "Viruses":
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
            if l[22] == "Viruses":
                virCnts += int(l[2])
            else:
                nonVirCnts += int(l[2])
        infh.close()
        outfh = open(output[0], 'w')
        outfh.write(f"Secondary NT viral:\t{virCnts}\n")
        outfh.write(f"Secondary NT non-viral:\t{nonVirCnts}\n")
        outfh.close()


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
    run:
        import plotly.graph_objects as go
        s0 = sum_counts(input[0])
        s1 = sum_counts(input[1])
        s2 = sum_counts(input[2])
        s3 = sum_counts(input[3])
        s4 = sum_counts(input[4])
        s5 = sum_counts(input[5])
        s6 = sum_counts(input[6])
        s7 = sum_counts(input[7])
        s7R1 = sum_counts(input[7], R1=True)
        s8 = sum_counts(input[8])
        s9 = sum_counts(input[9])
        s1disc = s0 - s1
        s2disc = s1 - s2
        s3disc = s2 - s3
        s4disc = s3 - s4
        s5disc = s4 - s5
        s6disc = s5 - s6
        s7disc = s6 - s7
        s8R2 = s7 - s7R1
        s8disc = s7R1 - s8
        s9rdnt = s8 - s9
        searches = {}
        for f in input[10:]:
            for l in stream_tsv(f):
                searches[l[0]] = int(l[1])
        paaV = searches['Primary AA viral:']
        paaU = searches['Primary AA non-viral:']
        saaV = searches['Secondary AA viral:']
        saaU = searches['Secondary AA non-viral:']
        pntV = searches['Primary NT viral:']
        pntU = searches['Primary NT non-viral:']
        sntV = searches['Secondary NT viral:']
        sntU = searches['Secondary NT non-viral:']
        labels = [
            "Raw reads",                # 0
            "5' primer",                # 1
            "3' read-through",          # 2
            "Primer-free adapter",      # 3
            "Adapter-free primer",      # 4
            "Vector removal",           # 5
            "Low-qual trim",            # 6
            "Host removal",             # 7
            "Duplicate removal",        # 8
            "Clustered seqs",           # 9
            "Redundant seqs",           # 10
            "Discard R2",               # 11
            "Discarded",                # 12
            "Seqtable with counts",     # 13
            "Primary AA viral",         # 14
            "Primary AA non-viral",     # 15
            "Secondary AA viral",       # 16
            "Secondary AA non-viral",   # 17
            "Primary NT viral",         # 18
            "Primary NT non-viral",     # 19
            "Secondary NT viral",       # 20
            "Secondary NT non-viral",   # 21
            "Viruses",                  # 22
            "Non-viral (virus-like)",   # 23
            "Non-viral"                 # 24
        ]
        source = [
            0, 0,       # s1
            1, 1,       # s2
            2, 2,       # s3
            3, 3,       # s4
            4, 4,       # s5
            5, 5,       # s6
            6, 6,       # s7
            7, 7, 7,    # s8.1
            11,         # s8.2
            8, 8,       # s9.1
            9, 10,      # s9.2
            13, 13,     # pAA
            14, 14,     # sAA
            16, 17,     # sAA
            15, 15,     # pNT
            19,         # pNT
            18, 18,     # sNT
            20, 21      # sNT
        ]
        target = [
            1, 12,      # s1
            2, 12,      # s2
            3, 12,      # s3
            4, 12,      # s4
            5, 12,      # s5
            6, 12,      # s6
            7, 12,      # s7
            8, 11, 12,  # s8.1
            12,         # s8.2
            9, 10,      # s9.1
            13, 13,     # s9.2
            14, 15,     # pAA
            16, 17,     # sAA
            22, 23,     # sAA
            18, 19,     # pNT
            25,         # pNT
            20, 21,     # sNT
            22, 23      # sNT
        ]
        values = [
            s1, s1disc,         # s1
            s2, s2disc,         # s2
            s3, s3disc,         # s3
            s4, s4disc,         # s4
            s5, s5disc,         # s5
            s6, s6disc,         # s6
            s7, s7disc,         # s7
            s8, s8R2, s8disc,   # s8.1
            s8R2,               # s8.2
            s9, s9rdnt,         # s9.1
            s9, s9rdnt,         # s9.2
            paaV, paaU,         # pAA
            saaV, saaU,         # sAA
            saaV, saaU,         # sAA
            pntV, pntU,         # pNT
            pntU,               # pNT
            sntV, sntU,         # sNT
            sntV, sntU          # sNT
        ]
        link = dict(source=source,target=target,value=value)
        data = go.Sankey(link=link)
        fig = go.Figure(data)
        fig.write_image(output[0])




