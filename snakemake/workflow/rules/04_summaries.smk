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
        outfh.write('sample\tstep\tR1\n')
        for sample in SAMPLES:
            R1c = file_len(os.path.join(QC, "CLUSTERED", f"{sample}_R1.deduped.out.fastq")) / 4
            outfh.write(f'{sample}\tStep_08\t{R1c}\n')

rule step09_counts:
    input:
        expand(os.path.join(QC, "CLUSTERED", "{sample}_R1_rep_seq.fasta"), sample=SAMPLES)
    output:
        report(os.path.join(SUMDIR, "Step09_counts.tsv"),
            caption = "../report/step09.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        outfh.write('sample\tstep\tR1_rep\n')
        for sample in SAMPLES:
            R1c = file_len(os.path.join(QC, "CLUSTERED", f"{sample}_R1_rep_seq.fasta")) / 2
            outfh.write(f'{sample}\tStep_09\t{R1c}\n')


















