rule tax_level_counts:
    """Generate a long table of hit counts at different taxon levels
    
    (excluding at the species level as that is essentially the bigtable.tsv file)
    """
    input:
        os.path.join(RESULTS, "bigtable.tsv")
    output:
        report(os.path.join(SUMDIR, "taxonLevelCounts.tsv"),
            caption = "../report/tax_level_counts.rst",
            category = "Step 1")
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    run:
        # the index position for each tax level in the bigtable.tsv - off by one because we capture the range as l[21:i]
        tsvIndex = {'kingdom': 22, 'phylum': 23, 'class': 24, 'order': 25, 'family': 26, 'genus': 27}
        idxStart = 21
        short = {'22':'k', '23':'p', '24':'c', '25':'o', '26':'f', '27':'g'}
        out = open(output[0], 'w')
        out.write('sampleID\ttaxonLevel\ttaxonPath\tcount\n')
        # re-read the file for each sample to keep the memory happy - this is probably not necessary
        for sample in SAMPLES:
            counts = {}         # counts[taxlevel][taxname] = int
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
                        except KeyError:
                            counts[t] = {}
                        taxPath = []
                        for o in range(idxStart,i):
                            taxPath.append(f'{short[str(o)]}:{l[o]}')
                        outPath = ','.join(taxPath)
                        try:
                            counts[t][outPath] += int(l[2])
                        except KeyError:
                            counts[t][outPath] = int(l[2])
            infh.close()
            for taxLevel in counts.keys():
                for taxName in counts[taxLevel].keys():
                    out.write(f'{sample}\t{taxLevel}\t{taxName}\t{counts[taxLevel][taxName]}\n')
        out.close()

