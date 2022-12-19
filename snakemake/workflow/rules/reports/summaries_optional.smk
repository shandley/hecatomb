rule rawReadCounts:
    output:
        report(os.path.join(dir.out.temp,"Step00_counts.tsv"),
            caption = "../report/step00.rst",
            category = "Preprocessing")
    run:
        collect_start_counts(samples, output[0])

rule hostRemovedCounts:
    input:
        expand(os.path.join(dir.out.temp, "p04", "{sample}_R1.all.fastq"),sample=samples.names)
    output:
        report(os.path.join(dir.out.temp, "Step01_counts.tsv"),
            caption = "../report/step01.rst",
            category = "Preprocessing")
    run:
        collect_counts(os.path.join(dir.out.temp, "p04"), ".all.fastq", "Host_removed", output[0])

rule hostRemovedCountsLongreads:
    input:
        expand(os.path.join(dir.out.temp, "p01", "{sample}_R1.all.fasta"), sample=samples.names)
    output:
        report(os.path.join(dir.out.temp, "Step01_counts.tsv"),
            caption = "../report/step01.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        for sample in samples.names:
            R1c = file_len(os.path.join(dir.out.temp, "p01", f"{sample}_R1.all.fasta")) / 2
            outfh.write(f'{sample}\tHost_removed\tR1\t{R1c}\n')
        outfh.close()

rule clusteredCounts:
    input:
        expand(os.path.join(dir.out.temp, "p05", "{sample}_R1_rep_seq.fasta"), sample=samples.names)
    output:
        report(os.path.join(dir.out.temp, "Step02_counts.tsv"),
            caption = "../report/step02.rst",
            category = "Preprocessing")
    run:
        outfh = open(output[0],'w')
        for sample in samples.names:
            R1c = file_len(os.path.join(dir.out.temp, "p05", f"{sample}_R1_rep_seq.fasta")) / 2
            outfh.write(f'{sample}\tR1_Clustered\tR1\t{R1c}\n')
        outfh.close()

# rule step01_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p01", "{sample}_{rn}.s1.out.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step01_counts.tsv"),
#             caption = "../report/step01.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p01"), ".s1.out.fastq", "Step_01", output[0])
#
# rule step02_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p02", "{sample}_{rn}.s2.out.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step02_counts.tsv"),
#             caption = "../report/step02.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p02"), ".s2.out.fastq", "Step_02", output[0])
#
# rule step03_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p03", "{sample}_{rn}.s3.out.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step03_counts.tsv"),
#             caption = "../report/step03.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p03"), ".s3.out.fastq", "Step_03", output[0])
#
# rule step04_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p04", "{sample}_{rn}.s4.out.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step04_counts.tsv"),
#             caption = "../report/step04.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p04"), ".s4.out.fastq", "Step_04", output[0])
#
# rule step05_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p05", "{sample}_{rn}.s5.out.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step05_counts.tsv"),
#             caption = "../report/step05.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p05"), ".s5.out.fastq", "Step_05", output[0])
#
# rule step06_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p06", "{sample}_{rn}.s6.out.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step06_counts.tsv"),
#             caption = "../report/step06.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p06"), ".s6.out.fastq", "Step_06", output[0])
#
# rule step07_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p07", "{sample}_{rn}.all.fastq"), sample=samples.names, rn=['R1', 'R2'])
#     output:
#         report(os.path.join(SUMDIR, "Step07_counts.tsv"),
#             caption = "../report/step07.rst",
#             category = "Preprocessing")
#     run:
#         collect_counts(os.path.join(dir.out.temp, "p07"), ".all.fastq", "Step_07", output[0])
#
# rule step08_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p08", "{sample}_R1.deduped.out.fastq"), sample=samples.names)
#     output:
#         report(os.path.join(SUMDIR, "Step08_counts.tsv"),
#             caption = "../report/step08.rst",
#             category = "Preprocessing")
#     run:
#         outfh = open(output[0],'w')
#         for sample in samples.names:
#             R1c = file_len(os.path.join(dir.out.temp, "p08", f"{sample}_R1.deduped.out.fastq")) / 4
#             outfh.write(f'{sample}\tStep_08\tR1\t{R1c}\n')
#         outfh.close()
#
# rule step09_counts:
#     input:
#         expand(os.path.join(dir.out.temp, "p09", "{sample}_R1_rep_seq.fasta"), sample=samples.names)
#     output:
#         report(os.path.join(SUMDIR, "Step09_counts.tsv"),
#             caption = "../report/step09.rst",
#             category = "Preprocessing")
#     run:
#         outfh = open(output[0],'w')
#         for sample in samples.names:
#             R1c = file_len(os.path.join(dir.out.temp, "p09", f"{sample}_R1_rep_seq.fasta")) / 2
#             outfh.write(f'{sample}\tStep_09\tR1\t{R1c}\n')
#         outfh.close()
#
#
### Collect rep seq counts following each search stage
rule step10_counts:
    input:
        classSeq = os.path.join(dir.out.primaryAA,"MMSEQS_AA_PRIMARY_classified.fasta"),
        unclassSeq = os.path.join(dir.out.primaryAA,"MMSEQS_AA_PRIMARY_unclassified.fasta")
    output:
        report(os.path.join(dir.out.temp, "Step10_counts.tsv"),
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
        tsv = os.path.join(dir.out.secondaryAA,"AA_bigtable.tsv")
    output:
        report(os.path.join(dir.out.temp, "Step11_counts.tsv"),
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
        classSeq = os.path.join(dir.out.primaryNT,"classified_seqs.fasta"),
        unclassSeq = os.path.join(dir.out.primaryNT,"unclassified_seqs.fasta")
    output:
        report(os.path.join(dir.out.temp,"Step12_counts.tsv"),
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
        os.path.join(dir.out.secondaryNT,"NT_bigtable.tsv")
    output:
        report(os.path.join(dir.out.temp,"Step13_counts.tsv"),
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


# rule sankey_diagram:
#     input:
#         os.path.join(SUMDIR,"Step00_counts.tsv"),
#         os.path.join(SUMDIR,"Step01_counts.tsv"),
#         os.path.join(SUMDIR,"Step02_counts.tsv"),
#         os.path.join(SUMDIR,"Step03_counts.tsv"),
#         os.path.join(SUMDIR,"Step04_counts.tsv"),
#         os.path.join(SUMDIR,"Step05_counts.tsv"),
#         os.path.join(SUMDIR,"Step06_counts.tsv"),
#         os.path.join(SUMDIR,"Step07_counts.tsv"),
#         os.path.join(SUMDIR,"Step08_counts.tsv"),
#         os.path.join(SUMDIR,"Step09_counts.tsv"),
#         os.path.join(SUMDIR,"Step10_counts.tsv"),
#         os.path.join(SUMDIR,"Step11_counts.tsv"),
#         os.path.join(SUMDIR,"Step12_counts.tsv"),
#         os.path.join(SUMDIR,"Step13_counts.tsv")
#     output:
#         report(os.path.join(SUMDIR, "Sankey.svg"),
#             caption = "../report/sankey.rst",
#             category = "Summary")
#     conda:
#         os.path.join('..', 'envs', 'plotly.yaml')
#     log:
#         os.path.join(dir.out.stderr, 'sankey_diagram.log')
#     script:
#         os.path.join('..', 'scripts', 'sankey.py')