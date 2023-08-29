
import yaml
import collections


rule rawReadCounts:
    output:
        report(os.path.join(dir["out"]["temp"],"raw_read_counts.yaml"),
            caption = os.path.join("..", "report", "raw_read_counts.rst"),
            category = "Preprocessing")
    params:
        samples = samples
    run:
        out_counts = dict()
        for sample in params.samples["names"]:
            out_counts[sample]["s1_raw_reads_R1"] = file_len(params.samples["reads"][sample]["R1"])
            try:
                out_counts[sample]["s1_raw_reads_R2"] = file_len(params.samples["reads"][sample]["R2"])
            except KeyError:
                pass
        with open(output[0],"w") as stream:
            yaml.dump(out_counts,stream)


rule trimmedHostRemovedCounts:
    input:
        expand(os.path.join(dir["out"]["temp"], "p04", "{sample}_R1.all.fastq"),sample=samples["names"])
    output:
        report(os.path.join(dir["out"]["temp"], "host_removed_counts.yaml"),
            caption = os.path.join("..", "report", "host_removed_counts.rst"),
            category = "Preprocessing")
    params:
        samples = samples
    run:
        out_counts = dict()
        for sample in params.samples["names"]:
            out_counts[sample]["s2_host_removed_R1"] = file_len(os.path.join(dir["out"]["temp"], "p04", f"{sample}_R1.all.fastq"))
            try:
                out_counts[sample]["s2_host_removed_R2"] = file_len(os.path.join(dir["out"]["temp"],"p04", f"{sample}_R2.all.fastq"))
            except FileNotFoundError:
                pass
        with open(output[0],"w") as stream:
            yaml.dump(out_counts,stream)


rule proteinAnnotations:
    input:
        vir = os.path.join(dir["out"]["secondaryAA"],"AA_bigtable.tsv"),
        nonvir = os.path.join(dir["out"]["secondaryAA"],"AA_bigtable.nonviral.tsv"),
    output:
        report(os.path.join(dir["out"]["temp"], "protein_annotations.yaml"),
            caption = os.path.join("..", "report", "protein_annotations.rst"),
            category = "Read annotation")
    params:
        samples=samples
    run:
        out_counts = dict()
        for sample in params.samples["names"]:
            out_counts[sample]["s3_aa_viral_R1"] = 0
            out_counts[sample]["s3_aa_nonviral_R1"] = 0
        for l in stream_tsv(input.vir, skip_header=True):
            out_counts[l[1]]["s3_aa_viral_R1"] += int(l[2])
        for l in stream_tsv(input.nonvir, skip_header=True):
            out_counts[l[1]]["s3_aa_nonviral_R1"] += int(l[2])
        with open(output[0],"w") as stream:
            yaml.dump(out_counts,stream)


rule nucleotideAnnotations:
    input:
        vir = os.path.join(dir["out"]["secondaryNT"], "NT_bigtable.tsv"),
        nonvir = os.path.join(dir["out"]["secondaryNT"], "NT_bigtable.nonviral.tsv"),
    output:
        report(os.path.join(dir["out"]["temp"], "nucleotide_annotations.yaml"),
            caption = os.path.join("..", "report", "nucleotide_annotations.rst"),
            category = "Read annotation")
    params:
        samples = samples
    run:
        out_counts = dict()
        for sample in params.samples["names"]:
            out_counts[sample]["s4_nt_viral_R1"] = 0
            out_counts[sample]["s4_nt_nonviral_R1"] = 0
        for l in stream_tsv(input.vir, skip_header=True):
            out_counts[l[1]]["s4_nt_viral_R1"] += int(l[2])
        for l in stream_tsv(input.nonvir, skip_header=True):
            out_counts[l[1]]["s4_nt_nonviral_R1"] += int(l[2])
        with open(output[0],"w") as stream:
            yaml.dump(out_counts,stream)


rule mappedCounts:
    input:
        os.path.join(dir["out"]["results"], "sample_coverage.tsv")
    output:
        report(os.path.join(dir["out"]["temp"],"mapped_counts.yaml"),
            caption = os.path.join("..", "report", "mapped_counts.rst"),
            category = "Assembly")
    params:
        samples = samples
    run:
        out_counts = dict()
        for sample in params.samples["names"]:
            out_counts[sample]["s5_mapped"] = 0
        for l in stream_tsv(input[0], skip_header=True):
            if l[0] != "Sample":
                out_counts[l[0]]["s5_mapped"] += int(l[2])
        with open(output[0],"w") as stream:
            yaml.dump(out_counts,stream)


rule unclassifiedSeqs:
    input:
        aaVir = os.path.join(dir["out"]["secondaryAA"],"AA_bigtable.tsv"),
        aaNonvir = os.path.join(dir["out"]["secondaryAA"],"AA_bigtable.nonviral.tsv"),
        ntVir = os.path.join(dir["out"]["secondaryNT"],"NT_bigtable.tsv"),
        ntNonvir = os.path.join(dir["out"]["secondaryNT"],"NT_bigtable.nonviral.tsv"),
        fa = os.path.join(dir["out"]["results"],"seqtable.fasta")
    output:
        report(os.path.join(dir["out"]["results"],"seqtable.unclassified.fasta"),
               caption = os.path.join("..", "report", "unclassified_seqtable.rst"),
            category="Read annotation")
    run:
        classSeq = {}
        for file in [input.aaVir, input.aaNonvir, input.ntVir, input.ntNonvir]:
            for line in stream_tsv(file, skip_header=True):
                classSeq[line[0]] = 1
        with open(output[0], "w") as out_fh:
            with open(input.fa, "r") as in_fh:
                for line in in_fh:
                    if line.startswith(">"):
                        id = line.strip().replace(">","")
                        seq = in_fh.readline().strip()
                        try:
                            classSeq[id]
                        except KeyError:
                            out_fh.write(f">{id}\n{seq}\n")
                    else:
                        sys.stderr.write(f"malformed {input.fa} file? expecting {line} to be fasta header, complain to Mike")
                        exit(1)


rule summaryTable:
    input:
        os.path.join(dir["out"]["temp"], "raw_read_counts.yaml"),
        os.path.join(dir["out"]["temp"], "host_removed_counts.yaml"),
        os.path.join(dir["out"]["temp"], "protein_annotations.yaml"),
        os.path.join(dir["out"]["temp"], "nucleotide_annotations.yaml"),
        os.path.join(dir["out"]["temp"], "mapped_counts.yaml"),
    output:
        report(os.path.join(dir["out"]["results"], "summary.tsv"),
               caption = os.path.join("..", "report", "summary_table.rst"),
            category="Read annotation")
    run:
        def recursive_merge(dict, merge_dict):
            def _update(d, u):
                for (key, value) in u.items():
                    if isinstance(value, collections.abc.Mapping):
                        d[key] = _update(d.get(key,{}),value)
                    else:
                        d[key] = value
                return d
            _update(dict, merge_dict)
        summary = {}
        for file in input:
            with open(file, "r") as stream:
                dict = yaml.safe_load(stream)
                recursive_merge(summary, dict)
        rows = list(summary.keys())
        cols = set()
        for row in rows:
            for col in summary[row]:
                cols.add(col)
        cols = list(cols)
        cols.sort()
        with open(output[0], "w") as out_fh:
            out_fh.write("\t".join(["sampleID"] + cols) + "\n")
            for row in rows:
                out_fh.write(row)
                for col in cols:
                    out_fh.write("\t" + str(summary[row][col]))
                out_fh.write("\n")
