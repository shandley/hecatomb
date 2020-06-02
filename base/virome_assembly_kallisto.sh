#!/bin/bash

# Prep output directories
mkdir -p ./assembly
mkdir -p ./assembly/stats
mkdir -p ./assembly/megahit_contigs
mkdir -p ./assembly/contig_dictionary

# Set variables
IN=./QC/step_7
OUT=./assembly
DIR=./assembly/contig_dictionary
MIN=1000 # Minimum contig length

# Set file names
for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

# Step 0: Tabulate some kmer statistics
bbcountunique.sh in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq \
	interval=2500 \
	out=$OUT/stats/"$F"_uniq_kmer_stats.txt \
	ow=t;

# Step 1: Digital Normalization
bbnorm.sh in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq extra=$IN/"$F"_singletons_R1.out.fastq,$IN/"$F"_singletons_R2.out.fastq \
	out=$OUT/"$F"_R1.norm.out.fastq out2=$OUT/"$F"_R2.norm.out.fastq outt=$OUT/"$F"_tossed.norm.fastq \
	target=20 mindepth=2 \
	hist=$OUT/"$F"_norm.hist \
	ow=t;

# Step 2: Assembly
megahit -1 $OUT/"$F"_R1.norm.out.fastq -2 $OUT/"$F"_R2.norm.out.fastq -r $IN/"$F"_singletons_R1.out.fastq,$IN/"$F"_singletons_R2.out.fastq -o $OUT/"$F"_megahit_out \
	--out-prefix "$F".mh;

done;

find . -type f -name "*.mh.contigs.fa" -exec cp -rfp {} $OUT/megahit_contigs \;

# Step 3: Contig assembly to create contig dictionary (flye)
cat $OUT/megahit_contigs/*.mh.contigs.fa > $OUT/contig_dictionary/all.mh.contigs.fa;

reformat.sh in=$OUT/contig_dictionary/all.mh.contigs.fa out=$OUT/contig_dictionary/all.mh.contigs_trunc.fa \
	ml=$MIN \
	ow=t;

rename.sh in=$OUT/contig_dictionary/all.mh.contigs_trunc.fa out=$OUT/contig_dictionary/all.mh.contigs_for_flye.fa

# Adjust path or bashrc accordingly. Was having some install issues with flye when developing the script
flye --subassemblies $OUT/contig_dictionary/all.mh.contigs_for_flye.fa -t 24 --meta --plasmids -o $OUT/contig_dictionary -g 1g

reformat.sh in=$OUT/contig_dictionary/assembly.fasta out=$OUT/contig_dictionary/contig_dictionary.fasta \
        ml=$MIN \
        ow=t;

# Step 3: Calculate contig coverage percent
# Set file names
for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

# Use bbmap to calculate verious stats that kallisto does not
bbmap.sh ref=$OUT/contig_dictionary/contig_dictionary.fasta in=$IN/"$F"_R1.s7.out.fastq in2=$IN/"$F"_R2.s7.out.fastq \
	out=$OUT/contig_dictionary/"$F".aln.sam.gz \
        kfilter=22 subfilter=15 maxindel=80 \
	ambiguous=random \
	physcov=t \
	covstats=$OUT/contig_dictionary/"$F".covstats \
	rpkm=$OUT/contig_dictionary/"$F".rpkm \
	statsfile=$OUT/contig_dictionary/"$F".statsfile \
	scafstats=$OUT/contig_dictionary/"$F".scafstats \
        ow=t;
done

# Calculate TPM with kallisto
# Create kallisto index
kallisto index -i $OUT/contig_dictionary/contig_dictionary.idx $OUT/contig_dictionary/contig_dictionary.fasta

for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

kallisto quant -i $OUT/contig_dictionary/contig_dictionary.idx -b 100 \
	$IN/"$F"_R1.s7.out.fastq $IN/"$F"_R2.s7.out.fastq \
	-o $DIR/"$F"_kallisto;

done

# Combine contig covereages from bbmap with kallisto results

# Extract contig IDs
grep ">" $DIR/assembly.fasta | sed 's/>//' > $DIR/contig.ids;

# Combine frags with coverage percent and calculate fragments per length (frags/contig length)
# Removes alignment statistic from mappings across < 90% of the length of the contig
for i in $DIR/*.covstats; do
	
	F=`basename $i .covstats`;

	tail -n+2 $DIR/"$F"_kallisto/abundance.tsv | cut -f5 > $DIR/"$F"_kallisto/"$F".temp;

	sed -i "1i$F" $DIR/"$F"_kallisto/"$F".temp;

	paste $DIR/"$F".covstats $DIR/"$F"_kallisto/"$F".temp > $DIR/"$F".temp;

	tail -n+2 $DIR/"$F".temp | awk '{ if ($5 >= 90) { print$1"\t"$11 } }' > $DIR/"$F".filt.tpm;

	sed -i "1iID\t$F" $DIR/"$F".filt.tpm;

done

rm $DIR/*.temp
