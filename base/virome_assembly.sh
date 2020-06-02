#!/bin/bash

# Prep output directories
mkdir -p ./assembly
mkdir -p ./assembly/stats
mkdir -p ./assembly/megahit_contigs
mkdir -p ./assembly/contig_dictionary

# Set variables
IN=./QC/step_7
OUT=./assembly
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

flye --subassemblies $OUT/contig_dictionary/all.mh.contigs_for_flye.fa -t 24 --meta --plasmids -o $OUT/contig_dictionary -g 1g

reformat.sh in=$OUT/contig_dictionary/assembly.fasta out=$OUT/contig_dictionary/contig_dictionary.fasta \
        ml=$MIN \
        ow=t;

# Step 3: Quantification by mapping
# Set file names
for i in $IN/*_R1.s7.out.fastq; do
        F=`basename $i _R1.s7.out.fastq`;

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

# Step 4: Create contig count table (covtable.all, this is similar to an otu_table or seqtable))
DIR=./assembly/contig_dictionary

# Extract contig IDs
grep ">" $DIR/assembly.fasta | sed 's/>//' > $DIR/contig.ids;

# Combine frags with coverage percent and calculate fragments per length (frags/contig length)
# Removes alignment statistic from mappings across < 90% of the length of the contig

# Output a table of counts
for i in $DIR/*.rpkm; do

        F=`basename $i .rpkm`;

        tail -n+5 $DIR/"$F".rpkm | cut -f7 > $DIR/"$F".frags;

        paste $DIR/"$F".covstats $DIR/"$F".frags > $DIR/"$F".tmp;

        tail -n+2 $DIR/"$F".tmp  | awk '{ if ($2 >= 90) { print$1"\t"$11 } }' > $DIR/"$F".frag.counts;

       sed -i "1iID\t$F" $DIR/"$F".frag.counts;
done

# Output a table of Frags / Length
for i in $DIR/*.rpkm; do

	F=`basename $i .rpkm`;

	tail -n+2 $DIR/"$F".tmp  | awk '{print$1"\t"$5"\t"($11 / $3)}' | awk '{ if ($2 >= 90) { print$1"\t"$3 } }' > $DIR/"$F".length.normalized;

	sed -i "1iID\t$F" $DIR/"$F".length.normalized;
done

# Remove temporary files
rm $DIR/*.tmp;
rm $DIR/*.frags;

# Output a table of FPKMs
for i in $DIR/*.rpkm; do

        F=`basename $i .rpkm`;

        tail -n+5 $DIR/"$F".rpkm | cut -f8 > $DIR/"$F".fpkm;

        paste $DIR/"$F".covstats $DIR/"$F".fpkm > $DIR/"$F".tmp;

       tail -n+2 $DIR/"$F".tmp  | awk '{ if ($5 >= 90) { print$1"\t"$11 } }'  > $DIR/"$F".fpkm.filt;

       sed -i "1iID\t$F" $DIR/"$F".fpkm.filt;
done

# Remove temporary files
rm $DIR/*.tmp;
rm $DIR/*.fpkm;

#### EXTRA CODE TO DELETE LATER LIES BELOW!####

#for i in $DIR/*.covstats; do
#        F=`basename $i .covstats`;

#        awk '{ if ($5 >= 90) { print } }' $DIR/"$F".covstats > $DIR/"$F"_filtered.covstats

# Extract the average coverage from individual sample mappings
#for i in $DIR/*.rpkm; do
#	F=`basename $i .rpkm`;
#
#	# Extract contig IDs and average fold coverage
#	tail -n+6 $DIR/"$F".rpkm | cut -f1,2,7 | awk '{print($3 / $2)}' > $DIR/"$F".length.norm;
#
#	sed -i "1i $F" $DIR/"$F".length.norm;
#done

# Combine individual average coverage files
#paste $DIR/*length.norm > $DIR/all.length.coverage;

# Extract contig IDs and join with average coverages
#grep ">" $DIR/assembly.fasta | sed 's/>//' > $DIR/contig.ids;
#sed -i '1iID' $DIR/contig.ids;
#paste $DIR/contig.ids $DIR/all.length.coverage > $DIR/covtable.all;

# Remove temporary files
#rm $DIR/*.length.norm;
#rm $DIR/contig.ids;

# Filter per sample low percent coverage contig mappings

#for i in $DIR/*.covstats; do
#	F=`basename $i .covstats`;

#	awk '{ if ($5 >= 90) { print } }' $DIR/"$F".covstats > $DIR/"$F"_filtered.covstats

