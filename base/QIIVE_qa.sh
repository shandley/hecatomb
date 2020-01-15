#!/bin/bash
# Set i/o paramaters
mkdir -p ./QC/stats/
IN=./QC
QC_DIRS=($IN/step_1 $IN/step_2 $IN/step_3 $IN/step_4 $IN/step_5 $IN/step_7 $IN/step_8)
OUT=./QC/stats

# Remove any zero byte files from QC folder
	find ./QC/ -maxdepth 2 -type f -size 0 > $OUT/zero_byte_files.txt;
	find ./QC/ -maxdepth 2 -type f -size 0 -exec rm -f {} \;-size 0 -exec rm -f {} \;

# Clumpify stats
	# Extract file size information
		ls -1 *.fastq.gz > ./clumped/clumped_filenames
		ls -l *.fastq.gz | awk -F ' ' '{ print$5 }' > ./clumped/base.fs
		ls -l ./clumped/*.fastq.clumped.gz | awk -F ' ' '{ print$5 }' > ./clumped/clumped.fs

	# Create combined file stats table
		paste -d '\t' ./clumped/clumped_filenames ./clumped/base.fs ./clumped/clumped.fs > ./clumped/clumping.stats
		sed -i '1isample_id\traw\tclumped' ./clumped/clumping.stats
		mv ./clumped/clumping.stats $OUT
		rm ./clumped/clumped_filenames ./clumped/base.fs ./clumped/clumped.fs

# Read count and length stats
	# Raw read stats
		seqkit stats -T -j 64 *R1* > $OUT/step_0_R1_readstats.txt;
		seqkit stats -T -j 64 *R2* > $OUT/step_0_R2_readstats.txt;

	# Read stats for bbduk.sh steps
	for i in ${QC_DIRS[*]}; do

		F=`basename $i`;	

		seqkit stats -T -j 64 "$i"/*_R1.s*.out.fastq > $OUT/"$F"_R1_readstats.txt;
	        seqkit stats -T -j 64 "$i"/*_R2.s*.out.fastq > $OUT/"$F"_R2_readstats.txt;

	done

	# Read stats for mapping steps (bbmap.sh) steps 6 and 9
	# Step 6
	seqkit stats -T -j 64 ./QC/step_6/*_R1.s6.out.fastq > $OUT/step_6_R1_readstats.txt;
       	seqkit stats -T -j 64 ./QC/step_6/*_R2.s6.out.fastq > $OUT/step_6_R2_readstats.txt;

	# Step 9
	seqkit stats -T -j 64 ./QC/step_9/*_bacterial.fastq > $OUT/step_9_bacterial_readstats.txt;
        seqkit stats -T -j 64 ./QC/step_9/*_viral_amb.fastq > $OUT/step_9_viral_amb_readstats.txt;

	# Read stats from dereplication (cluster_count.sh)
	seqkit stats -T -j 64 ./QC/step_9/clustered/*.deduped.out.fastq > $OUT/deduped_readstats.txt;
	seqkit stats -T	-j 64 ./QC/step_9/clustered/*_reformated.fasta > $OUT/dereped_readstats.txt;

	# Compile read statistics

	# Read 1 counts & lengths
	for i in `(cd $OUT && ls *_R1_readstats.txt)`; do
		awk '{print$4}' $OUT/$i > $OUT/"$i"_R1_counts;

	        awk '{print$7}' $OUT/$i > $OUT/"$i"_R1_avg_length;

		awk '{print$6}' $OUT/$i > $OUT/"$i"_R1_min_length;

		awk '{print$8}'	$OUT/$i	> $OUT/"$i"_R1_max_length;

	done

	# Read 2 counts and lengths
	for i in `(cd $OUT && ls *_R2_readstats.txt)`; do
        	awk '{print$4}' $OUT/$i > $OUT/"$i"_R2_counts;

	        awk '{print$7}' $OUT/$i > $OUT/"$i"_R2_avg_length;

                awk '{print$6}'	$OUT/$i	> $OUT/"$i"_R2_min_length;

                awk '{print$8}' $OUT/$i > $OUT/"$i"_R2_max_length;

	done

	# Step 6 and 9 counts and lengths
		awk '{print$4}' $OUT/step_9_viral_amb_readstats.txt > $OUT/step_9_viral_amb_readstats.txt_counts;
		awk '{print$4}' $OUT/step_9_bacterial_readstats.txt > $OUT/step_9_bacterial_readstats.txt_counts;

                awk '{print$7}' $OUT/step_9_viral_amb_readstats.txt > $OUT/step_9_viral_amb_readstats.txt_avg_length;
                awk '{print$7}' $OUT/step_9_bacterial_readstats.txt > $OUT/step_9_bacterial_readstats.txt_avg_length;

                awk '{print$6}' $OUT/step_9_viral_amb_readstats.txt > $OUT/step_9_viral_amb_readstats.txt_min_length;
                awk '{print$6}' $OUT/step_9_bacterial_readstats.txt > $OUT/step_9_bacterial_readstats.txt_min_length;

       	       	awk '{print$8}' $OUT/step_9_viral_amb_readstats.txt > $OUT/step_9_viral_amb_readstats.txt_max_length;
                awk '{print$8}' $OUT/step_9_bacterial_readstats.txt > $OUT/step_9_bacterial_readstats.txt_max_length;

	# Dereplication (cluster_count.sh) counts and lengths
		awk '{print$4}' $OUT/deduped_readstats.txt > $OUT/deduped.txt_counts;
		awk '{print$7}' $OUT/deduped_readstats.txt > $OUT/deduped.txt_avg_length;
		awk '{print$6}' $OUT/deduped_readstats.txt > $OUT/deduped.txt_min_length;
		awk '{print$8}' $OUT/deduped_readstats.txt > $OUT/deduped.txt_max_length;

		awk '{print$4}' $OUT/dereped_readstats.txt > $OUT/dereped.txt_counts;
                awk '{print$7}' $OUT/dereped_readstats.txt > $OUT/dereped.txt_avg_length;
               	awk '{print$6}' $OUT/dereped_readstats.txt > $OUT/dereped.txt_min_length;
               	awk '{print$8}' $OUT/dereped_readstats.txt > $OUT/dereped.txt_max_length;

	# Create read_count output table
		# Extract filenames
		cut -f1 $OUT/step_0_R1_readstats.txt > $OUT/filenames;

		# Tabuluate Read 1 counts	
		paste -d "\t" $OUT/*_R1_counts $OUT/step_9_viral_amb_readstats.txt_counts $OUT/step_9_bacterial_readstats.txt_counts $OUT/deduped.txt_counts $OUT/dereped.txt_counts | tail -n +2 > $OUT/R1_counts;
		sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R1_counts;
	        paste $OUT/filenames $OUT/R1_counts > $OUT/R1_readcount.table;

		# Tabluate Read 2 counts
		paste -d "\t" $OUT/*_R2_counts $OUT/step_9_viral_amb_readstats.txt_counts $OUT/step_9_bacterial_readstats.txt_counts $OUT/deduped.txt_counts $OUT/dereped.txt_counts | tail -n +2 > $OUT/R2_counts;
		sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R2_counts;
	        paste $OUT/filenames $OUT/R2_counts > $OUT/R2_readcount.table;

	# Create average read_length output table
		# Tabuluate Read 1 average length
        	        paste -d "\t" $OUT/*_R1_avg_length $OUT/step_9_viral_amb_readstats.txt_avg_length $OUT/step_9_bacterial_readstats.txt_avg_length $OUT/deduped.txt_avg_length $OUT/dereped.txt_avg_length | tail -n +2 > $OUT/R1_avg_length;
                	sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R1_avg_length;
	                paste $OUT/filenames $OUT/R1_avg_length > $OUT/R1_avg_length.table;

	       	# Tabuluate Read 2 average length
        	       	paste -d "\t" $OUT/*_R2_avg_length $OUT/step_9_viral_amb_readstats.txt_avg_length $OUT/step_9_bacterial_readstats.txt_avg_length $OUT/deduped.txt_avg_length $OUT/dereped.txt_avg_length | tail -n +2 > $OUT/R2_avg_length;
                	sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R2_avg_length;
                	paste $OUT/filenames $OUT/R2_avg_length > $OUT/R2_avg_length.table;

        # Create minimum read_length output table
		# Tabuluate Read 1 min lnegth
                        paste -d "\t" $OUT/*_R1_min_length $OUT/step_9_viral_amb_readstats.txt_min_length $OUT/step_9_bacterial_readstats.txt_min_length $OUT/deduped.txt_min_length $OUT/dereped.txt_min_length | tail -n +2 > $OUT/R1_min_length;
                        sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R1_min_length;
                        paste $OUT/filenames $OUT/R1_min_length > $OUT/R1_min_length.table;

                # Tabuluate Read 2 min length
                        paste -d "\t" $OUT/*_R2_min_length $OUT/step_9_viral_amb_readstats.txt_min_length $OUT/step_9_bacterial_readstats.txt_min_length $OUT/deduped.txt_min_length $OUT/dereped.txt_min_length | tail -n +2 > $OUT/R2_min_length;
                        sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R2_min_length;
                        paste $OUT/filenames $OUT/R2_min_length > $OUT/R2_min_length.table;

	# Create maximum read_length output table
		# Tabuluate Read 1 max length
                        paste -d "\t" $OUT/*_R1_max_length $OUT/step_9_viral_amb_readstats.txt_max_length $OUT/step_9_bacterial_readstats.txt_max_length $OUT/deduped.txt_max_length $OUT/dereped.txt_max_length | tail -n +2 > $OUT/R1_max_length;
                        sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R1_max_length;
                        paste $OUT/filenames $OUT/R1_max_length > $OUT/R1_max_length.table;

                # Tabuluate Read 2 max length
                        paste -d "\t" $OUT/*_R2_max_length $OUT/step_9_viral_amb_readstats.txt_max_length $OUT/step_9_bacterial_readstats.txt_max_length $OUT/deduped.txt_max_length $OUT/dereped.txt_max_length | tail -n +2 > $OUT/R2_max_length;
                        sed -i '1istep_0\tstep_1\tstep_2\tstep_3\tstep_4\tstep_5\tstep_6\tstep_7\tstep_8\tstep_9_viral_amb\tstep_9_bacterial\tdeduplication\tdereplication' $OUT/R2_max_length;
                        paste $OUT/filenames $OUT/R2_max_length > $OUT/R2_max_length.table;

	# Clean up

	rm $OUT/*_counts;
	rm $OUT/*_avg_length;
	rm $OUT/*_min_length;
	rm $OUT/*_max_length;


# Primer B evaluation
	for i in ./QC/step_1/*.stats; do

		NAME=`basename $i .stats`;

		grep -v "#" ./QC/step_1/"$NAME".stats | cut -f1,2 | sed "1iprimerB\t$NAME" > ./QC/step_1/"$NAME"_primerB_num_reads.out;

		grep -v "#" ./QC/step_1/"$NAME".stats | cut -f1,3 | sed "1iprimerB\t$NAME" | sed 's/%//' > ./QC/step_1/"$NAME"_primerB_perc_reads.out;

	done

	# Remove any zero byte files
		find ./QC/step_1/ -type f -size 0 -exec rm -f {} \;
	
	# Join primerB filed
	primerB_QC.R

# Contaminant (vector) stats
contaminant_stats.sh

# Bacterial contaminants
bacterial_contaminants.sh

