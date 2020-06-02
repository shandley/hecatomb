#!/bin/bash

# One-liner to rename fastq files

for f in *.gz; do mv "$f" $(echo $f | sed 's/_001.fastq.gz/.fastq.gz/g'); done

